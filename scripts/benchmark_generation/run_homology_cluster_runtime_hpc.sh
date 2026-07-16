#!/usr/bin/env bash
# Scratch-first runtime orchestration shared by the pilot and six-task HPC entrypoints.

set -euo pipefail

git_in_dir() {
    local directory="$1"
    shift
    (cd "$directory" && git "$@")
}

RUNTIME_KIND="${HOMOLOGY_RUNTIME_KIND:-}"
case "$RUNTIME_KIND" in
    pilot|array) ;;
    *) echo "HOMOLOGY_RUNTIME_KIND must be pilot or array" >&2; exit 2 ;;
esac

TASK_ID="${SGE_TASK_ID:-}"
case "$TASK_ID" in
    1) IDENTITY=30 ;;
    2) IDENTITY=25 ;;
    3) IDENTITY=20 ;;
    4) IDENTITY=15 ;;
    5) IDENTITY=10 ;;
    6) IDENTITY=5 ;;
    *) echo "SGE_TASK_ID must be one integer from 1 through 6" >&2; exit 2 ;;
esac
if [[ "$RUNTIME_KIND" == "pilot" && "$TASK_ID" != "1" ]]; then
    echo "The pilot entrypoint is locked to task 1 / 30% identity" >&2
    exit 2
fi

UNIPROT_SOURCE_SCOPE="${UNIPROT_SOURCE_SCOPE:-sprot-and-trembl}"
case "$UNIPROT_SOURCE_SCOPE" in
    sprot-only|trembl-only|sprot-and-trembl) ;;
    *) echo "Unsupported UNIPROT_SOURCE_SCOPE=$UNIPROT_SOURCE_SCOPE" >&2; exit 2 ;;
esac

SPLIT_POLICY="${SPLIT_POLICY:-sequence-balanced}"
TRAINING_POPULATION="${TRAINING_POPULATION:-annotated-only}"
SEED="${SEED:-0}"
MIN_COUNT="${MIN_COUNT:-50}"
[[ "$SPLIT_POLICY" == "sequence-balanced" || "$SPLIT_POLICY" == "cluster-count-random" ]] || {
    echo "Unsupported SPLIT_POLICY=$SPLIT_POLICY" >&2
    exit 2
}
[[ "$TRAINING_POPULATION" == "annotated-only" ]] || {
    echo "Only TRAINING_POPULATION=annotated-only is implemented" >&2
    exit 2
}
[[ "$SEED" =~ ^[0-9]+$ ]] || { echo "SEED must be a non-negative integer" >&2; exit 2; }
[[ "$MIN_COUNT" =~ ^[1-9][0-9]*$ ]] || { echo "MIN_COUNT must be positive" >&2; exit 2; }

TEST_MODE="${HOMOLOGY_RUNTIME_TEST_MODE:-0}"
case "$TEST_MODE" in 0|1) ;; *) echo "HOMOLOGY_RUNTIME_TEST_MODE must be 0 or 1" >&2; exit 2 ;; esac
if [[ "$TEST_MODE" != "1" ]]; then
    for name in HOMOLOGY_RUNTIME_TEST_BUILD_COMMAND HOMOLOGY_RUNTIME_TEST_COPY_COMMAND; do
        [[ -z "${!name:-}" ]] || {
            echo "$name is test-only; HOMOLOGY_RUNTIME_TEST_MODE is not enabled" >&2
            exit 2
        }
    done
fi

SUBMISSION_ROOT="${FRAMEWORK_SOURCE_ROOT:-${SGE_O_WORKDIR:-$PWD}}"
[[ "$SUBMISSION_ROOT" == /* && -d "$SUBMISSION_ROOT" ]] || {
    echo "FRAMEWORK_SOURCE_ROOT/SGE_O_WORKDIR must identify an absolute framework checkout" >&2
    exit 2
}
SUBMISSION_ROOT="$(cd "$SUBMISSION_ROOT" && pwd -P)"
FRAMEWORK_REPO_URL="${FRAMEWORK_REPO_URL:-https://github.com/nadroj0-0/Protein-Benchmark-Framework-Dissertation.git}"
FRAMEWORK_REVISION="${FRAMEWORK_REVISION:-}"
if [[ "$TEST_MODE" != "1" && -z "$FRAMEWORK_REVISION" ]]; then
    FRAMEWORK_REVISION="$(git_in_dir "$SUBMISSION_ROOT" rev-parse HEAD)"
fi
if [[ "$TEST_MODE" == "1" && -z "$FRAMEWORK_REVISION" ]]; then
    FRAMEWORK_REVISION="$(printf 'a%.0s' {1..40})"
fi
[[ "$FRAMEWORK_REVISION" =~ ^[0-9a-f]{40}$ ]] || {
    echo "FRAMEWORK_REVISION must be exactly 40 lowercase hexadecimal characters" >&2
    exit 2
}

JOB_KEY="${JOB_ID:-local}"
RUN_ID="${RUN_ID:-runtime-${JOB_KEY}}"
for name in JOB_KEY RUN_ID; do
    value="${!name}"
    [[ "$value" =~ ^[A-Za-z0-9._-]+$ && "$value" =~ [A-Za-z0-9] ]] || {
        echo "$name contains unsafe path characters" >&2
        exit 2
    }
done

WORK_BASE="${WORK_BASE:-/scratch0}"
[[ "$WORK_BASE" == /* && "$WORK_BASE" != "/" && -d "$WORK_BASE" ]] || {
    echo "WORK_BASE must be an existing absolute non-root directory" >&2
    exit 2
}
WORK_BASE="$(cd "$WORK_BASE" && pwd -P)"
[[ "$WORK_BASE" != "/" ]] || { echo "Resolved WORK_BASE must not be /" >&2; exit 2; }
WORK="$WORK_BASE/homology_runtime_${JOB_KEY}_${TASK_ID}_${IDENTITY}_${RUN_ID}"
[[ ! -e "$WORK" ]] || { echo "Refusing to reuse scratch path: $WORK" >&2; exit 2; }
mkdir -p "$WORK/artifacts/logs" "$WORK/inputs" "$WORK/tools" "$WORK/tmp"
touch "$WORK/.homology-runtime-owned"
export TMPDIR="$WORK/tmp"
export TMP="$WORK/tmp"
export TEMP="$WORK/tmp"

ARTIFACTS="$WORK/artifacts"
INPUT_ROOT="$WORK/inputs"
TOOLS_ROOT="$WORK/tools"
FRAMEWORK_DIR="$WORK/Protein-Benchmark-Framework-Dissertation"
RESULTS_ROOT="${RESULTS_ROOT:-$HOME/homology_cluster_benchmark_results}"
[[ "$RESULTS_ROOT" == /* && "$RESULTS_ROOT" != "/" ]] || {
    echo "RESULTS_ROOT must be an absolute non-root path" >&2
    exit 2
}
REVISION_TAG="${FRAMEWORK_REVISION:0:12}"
FINAL_ROOT="$RESULTS_ROOT/runtime_${RUNTIME_KIND}/source_${UNIPROT_SOURCE_SCOPE}/framework_${REVISION_TAG}/run_${RUN_ID}/job_${JOB_KEY}/task_${TASK_ID}_identity_${IDENTITY}"
PARTIAL_ROOT="${FINAL_ROOT}.partial-${JOB_KEY}-${TASK_ID}"
LOG_FILE="$ARTIFACTS/logs/runtime.log"
LOG_PIPE="$WORK/.runtime-log.pipe"
TEE_PID=""
DISK_MONITOR_PID=""
DISK_MONITOR_INTERVAL_SECONDS="${DISK_MONITOR_INTERVAL_SECONDS:-120}"
exec 3>&1 4>&2

[[ "$DISK_MONITOR_INTERVAL_SECONDS" =~ ^[1-9][0-9]*$ ]] || {
    echo "DISK_MONITOR_INTERVAL_SECONDS must be a positive integer" >&2
    exit 2
}

DISK_USAGE_LOG="$ARTIFACTS/logs/disk_usage.tsv"
DISK_BREAKDOWN_LOG="$ARTIFACTS/logs/disk_usage_by_path.tsv"

path_kib() {
    local path="$1"
    if [[ -e "$path" ]]; then
        du -sk "$path" 2>/dev/null | awk 'NR==1 {print $1}'
    else
        printf '0\n'
    fi
}

record_disk_usage() {
    local stage="$1"
    local work_kib free_kib used_kib
    work_kib="$(path_kib "$WORK")"
    read -r used_kib free_kib < <(df -Pk "$WORK" | awk 'NR==2 {print $3, $4}')
    printf '%s\t%s\t%s\t%s\t%s\n' \
        "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$stage" \
        "$((work_kib * 1024))" "$((used_kib * 1024))" "$((free_kib * 1024))" \
        >> "$DISK_USAGE_LOG"
}

record_disk_breakdown() {
    local stage="$1"
    local label path kib
    for label in inputs tools tmp artifacts framework; do
        case "$label" in
            inputs) path="$INPUT_ROOT" ;;
            tools) path="$TOOLS_ROOT" ;;
            tmp) path="$WORK/tmp" ;;
            artifacts) path="$ARTIFACTS" ;;
            framework) path="$FRAMEWORK_DIR" ;;
        esac
        kib="$(path_kib "$path")"
        printf '%s\t%s\t%s\t%s\n' \
            "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$stage" "$label" "$((kib * 1024))" \
            >> "$DISK_BREAKDOWN_LOG"
    done
}

checkpoint_disk_usage() {
    local stage="$1"
    record_disk_usage "$stage"
    record_disk_breakdown "$stage"
}

start_disk_monitor() {
    (
        sleep_pid=""
        trap '[[ -z "$sleep_pid" ]] || kill "$sleep_pid" 2>/dev/null; exit 0' TERM INT
        while true; do
            record_disk_usage periodic
            sleep "$DISK_MONITOR_INTERVAL_SECONDS" &
            sleep_pid=$!
            wait "$sleep_pid"
            sleep_pid=""
        done
    ) &
    DISK_MONITOR_PID=$!
}

stop_disk_monitor() {
    if [[ -n "$DISK_MONITOR_PID" ]]; then
        kill "$DISK_MONITOR_PID" 2>/dev/null || true
        wait "$DISK_MONITOR_PID" 2>/dev/null || true
        DISK_MONITOR_PID=""
    fi
}

summarize_disk_usage() {
    local summary="$ARTIFACTS/logs/disk_usage_summary.tsv"
    awk -F '\t' -v interval="$DISK_MONITOR_INTERVAL_SECONDS" '
        NR == 1 { next }
        {
            samples++
            if ($3 > peak || samples == 1) { peak=$3; peak_time=$1; peak_stage=$2 }
            if ($5 < min_free || samples == 1) { min_free=$5 }
            if (samples == 1) { first=$3 }
            last=$3
        }
        END {
            print "metric\tvalue"
            print "sample_count\t" samples
            print "monitor_interval_seconds\t" interval
            print "first_work_bytes\t" first
            print "final_work_bytes\t" last
            print "peak_work_bytes\t" peak
            print "peak_work_gib\t" peak / (1024 * 1024 * 1024)
            print "peak_timestamp_utc\t" peak_time
            print "peak_sample_stage\t" peak_stage
            print "minimum_filesystem_free_bytes\t" min_free
            print "measurement_scope\tjob-owned_work_directory"
            print "measurement_method\tdu_-sk_periodic_plus_stage_checkpoints"
        }
    ' "$DISK_USAGE_LOG" > "$summary"
}

copy_back() {
    local copy_status=0
    mkdir -p "$(dirname "$FINAL_ROOT")" || return 1
    if [[ -e "$FINAL_ROOT" || -e "$PARTIAL_ROOT" ]]; then
        echo "Refusing to overwrite an existing final or partial result: $FINAL_ROOT" >&2
        return 1
    fi
    mkdir -p "$PARTIAL_ROOT" || return 1
    if [[ "$TEST_MODE" == "1" && -n "${HOMOLOGY_RUNTIME_TEST_COPY_COMMAND:-}" ]]; then
        "${HOMOLOGY_RUNTIME_TEST_COPY_COMMAND}" "$ARTIFACTS" "$PARTIAL_ROOT" || copy_status=$?
    else
        cp -a "$ARTIFACTS/." "$PARTIAL_ROOT/" || copy_status=$?
    fi
    if [[ "$copy_status" != "0" ]]; then
        echo "Copy-back failed with status $copy_status; removing incomplete home copy" >&2
        rm -rf -- "$PARTIAL_ROOT"
        return "$copy_status"
    fi
    printf '%s\n' "$FINAL_ROOT" > "$PARTIAL_ROOT/FINAL_RESULT_PATH.txt" || {
        rm -rf -- "$PARTIAL_ROOT"
        return 1
    }
    mv "$PARTIAL_ROOT" "$FINAL_ROOT" || {
        rm -rf -- "$PARTIAL_ROOT"
        return 1
    }
    return 0
}

cleanup() {
    local status=$?
    local copy_status=0
    trap - EXIT INT TERM
    set +e
    checkpoint_disk_usage cleanup-before-copy
    stop_disk_monitor
    summarize_disk_usage
    if [[ -n "$TEE_PID" ]]; then
        exec 1>&3 2>&4
        wait "$TEE_PID"
        rm -f -- "$LOG_PIPE"
        TEE_PID=""
    fi
    echo
    echo "Final job status before copy-back: $status"
    echo "Copying important outputs and logs to: $FINAL_ROOT"
    copy_back || copy_status=$?
    if [[ "$copy_status" != "0" ]]; then
        echo "Result copy-back failed. The Grid Engine output log is the retained diagnostic." >&2
        [[ "$status" != "0" ]] || status=74
    fi
    if [[ "$WORK" == "$WORK_BASE"/homology_runtime_* && -f "$WORK/.homology-runtime-owned" ]]; then
        echo "Removing job-owned scratch directory: $WORK"
        rm -rf -- "$WORK"
    else
        echo "Refusing unsafe scratch cleanup for $WORK" >&2
        [[ "$status" != "0" ]] || status=75
    fi
    if [[ "$copy_status" == "0" ]]; then
        echo "Results copied to: $FINAL_ROOT"
    fi
    exit "$status"
}
trap cleanup EXIT
trap 'exit 130' INT
trap 'exit 143' TERM

mkfifo "$LOG_PIPE"
tee -a "$LOG_FILE" < "$LOG_PIPE" >&3 &
TEE_PID=$!
exec > "$LOG_PIPE" 2>&1

printf 'timestamp_utc\tstage\twork_bytes\tfilesystem_used_bytes\tfilesystem_free_bytes\n' \
    > "$DISK_USAGE_LOG"
printf 'timestamp_utc\tstage\tpath_role\tallocated_bytes\n' > "$DISK_BREAKDOWN_LOG"
checkpoint_disk_usage scratch-created
start_disk_monitor

echo "Host             : $(hostname)"
echo "Job ID           : $JOB_KEY"
echo "Task ID          : $TASK_ID"
echo "Run kind         : $RUNTIME_KIND"
echo "Identity         : $IDENTITY%"
echo "UniProt scope    : $UNIPROT_SOURCE_SCOPE"
echo "Framework commit : $FRAMEWORK_REVISION"
echo "Scratch          : $WORK"
echo "Final output     : $FINAL_ROOT"
echo "Pilot prerequisite for array: no"

if [[ "$TEST_MODE" == "1" ]]; then
    [[ -n "${HOMOLOGY_RUNTIME_TEST_BUILD_COMMAND:-}" ]] || {
        echo "Test mode requires HOMOLOGY_RUNTIME_TEST_BUILD_COMMAND" >&2
        exit 2
    }
    "${HOMOLOGY_RUNTIME_TEST_BUILD_COMMAND}" "$ARTIFACTS"
    checkpoint_disk_usage test-build-complete
    exit 0
fi

for command in git wget tar gzip awk sha256sum md5sum; do
    command -v "$command" >/dev/null 2>&1 || { echo "Missing required command: $command" >&2; exit 1; }
done

CONDA_EXE="${CONDA_EXE:-/share/apps/miniforge3_mamba/bin/conda}"
MMFP_ENV_DIR="${MMFP_ENV_DIR:-$HOME/.conda/envs/mmfp}"
[[ -x "$CONDA_EXE" && -d "$MMFP_ENV_DIR" ]] || {
    echo "The existing mmfp Conda environment is unavailable" >&2
    exit 1
}
eval "$("$CONDA_EXE" shell.bash hook)"
conda activate "$MMFP_ENV_DIR"
PYTHON_BIN="$MMFP_ENV_DIR/bin/python"
[[ -x "$PYTHON_BIN" ]] || { echo "mmfp Python is unavailable: $PYTHON_BIN" >&2; exit 1; }
echo "Using existing Conda environment: $CONDA_PREFIX"

echo "Cloning the pinned framework revision into scratch"
git clone --quiet "$FRAMEWORK_REPO_URL" "$FRAMEWORK_DIR"
git_in_dir "$FRAMEWORK_DIR" checkout --quiet --detach "$FRAMEWORK_REVISION"
[[ "$(git_in_dir "$FRAMEWORK_DIR" rev-parse HEAD)" == "$FRAMEWORK_REVISION" ]] || {
    echo "Scratch checkout does not match FRAMEWORK_REVISION" >&2
    exit 1
}
[[ -z "$(git_in_dir "$FRAMEWORK_DIR" status --porcelain)" ]] || {
    echo "Scratch framework checkout is unexpectedly dirty" >&2
    exit 1
}

# Python runs inside the minimal MMFP Singularity image, which intentionally has
# no Git executable. Export the state already verified above so the builder can
# retain its production provenance gate without repeating Git inside the image.
export HOMOLOGY_HOST_GIT_VERIFIED_COMMIT="$FRAMEWORK_REVISION"
export HOMOLOGY_HOST_GIT_VERIFIED_CLEAN=1
export HOMOLOGY_HOST_GIT_VERIFIED_REPOSITORY="$FRAMEWORK_DIR"
export SINGULARITYENV_HOMOLOGY_HOST_GIT_VERIFIED_COMMIT="$FRAMEWORK_REVISION"
export SINGULARITYENV_HOMOLOGY_HOST_GIT_VERIFIED_CLEAN=1
export SINGULARITYENV_HOMOLOGY_HOST_GIT_VERIFIED_REPOSITORY="$FRAMEWORK_DIR"
checkpoint_disk_usage framework-cloned
# shellcheck source=../reproduction_common.sh
source "$FRAMEWORK_DIR/scripts/reproduction_common.sh"
validate_mmfp_env "$PYTHON_BIN"

MINIMUM_SCRATCH_GB="${MINIMUM_SCRATCH_GB:-275}"
[[ "$MINIMUM_SCRATCH_GB" =~ ^[1-9][0-9]*$ ]] || {
    echo "MINIMUM_SCRATCH_GB must be a positive integer" >&2
    exit 2
}
free_kb="$(df -Pk "$WORK" | awk 'NR==2 {print $4}')"
required_kb="$((MINIMUM_SCRATCH_GB * 1024 * 1024))"
if (( free_kb < required_kb )); then
    echo "Scratch has ${free_kb} KiB free; ${required_kb} KiB is required before downloads" >&2
    exit 1
fi

UNIREF90_URL="https://ftp.uniprot.org/pub/databases/uniprot/current_release/uniref/uniref90/uniref90.fasta.gz"
IDMAPPING_URL="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz"
SPROT_URL="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz"
TREMBL_URL="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.dat.gz"
GOA_URL="https://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz"
GOA_MD5_URL="https://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz.md5"
GO_OBO_URL="https://release.geneontology.org/2026-06-19/ontology/go-basic.obo"
UNIPROT_RELNOTES_URL="https://ftp.uniprot.org/pub/databases/uniprot/current_release/relnotes.txt"
GOA_RELEASES_URL="https://ftp.ebi.ac.uk/pub/databases/GO/goa/current_release_numbers.txt"

download_file() {
    local url="$1"
    local destination="$2"
    local partial="${destination}.part"
    echo "Downloading $url"
    wget --continue --tries=5 --timeout=60 --progress=dot:giga -O "$partial" "$url"
    [[ -s "$partial" ]] || { echo "Download is empty: $url" >&2; return 1; }
    mv "$partial" "$destination"
}

stage_or_download() {
    local role="$1"
    local supplied_path="$2"
    local url="$3"
    local destination="$4"
    local expected_sha="$5"
    local acquisition
    mkdir -p "$(dirname "$destination")"
    if [[ -n "$supplied_path" ]]; then
        [[ -f "$supplied_path" && -s "$supplied_path" ]] || {
            echo "Provided $role path is missing or empty: $supplied_path" >&2
            return 1
        }
        echo "Staging provided $role into scratch: $supplied_path"
        cp -p "$supplied_path" "$destination"
        acquisition="provided-path-staged-to-scratch"
    else
        download_file "$url" "$destination"
        acquisition="downloaded-to-scratch"
    fi
    [[ -s "$destination" ]] || { echo "Staged $role is empty: $destination" >&2; return 1; }
    local observed_sha
    observed_sha="$(sha256sum "$destination" | awk '{print $1}')"
    if [[ -n "$expected_sha" && "$observed_sha" != "$expected_sha" ]]; then
        echo "$role SHA-256 mismatch: expected=$expected_sha observed=$observed_sha" >&2
        return 1
    fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$role" "$destination" "$url" "$acquisition" "$(stat -c '%s' "$destination")" "$observed_sha" \
        >> "$ARTIFACTS/logs/runtime_input_staging.tsv"
    checkpoint_disk_usage "input-${role}-staged"
}

printf 'role\tpath\tofficial_url\tacquisition\tsize_bytes\tsha256\n' \
    > "$ARTIFACTS/logs/runtime_input_staging.tsv"

needs_uniprot_download=0
for name in UNIREF90_FASTA IDMAPPING; do
    [[ -n "${!name:-}" ]] || needs_uniprot_download=1
done
if [[ "$UNIPROT_SOURCE_SCOPE" != "trembl-only" && -z "${UNIPROT_SPROT_SEQUENCES:-}" ]]; then
    needs_uniprot_download=1
fi
if [[ "$UNIPROT_SOURCE_SCOPE" != "sprot-only" && -z "${UNIPROT_TREMBL_SEQUENCES:-}" ]]; then
    needs_uniprot_download=1
fi
if [[ "$needs_uniprot_download" == "1" ]]; then
    download_file "$UNIPROT_RELNOTES_URL" "$ARTIFACTS/logs/uniprot_relnotes.txt"
    grep -Eq 'UniProt Release[[:space:]]+2026_02([^0-9]|$)' "$ARTIFACTS/logs/uniprot_relnotes.txt" || {
        echo "The UniProt current_release endpoint is no longer release 2026_02" >&2
        exit 1
    }
fi
if [[ -z "${GOA:-}" ]]; then
    download_file "$GOA_RELEASES_URL" "$ARTIFACTS/logs/goa_current_release_numbers.txt"
    awk '$1 == "uniprot" && $2 == "234" && $3 == "2026-06-17" {found=1} END {exit !found}' \
        "$ARTIFACTS/logs/goa_current_release_numbers.txt" || {
        echo "The GOA current endpoint is no longer UniProt-GOA release 234 (2026-06-17)" >&2
        exit 1
    }
fi

stage_or_download uniref90_fasta "${UNIREF90_FASTA:-}" "$UNIREF90_URL" \
    "$INPUT_ROOT/uniref90.fasta.gz" "${UNIREF90_FASTA_SHA256:-}"
stage_or_download idmapping "${IDMAPPING:-}" "$IDMAPPING_URL" \
    "$INPUT_ROOT/idmapping_selected.tab.gz" "${IDMAPPING_SHA256:-}"
if [[ "$UNIPROT_SOURCE_SCOPE" != "trembl-only" ]]; then
    stage_or_download uniprot_sprot_sequences "${UNIPROT_SPROT_SEQUENCES:-}" "$SPROT_URL" \
        "$INPUT_ROOT/uniprot_sprot.dat.gz" "${UNIPROT_SPROT_SEQUENCES_SHA256:-}"
fi
if [[ "$UNIPROT_SOURCE_SCOPE" != "sprot-only" ]]; then
    stage_or_download uniprot_trembl_sequences "${UNIPROT_TREMBL_SEQUENCES:-}" "$TREMBL_URL" \
        "$INPUT_ROOT/uniprot_trembl.dat.gz" "${UNIPROT_TREMBL_SEQUENCES_SHA256:-}"
fi
stage_or_download goa "${GOA:-}" "$GOA_URL" \
    "$INPUT_ROOT/goa_uniprot_all.gaf.234.gz" "${GOA_SHA256:-}"
stage_or_download go_obo "${GO_OBO:-}" "$GO_OBO_URL" \
    "$INPUT_ROOT/go-basic.obo" "${GO_OBO_SHA256:-}"

if [[ "$needs_uniprot_download" == "1" ]]; then
    download_file "$UNIPROT_RELNOTES_URL" "$ARTIFACTS/logs/uniprot_relnotes_after_download.txt"
    [[ "$(sha256sum "$ARTIFACTS/logs/uniprot_relnotes.txt" | awk '{print $1}')" == \
       "$(sha256sum "$ARTIFACTS/logs/uniprot_relnotes_after_download.txt" | awk '{print $1}')" ]] || {
        echo "UniProt current_release changed while this task was downloading inputs" >&2
        exit 1
    }
fi
if [[ -z "${GOA:-}" ]]; then
    download_file "$GOA_RELEASES_URL" "$ARTIFACTS/logs/goa_current_release_numbers_after_download.txt"
    [[ "$(awk '$1 == "uniprot" {print}' "$ARTIFACTS/logs/goa_current_release_numbers.txt")" == \
       "$(awk '$1 == "uniprot" {print}' "$ARTIFACTS/logs/goa_current_release_numbers_after_download.txt")" ]] || {
        echo "GOA current release metadata changed while this task was downloading inputs" >&2
        exit 1
    }
fi

if [[ -z "${GOA:-}" ]]; then
    download_file "$GOA_MD5_URL" "$ARTIFACTS/logs/goa_uniprot_all.gaf.gz.md5"
    expected_md5="$(awk 'NR==1 {print $1}' "$ARTIFACTS/logs/goa_uniprot_all.gaf.gz.md5")"
    observed_md5="$(md5sum "$INPUT_ROOT/goa_uniprot_all.gaf.234.gz" | awk '{print $1}')"
    [[ "$expected_md5" =~ ^[0-9a-fA-F]{32}$ && "${expected_md5,,}" == "$observed_md5" ]] || {
        echo "Downloaded GOA file does not match the official MD5 sidecar" >&2
        exit 1
    }
fi

MMSEQS_RELEASE_TAG="${MMSEQS_RELEASE_TAG:-18-8cc5c}"
EXPECTED_MMSEQS_VERSION="${EXPECTED_MMSEQS_VERSION:-$MMSEQS_RELEASE_TAG}"
EXPECTED_MMSEQS_BINARY_VERSION="${EXPECTED_MMSEQS_BINARY_VERSION:-8cc5ce367b5638c4306c2d7cfc652dd099a4643f}"
if [[ -n "${MMSEQS_BIN:-}" ]]; then
    [[ -x "$MMSEQS_BIN" ]] || { echo "MMSEQS_BIN is not executable: $MMSEQS_BIN" >&2; exit 1; }
else
    [[ "$EXPECTED_MMSEQS_VERSION" == "$MMSEQS_RELEASE_TAG" ]] || {
        echo "Downloaded MMseqs2 release tag must equal EXPECTED_MMSEQS_VERSION" >&2
        exit 1
    }
    if grep -qw avx2 /proc/cpuinfo; then
        MMSEQS_ARCHIVE_NAME="mmseqs-linux-avx2.tar.gz"
    else
        MMSEQS_ARCHIVE_NAME="mmseqs-linux-sse41.tar.gz"
    fi
    MMSEQS_URL="https://github.com/soedinglab/MMseqs2/releases/download/$MMSEQS_RELEASE_TAG/$MMSEQS_ARCHIVE_NAME"
    download_file "$MMSEQS_URL" "$TOOLS_ROOT/$MMSEQS_ARCHIVE_NAME"
    tar -xzf "$TOOLS_ROOT/$MMSEQS_ARCHIVE_NAME" -C "$TOOLS_ROOT"
    MMSEQS_BIN="$TOOLS_ROOT/mmseqs/bin/mmseqs"
    [[ -x "$MMSEQS_BIN" ]] || { echo "Downloaded MMseqs2 executable is missing" >&2; exit 1; }
fi
observed_mmseqs_version="$("$MMSEQS_BIN" version | tr -d '[:space:]')"
[[ "$observed_mmseqs_version" == "$EXPECTED_MMSEQS_BINARY_VERSION" ]] || {
    echo "MMseqs2 binary mismatch: expected_binary_identity=$EXPECTED_MMSEQS_BINARY_VERSION observed_binary_identity=$observed_mmseqs_version" >&2
    exit 1
}
{
    echo "release_tag=$MMSEQS_RELEASE_TAG"
    echo "expected_version=$EXPECTED_MMSEQS_VERSION"
    echo "expected_binary_version=$EXPECTED_MMSEQS_BINARY_VERSION"
    echo "observed_version=$observed_mmseqs_version"
    echo "executable=$MMSEQS_BIN"
    echo "executable_sha256=$(sha256sum "$MMSEQS_BIN" | awk '{print $1}')"
} > "$ARTIFACTS/logs/mmseqs_runtime.txt"
checkpoint_disk_usage mmseqs-installed

MANIFEST="$ARTIFACTS/contracts/frozen_input_manifest.json"
ATTRITION_POLICY="$ARTIFACTS/contracts/runtime_attrition_policy.json"
mkdir -p "$(dirname "$MANIFEST")"
contract_command=(
    "$PYTHON_BIN" -m homology_cluster_benchmark.runtime_contract prepare
    --manifest-out "$MANIFEST"
    --policy-out "$ATTRITION_POLICY"
    --source-scope "$UNIPROT_SOURCE_SCOPE"
    --framework-revision "$FRAMEWORK_REVISION"
    --uniref90-fasta "$INPUT_ROOT/uniref90.fasta.gz"
    --uniref90-fasta-url "$UNIREF90_URL"
    --uniref90-fasta-acquisition "$(awk -F '\t' '$1=="uniref90_fasta" {print $4}' "$ARTIFACTS/logs/runtime_input_staging.tsv")"
    --idmapping "$INPUT_ROOT/idmapping_selected.tab.gz"
    --idmapping-url "$IDMAPPING_URL"
    --idmapping-acquisition "$(awk -F '\t' '$1=="idmapping" {print $4}' "$ARTIFACTS/logs/runtime_input_staging.tsv")"
    --goa "$INPUT_ROOT/goa_uniprot_all.gaf.234.gz"
    --goa-url "$GOA_URL"
    --goa-acquisition "$(awk -F '\t' '$1=="goa" {print $4}' "$ARTIFACTS/logs/runtime_input_staging.tsv")"
    --go-obo "$INPUT_ROOT/go-basic.obo"
    --go-obo-url "$GO_OBO_URL"
    --go-obo-acquisition "$(awk -F '\t' '$1=="go_obo" {print $4}' "$ARTIFACTS/logs/runtime_input_staging.tsv")"
)
if [[ "$UNIPROT_SOURCE_SCOPE" != "trembl-only" ]]; then
    contract_command+=(
        --uniprot-sprot-sequences "$INPUT_ROOT/uniprot_sprot.dat.gz"
        --uniprot-sprot-sequences-url "$SPROT_URL"
        --uniprot-sprot-sequences-acquisition "$(awk -F '\t' '$1=="uniprot_sprot_sequences" {print $4}' "$ARTIFACTS/logs/runtime_input_staging.tsv")"
    )
fi
if [[ "$UNIPROT_SOURCE_SCOPE" != "sprot-only" ]]; then
    contract_command+=(
        --uniprot-trembl-sequences "$INPUT_ROOT/uniprot_trembl.dat.gz"
        --uniprot-trembl-sequences-url "$TREMBL_URL"
        --uniprot-trembl-sequences-acquisition "$(awk -F '\t' '$1=="uniprot_trembl_sequences" {print $4}' "$ARTIFACTS/logs/runtime_input_staging.tsv")"
    )
fi
export PYTHONPATH="$FRAMEWORK_DIR/benchmark_builders/homology_cluster/src${PYTHONPATH:+:$PYTHONPATH}"
"${contract_command[@]}" | tee "$ARTIFACTS/logs/runtime_contract.json"
checkpoint_disk_usage runtime-contract-created

input_sha() {
    awk -F '\t' -v role="$1" '$1==role {print $6}' "$ARTIFACTS/logs/runtime_input_staging.tsv"
}

echo "Running the existing homology benchmark builder"
SCRATCH_SAFETY_MULTIPLIER="${SCRATCH_SAFETY_MULTIPLIER:-1}"
MMSEQS_WORK_MULTIPLIER="${MMSEQS_WORK_MULTIPLIER:-1}"
PUBLICATION_SAFETY_MULTIPLIER="${PUBLICATION_SAFETY_MULTIPLIER:-1}"
builder_environment=(
    PYTHON_BIN="$PYTHON_BIN"
    IDENTITY="$IDENTITY"
    SPLIT_POLICY="$SPLIT_POLICY"
    TRAINING_POPULATION="$TRAINING_POPULATION"
    UNIPROT_SOURCE_SCOPE="$UNIPROT_SOURCE_SCOPE"
    OUTPUT_ROOT="$ARTIFACTS/benchmark"
    TEMP_DIR="$WORK/tmp"
    THREADS="${NSLOTS:-8}"
    SEED="$SEED"
    MIN_COUNT="$MIN_COUNT"
    MMSEQS_BIN="$MMSEQS_BIN"
    EXPECTED_MMSEQS_VERSION="$EXPECTED_MMSEQS_VERSION"
    FROZEN_INPUT_MANIFEST="$MANIFEST"
    ATTRITION_POLICY="$ATTRITION_POLICY"
    FRAMEWORK_REVISION="$FRAMEWORK_REVISION"
    DIAGNOSTIC_PILOT="$([[ "$RUNTIME_KIND" == "pilot" ]] && echo 1 || echo 0)"
    RUN_ID="$RUN_ID"
    REQUESTED_SLOTS="${NSLOTS:-8}"
    ALLOCATED_SLOTS="${NSLOTS:-8}"
    UNIPROT_RELEASE=2026_02
    GOA_RELEASE=234
    ONTOLOGY_RELEASE=releases/2026-06-15
    SCRATCH_SAFETY_MULTIPLIER="$SCRATCH_SAFETY_MULTIPLIER"
    MMSEQS_WORK_MULTIPLIER="$MMSEQS_WORK_MULTIPLIER"
    PUBLICATION_SAFETY_MULTIPLIER="$PUBLICATION_SAFETY_MULTIPLIER"
    NO_DOWNLOADS=1
    UNIREF90_FASTA="$INPUT_ROOT/uniref90.fasta.gz"
    UNIREF90_FASTA_URL="$UNIREF90_URL"
    UNIREF90_FASTA_SHA256="$(input_sha uniref90_fasta)"
    IDMAPPING="$INPUT_ROOT/idmapping_selected.tab.gz"
    IDMAPPING_URL="$IDMAPPING_URL"
    IDMAPPING_SHA256="$(input_sha idmapping)"
    GOA="$INPUT_ROOT/goa_uniprot_all.gaf.234.gz"
    GOA_URL="$GOA_URL"
    GOA_SHA256="$(input_sha goa)"
    GO_OBO="$INPUT_ROOT/go-basic.obo"
    GO_OBO_URL="$GO_OBO_URL"
    GO_OBO_SHA256="$(input_sha go_obo)"
    LOG_FILE="$ARTIFACTS/logs/builder.log"
)
if [[ "$UNIPROT_SOURCE_SCOPE" != "trembl-only" ]]; then
    builder_environment+=(
        UNIPROT_SPROT_SEQUENCES="$INPUT_ROOT/uniprot_sprot.dat.gz"
        UNIPROT_SPROT_SEQUENCES_URL="$SPROT_URL"
        UNIPROT_SPROT_SEQUENCES_SHA256="$(input_sha uniprot_sprot_sequences)"
    )
fi
if [[ "$UNIPROT_SOURCE_SCOPE" != "sprot-only" ]]; then
    builder_environment+=(
        UNIPROT_TREMBL_SEQUENCES="$INPUT_ROOT/uniprot_trembl.dat.gz"
        UNIPROT_TREMBL_SEQUENCES_URL="$TREMBL_URL"
        UNIPROT_TREMBL_SEQUENCES_SHA256="$(input_sha uniprot_trembl_sequences)"
    )
fi
checkpoint_disk_usage builder-starting
env -u PERSISTENT_RESULTS_ROOT "${builder_environment[@]}" \
    bash "$FRAMEWORK_DIR/scripts/benchmark_generation/run_homology_cluster_benchmark.sh"
checkpoint_disk_usage builder-complete

COMPLETION_LIST="$ARTIFACTS/logs/completion_markers.txt"
find "$ARTIFACTS/benchmark" -type f -name RUN_COMPLETE.json -print > "$COMPLETION_LIST"
completion_count="$(awk 'END {print NR}' "$COMPLETION_LIST")"
[[ "$completion_count" == "1" ]] || {
    echo "Expected exactly one completed publication; found $completion_count" >&2
    exit 1
}
RUN_DIR="$(dirname "$(cat "$COMPLETION_LIST")")"
"$PYTHON_BIN" -m homology_cluster_benchmark validate --run-dir "$RUN_DIR"
"$PYTHON_BIN" -m homology_cluster_benchmark.runtime_contract review \
    --run-dir "$RUN_DIR" \
    --output-dir "$ARTIFACTS/review" \
    --run-kind "$RUNTIME_KIND"
checkpoint_disk_usage validation-review-complete

cat > "$ARTIFACTS/TASK_SUMMARY.txt" <<EOF
run_kind=$RUNTIME_KIND
job_id=$JOB_KEY
task_id=$TASK_ID
identity_percent=$IDENTITY
uniprot_source_scope=$UNIPROT_SOURCE_SCOPE
framework_revision=$FRAMEWORK_REVISION
pilot_required_for_array=false
publication_directory=$RUN_DIR
EOF

echo "Homology runtime task completed and passed automatic review"
