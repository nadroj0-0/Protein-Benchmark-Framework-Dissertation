#!/usr/bin/env bash
# UCL Grid Engine wrapper. Resource directives are provisional until a measured smoke run.

#$ -l tmem=64G
#$ -l tscratch=200G
#$ -l scratch0free=200G
#$ -l h_rt=72:0:0
#$ -pe smp 8
#$ -j y
#$ -N homology_cluster

set -euo pipefail

portable_sha256() {
    if command -v sha256sum >/dev/null 2>&1; then
        sha256sum "$1" | awk '{print $1}'
    elif command -v shasum >/dev/null 2>&1; then
        shasum -a 256 "$1" | awk '{print $1}'
    else
        echo "Neither sha256sum nor shasum is available for queued-evidence verification" >&2
        return 1
    fi
}

verify_queued_hash() {
    local path_variable="$1"
    local hash_variable="$2"
    local path_value="${!path_variable:-}"
    local expected_value="${!hash_variable:-}"
    [[ -f "$path_value" && "$expected_value" =~ ^[0-9a-f]{64}$ ]] || {
        echo "Queued evidence binding is missing for $path_variable/$hash_variable" >&2
        return 1
    }
    [[ "$(portable_sha256 "$path_value")" == "$expected_value" ]] || {
        echo "Queued evidence changed after launcher review: $path_variable" >&2
        return 1
    }
}

TEST_MODE="${HOMOLOGY_WRAPPER_TEST_MODE:-0}"
TASK_ID="${SGE_TASK_ID:-}"
case "$TASK_ID" in
    1) MAPPED_IDENTITY=30 ;;
    2) MAPPED_IDENTITY=25 ;;
    3) MAPPED_IDENTITY=20 ;;
    4) MAPPED_IDENTITY=15 ;;
    5) MAPPED_IDENTITY=10 ;;
    6) MAPPED_IDENTITY=5 ;;
    *) echo "SGE_TASK_ID must be one integer from 1 through 6" >&2; exit 2 ;;
esac
if [[ -n "${IDENTITY:-}" && "$IDENTITY" != "$MAPPED_IDENTITY" ]]; then
    echo "IDENTITY=$IDENTITY conflicts with locked SGE_TASK_ID=$TASK_ID mapping to $MAPPED_IDENTITY" >&2
    exit 2
fi
IDENTITY="$MAPPED_IDENTITY"
TRAINING_POPULATION="${TRAINING_POPULATION:-annotated-only}"
if [[ "$TRAINING_POPULATION" != "annotated-only" ]]; then
    echo "TRAINING_POPULATION=$TRAINING_POPULATION is unsupported; only annotated-only is implemented." >&2
    echo "No zero-negative, homology-transfer, or representative-label policy is authorized." >&2
    exit 2
fi
FIXTURE_MODE_VALUE="${FIXTURE_MODE:-0}"
DIAGNOSTIC_PILOT_VALUE="${DIAGNOSTIC_PILOT:-0}"
case "$FIXTURE_MODE_VALUE" in 0|1) ;; *) echo "FIXTURE_MODE must be 0 or 1" >&2; exit 2 ;; esac
case "$DIAGNOSTIC_PILOT_VALUE" in 0|1) ;; *) echo "DIAGNOSTIC_PILOT must be 0 or 1" >&2; exit 2 ;; esac
UNIPROT_SOURCE_SCOPE="${UNIPROT_SOURCE_SCOPE:-}"
if [[ -z "$UNIPROT_SOURCE_SCOPE" && "$FIXTURE_MODE_VALUE" == "1" ]]; then
    UNIPROT_SOURCE_SCOPE="sprot-only"
fi
case "$UNIPROT_SOURCE_SCOPE" in
    sprot-only|trembl-only|sprot-and-trembl) ;;
    *) echo "Production array tasks require explicit UNIPROT_SOURCE_SCOPE" >&2; exit 2 ;;
esac
if [[ "$DIAGNOSTIC_PILOT_VALUE" == "1" && "$TASK_ID" != "1" ]]; then
    echo "DIAGNOSTIC_PILOT is locked to SGE_TASK_ID=1 / identity 30%" >&2
    exit 2
fi
if [[ "$TEST_MODE" == "1" && "$FIXTURE_MODE_VALUE" != "1" ]]; then
    echo "HOMOLOGY_WRAPPER_TEST_MODE is permitted only with FIXTURE_MODE=1" >&2
    exit 2
fi
if [[ "$FIXTURE_MODE_VALUE" != "1" ]]; then
    [[ -n "${FROZEN_INPUT_MANIFEST:-}" && -f "$FROZEN_INPUT_MANIFEST" ]] || {
        echo "Production HPC run requires FROZEN_INPUT_MANIFEST before staging or checkout" >&2
        exit 2
    }
    [[ -n "${EXPECTED_MMSEQS_VERSION:-}" ]] || {
        echo "Production HPC run requires exact EXPECTED_MMSEQS_VERSION" >&2
        exit 2
    }
    [[ "${FRAMEWORK_REVISION:-}" =~ ^[0-9a-f]{40}$ ]] || {
        echo "Production HPC run requires FRAMEWORK_REVISION as exactly 40 lowercase hex characters" >&2
        exit 2
    }
    [[ "${NO_DOWNLOADS:-}" == "1" ]] || {
        echo "Production array tasks require NO_DOWNLOADS=1 and shared local frozen inputs" >&2
        exit 2
    }
    verify_queued_hash FROZEN_INPUT_MANIFEST EXPECTED_FROZEN_INPUT_MANIFEST_SHA256
    if [[ "$DIAGNOSTIC_PILOT_VALUE" != "1" ]]; then
        [[ -n "${ATTRITION_POLICY:-}" && -f "$ATTRITION_POLICY" ]] || {
            echo "Production array run requires reviewed ATTRITION_POLICY" >&2
            exit 2
        }
        for binding in \
            "ATTRITION_POLICY EXPECTED_ATTRITION_POLICY_SHA256" \
            "PILOT_APPROVAL EXPECTED_PILOT_APPROVAL_SHA256" \
            "PILOT_COMPLETION_MARKER EXPECTED_PILOT_COMPLETION_MARKER_SHA256" \
            "PILOT_ATTRITION_REPORT EXPECTED_PILOT_ATTRITION_REPORT_SHA256" \
            "PILOT_TASK_CONTEXT EXPECTED_PILOT_TASK_CONTEXT_SHA256" \
            "PILOT_MEASUREMENT_EVIDENCE EXPECTED_PILOT_MEASUREMENT_EVIDENCE_SHA256"
        do
            read -r path_variable hash_variable <<< "$binding"
            verify_queued_hash "$path_variable" "$hash_variable"
        done
        [[ -n "${PILOT_RUN_DIR:-}" && -d "$PILOT_RUN_DIR" ]] || {
            echo "Production array run requires accessible PILOT_RUN_DIR" >&2
            exit 2
        }
    fi
    for required_path_variable in UNIREF90_FASTA IDMAPPING GOA GO_OBO; do
        [[ -n "${!required_path_variable:-}" && -f "${!required_path_variable}" ]] || {
            echo "Production array task requires shared local $required_path_variable" >&2
            exit 2
        }
    done
    if [[ "$UNIPROT_SOURCE_SCOPE" != "trembl-only" ]]; then
        [[ -n "${UNIPROT_SPROT_SEQUENCES:-}" && -f "$UNIPROT_SPROT_SEQUENCES" ]] || {
            echo "Selected scope requires shared local UNIPROT_SPROT_SEQUENCES" >&2
            exit 2
        }
    elif [[ -n "${UNIPROT_SPROT_SEQUENCES:-}" ]]; then
        echo "trembl-only forbids UNIPROT_SPROT_SEQUENCES" >&2
        exit 2
    fi
    if [[ "$UNIPROT_SOURCE_SCOPE" != "sprot-only" ]]; then
        [[ -n "${UNIPROT_TREMBL_SEQUENCES:-}" && -f "$UNIPROT_TREMBL_SEQUENCES" ]] || {
            echo "Selected scope requires shared local UNIPROT_TREMBL_SEQUENCES" >&2
            exit 2
        }
    elif [[ -n "${UNIPROT_TREMBL_SEQUENCES:-}" ]]; then
        echo "sprot-only forbids UNIPROT_TREMBL_SEQUENCES" >&2
        exit 2
    fi
fi

if [[ "$TEST_MODE" != "1" ]]; then
    for test_variable in \
        HOMOLOGY_BUILDER_COMMAND HOMOLOGY_VALIDATE_COMMAND HOMOLOGY_COPY_COMMAND \
        HOMOLOGY_FRAMEWORK_DIR HOMOLOGY_SKIP_CONDA HOMOLOGY_SCRATCH_FREE_BYTES_OVERRIDE \
        HOMOLOGY_PERSISTENT_FREE_BYTES_OVERRIDE HOMOLOGY_PRECOPY_FREE_BYTES_OVERRIDE \
        HOMOLOGY_SIGNAL_GRACE_SECONDS CAPACITY_PYTHON
    do
        if [[ -n "${!test_variable:-}" ]]; then
            echo "$test_variable is a test-only override; set no production override" >&2
            exit 2
        fi
    done
fi

SPLIT_POLICY="${SPLIT_POLICY:-sequence-balanced}"
SEED="${SEED:-0}"
MIN_COUNT="${MIN_COUNT:-50}"
case "$SPLIT_POLICY" in
    cluster-count-random|sequence-balanced) ;;
    *) echo "SPLIT_POLICY has an unsupported value: $SPLIT_POLICY" >&2; exit 2 ;;
esac
[[ "$SEED" =~ ^[0-9]+$ ]] || { echo "SEED must be a non-negative integer" >&2; exit 2; }
[[ "$MIN_COUNT" =~ ^[1-9][0-9]*$ ]] || { echo "MIN_COUNT must be a positive integer" >&2; exit 2; }
if [[ "$FIXTURE_MODE_VALUE" == "1" ]]; then
    JOB_KEY="${JOB_ID:-fixture-job}"
else
    [[ -n "${JOB_ID:-}" ]] || { echo "Production array task requires JOB_ID" >&2; exit 2; }
    JOB_KEY="$JOB_ID"
fi
WORK_BASE="${WORK_BASE:-/scratch0}"
[[ "$WORK_BASE" == /* && "$WORK_BASE" != "/" && -d "$WORK_BASE" && ! -L "$WORK_BASE" ]] || {
    echo "WORK_BASE must be an existing absolute non-root directory, not a symlink" >&2
    exit 2
}
WORK_BASE="$(cd "$WORK_BASE" && pwd -P)"
RUN_ID="${RUN_ID:-$(date -u +%Y%m%dT%H%M%SZ)-$$-${RANDOM}}"
for component_name in JOB_KEY TASK_ID IDENTITY UNIPROT_SOURCE_SCOPE RUN_ID; do
    component_value="${!component_name}"
    [[ "$component_value" =~ ^[A-Za-z0-9._-]+$ && "$component_value" =~ [A-Za-z0-9] ]] || {
        echo "$component_name contains unsafe path characters" >&2
        exit 2
    }
done
WORK="$WORK_BASE/homology_cluster_${JOB_KEY}_${TASK_ID}_${IDENTITY}_${UNIPROT_SOURCE_SCOPE}_${RUN_ID}"
FRAMEWORK_REPO_URL="${FRAMEWORK_REPO_URL:-https://github.com/nadroj0-0/Protein-Benchmark-Framework-Dissertation.git}"
FRAMEWORK_REVISION="${FRAMEWORK_REVISION:-}"
FRAMEWORK_DIR="${HOMOLOGY_FRAMEWORK_DIR:-$WORK/Protein-Benchmark-Framework-Dissertation}"
SCRATCH_INPUTS="$WORK/inputs"
SCRATCH_OUTPUTS="$WORK/run"
SCRATCH_TEMP="$WORK/tmp"
RESULTS_ROOT="${RESULTS_ROOT:-$HOME/homology_cluster_benchmark_results}"
PERSISTENT_RESULTS_ROOT="${PERSISTENT_RESULTS_ROOT:-$RESULTS_ROOT}"
if [[ "$FIXTURE_MODE_VALUE" == "1" ]]; then
    REVISION_TAG="fixture"
else
    REVISION_TAG="${FRAMEWORK_REVISION:0:12}"
fi
RUN_TAG="source_${UNIPROT_SOURCE_SCOPE}/framework_${REVISION_TAG}/${RUN_ID}/job_${JOB_KEY}/task_${TASK_ID}/identity_${IDENTITY}_${SPLIT_POLICY}_${TRAINING_POPULATION}_seed_${SEED}_min_count_${MIN_COUNT}"
FINAL_RUN_ROOT="$PERSISTENT_RESULTS_ROOT/$RUN_TAG"
PARTIAL_RUN_ROOT="${FINAL_RUN_ROOT}.partial-${JOB_KEY}"
FAILED_RUN_ROOT="${FINAL_RUN_ROOT}.failed"
FAILED_PARTIAL_ROOT="${FAILED_RUN_ROOT}.partial-${JOB_KEY}"
CLAIM_ROOT="${FINAL_RUN_ROOT}.claim"
CONDA_EXE="${CONDA_EXE:-/share/apps/miniforge3_mamba/bin/conda}"
MMFP_ENV_DIR="${MMFP_ENV_DIR:-$HOME/.conda/envs/mmfp}"
SCRATCH_RESERVE_BYTES="${SCRATCH_RESERVE_BYTES:-0}"
PERSISTENT_RESERVE_BYTES="${PERSISTENT_RESERVE_BYTES:-0}"
MMSEQS_WORK_MULTIPLIER="${MMSEQS_WORK_MULTIPLIER:-8}"
PUBLICATION_SAFETY_MULTIPLIER="${PUBLICATION_SAFETY_MULTIPLIER:-2}"
COPY_SAFETY_MULTIPLIER="${COPY_SAFETY_MULTIPLIER:-1.10}"
ACTIVE_CHILD_PID=""
SIGNAL_NAME=""
SIGNAL_STATUS=0
JOB_SUCCEEDED=0
COPY_FAILED=0
FAILURE_STAGE="initialization"
SIGNAL_WATCHDOG_PID=""
SIGNAL_GRACE_SECONDS="${HOMOLOGY_SIGNAL_GRACE_SECONDS:-10}"
WORK_OWNED=0
PERSISTENT_CLAIM_OWNED=0
PARTIAL_OWNED=0
FINAL_OWNED=0

ALLOCATED_SLOTS="${NSLOTS:-}"
if [[ "$FIXTURE_MODE_VALUE" == "1" && -z "$ALLOCATED_SLOTS" ]]; then
    ALLOCATED_SLOTS=1
fi
REQUESTED_THREADS="${THREADS:-$ALLOCATED_SLOTS}"
[[ "$ALLOCATED_SLOTS" =~ ^[1-9][0-9]*$ && "$REQUESTED_THREADS" =~ ^[1-9][0-9]*$ ]] || {
    echo "NSLOTS and THREADS must be positive integers" >&2
    exit 1
}
if [[ "$REQUESTED_THREADS" != "$ALLOCATED_SLOTS" ]]; then
    echo "THREADS=$REQUESTED_THREADS must equal scheduler-provided NSLOTS=$ALLOCATED_SLOTS" >&2
    exit 1
fi
if [[ -n "${REQUESTED_SLOTS:-}" && "$REQUESTED_SLOTS" != "8" ]]; then
    echo "REQUESTED_SLOTS must equal the locked smp request of 8" >&2
    exit 1
fi
if [[ "$FIXTURE_MODE_VALUE" != "1" && "$ALLOCATED_SLOTS" != "8" ]]; then
    echo "Production array task requires the requested smp allocation of NSLOTS=8" >&2
    exit 1
fi

identity_directory() {
    if [[ "$IDENTITY" == "5" ]]; then
        printf 'identity_05'
    else
        printf 'identity_%s' "$IDENTITY"
    fi
}

RUN_RELATIVE_PATH="source_${UNIPROT_SOURCE_SCOPE}/framework_${REVISION_TAG}/$(identity_directory)/$SPLIT_POLICY/$TRAINING_POPULATION/seed_$SEED/min_count_$MIN_COUNT"

file_size() {
    LC_ALL=C wc -c < "$1" | awk '{print $1}'
}

tree_size() {
    local kilobytes
    kilobytes="$(du -sk "$1" | awk '{print $1}')"
    printf '%s' "$((kilobytes * 1024))"
}

free_bytes() {
    local path="$1"
    local override="$2"
    if [[ -n "$override" ]]; then
        printf '%s' "$override"
        return
    fi
    df -Pk "$path" | awk 'NR==2 {printf "%.0f", $4 * 1024}'
}

scaled_bytes() {
    awk -v bytes="$1" -v multiplier="$2" 'BEGIN {printf "%.0f", bytes * multiplier}'
}

require_capacity() {
    local label="$1"
    local available="$2"
    local required="$3"
    if (( available < required )); then
        echo "$label capacity preflight failed: free=$available required=$required bytes" >&2
        return 1
    fi
}

require_nonnegative_integer() {
    local label="$1"
    local value="$2"
    [[ "$value" =~ ^[0-9]+$ ]] || {
        echo "$label must be a non-negative integer; observed=$value" >&2
        exit 2
    }
}

require_safety_multiplier() {
    local label="$1"
    local value="$2"
    [[ "$value" =~ ^[0-9]+([.][0-9]+)?$ ]] || {
        echo "$label must be a finite numeric value of at least 1; observed=$value" >&2
        exit 2
    }
    awk -v value="$value" 'BEGIN {exit !(value >= 1)}' || {
        echo "$label must be at least 1; observed=$value" >&2
        exit 2
    }
}

check_signal() {
    if [[ "$SIGNAL_STATUS" != "0" ]]; then
        exit "$SIGNAL_STATUS"
    fi
}

wait_for_active_child() {
    local child_status=0
    while true; do
        wait "$ACTIVE_CHILD_PID"
        child_status=$?
        if ! kill -0 "$ACTIVE_CHILD_PID" 2>/dev/null; then
            break
        fi
    done
    if [[ -n "$SIGNAL_WATCHDOG_PID" ]]; then
        kill "$SIGNAL_WATCHDOG_PID" 2>/dev/null || true
        wait "$SIGNAL_WATCHDOG_PID" 2>/dev/null || true
        SIGNAL_WATCHDOG_PID=""
    fi
    ACTIVE_CHILD_PID=""
    if [[ "$SIGNAL_STATUS" != "0" ]]; then
        return "$SIGNAL_STATUS"
    fi
    return "$child_status"
}

local_input_bytes() {
    local total=0
    local variable value
    for variable in \
        UNIREF90_FASTA IDMAPPING UNIPROT_SPROT_SEQUENCES UNIPROT_TREMBL_SEQUENCES \
        GOA GO_OBO FROZEN_INPUT_MANIFEST ATTRITION_POLICY ATTRITION_OVERRIDE \
        CLUSTER_ASSIGNMENTS
    do
        value="${!variable:-}"
        if [[ -n "$value" && -f "$value" ]]; then
            total=$((total + $(file_size "$value")))
        fi
    done
    printf '%s' "$total"
}

manifest_input_bytes() {
    if [[ -z "${FROZEN_INPUT_MANIFEST:-}" || ! -f "$FROZEN_INPUT_MANIFEST" ]]; then
        printf '0'
        return
    fi
    "$PYTHON_BIN" -c \
        'import json,sys; print(sum(int(x["size_bytes"]) for x in json.load(open(sys.argv[1]))["inputs"]))' \
        "$FROZEN_INPUT_MANIFEST"
}

manifest_uniref_bytes() {
    if [[ -n "${UNIREF90_FASTA:-}" && -f "$UNIREF90_FASTA" ]]; then
        file_size "$UNIREF90_FASTA"
        return
    fi
    "$PYTHON_BIN" -c \
        'import json,sys; print(next(int(x["size_bytes"]) for x in json.load(open(sys.argv[1]))["inputs"] if x["name"]=="uniref90_fasta"))' \
        "$FROZEN_INPUT_MANIFEST"
}

copy_tree() {
    local source="$1"
    local destination="$2"
    mkdir -p "$destination"
    if [[ -n "${HOMOLOGY_COPY_COMMAND:-}" ]]; then
        "$HOMOLOGY_COPY_COMMAND" "$source" "$destination"
    else
        cp -a "$source/." "$destination/"
    fi
}

validate_run() {
    local run_dir="$1"
    if [[ -n "${HOMOLOGY_VALIDATE_COMMAND:-}" ]]; then
        "$HOMOLOGY_VALIDATE_COMMAND" "$run_dir"
    else
        PYTHONPATH="$FRAMEWORK_DIR/benchmark_builders/homology_cluster/src${PYTHONPATH:+:$PYTHONPATH}" \
            "$PYTHON_BIN" -m homology_cluster_benchmark validate --run-dir "$run_dir"
    fi
}

forward_signal() {
    SIGNAL_NAME="$1"
    SIGNAL_STATUS="$2"
    FAILURE_STAGE="signal-$1"
    echo "Received $1; forwarding it to active child" >&2
    if [[ -n "$ACTIVE_CHILD_PID" ]]; then
        kill -s "$1" "$ACTIVE_CHILD_PID" 2>/dev/null || true
        # Non-interactive Bash starts asynchronous children with SIGINT ignored. Deliver TERM as
        # the portable termination mechanism while retaining the wrapper's conventional 130 exit.
        if [[ "$1" == "INT" ]]; then
            kill -s TERM "$ACTIVE_CHILD_PID" 2>/dev/null || true
        fi
        local child_pid="$ACTIVE_CHILD_PID"
        (
            sleep "$SIGNAL_GRACE_SECONDS"
            kill -s KILL "$child_pid" 2>/dev/null || true
        ) &
        SIGNAL_WATCHDOG_PID=$!
    fi
}
trap 'forward_signal INT 130' INT
trap 'forward_signal TERM 143' TERM

json_escape() {
    local value="$1"
    value="${value//\\/\\\\}"
    value="${value//\"/\\\"}"
    value="${value//$'\n'/\\n}"
    printf '%s' "$value"
}

write_failure_metadata() {
    local destination="$1"
    local exit_status="$2"
    mkdir -p "$destination"
    printf '{\n  "complete": false,\n  "exit_code": %s,\n  "signal": "%s",\n  "failure_stage": "%s",\n  "scratch_path": "%s",\n  "final_published": false\n}\n' \
        "$exit_status" "$(json_escape "$SIGNAL_NAME")" \
        "$(json_escape "$FAILURE_STAGE")" "$(json_escape "$WORK")" \
        > "$destination/FAILURE.json"
}

publish_failure() {
    local exit_status="$1"
    [[ "$PERSISTENT_CLAIM_OWNED" == "1" ]] || {
        echo "No task-owned persistent claim; refusing to publish failure diagnostics" >&2
        return 1
    }
    [[ ! -e "$FAILED_RUN_ROOT" && ! -e "$FAILED_PARTIAL_ROOT" ]] || {
        echo "Failure destination already exists; diagnostics could not be published" >&2
        return 1
    }
    write_failure_metadata "$SCRATCH_OUTPUTS/logs" "$exit_status" || true
    if [[ "$FINAL_OWNED" == "1" && -d "$FINAL_RUN_ROOT" ]]; then
        find "$FINAL_RUN_ROOT" -name RUN_COMPLETE.json -type f -exec rm -f {} \;
        mv "$FINAL_RUN_ROOT" "$FAILED_PARTIAL_ROOT" || return 1
        FINAL_OWNED=0
    elif [[ "$PARTIAL_OWNED" == "1" && -d "$PARTIAL_RUN_ROOT" ]]; then
        find "$PARTIAL_RUN_ROOT" -name RUN_COMPLETE.json -type f -exec rm -f {} \;
        mv "$PARTIAL_RUN_ROOT" "$FAILED_PARTIAL_ROOT" || return 1
        PARTIAL_OWNED=0
    else
        mkdir -p "$(dirname "$FAILED_PARTIAL_ROOT")"
        mkdir "$FAILED_PARTIAL_ROOT" || return 1
        if [[ -d "$SCRATCH_OUTPUTS" ]]; then
            copy_tree "$SCRATCH_OUTPUTS" "$FAILED_PARTIAL_ROOT" || return 1
        fi
    fi
    find "$FAILED_PARTIAL_ROOT" -name RUN_COMPLETE.json -type f -exec rm -f {} \;
    write_failure_metadata "$FAILED_PARTIAL_ROOT" "$exit_status" || return 1
    mv "$FAILED_PARTIAL_ROOT" "$FAILED_RUN_ROOT" || return 1
    echo "Failure diagnostics published atomically: $FAILED_RUN_ROOT" >&2
    return 0
}

safe_cleanup_work() {
    [[ "$WORK_OWNED" == "1" ]] || return 0
    [[ -e "$WORK" ]] || return 0
    [[ -n "$WORK" && "$WORK" != "/" && "$WORK" != "$WORK_BASE" ]] || {
        echo "Refusing unsafe scratch cleanup path: $WORK" >&2
        return 1
    }
    local resolved_base resolved_work expected_owner observed_owner
    resolved_base="$(cd "$WORK_BASE" 2>/dev/null && pwd -P)" || return 1
    resolved_work="$(cd "$WORK" 2>/dev/null && pwd -P)" || return 1
    case "$resolved_work" in
        "$resolved_base"/homology_cluster_*) ;;
        *) echo "Refusing scratch cleanup outside owned work base: $resolved_work" >&2; return 1 ;;
    esac
    expected_owner="$JOB_KEY:$TASK_ID:$IDENTITY:$UNIPROT_SOURCE_SCOPE:$RUN_ID"
    [[ -f "$resolved_work/.homology-task-owner" ]] || {
        echo "Refusing scratch cleanup without ownership marker: $resolved_work" >&2
        return 1
    }
    observed_owner="$(<"$resolved_work/.homology-task-owner")"
    [[ "$observed_owner" == "$expected_owner" ]] || {
        echo "Refusing scratch cleanup with mismatched ownership marker" >&2
        return 1
    }
    cd "$HOME" 2>/dev/null || cd "$resolved_base"
    rm -rf -- "$resolved_work" || return 1
    [[ ! -e "$resolved_work" && ! -L "$resolved_work" ]] || return 1
    WORK_OWNED=0
}

safe_cleanup_claim() {
    [[ "$PERSISTENT_CLAIM_OWNED" == "1" ]] || return 0
    [[ -d "$CLAIM_ROOT" && ! -L "$CLAIM_ROOT" ]] || {
        echo "Refusing to remove a missing or unsafe persistent claim" >&2
        return 1
    }
    local expected_owner observed_owner
    expected_owner="$JOB_KEY:$TASK_ID:$IDENTITY:$UNIPROT_SOURCE_SCOPE:$RUN_ID"
    [[ -f "$CLAIM_ROOT/.homology-task-owner" ]] || {
        echo "Refusing to remove persistent claim without ownership marker" >&2
        return 1
    }
    observed_owner="$(<"$CLAIM_ROOT/.homology-task-owner")"
    [[ "$observed_owner" == "$expected_owner" ]] || {
        echo "Refusing to remove persistent claim with mismatched ownership marker" >&2
        return 1
    }
    rm "$CLAIM_ROOT/.homology-task-owner" || return 1
    rmdir "$CLAIM_ROOT" || return 1
    [[ ! -e "$CLAIM_ROOT" && ! -L "$CLAIM_ROOT" ]] || return 1
    PERSISTENT_CLAIM_OWNED=0
}

cleanup() {
    local status=$?
    trap - EXIT
    set +e
    if [[ "$SIGNAL_STATUS" != "0" ]]; then
        JOB_SUCCEEDED=0
    fi
    local failure_published=0
    if [[ "$JOB_SUCCEEDED" != "1" && "$PERSISTENT_CLAIM_OWNED" == "1" ]]; then
        publish_failure "$status"
        failure_published=$?
        if [[ "$failure_published" -ne 0 && "$status" -eq 0 ]]; then
            status=1
        fi
    fi
    safe_cleanup_work
    local cleanup_status=$?
    if [[ "$cleanup_status" -ne 0 && "$status" -eq 0 ]]; then
        status=1
    fi
    safe_cleanup_claim
    local claim_cleanup_status=$?
    if [[ "$claim_cleanup_status" -ne 0 && "$status" -eq 0 ]]; then
        status=1
    fi
    if [[ "$COPY_FAILED" == "1" ]]; then
        echo "Persistent copy failed; job-owned scratch cleanup was still attempted" >&2
    fi
    if [[ "$SIGNAL_STATUS" != "0" ]]; then
        status="$SIGNAL_STATUS"
    fi
    exit "$status"
}
trap cleanup EXIT

stage_if_local() {
    local variable_name="$1"
    local value="${!variable_name:-}"
    [[ -n "$value" ]] || return 0
    [[ -f "$value" ]] || { echo "$variable_name does not exist: $value" >&2; exit 1; }
    if [[ "$FIXTURE_MODE_VALUE" != "1" ]]; then
        # Production arrays consume the reviewed shared frozen collection read-only. Each task's
        # Python builder rechecks hashes; copying multi-gigabyte sources six times is prohibited.
        export "$variable_name"
        return 0
    fi
    local destination_dir="$SCRATCH_INPUTS/$variable_name"
    mkdir -p "$destination_dir"
    local destination="$destination_dir/$(basename "$value")"
    cp -Lp "$value" "$destination"
    printf -v "$variable_name" '%s' "$destination"
    export "$variable_name"
}

echo "Host        : $(hostname)"
echo "Job ID      : $JOB_KEY"
echo "Array task  : $TASK_ID"
echo "Identity    : $IDENTITY%"
echo "Source scope: $UNIPROT_SOURCE_SCOPE"
echo "Run ID      : $RUN_ID"
echo "Split policy: $SPLIT_POLICY"
echo "Scratch     : $WORK"
echo "Final output: $FINAL_RUN_ROOT"
echo "Resources   : provisional 64G tmem / 200G scratch / 72h"

[[ ! -e "$WORK" && ! -L "$WORK" ]] || {
    echo "Refusing to reuse pre-existing task scratch path: $WORK" >&2
    exit 2
}
mkdir "$WORK"
WORK_OWNED=1
printf '%s\n' "$JOB_KEY:$TASK_ID:$IDENTITY:$UNIPROT_SOURCE_SCOPE:$RUN_ID" > "$WORK/.homology-task-owner"
mkdir -p "$SCRATCH_INPUTS" "$SCRATCH_OUTPUTS/logs" "$SCRATCH_TEMP" "$RESULTS_ROOT" "$PERSISTENT_RESULTS_ROOT"

RESOLVED_RESULTS_ROOT="$(cd "$RESULTS_ROOT" && pwd -P)"
RESOLVED_PERSISTENT_RESULTS_ROOT="$(cd "$PERSISTENT_RESULTS_ROOT" && pwd -P)"
[[ "$RESOLVED_RESULTS_ROOT" == "$RESOLVED_PERSISTENT_RESULTS_ROOT" ]] || {
    echo "RESULTS_ROOT and PERSISTENT_RESULTS_ROOT must resolve to the same destination" >&2
    exit 2
}
RESULTS_ROOT="$RESOLVED_RESULTS_ROOT"
PERSISTENT_RESULTS_ROOT="$RESOLVED_PERSISTENT_RESULTS_ROOT"
FINAL_RUN_ROOT="$PERSISTENT_RESULTS_ROOT/$RUN_TAG"
PARTIAL_RUN_ROOT="${FINAL_RUN_ROOT}.partial-${JOB_KEY}"
FAILED_RUN_ROOT="${FINAL_RUN_ROOT}.failed"
FAILED_PARTIAL_ROOT="${FAILED_RUN_ROOT}.partial-${JOB_KEY}"
CLAIM_ROOT="${FINAL_RUN_ROOT}.claim"
mkdir -p "$(dirname "$FINAL_RUN_ROOT")"
for reserved_path in "$FINAL_RUN_ROOT" "$PARTIAL_RUN_ROOT" "$FAILED_RUN_ROOT" "$FAILED_PARTIAL_ROOT"; do
    [[ ! -e "$reserved_path" && ! -L "$reserved_path" ]] || {
        echo "Refusing to overwrite pre-existing persistent task path: $reserved_path" >&2
        exit 1
    }
done
[[ ! -e "$CLAIM_ROOT" && ! -L "$CLAIM_ROOT" ]] || {
    echo "Persistent task claim already exists: $CLAIM_ROOT" >&2
    exit 1
}
mkdir "$CLAIM_ROOT" || {
    echo "Could not atomically reserve persistent task path: $CLAIM_ROOT" >&2
    exit 1
}
printf '%s\n' "$JOB_KEY:$TASK_ID:$IDENTITY:$UNIPROT_SOURCE_SCOPE:$RUN_ID" \
    > "$CLAIM_ROOT/.homology-task-owner"
PERSISTENT_CLAIM_OWNED=1

require_nonnegative_integer "SEED" "$SEED"
require_nonnegative_integer "MIN_COUNT" "$MIN_COUNT"
(( MIN_COUNT >= 1 )) || { echo "MIN_COUNT must be at least 1" >&2; exit 2; }
require_nonnegative_integer "SCRATCH_RESERVE_BYTES" "$SCRATCH_RESERVE_BYTES"
require_nonnegative_integer "PERSISTENT_RESERVE_BYTES" "$PERSISTENT_RESERVE_BYTES"
require_nonnegative_integer "SIGNAL_GRACE_SECONDS" "$SIGNAL_GRACE_SECONDS"
(( SIGNAL_GRACE_SECONDS >= 1 )) || { echo "SIGNAL_GRACE_SECONDS must be at least 1" >&2; exit 2; }
require_safety_multiplier "MMSEQS_WORK_MULTIPLIER" "$MMSEQS_WORK_MULTIPLIER"
require_safety_multiplier "PUBLICATION_SAFETY_MULTIPLIER" "$PUBLICATION_SAFETY_MULTIPLIER"
require_safety_multiplier "COPY_SAFETY_MULTIPLIER" "$COPY_SAFETY_MULTIPLIER"
for override_name in \
    HOMOLOGY_SCRATCH_FREE_BYTES_OVERRIDE HOMOLOGY_PERSISTENT_FREE_BYTES_OVERRIDE \
    HOMOLOGY_PRECOPY_FREE_BYTES_OVERRIDE
do
    if [[ -n "${!override_name:-}" ]]; then
        require_nonnegative_integer "$override_name" "${!override_name}"
    fi
done

FAILURE_STAGE="environment-activation"
if [[ "$TEST_MODE" == "1" && "${HOMOLOGY_SKIP_CONDA:-0}" == "1" ]]; then
    export CONDA_PREFIX="${CONDA_PREFIX:-$MMFP_ENV_DIR}"
    export PYTHON_BIN="${PYTHON_BIN:-python3}"
else
    [[ -x "$CONDA_EXE" ]] || { echo "Missing configured Conda executable: $CONDA_EXE" >&2; exit 1; }
    [[ -d "$MMFP_ENV_DIR" ]] || {
        echo "Required existing mmfp environment is absent: $MMFP_ENV_DIR" >&2
        echo "This wrapper will not create environments or install packages." >&2
        exit 1
    }
    eval "$("$CONDA_EXE" shell.bash hook)"
    conda activate "$MMFP_ENV_DIR"
    export PYTHON_BIN="$(command -v python)"
fi
EXPECTED_PREFIX="$(cd "$MMFP_ENV_DIR" && pwd -P)"
[[ -n "${CONDA_PREFIX:-}" && -d "$CONDA_PREFIX" ]] || {
    echo "Activated environment did not set a valid CONDA_PREFIX" >&2; exit 1;
}
OBSERVED_PREFIX="$(cd "$CONDA_PREFIX" && pwd -P)"
[[ "$EXPECTED_PREFIX" == "$OBSERVED_PREFIX" ]] || {
    echo "Activated CONDA_PREFIX mismatch: expected=$EXPECTED_PREFIX observed=$OBSERVED_PREFIX" >&2
    exit 1
}
[[ -n "$PYTHON_BIN" ]] || { echo "Activated environment has no Python executable" >&2; exit 1; }
check_signal

FAILURE_STAGE="capacity-before-input-staging"
LOCAL_INPUT_BYTES="$(local_input_bytes)"
MANIFEST_INPUT_BYTES="$(manifest_input_bytes)"
if (( MANIFEST_INPUT_BYTES > LOCAL_INPUT_BYTES )); then
    INPUT_BUDGET_BYTES="$MANIFEST_INPUT_BYTES"
else
    INPUT_BUDGET_BYTES="$LOCAL_INPUT_BYTES"
fi
UNIREF_BYTES="$(manifest_uniref_bytes)"
MMSEQS_ESTIMATE_BYTES="$(scaled_bytes "$UNIREF_BYTES" "$MMSEQS_WORK_MULTIPLIER")"
SCRATCH_REQUIRED_BYTES=$((INPUT_BUDGET_BYTES + MMSEQS_ESTIMATE_BYTES + SCRATCH_RESERVE_BYTES))
PUBLICATION_ESTIMATE_BYTES="$(scaled_bytes "$INPUT_BUDGET_BYTES" "$PUBLICATION_SAFETY_MULTIPLIER")"
PERSISTENT_REQUIRED_BYTES=$((PUBLICATION_ESTIMATE_BYTES + PERSISTENT_RESERVE_BYTES))
SCRATCH_FREE_BYTES="$(free_bytes "$WORK" "${HOMOLOGY_SCRATCH_FREE_BYTES_OVERRIDE:-}")"
PERSISTENT_FREE_BYTES="$(free_bytes "$PERSISTENT_RESULTS_ROOT" "${HOMOLOGY_PERSISTENT_FREE_BYTES_OVERRIDE:-}")"
printf '{"local_input_bytes":%s,"manifest_input_bytes":%s,"input_budget_bytes":%s,"uniref90_bytes":%s,"mmseqs_work_multiplier":%s,"mmseqs_estimate_bytes":%s,"scratch_reserve_bytes":%s,"scratch_required_bytes":%s,"scratch_free_bytes":%s,"publication_safety_multiplier":%s,"publication_estimate_bytes":%s,"persistent_reserve_bytes":%s,"persistent_required_bytes":%s,"persistent_free_bytes":%s,"estimates_exact":false}\n' \
    "$LOCAL_INPUT_BYTES" "$MANIFEST_INPUT_BYTES" "$INPUT_BUDGET_BYTES" "$UNIREF_BYTES" \
    "$MMSEQS_WORK_MULTIPLIER" "$MMSEQS_ESTIMATE_BYTES" "$SCRATCH_RESERVE_BYTES" \
    "$SCRATCH_REQUIRED_BYTES" "$SCRATCH_FREE_BYTES" "$PUBLICATION_SAFETY_MULTIPLIER" \
    "$PUBLICATION_ESTIMATE_BYTES" "$PERSISTENT_RESERVE_BYTES" "$PERSISTENT_REQUIRED_BYTES" \
    "$PERSISTENT_FREE_BYTES" > "$SCRATCH_OUTPUTS/logs/hpc_capacity_preflight.json"
require_capacity "Scratch before input staging/MMseqs" "$SCRATCH_FREE_BYTES" "$SCRATCH_REQUIRED_BYTES"
require_capacity "Persistent results before run" "$PERSISTENT_FREE_BYTES" "$PERSISTENT_REQUIRED_BYTES"
check_signal

FAILURE_STAGE="framework-checkout"
if [[ "$TEST_MODE" == "1" && -n "${HOMOLOGY_FRAMEWORK_DIR:-}" ]]; then
    [[ -d "$FRAMEWORK_DIR" ]] || { echo "Test framework directory is missing: $FRAMEWORK_DIR" >&2; exit 1; }
else
    git clone "$FRAMEWORK_REPO_URL" "$FRAMEWORK_DIR"
fi
if [[ "$FIXTURE_MODE_VALUE" != "1" ]]; then
    git -C "$FRAMEWORK_DIR" checkout --detach "$FRAMEWORK_REVISION"
    OBSERVED_FRAMEWORK_REVISION="$(git -C "$FRAMEWORK_DIR" rev-parse HEAD)"
    [[ "$OBSERVED_FRAMEWORK_REVISION" == "$FRAMEWORK_REVISION" ]] || {
        echo "Checked-out HEAD differs from FRAMEWORK_REVISION" >&2
        exit 1
    }
    if git -C "$FRAMEWORK_DIR" symbolic-ref -q HEAD >/dev/null 2>&1; then
        echo "Production framework checkout is not detached" >&2
        exit 1
    fi
    [[ -z "$(git -C "$FRAMEWORK_DIR" status --porcelain)" ]] || {
        echo "Production framework checkout is dirty" >&2
        exit 1
    }

    # shellcheck source=../../scripts/reproduction_common.sh
    source "$FRAMEWORK_DIR/scripts/reproduction_common.sh"
    validate_mmfp_env "$PYTHON_BIN"
fi

if [[ "$FIXTURE_MODE_VALUE" != "1" && "$DIAGNOSTIC_PILOT_VALUE" != "1" ]]; then
    for binding in \
        "FROZEN_INPUT_MANIFEST EXPECTED_FROZEN_INPUT_MANIFEST_SHA256" \
        "ATTRITION_POLICY EXPECTED_ATTRITION_POLICY_SHA256" \
        "PILOT_APPROVAL EXPECTED_PILOT_APPROVAL_SHA256" \
        "PILOT_COMPLETION_MARKER EXPECTED_PILOT_COMPLETION_MARKER_SHA256" \
        "PILOT_ATTRITION_REPORT EXPECTED_PILOT_ATTRITION_REPORT_SHA256" \
        "PILOT_TASK_CONTEXT EXPECTED_PILOT_TASK_CONTEXT_SHA256" \
        "PILOT_MEASUREMENT_EVIDENCE EXPECTED_PILOT_MEASUREMENT_EVIDENCE_SHA256"
    do
        read -r path_variable hash_variable <<< "$binding"
        verify_queued_hash "$path_variable" "$hash_variable"
    done
    PYTHONPATH="$FRAMEWORK_DIR/benchmark_builders/homology_cluster/src${PYTHONPATH:+:$PYTHONPATH}" \
        "$PYTHON_BIN" -m homology_cluster_benchmark authorize-array \
        --attrition-policy "$ATTRITION_POLICY" \
        --pilot-approval "$PILOT_APPROVAL" \
        --pilot-completion-marker "$PILOT_COMPLETION_MARKER" \
        --pilot-attrition-report "$PILOT_ATTRITION_REPORT" \
        --pilot-run-dir "$PILOT_RUN_DIR" \
        --pilot-task-context "$PILOT_TASK_CONTEXT" \
        --pilot-measurement-evidence "$PILOT_MEASUREMENT_EVIDENCE" \
        --frozen-input-manifest "$FROZEN_INPUT_MANIFEST" \
        --framework-revision "$FRAMEWORK_REVISION" \
        --uniprot-source-scope "$UNIPROT_SOURCE_SCOPE" \
        --split-policy "$SPLIT_POLICY" \
        --training-population "$TRAINING_POPULATION" \
        --expected-mmseqs-version "$EXPECTED_MMSEQS_VERSION" \
        --uniprot-release "${UNIPROT_RELEASE:-2026_02}" \
        --goa-release "${GOA_RELEASE:-234}" \
        --ontology-release "${ONTOLOGY_RELEASE:-releases/2026-06-15}"
fi

FAILURE_STAGE="input-staging"
stage_if_local UNIREF90_FASTA
stage_if_local IDMAPPING
stage_if_local UNIPROT_SPROT_SEQUENCES
stage_if_local UNIPROT_TREMBL_SEQUENCES
stage_if_local GOA
stage_if_local GO_OBO
stage_if_local FROZEN_INPUT_MANIFEST
stage_if_local ATTRITION_POLICY
stage_if_local ATTRITION_OVERRIDE
stage_if_local CLUSTER_ASSIGNMENTS
check_signal

MMSEQS_BIN="${MMSEQS_BIN:-mmseqs}"
if [[ -z "${CLUSTER_ASSIGNMENTS:-}" ]]; then
    if [[ "$MMSEQS_BIN" == */* ]]; then
        [[ -x "$MMSEQS_BIN" ]] || { echo "MMseqs2 is not executable: $MMSEQS_BIN" >&2; exit 1; }
    else
        MMSEQS_BIN="$(command -v "$MMSEQS_BIN" || true)"
        [[ -n "$MMSEQS_BIN" ]] || { echo "MMseqs2 is unavailable in the activated environment" >&2; exit 1; }
    fi
fi

export IDENTITY SPLIT_POLICY TRAINING_POPULATION SEED MIN_COUNT MMSEQS_BIN
export OUTPUT_ROOT="$SCRATCH_OUTPUTS/benchmark"
export TEMP_DIR="$SCRATCH_TEMP"
export PERSISTENT_RESULTS_ROOT
export THREADS="$ALLOCATED_SLOTS"
export REQUESTED_SLOTS=8 ALLOCATED_SLOTS
export UNIPROT_SOURCE_SCOPE FRAMEWORK_REVISION RUN_ID
export DIAGNOSTIC_PILOT="$DIAGNOSTIC_PILOT_VALUE"
export LOG_FILE="$SCRATCH_OUTPUTS/logs/homology_cluster_builder.log"
printf '{"job_id":"%s","sge_task_id":%s,"identity_percent":%s,"uniprot_source_scope":"%s","run_id":"%s","framework_revision":"%s","requested_smp_slots":8,"nslots":%s,"mmseqs_threads":%s,"maximum_array_cpu_slots_if_all_six_run":48}\n' \
    "$JOB_KEY" "$TASK_ID" "$IDENTITY" "$UNIPROT_SOURCE_SCOPE" "$RUN_ID" \
    "$FRAMEWORK_REVISION" "$ALLOCATED_SLOTS" "$THREADS" \
    > "$SCRATCH_OUTPUTS/logs/hpc_task_context.json"

FAILURE_STAGE="benchmark-builder"
cd "$FRAMEWORK_DIR"
if [[ -n "${HOMOLOGY_BUILDER_COMMAND:-}" ]]; then
    "$HOMOLOGY_BUILDER_COMMAND" &
else
    bash scripts/benchmark_generation/run_homology_cluster_benchmark.sh &
fi
ACTIVE_CHILD_PID=$!
set +e
wait_for_active_child
BUILDER_STATUS=$?
set -e
check_signal
if [[ "$BUILDER_STATUS" != "0" ]]; then
    echo "Benchmark builder failed with status $BUILDER_STATUS" >&2
    exit "$BUILDER_STATUS"
fi

FAILURE_STAGE="scratch-validation"
SCRATCH_RUN="$SCRATCH_OUTPUTS/benchmark/$RUN_RELATIVE_PATH"
validate_run "$SCRATCH_RUN"
check_signal

FAILURE_STAGE="persistent-capacity-before-copy"
OBSERVED_SCRATCH_USAGE_BYTES="$(tree_size "$SCRATCH_OUTPUTS")"
COPY_ESTIMATE_BYTES="$(scaled_bytes "$OBSERVED_SCRATCH_USAGE_BYTES" "$COPY_SAFETY_MULTIPLIER")"
PRECOPY_REQUIRED_BYTES=$((COPY_ESTIMATE_BYTES + PERSISTENT_RESERVE_BYTES))
PRECOPY_FREE_BYTES="$(free_bytes "$PERSISTENT_RESULTS_ROOT" "${HOMOLOGY_PRECOPY_FREE_BYTES_OVERRIDE:-}")"
printf '{"observed_scratch_usage_bytes":%s,"usage_measurement":"du_allocated_kib_rounded","copy_safety_multiplier":%s,"copy_estimate_bytes":%s,"persistent_reserve_bytes":%s,"persistent_required_bytes":%s,"persistent_free_bytes":%s,"estimates_exact":false}\n' \
    "$OBSERVED_SCRATCH_USAGE_BYTES" "$COPY_SAFETY_MULTIPLIER" "$COPY_ESTIMATE_BYTES" \
    "$PERSISTENT_RESERVE_BYTES" "$PRECOPY_REQUIRED_BYTES" "$PRECOPY_FREE_BYTES" \
    > "$SCRATCH_OUTPUTS/logs/hpc_capacity_precopy.json"
require_capacity "Persistent results before copy" "$PRECOPY_FREE_BYTES" "$PRECOPY_REQUIRED_BYTES"
check_signal

FAILURE_STAGE="copy-to-persistent-partial"
mkdir "$PARTIAL_RUN_ROOT"
printf '%s\n' "$JOB_KEY:$TASK_ID:$IDENTITY:$UNIPROT_SOURCE_SCOPE:$RUN_ID" \
    > "$PARTIAL_RUN_ROOT/.homology-task-owner"
PARTIAL_OWNED=1
set +e
copy_tree "$SCRATCH_OUTPUTS" "$PARTIAL_RUN_ROOT"
COPY_STATUS=$?
set -e
if [[ "$COPY_STATUS" != "0" ]]; then
    COPY_FAILED=1
    echo "Persistent result copy failed; diagnostics and mandatory scratch cleanup will be attempted" >&2
    exit "$COPY_STATUS"
fi
check_signal

FAILURE_STAGE="copied-validation"
COPIED_RUN="$PARTIAL_RUN_ROOT/benchmark/$RUN_RELATIVE_PATH"
validate_run "$COPIED_RUN"
check_signal

FAILURE_STAGE="atomic-final-rename"
check_signal
mv "$PARTIAL_RUN_ROOT" "$FINAL_RUN_ROOT"
PARTIAL_OWNED=0
FINAL_OWNED=1
check_signal
JOB_SUCCEEDED=1
check_signal
echo "Finished successfully: $FINAL_RUN_ROOT"
