#!/usr/bin/env bash
# UCL Grid Engine wrapper. Resource directives are provisional until a measured smoke run.

#$ -l tmem=64G
#$ -l tscratch=200G
#$ -l scratch0free=200G
#$ -l h_rt=72:0:0
#$ -j y
#$ -N homology_cluster
#$ -V

set -euo pipefail

IDENTITY="${IDENTITY:-}"
case "$IDENTITY" in
    30|25|20|15|10|5) ;;
    *) echo "IDENTITY must be exactly one of 30,25,20,15,10,5" >&2; exit 2 ;;
esac
TRAINING_POPULATION="${TRAINING_POPULATION:-annotated-only}"
if [[ "$TRAINING_POPULATION" != "annotated-only" ]]; then
    echo "TRAINING_POPULATION=$TRAINING_POPULATION is unsupported; only annotated-only is implemented." >&2
    echo "No zero-negative, homology-transfer, or representative-label policy is authorized." >&2
    exit 2
fi
FIXTURE_MODE_VALUE="${FIXTURE_MODE:-0}"
TEST_MODE="${HOMOLOGY_WRAPPER_TEST_MODE:-0}"
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
JOB_KEY="${JOB_ID:-manual}"
WORK_BASE="${WORK_BASE:-/scratch0}"
WORK="$WORK_BASE/homology_cluster_${JOB_KEY}_${IDENTITY}"
FRAMEWORK_REPO_URL="${FRAMEWORK_REPO_URL:-https://github.com/nadroj0-0/Protein-Benchmark-Framework-Dissertation.git}"
FRAMEWORK_REVISION="${FRAMEWORK_REVISION:-}"
FRAMEWORK_DIR="${HOMOLOGY_FRAMEWORK_DIR:-$WORK/Protein-Benchmark-Framework-Dissertation}"
SCRATCH_INPUTS="$WORK/inputs"
SCRATCH_OUTPUTS="$WORK/run"
SCRATCH_TEMP="$WORK/tmp"
RESULTS_ROOT="${RESULTS_ROOT:-$HOME/homology_cluster_benchmark_results}"
PERSISTENT_RESULTS_ROOT="${PERSISTENT_RESULTS_ROOT:-$RESULTS_ROOT}"
RUN_TAG="identity_${IDENTITY}_${SPLIT_POLICY}_${TRAINING_POPULATION}_seed_${SEED}_min_count_${MIN_COUNT}_${JOB_KEY}_$(date +%Y%m%d_%H%M%S)"
FINAL_RUN_ROOT="$PERSISTENT_RESULTS_ROOT/$RUN_TAG"
PARTIAL_RUN_ROOT="${FINAL_RUN_ROOT}.partial-${JOB_KEY}"
FAILED_RUN_ROOT="${FINAL_RUN_ROOT}.failed"
FAILED_PARTIAL_ROOT="${FAILED_RUN_ROOT}.partial-${JOB_KEY}"
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

identity_directory() {
    if [[ "$IDENTITY" == "5" ]]; then
        printf 'identity_05'
    else
        printf 'identity_%s' "$IDENTITY"
    fi
}

RUN_RELATIVE_PATH="$(identity_directory)/$SPLIT_POLICY/$TRAINING_POPULATION/seed_$SEED/min_count_$MIN_COUNT"

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
        UNIREF90_FASTA IDMAPPING UNIPROT_SEQUENCES GOA GO_OBO \
        FROZEN_INPUT_MANIFEST CLUSTER_ASSIGNMENTS
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
    [[ ! -e "$FAILED_RUN_ROOT" && ! -e "$FAILED_PARTIAL_ROOT" ]] || {
        echo "Failure destination already exists; preserving scratch at: $WORK" >&2
        return 1
    }
    write_failure_metadata "$SCRATCH_OUTPUTS/logs" "$exit_status" || true
    if [[ -e "$FINAL_RUN_ROOT" ]]; then
        find "$FINAL_RUN_ROOT" -name RUN_COMPLETE.json -type f -exec rm -f {} \;
        mv "$FINAL_RUN_ROOT" "$FAILED_PARTIAL_ROOT" || return 1
    elif [[ -e "$PARTIAL_RUN_ROOT" ]]; then
        find "$PARTIAL_RUN_ROOT" -name RUN_COMPLETE.json -type f -exec rm -f {} \;
        mv "$PARTIAL_RUN_ROOT" "$FAILED_PARTIAL_ROOT" || return 1
    else
        mkdir -p "$FAILED_PARTIAL_ROOT"
        if [[ -d "$SCRATCH_OUTPUTS" ]]; then
            copy_tree "$SCRATCH_OUTPUTS" "$FAILED_PARTIAL_ROOT" || return 1
        fi
    fi
    find "$FAILED_PARTIAL_ROOT" -name RUN_COMPLETE.json -type f -exec rm -f {} \;
    write_failure_metadata "$FAILED_PARTIAL_ROOT" "$exit_status" || return 1
    mv "$FAILED_PARTIAL_ROOT" "$FAILED_RUN_ROOT" || return 1
    echo "Failure diagnostics published atomically: $FAILED_RUN_ROOT" >&2
    if [[ "$COPY_FAILED" == "1" ]]; then
        echo "Result copy failed; preserving recoverable scratch exactly at: $WORK" >&2
    fi
    return 0
}

cleanup() {
    local status=$?
    trap - EXIT
    set +e
    if [[ "$SIGNAL_STATUS" != "0" ]]; then
        JOB_SUCCEEDED=0
    fi
    if [[ "$JOB_SUCCEEDED" == "1" ]]; then
        cd "$HOME"
        rm -rf "$WORK"
    else
        publish_failure "$status"
        local failure_published=$?
        if [[ "$COPY_FAILED" == "1" ]]; then
            echo "Scratch retained for recovery: $WORK" >&2
        elif [[ "$failure_published" -eq 0 ]]; then
            cd "$HOME"
            rm -rf "$WORK"
        else
            echo "Scratch retained for recovery: $WORK" >&2
            if [[ "$status" -eq 0 ]]; then
                status=1
            fi
        fi
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
    local destination_dir="$SCRATCH_INPUTS/$variable_name"
    mkdir -p "$destination_dir"
    local destination="$destination_dir/$(basename "$value")"
    cp -Lp "$value" "$destination"
    printf -v "$variable_name" '%s' "$destination"
    export "$variable_name"
}

echo "Host        : $(hostname)"
echo "Job ID      : $JOB_KEY"
echo "Identity    : $IDENTITY%"
echo "Split policy: $SPLIT_POLICY"
echo "Scratch     : $WORK"
echo "Final output: $FINAL_RUN_ROOT"
echo "Resources   : provisional 64G tmem / 200G scratch / 72h"

mkdir -p "$WORK" "$SCRATCH_INPUTS" "$SCRATCH_OUTPUTS/logs" "$SCRATCH_TEMP" "$RESULTS_ROOT" "$PERSISTENT_RESULTS_ROOT"

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
    if [[ -n "$FRAMEWORK_REVISION" ]]; then
        git -C "$FRAMEWORK_DIR" checkout --detach "$FRAMEWORK_REVISION"
    fi
fi

FAILURE_STAGE="input-staging"
stage_if_local UNIREF90_FASTA
stage_if_local IDMAPPING
stage_if_local UNIPROT_SEQUENCES
stage_if_local GOA
stage_if_local GO_OBO
stage_if_local FROZEN_INPUT_MANIFEST
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
ALLOCATED_SLOTS="${NSLOTS:-1}"
REQUESTED_THREADS="${THREADS:-$ALLOCATED_SLOTS}"
[[ "$ALLOCATED_SLOTS" =~ ^[1-9][0-9]*$ && "$REQUESTED_THREADS" =~ ^[1-9][0-9]*$ ]] || {
    echo "NSLOTS and THREADS must be positive integers" >&2; exit 1;
}
if (( REQUESTED_THREADS > ALLOCATED_SLOTS )); then
    echo "THREADS=$REQUESTED_THREADS exceeds allocated NSLOTS=$ALLOCATED_SLOTS" >&2
    exit 1
fi
export THREADS="$REQUESTED_THREADS"
export LOG_FILE="$SCRATCH_OUTPUTS/logs/homology_cluster_builder.log"

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
[[ ! -e "$FINAL_RUN_ROOT" && ! -e "$PARTIAL_RUN_ROOT" && ! -e "$FAILED_RUN_ROOT" ]] || {
    echo "Refusing to overwrite persistent final/partial/failed output" >&2; exit 1;
}
set +e
copy_tree "$SCRATCH_OUTPUTS" "$PARTIAL_RUN_ROOT"
COPY_STATUS=$?
set -e
if [[ "$COPY_STATUS" != "0" ]]; then
    COPY_FAILED=1
    echo "Persistent result copy failed; scratch will be preserved" >&2
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
check_signal
JOB_SUCCEEDED=1
check_signal
echo "Finished successfully: $FINAL_RUN_ROOT"
