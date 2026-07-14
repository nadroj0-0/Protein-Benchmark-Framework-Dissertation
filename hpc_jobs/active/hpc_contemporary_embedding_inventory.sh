#!/usr/bin/env bash
# UCL Grid Engine wrapper for a real contemporary-benchmark embedding inventory.
# Published archives are downloaded/extracted only in scratch; reports are copied home.

#$ -l tmem=32G
#$ -l tscratch=20G
#$ -l scratch0free=20G
#$ -l h_rt=24:0:0
#$ -j y
#$ -N cont_emb_inv
#$ -V

set -euo pipefail

JOB_TOKEN="${JOB_ID:-manual_$$}"
RUN_TAG="${JOB_TOKEN}_$(date +%Y%m%d_%H%M%S)"
WORK="/scratch0/contemporary_embedding_inventory_${JOB_TOKEN}"
RESULTS_ROOT="${RESULTS_ROOT:-$HOME/contemporary_embedding_inventory_results}"
FINAL_RUN_ROOT="${RESULTS_ROOT}/${RUN_TAG}"
FRAMEWORK_REPO_URL="${FRAMEWORK_REPO_URL:-https://github.com/nadroj0-0/Protein-Benchmark-Framework-Dissertation.git}"
FRAMEWORK_DIR="${WORK}/Protein-Benchmark-Framework-Dissertation"
WORKFLOW_WORK_DIR="${WORK}/inventory_work"
SCRATCH_RESULT_ROOT="${WORK}/result"
WORKFLOW_LOG="${WORK}/workflow.log"

WORK_OWNED=0
RESULTS_COPIED=0

die() {
  echo "ERROR: $*" >&2
  exit 2
}

copy_results() {
  local copy_status=0
  if [[ "$RESULTS_COPIED" == "1" ]]; then
    return 0
  fi
  if [[ -e "$FINAL_RUN_ROOT" ]]; then
    echo "Refusing to overwrite result directory: $FINAL_RUN_ROOT" >&2
    return 1
  fi
  mkdir -p "$FINAL_RUN_ROOT/logs" || return 1
  if [[ -d "$SCRATCH_RESULT_ROOT" ]]; then
    cp -a "$SCRATCH_RESULT_ROOT/." "$FINAL_RUN_ROOT/" || copy_status=$?
  fi
  if [[ -f "$WORKFLOW_LOG" ]]; then
    cp -p "$WORKFLOW_LOG" "$FINAL_RUN_ROOT/logs/workflow.log" || copy_status=$?
  fi
  if [[ "$copy_status" == "0" ]]; then
    RESULTS_COPIED=1
  fi
  return "$copy_status"
}

cleanup() {
  local status=$?
  local copy_status=0
  set +e
  echo
  echo "==> Final workflow status: $status"
  echo "==> Copying reports and lists to: $FINAL_RUN_ROOT"
  copy_results || copy_status=$?

  if [[ "$WORK_OWNED" == "1" ]]; then
    if [[ "$WORK" == /scratch0/contemporary_embedding_inventory_* && ! -L "$WORK" ]]; then
      echo "==> Removing job-owned scratch directory: $WORK"
      cd "$HOME"
      rm -rf "$WORK"
    else
      echo "Refusing unsafe scratch cleanup path: $WORK" >&2
      [[ "$status" != "0" ]] || status=1
    fi
  fi

  if [[ "$status" == "0" && "$copy_status" != "0" ]]; then
    status="$copy_status"
  fi
  exit "$status"
}
trap cleanup EXIT
trap 'echo "Received termination signal"; exit 130' INT TERM

if [[ -e "$WORK" ]]; then
  die "Scratch path already exists; refusing to reuse it: $WORK"
fi
mkdir "$WORK"
WORK_OWNED=1
mkdir -p "$RESULTS_ROOT"

echo "Host          : $(hostname)"
echo "Job ID        : ${JOB_ID:-manual}"
echo "Scratch       : $WORK"
echo "Final output  : $FINAL_RUN_ROOT"
echo "Benchmark dir : ${BENCHMARK_DIR:-individual CSV overrides}"
echo "Policy        : ${INVENTORY_POLICY:-maximize-coverage}"
echo "Report level  : ${REPORT_LEVEL:-compact}"
echo "Started       : $(date)"
echo

echo "==> Cloning the dissertation framework into scratch"
git clone "$FRAMEWORK_REPO_URL" "$FRAMEWORK_DIR"
cd "$FRAMEWORK_DIR"

# Use the normal shared environment. On the configured UCL cluster this reuses
# the existing mmfp environment rather than creating an isolated environment.
source scripts/reproduction_common.sh
load_framework_paths "$FRAMEWORK_DIR"
activate_or_create_mmfp_env
PYTHON_BIN="$(command -v python)"

COMMAND=(
  bash "${FRAMEWORK_DIR}/scripts/verification/run_contemporary_embedding_inventory.sh"
  --work-dir "$WORKFLOW_WORK_DIR"
  --output-dir "$SCRATCH_RESULT_ROOT"
  --python-bin "$PYTHON_BIN"
  --config "${INVENTORY_CONFIG:-${FRAMEWORK_DIR}/configs/embedding_inventory.contemporary.json}"
  --policy "${INVENTORY_POLICY:-maximize-coverage}"
  --report-level "${REPORT_LEVEL:-compact}"
)

if [[ -n "${BENCHMARK_DIR:-}" ]]; then
  COMMAND+=(--benchmark-dir "$BENCHMARK_DIR")
fi

add_csv_override() {
  local variable_name="$1"
  local csv_name="$2"
  local value="${!variable_name:-}"
  if [[ -n "$value" ]]; then
    COMMAND+=(--benchmark-csv "${csv_name}=${value}")
  fi
}

add_csv_override BP_TRAINING_CSV bp-training.csv
add_csv_override BP_VALIDATION_CSV bp-validation.csv
add_csv_override BP_TEST_CSV bp-test.csv
add_csv_override CC_TRAINING_CSV cc-training.csv
add_csv_override CC_VALIDATION_CSV cc-validation.csv
add_csv_override CC_TEST_CSV cc-test.csv
add_csv_override MF_TRAINING_CSV mf-training.csv
add_csv_override MF_VALIDATION_CSV mf-validation.csv
add_csv_override MF_TEST_CSV mf-test.csv

if [[ -n "${SOURCE_BENCHMARK_DIR:-}" ]]; then
  COMMAND+=(--source-benchmark-dir "$SOURCE_BENCHMARK_DIR")
fi
if [[ -n "${PUBLISHED_EMBEDDING_ARCHIVE_DIR:-}" ]]; then
  COMMAND+=(--embedding-archive-dir "$PUBLISHED_EMBEDDING_ARCHIVE_DIR")
fi
if [[ -n "${ALIASES_FILE:-}" ]]; then
  COMMAND+=(--aliases "$ALIASES_FILE")
fi

if [[ -z "${BENCHMARK_DIR:-}" ]]; then
  for variable_name in \
    BP_TRAINING_CSV BP_VALIDATION_CSV BP_TEST_CSV \
    CC_TRAINING_CSV CC_VALIDATION_CSV CC_TEST_CSV \
    MF_TRAINING_CSV MF_VALIDATION_CSV MF_TEST_CSV; do
    [[ -n "${!variable_name:-}" ]] || \
      die "Supply BENCHMARK_DIR or all nine per-CSV environment variables (missing $variable_name)"
  done
fi

echo "==> Running the reusable inventory workflow"
printf 'Command:'
printf ' %q' "${COMMAND[@]}"
printf '\n\n'

set +e
"${COMMAND[@]}" 2>&1 | tee "$WORKFLOW_LOG"
WORKFLOW_STATUS=${PIPESTATUS[0]}
set -e

if [[ "$WORKFLOW_STATUS" != "0" ]]; then
  exit "$WORKFLOW_STATUS"
fi

copy_results
echo
echo "Finished successfully: $(date)"
echo "Read first: ${FINAL_RUN_ROOT}/job_summary.md"
