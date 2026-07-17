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

CLI_ARTIFACT_CATALOG="${ARTIFACT_CATALOG:-}"
if [[ "${1:-}" == "--artifact-catalog" ]]; then
  [[ $# -ge 2 ]] || { echo "--artifact-catalog requires a path" >&2; exit 2; }
  CLI_ARTIFACT_CATALOG="$2"
  shift 2
fi
[[ $# -eq 0 ]] || { echo "Unknown argument: $1" >&2; exit 2; }

JOB_TOKEN="${JOB_ID:-manual_$$}"
RUN_TAG="${JOB_TOKEN}_$(date +%Y%m%d_%H%M%S)"
WORK="/scratch0/contemporary_embedding_inventory_${JOB_TOKEN}"
RESULTS_ROOT="${RESULTS_ROOT:-$HOME/contemporary_embedding_inventory_results}"
FINAL_RUN_ROOT="${RESULTS_ROOT}/${RUN_TAG}"
FAILED_RUN_ROOT="${FINAL_RUN_ROOT}.failed"
FRAMEWORK_REPO_URL="${FRAMEWORK_REPO_URL:-https://github.com/nadroj0-0/Protein-Benchmark-Framework-Dissertation.git}"
SUBMISSION_DIR="${SGE_O_WORKDIR:-$PWD}"
FRAMEWORK_COMMIT="${FRAMEWORK_COMMIT:-}"
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
  local workflow_status="$1"
  local copy_status=0
  local destination="$FINAL_RUN_ROOT"
  local staging="${FINAL_RUN_ROOT}.staging-${JOB_TOKEN}"
  if [[ "$RESULTS_COPIED" == "1" ]]; then
    return 0
  fi
  if [[ "$workflow_status" != "0" ]]; then
    destination="$FAILED_RUN_ROOT"
    staging="${FAILED_RUN_ROOT}.staging-${JOB_TOKEN}"
  fi
  if [[ -e "$destination" || -e "$staging" ]]; then
    echo "Refusing to overwrite result or staging directory: $destination / $staging" >&2
    return 1
  fi
  mkdir -p "$staging/logs" || return 1
  if [[ -d "$SCRATCH_RESULT_ROOT" ]]; then
    cp -a "$SCRATCH_RESULT_ROOT/." "$staging/" || copy_status=$?
  fi
  if [[ -f "$WORKFLOW_LOG" ]]; then
    cp -p "$WORKFLOW_LOG" "$staging/logs/workflow.log" || copy_status=$?
  fi
  if [[ "$copy_status" == "0" && "$workflow_status" == "0" ]]; then
    [[ -f "$staging/WORKFLOW_COMPLETE.json" ]] || copy_status=1
    [[ -f "$staging/inventory/RUN_COMPLETE.json" ]] || copy_status=1
    [[ -f "$staging/inventory/output_manifest.json" ]] || copy_status=1
  elif [[ "$workflow_status" != "0" ]]; then
    rm -f "$staging/WORKFLOW_COMPLETE.json"
    printf '{"complete":false,"workflow_exit_status":%s}\n' "$workflow_status" \
      > "$staging/WORKFLOW_FAILED.json" || copy_status=$?
  fi
  if [[ "$copy_status" == "0" ]]; then
    mv "$staging" "$destination" || copy_status=$?
  fi
  if [[ "$copy_status" == "0" ]]; then
    RESULTS_COPIED=1
    echo "==> Published result directory atomically: $destination"
  elif [[ -d "$staging" && ! -L "$staging" ]]; then
    rm -rf "$staging"
  fi
  return "$copy_status"
}

cleanup() {
  local status=$?
  local copy_status=0
  set +e
  echo
  echo "==> Final workflow status: $status"
  echo "==> Publishing reports and lists"
  copy_results "$status" || copy_status=$?

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

if [[ -z "$FRAMEWORK_COMMIT" ]]; then
  [[ -d "$SUBMISSION_DIR/.git" ]] || \
    die "Submit from a clean framework Git checkout or pass FRAMEWORK_COMMIT"
  [[ -z "$(git -C "$SUBMISSION_DIR" status --porcelain)" ]] || \
    die "Submission checkout has uncommitted changes; commit them before qsub"
  FRAMEWORK_COMMIT="$(git -C "$SUBMISSION_DIR" rev-parse HEAD)"
fi
[[ "$FRAMEWORK_COMMIT" =~ ^[0-9a-fA-F]{40}$ ]] || \
  die "FRAMEWORK_COMMIT must be a complete 40-character Git commit"

echo "Framework Git: $FRAMEWORK_COMMIT"
echo "==> Cloning the pinned dissertation framework into scratch"
git clone --no-checkout "$FRAMEWORK_REPO_URL" "$FRAMEWORK_DIR"
git -C "$FRAMEWORK_DIR" checkout --detach "$FRAMEWORK_COMMIT"
[[ "$(git -C "$FRAMEWORK_DIR" rev-parse HEAD)" == "$FRAMEWORK_COMMIT" ]] || \
  die "Scratch checkout does not match FRAMEWORK_COMMIT"
cd "$FRAMEWORK_DIR"

# Use the normal shared environment. On the configured UCL cluster this reuses
# the existing mmfp environment rather than creating an isolated environment.
source scripts/reproduction_common.sh
export ARTIFACT_CATALOG="$CLI_ARTIFACT_CATALOG"
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
if [[ -n "${ARTIFACT_CATALOG:-}" ]]; then
  COMMAND+=(--artifact-catalog "$ARTIFACT_CATALOG")
fi

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
if [[ -n "${PFP_REFERENCE_DIR:-}" ]]; then
  COMMAND+=(--pfp-reference-dir "$PFP_REFERENCE_DIR")
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

copy_results 0
echo
echo "Finished successfully: $(date)"
echo "Read first: ${FINAL_RUN_ROOT}/job_summary.md"
