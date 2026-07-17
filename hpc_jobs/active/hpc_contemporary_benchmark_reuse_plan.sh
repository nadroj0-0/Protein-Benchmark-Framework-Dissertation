#!/usr/bin/env bash
# UCL Grid Engine wrapper for comparing the completed contemporary CSVs with
# Zijian's canonical CAFA3 benchmark scope. Inputs are staged in scratch.

#$ -l tmem=32G
#$ -l tscratch=8G
#$ -l scratch0free=8G
#$ -l h_rt=12:0:0
#$ -j y
#$ -N cont_reuse
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
WORK="/scratch0/contemporary_benchmark_reuse_${JOB_TOKEN}"
RESULTS_ROOT="${RESULTS_ROOT:-$HOME/contemporary_benchmark_reuse_results}"
FINAL_RUN_ROOT="${RESULTS_ROOT}/${RUN_TAG}"
FAILED_RUN_ROOT="${FINAL_RUN_ROOT}.failed"
FRAMEWORK_REPO_URL="${FRAMEWORK_REPO_URL:-https://github.com/nadroj0-0/Protein-Benchmark-Framework-Dissertation.git}"
SUBMISSION_DIR="${SGE_O_WORKDIR:-$PWD}"
FRAMEWORK_COMMIT="${FRAMEWORK_COMMIT:-}"
FRAMEWORK_DIR="${WORK}/Protein-Benchmark-Framework-Dissertation"
WORKFLOW_WORK_DIR="${WORK}/reuse_work"
SCRATCH_RESULT_ROOT="${WORK}/result"
WORKFLOW_LOG="${WORK}/workflow.log"

WORK_OWNED=0
RESULTS_COPIED=0

die() {
  echo "ERROR: $*" >&2
  exit 2
}

git_in_dir() {
  local directory="$1"
  shift
  (cd "$directory" && git "$@")
}

copy_results() {
  local workflow_status="$1"
  local copy_status=0
  local destination="$FINAL_RUN_ROOT"
  local staging="${FINAL_RUN_ROOT}.staging-${JOB_TOKEN}"
  [[ "$RESULTS_COPIED" == "0" ]] || return 0
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
    [[ -f "$staging/plan/RUN_COMPLETE.json" ]] || copy_status=1
    [[ -f "$staging/plan/output_manifest.json" ]] || copy_status=1
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
  copy_results "$status" || copy_status=$?
  if [[ "$WORK_OWNED" == "1" ]]; then
    if [[ "$WORK" == /scratch0/contemporary_benchmark_reuse_* && ! -L "$WORK" ]]; then
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

[[ -n "${TARGET_BENCHMARK_DIR:-}" ]] || \
  die "Pass TARGET_BENCHMARK_DIR=/path/to/completed/contemporary/outputs with qsub -v"
[[ -d "$TARGET_BENCHMARK_DIR" ]] || \
  die "TARGET_BENCHMARK_DIR does not exist: $TARGET_BENCHMARK_DIR"
[[ ! -e "$WORK" ]] || die "Scratch path already exists: $WORK"
mkdir "$WORK"
WORK_OWNED=1
mkdir -p "$RESULTS_ROOT"

echo "Host          : $(hostname)"
echo "Job ID        : ${JOB_ID:-manual}"
echo "Scratch       : $WORK"
echo "Final output  : $FINAL_RUN_ROOT"
echo "Target CSVs   : $TARGET_BENCHMARK_DIR"
echo "Embedded CSVs : ${EMBEDDED_BENCHMARK_DIR:-Zenodo 7409660 download}"
echo "Started       : $(date)"
echo

if [[ -z "$FRAMEWORK_COMMIT" ]]; then
  [[ -d "$SUBMISSION_DIR/.git" ]] || \
    die "Submit from a clean framework Git checkout or pass FRAMEWORK_COMMIT"
  [[ -z "$(git_in_dir "$SUBMISSION_DIR" status --porcelain)" ]] || \
    die "Submission checkout has uncommitted changes; commit them before qsub"
  FRAMEWORK_COMMIT="$(git_in_dir "$SUBMISSION_DIR" rev-parse HEAD)"
fi
[[ "$FRAMEWORK_COMMIT" =~ ^[0-9a-fA-F]{40}$ ]] || \
  die "FRAMEWORK_COMMIT must be a complete 40-character Git commit"

echo "Framework Git: $FRAMEWORK_COMMIT"
git clone --no-checkout "$FRAMEWORK_REPO_URL" "$FRAMEWORK_DIR"
git_in_dir "$FRAMEWORK_DIR" checkout --detach "$FRAMEWORK_COMMIT"
[[ "$(git_in_dir "$FRAMEWORK_DIR" rev-parse HEAD)" == "$FRAMEWORK_COMMIT" ]] || \
  die "Scratch checkout does not match FRAMEWORK_COMMIT"
cd "$FRAMEWORK_DIR"

source scripts/reproduction_common.sh
export ARTIFACT_CATALOG="$CLI_ARTIFACT_CATALOG"
load_framework_paths "$FRAMEWORK_DIR"
activate_or_create_mmfp_env
PYTHON_BIN="$(command -v python)"

COMMAND=(
  bash "${FRAMEWORK_DIR}/scripts/verification/run_benchmark_reuse_plan.sh"
  --target-benchmark "${TARGET_BENCHMARK_NAME:-contemporary_2025_2026}=${TARGET_BENCHMARK_DIR}"
  --work-dir "$WORKFLOW_WORK_DIR"
  --output-dir "$SCRATCH_RESULT_ROOT"
  --python-bin "$PYTHON_BIN"
)
if [[ -n "${ARTIFACT_CATALOG:-}" ]]; then
  COMMAND+=(--artifact-catalog "$ARTIFACT_CATALOG")
fi
if [[ -n "${EMBEDDED_BENCHMARK_DIR:-}" ]]; then
  COMMAND+=(
    --embedded-benchmark
    "${EMBEDDED_BENCHMARK_NAME:-cafa3_zijian}=${EMBEDDED_BENCHMARK_DIR}"
  )
fi

echo "==> Running benchmark-level reuse comparison"
printf 'Command:'
printf ' %q' "${COMMAND[@]}"
printf '\n\n'

set +e
"${COMMAND[@]}" 2>&1 | tee "$WORKFLOW_LOG"
WORKFLOW_STATUS=${PIPESTATUS[0]}
set -e
[[ "$WORKFLOW_STATUS" == "0" ]] || exit "$WORKFLOW_STATUS"

copy_results 0
echo
echo "Finished successfully: $(date)"
echo "Read first: ${FINAL_RUN_ROOT}/plan/summary.md"
