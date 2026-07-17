#!/usr/bin/env bash
# UCL Grid Engine wrapper for contemporary embedding generation and assembly.

#$ -l tmem=24G
#$ -l tscratch=100G
#$ -l scratch0free=300G
#$ -l h_rt=96:0:0
# Optional Zeus-only restriction (currently disabled because both queues are offline):
# #$ -q gpu.q@zeus1.local,gpu.q@zeus2.local
#$ -l gpu=true
#$ -pe gpu 3
#$ -j y
#$ -N cont_embeddings
#$ -V

set -euo pipefail

usage() {
  cat <<'EOF'
Usage: qsub hpc_jobs/active/hpc_contemporary_embedding_generation.sh \
  --target-benchmark-dir /absolute/path/to/outputs \
  --reuse-plan-dir /absolute/path/to/plan \
  [--results-root /absolute/output/path] \
  [--text-cutoff-date YYYY-MM-DD]
EOF
}

die() {
  echo "ERROR: $*" >&2
  exit 2
}

TARGET_BENCHMARK_DIR=""
REUSE_PLAN_DIR=""
CLI_RESULTS_ROOT=""
CLI_TEXT_CUTOFF_DATE=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --target-benchmark-dir)
      [[ $# -ge 2 ]] || die "--target-benchmark-dir requires a path"
      TARGET_BENCHMARK_DIR="$2"
      shift 2
      ;;
    --reuse-plan-dir)
      [[ $# -ge 2 ]] || die "--reuse-plan-dir requires a path"
      REUSE_PLAN_DIR="$2"
      shift 2
      ;;
    --results-root)
      [[ $# -ge 2 ]] || die "--results-root requires a path"
      CLI_RESULTS_ROOT="$2"
      shift 2
      ;;
    --text-cutoff-date)
      [[ $# -ge 2 ]] || die "--text-cutoff-date requires YYYY-MM-DD"
      CLI_TEXT_CUTOFF_DATE="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      usage >&2
      die "Unknown argument: $1"
      ;;
  esac
done

[[ -n "$TARGET_BENCHMARK_DIR" ]] || { usage >&2; die "--target-benchmark-dir is required"; }
[[ -n "$REUSE_PLAN_DIR" ]] || { usage >&2; die "--reuse-plan-dir is required"; }

JOB_TOKEN="${JOB_ID:-manual_$$}"
RUN_TAG="${JOB_TOKEN}_$(date +%Y%m%d_%H%M%S)"
WORK="/scratch0/contemporary_embedding_generation_${JOB_TOKEN}"
RESULTS_ROOT="${CLI_RESULTS_ROOT:-${RESULTS_ROOT:-$HOME/contemporary_embedding_generation_results}}"
FINAL_RUN_ROOT="${RESULTS_ROOT}/${RUN_TAG}"
FAILED_RUN_ROOT="${FINAL_RUN_ROOT}.failed"
FRAMEWORK_REPO_URL="${FRAMEWORK_REPO_URL:-https://github.com/nadroj0-0/Protein-Benchmark-Framework-Dissertation.git}"
FRAMEWORK_COMMIT="${FRAMEWORK_COMMIT:-}"
FRAMEWORK_DIR="$WORK/Protein-Benchmark-Framework-Dissertation"
PFP_REPO_URL="${PFP_REPO_URL:-https://github.com/psipred/PFP.git}"
PFP_COMMIT="${PFP_COMMIT:-1e04fd6d6d3c40458fd41ec1a881ed6e24de768e}"
PFP_DIR="$WORK/PFP"
WORKFLOW_WORK_DIR="$WORK/workflow"
SCRATCH_RESULT_ROOT="$WORK/result"
WORKFLOW_LOG="$WORK/workflow.log"
SUBMISSION_DIR="${SGE_O_WORKDIR:-$PWD}"

WORK_OWNED=0
RESULTS_COPIED=0

git_in_dir() {
  local directory="$1"
  shift
  (cd "$directory" && git "$@")
}

copy_results() {
  local workflow_status="$1"
  local destination="$FINAL_RUN_ROOT"
  local staging="${FINAL_RUN_ROOT}.staging-${JOB_TOKEN}"
  local copy_status=0
  [[ "$RESULTS_COPIED" == "0" ]] || return 0
  if [[ "$workflow_status" != "0" ]]; then
    destination="$FAILED_RUN_ROOT"
    staging="${FAILED_RUN_ROOT}.staging-${JOB_TOKEN}"
  fi
  [[ ! -e "$destination" && ! -e "$staging" ]] || return 1
  mkdir -p "$staging/logs" || return 1
  if [[ -d "$SCRATCH_RESULT_ROOT" ]]; then
    cp -a "$SCRATCH_RESULT_ROOT/." "$staging/" || copy_status=$?
  fi
  if [[ -f "$WORKFLOW_LOG" ]]; then
    cp -p "$WORKFLOW_LOG" "$staging/logs/workflow.log" || copy_status=$?
  fi
  if [[ "$workflow_status" == "0" ]]; then
    [[ -f "$staging/WORKFLOW_COMPLETE.json" ]] || copy_status=1
    [[ -f "$staging/archives/contemporary_embedding_cache.tar.gz" ]] || copy_status=1
    [[ -f "$staging/reports/assembly/assembly_summary.json" ]] || copy_status=1
  else
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
  trap - EXIT
  set +e
  echo
  echo "==> Final workflow status: $status"
  copy_results "$status" || copy_status=$?
  if [[ "$WORK_OWNED" == "1" ]]; then
    if [[ "$WORK" == /scratch0/contemporary_embedding_generation_* && ! -L "$WORK" ]]; then
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

[[ -d "$TARGET_BENCHMARK_DIR" ]] || die "Target benchmark does not exist: $TARGET_BENCHMARK_DIR"
[[ -d "$REUSE_PLAN_DIR" ]] || die "Reuse plan does not exist: $REUSE_PLAN_DIR"
[[ ! -e "$WORK" ]] || die "Scratch path already exists: $WORK"
mkdir -p "$WORK/tmp"
WORK_OWNED=1
export TMPDIR="$WORK/tmp"
export TMP="$WORK/tmp"
export TEMP="$WORK/tmp"
mkdir -p "$RESULTS_ROOT"

echo "Host              : $(hostname)"
echo "Job ID            : ${JOB_ID:-manual}"
echo "Scratch           : $WORK"
echo "Final output      : $FINAL_RUN_ROOT"
echo "Target benchmark  : $TARGET_BENCHMARK_DIR"
echo "Reuse plan        : $REUSE_PLAN_DIR"
echo "PFP commit        : $PFP_COMMIT"
echo "Preflight/split   : ${PREFLIGHT_PER_SPLIT:-2}"
echo "Started           : $(date)"
echo

if [[ -z "$FRAMEWORK_COMMIT" ]]; then
  [[ -d "$SUBMISSION_DIR/.git" ]] || \
    die "Submit from a clean framework checkout or pass FRAMEWORK_COMMIT"
  [[ -z "$(git_in_dir "$SUBMISSION_DIR" status --porcelain)" ]] || \
    die "Submission checkout has uncommitted changes"
  FRAMEWORK_COMMIT="$(git_in_dir "$SUBMISSION_DIR" rev-parse HEAD)"
fi
[[ "$FRAMEWORK_COMMIT" =~ ^[0-9a-fA-F]{40}$ ]] || die "FRAMEWORK_COMMIT must be a full commit"
[[ "$PFP_COMMIT" =~ ^[0-9a-fA-F]{40}$ ]] || die "PFP_COMMIT must be a full commit"

echo "==> Cloning pinned framework and PFP checkouts into scratch"
git clone --no-checkout "$FRAMEWORK_REPO_URL" "$FRAMEWORK_DIR"
git_in_dir "$FRAMEWORK_DIR" checkout --detach "$FRAMEWORK_COMMIT"
git clone --no-checkout "$PFP_REPO_URL" "$PFP_DIR"
git_in_dir "$PFP_DIR" checkout --detach "$PFP_COMMIT"
[[ "$(git_in_dir "$FRAMEWORK_DIR" rev-parse HEAD)" == "$FRAMEWORK_COMMIT" ]] || \
  die "Framework checkout mismatch"
[[ "$(git_in_dir "$PFP_DIR" rev-parse HEAD)" == "$PFP_COMMIT" ]] || \
  die "PFP checkout mismatch"

cd "$FRAMEWORK_DIR"
source scripts/reproduction_common.sh
load_framework_paths "$FRAMEWORK_DIR"
activate_or_create_mmfp_env
PYTHON_BIN="$(command -v python)"

TEXT_CUTOFF_DATE="${CLI_TEXT_CUTOFF_DATE:-${TEXT_CUTOFF_DATE:-}}"
if [[ -z "$TEXT_CUTOFF_DATE" ]]; then
  BUILD_MANIFEST="$(dirname "$TARGET_BENCHMARK_DIR")/reports/build_manifest.json"
  [[ -f "$BUILD_MANIFEST" ]] || \
    die "Cannot derive text cutoff; pass TEXT_CUTOFF_DATE or provide the completed run layout"
  TEXT_CUTOFF_DATE="$($PYTHON_BIN - "$BUILD_MANIFEST" <<'PY'
import json
import sys
value = str(json.load(open(sys.argv[1]))["t0_cutoff"])
if len(value) != 8 or not value.isdigit():
    raise SystemExit("Invalid t0_cutoff in benchmark manifest: " + value)
print(value[:4] + "-" + value[4:6] + "-" + value[6:])
PY
)"
fi

COMMAND=(
  bash "$FRAMEWORK_DIR/scripts/embeddings/run_contemporary_embedding_generation.sh"
  --target-benchmark-dir "$TARGET_BENCHMARK_DIR"
  --reuse-plan-dir "$REUSE_PLAN_DIR"
  --pfp-root "$PFP_DIR"
  --work-dir "$WORKFLOW_WORK_DIR"
  --output-dir "$SCRATCH_RESULT_ROOT"
  --text-cutoff-date "$TEXT_CUTOFF_DATE"
)

echo "==> Running contemporary embedding generation"
printf 'Command:'
printf ' %q' "${COMMAND[@]}"
printf '\n\n'

set +e
PYTHON_BIN="$PYTHON_BIN" \
PFP_COMMIT="$PFP_COMMIT" \
FRAMEWORK_COMMIT="$FRAMEWORK_COMMIT" \
PREFLIGHT_PER_SPLIT="${PREFLIGHT_PER_SPLIT:-2}" \
"${COMMAND[@]}" 2>&1 | tee "$WORKFLOW_LOG"
WORKFLOW_STATUS=${PIPESTATUS[0]}
set -e
[[ "$WORKFLOW_STATUS" == "0" ]] || exit "$WORKFLOW_STATUS"

copy_results 0
echo
echo "Finished successfully: $(date)"
echo "Final archive: $FINAL_RUN_ROOT/archives/contemporary_embedding_cache.tar.gz"
