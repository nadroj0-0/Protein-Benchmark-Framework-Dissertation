#!/usr/bin/env bash
# UCL Grid Engine wrapper for the complete CAFA3 regeneration/train/eval audit.

#$ -l tmem=24G
#$ -l tscratch=100G
#$ -l scratch0free=300G
#$ -l h_rt=96:0:0
#$ -l gpu=true
#$ -pe gpu 3
#$ -j y
#$ -N cafa3_full_repro
#$ -V

set -euo pipefail

usage() {
  cat <<'EOF'
Usage: qsub hpc_jobs/active/hpc_cafa3_full_from_scratch_reproduction.sh \
  [--results-root /absolute/path] \
  [--embedding-state-root /SAN/bioinf/bmpfp/embedding_states/cafa3_full_reproduction] \
  [--embedding-mode initial|resume] \
  [--text-cutoff-date YYYY-MM-DD]

The job downloads every external input into node-local scratch. No persistent
database, CSV, PFP, or published-embedding input path is required. Validated
regenerated arrays and authenticated AlphaFold PDBs are accumulated only under
the explicit persistent embedding-state root.
EOF
}

die() {
  echo "ERROR: $*" >&2
  exit 2
}

CLI_RESULTS_ROOT=""
TEXT_CUTOFF_DATE="2016-02-17"
EMBEDDING_STATE_ROOT="/SAN/bioinf/bmpfp/embedding_states/cafa3_full_reproduction"
EMBEDDING_MODE="initial"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --results-root)
      [[ $# -ge 2 ]] || die "--results-root requires a path"
      CLI_RESULTS_ROOT="$2"
      shift 2
      ;;
    --text-cutoff-date)
      [[ $# -ge 2 ]] || die "--text-cutoff-date requires YYYY-MM-DD"
      TEXT_CUTOFF_DATE="$2"
      shift 2
      ;;
    --embedding-state-root)
      [[ $# -ge 2 ]] || die "--embedding-state-root requires a path"
      EMBEDDING_STATE_ROOT="$2"
      shift 2
      ;;
    --embedding-mode)
      [[ $# -ge 2 ]] || die "--embedding-mode requires initial or resume"
      EMBEDDING_MODE="$2"
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

[[ "$TEXT_CUTOFF_DATE" =~ ^[0-9]{4}-[0-9]{2}-[0-9]{2}$ ]] || \
  die "--text-cutoff-date must use YYYY-MM-DD"
[[ "$EMBEDDING_MODE" == "initial" || "$EMBEDDING_MODE" == "resume" ]] || \
  die "--embedding-mode must be initial or resume"
mkdir -p "$EMBEDDING_STATE_ROOT"
[[ -d "$EMBEDDING_STATE_ROOT" ]] || die "Cannot access embedding state: $EMBEDDING_STATE_ROOT"

JOB_TOKEN="${JOB_ID:-manual_$$}"
RUN_TAG="${JOB_TOKEN}_$(date +%Y%m%d_%H%M%S)"
WORK="/scratch0/cafa3_full_reproduction_${JOB_TOKEN}"
RESULTS_ROOT="${CLI_RESULTS_ROOT:-${RESULTS_ROOT:-$HOME/cafa3_full_reproduction_results}}"
FINAL_RUN_ROOT="${RESULTS_ROOT}/${RUN_TAG}"
FAILED_RUN_ROOT="${FINAL_RUN_ROOT}.failed"
INCOMPLETE_RUN_ROOT="${FINAL_RUN_ROOT}.incomplete"
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
  if [[ "$workflow_status" == "0" && -f "$SCRATCH_RESULT_ROOT/GENERATION_INCOMPLETE.json" ]]; then
    destination="$INCOMPLETE_RUN_ROOT"
    staging="${INCOMPLETE_RUN_ROOT}.staging-${JOB_TOKEN}"
  elif [[ "$workflow_status" != "0" ]]; then
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
    if [[ -f "$staging/GENERATION_INCOMPLETE.json" ]]; then
      [[ -f "$staging/reports/embedding_state/needs_retry.tsv" ]] || copy_status=1
      [[ -f "$staging/reports/embedding_state/coverage.json" ]] || copy_status=1
      [[ ! -f "$staging/WORKFLOW_COMPLETE.json" ]] || copy_status=1
    else
      [[ -f "$staging/WORKFLOW_COMPLETE.json" ]] || copy_status=1
      [[ -f "$staging/cafa3_full_reproduction_report.md" ]] || copy_status=1
      [[ -f "$staging/cafa3_full_reproduction_report.json" ]] || copy_status=1
      [[ -f "$staging/reports/embedding_comparison_summary.json" ]] || copy_status=1
      [[ -f "$staging/reports/evaluation/reproduction_summary.json" ]] || copy_status=1
    fi
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
    if [[ "$WORK" == /scratch0/cafa3_full_reproduction_* && ! -L "$WORK" ]]; then
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
echo "Embedding state   : $EMBEDDING_STATE_ROOT"
echo "Embedding mode    : $EMBEDDING_MODE"
echo "PFP commit        : $PFP_COMMIT"
echo "Text cutoff       : $TEXT_CUTOFF_DATE"
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
[[ "$FRAMEWORK_COMMIT" =~ ^[0-9a-fA-F]{40}$ ]] || \
  die "FRAMEWORK_COMMIT must be a full commit"
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

COMMAND=(
  bash "$FRAMEWORK_DIR/scripts/reproduction/run_cafa3_full_from_scratch_reproduction.sh"
  --pfp-root "$PFP_DIR"
  --work-dir "$WORKFLOW_WORK_DIR"
  --output-dir "$SCRATCH_RESULT_ROOT"
  --embedding-state-root "$EMBEDDING_STATE_ROOT"
  --embedding-mode "$EMBEDDING_MODE"
  --text-cutoff-date "$TEXT_CUTOFF_DATE"
)

echo "==> Running complete CAFA3 reproduction"
printf 'Command:'
printf ' %q' "${COMMAND[@]}"
printf '\n\n'

set +e
PYTHON_BIN="$PYTHON_BIN" \
PREFLIGHT_PER_SPLIT="${PREFLIGHT_PER_SPLIT:-2}" \
"${COMMAND[@]}" 2>&1 | tee "$WORKFLOW_LOG"
WORKFLOW_STATUS=${PIPESTATUS[0]}
set -e
[[ "$WORKFLOW_STATUS" == "0" ]] || exit "$WORKFLOW_STATUS"

copy_results 0
echo
if [[ -f "$SCRATCH_RESULT_ROOT/GENERATION_INCOMPLETE.json" ]]; then
  echo "Finished with durable incomplete embedding state: $(date)"
  echo "Retry pairs: $EMBEDDING_STATE_ROOT/needs_retry.tsv"
else
  echo "Finished successfully: $(date)"
  echo "Report: $FINAL_RUN_ROOT/cafa3_full_reproduction_report.md"
fi
