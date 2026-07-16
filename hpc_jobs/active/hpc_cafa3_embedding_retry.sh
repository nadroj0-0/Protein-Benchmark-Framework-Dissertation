#!/usr/bin/env bash
# UCL Grid Engine wrapper for one modality-level CAFA3 embedding retry.

#$ -l tmem=24G
#$ -l tscratch=100G
#$ -l scratch0free=300G
#$ -l h_rt=96:0:0
#$ -l gpu=true
#$ -pe gpu 1
#$ -j y
#$ -N cafa3_emb_retry
#$ -V

set -euo pipefail

usage() {
  cat <<'EOF'
Usage: qsub hpc_jobs/active/hpc_cafa3_embedding_retry.sh \
  --modality sequence|text|structure|ppi \
  [--embedding-state-root /SAN/bioinf/bmpfp/embedding_states/cafa3_full_reproduction] \
  [--results-root /absolute/path] \
  [--text-cutoff-date YYYY-MM-DD]

This wrapper never submits another job. It retries only missing pairs for one
modality, publishes valid arrays into the single persistent state, copies only
compact reports home, and always removes job-owned scratch.
EOF
}

die() { echo "ERROR: $*" >&2; exit 2; }

MODALITY=""
EMBEDDING_STATE_ROOT="/SAN/bioinf/bmpfp/embedding_states/cafa3_full_reproduction"
CLI_RESULTS_ROOT=""
TEXT_CUTOFF_DATE="2016-02-17"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --modality) MODALITY="$2"; shift 2 ;;
    --embedding-state-root) EMBEDDING_STATE_ROOT="$2"; shift 2 ;;
    --results-root) CLI_RESULTS_ROOT="$2"; shift 2 ;;
    --text-cutoff-date) TEXT_CUTOFF_DATE="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) usage >&2; die "Unknown argument: $1" ;;
  esac
done
case "$MODALITY" in sequence|text|structure|ppi) ;; *) usage >&2; die "--modality is required" ;; esac
[[ -d "$EMBEDDING_STATE_ROOT" ]] || die "Embedding state does not exist: $EMBEDDING_STATE_ROOT"

JOB_TOKEN="${JOB_ID:-manual_$$}"
RUN_TAG="${JOB_TOKEN}_$(date +%Y%m%d_%H%M%S)_${MODALITY}"
WORK="/scratch0/cafa3_embedding_retry_${JOB_TOKEN}"
RESULTS_ROOT="${CLI_RESULTS_ROOT:-${RESULTS_ROOT:-$HOME/cafa3_embedding_retry_results}}"
FINAL_RUN_ROOT="$RESULTS_ROOT/$RUN_TAG"
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

git_in_dir() { local directory="$1"; shift; (cd "$directory" && git "$@"); }

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
  [[ ! -d "$SCRATCH_RESULT_ROOT" ]] || cp -a "$SCRATCH_RESULT_ROOT/." "$staging/" || copy_status=$?
  [[ ! -f "$WORKFLOW_LOG" ]] || cp -p "$WORKFLOW_LOG" "$staging/logs/workflow.log" || copy_status=$?
  if [[ "$workflow_status" == "0" ]]; then
    [[ -f "$staging/RETRY_COMPLETE.json" ]] || copy_status=1
    if [[ ! -f "$staging/EMBEDDING_GATE_PASSED.json" && \
          ! -f "$staging/GENERATION_INCOMPLETE.json" && \
          ! -f "$staging/reports/embedding_state/EMBEDDING_GATE_PASSED.json" && \
          ! -f "$staging/reports/embedding_state/GENERATION_INCOMPLETE.json" ]]; then
      copy_status=1
    fi
  else
    rm -f "$staging/RETRY_COMPLETE.json" "$staging/EMBEDDING_GATE_PASSED.json"
    printf '{"complete":false,"workflow_exit_status":%s}\n' "$workflow_status" \
      > "$staging/WORKFLOW_FAILED.json" || copy_status=$?
  fi
  if [[ "$copy_status" == "0" ]]; then
    mv "$staging" "$destination" || copy_status=$?
  fi
  if [[ "$copy_status" == "0" ]]; then
    RESULTS_COPIED=1
    echo "==> Published compact retry report: $destination"
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
  copy_results "$status" || copy_status=$?
  if [[ "$WORK_OWNED" == "1" ]]; then
    if [[ "$WORK" == /scratch0/cafa3_embedding_retry_* && ! -L "$WORK" ]]; then
      cd "$HOME"
      rm -rf "$WORK"
    else
      echo "Refusing unsafe scratch cleanup path: $WORK" >&2
      [[ "$status" != "0" ]] || status=1
    fi
  fi
  if [[ "$status" == "0" && "$copy_status" != "0" ]]; then status="$copy_status"; fi
  exit "$status"
}
trap cleanup EXIT
trap 'echo "Received termination signal"; exit 130' INT TERM

[[ ! -e "$WORK" ]] || die "Scratch path already exists: $WORK"
mkdir -p "$WORK/tmp" "$RESULTS_ROOT"
WORK_OWNED=1
export TMPDIR="$WORK/tmp" TMP="$WORK/tmp" TEMP="$WORK/tmp"

if [[ -z "$FRAMEWORK_COMMIT" ]]; then
  [[ -d "$SUBMISSION_DIR/.git" ]] || die "Submit from a clean framework checkout"
  [[ -z "$(git_in_dir "$SUBMISSION_DIR" status --porcelain)" ]] || die "Submission checkout is dirty"
  FRAMEWORK_COMMIT="$(git_in_dir "$SUBMISSION_DIR" rev-parse HEAD)"
fi
[[ "$FRAMEWORK_COMMIT" =~ ^[0-9a-fA-F]{40}$ ]] || die "FRAMEWORK_COMMIT must be complete"
[[ "$PFP_COMMIT" =~ ^[0-9a-fA-F]{40}$ ]] || die "PFP_COMMIT must be complete"

echo "Host              : $(hostname)"
echo "Job ID            : ${JOB_ID:-manual}"
echo "Modality          : $MODALITY"
echo "Persistent state  : $EMBEDDING_STATE_ROOT"
echo "Scratch           : $WORK"
echo "Final report      : $FINAL_RUN_ROOT"

git clone --no-checkout "$FRAMEWORK_REPO_URL" "$FRAMEWORK_DIR"
git_in_dir "$FRAMEWORK_DIR" checkout --detach "$FRAMEWORK_COMMIT"
git clone --no-checkout "$PFP_REPO_URL" "$PFP_DIR"
git_in_dir "$PFP_DIR" checkout --detach "$PFP_COMMIT"

cd "$FRAMEWORK_DIR"
source scripts/reproduction_common.sh
load_framework_paths "$FRAMEWORK_DIR"
activate_or_create_mmfp_env
PYTHON_BIN="$(command -v python)"

command=(
  bash "$FRAMEWORK_DIR/scripts/reproduction/run_cafa3_embedding_retry.sh"
  --pfp-root "$PFP_DIR"
  --work-dir "$WORKFLOW_WORK_DIR"
  --output-dir "$SCRATCH_RESULT_ROOT"
  --embedding-state-root "$EMBEDDING_STATE_ROOT"
  --modality "$MODALITY"
  --text-cutoff-date "$TEXT_CUTOFF_DATE"
)
printf 'Command:'; printf ' %q' "${command[@]}"; printf '\n'
set +e
PYTHON_BIN="$PYTHON_BIN" "${command[@]}" 2>&1 | tee "$WORKFLOW_LOG"
WORKFLOW_STATUS=${PIPESTATUS[0]}
set -e
[[ "$WORKFLOW_STATUS" == "0" ]] || exit "$WORKFLOW_STATUS"
copy_results 0
echo "Finished: $(date)"
