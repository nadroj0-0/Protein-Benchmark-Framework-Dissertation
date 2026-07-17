#!/usr/bin/env bash
# UCL Grid Engine wrapper for one contemporary embedding-modality retry.

#$ -l tmem=24G
#$ -l tscratch=120G
#$ -l scratch0free=300G
#$ -l h_rt=96:0:0
#$ -l gpu=true
#$ -pe gpu 1
#$ -j y
#$ -N cont_emb_retry
#$ -V

set -euo pipefail

usage() {
  cat <<'EOF'
Usage: qsub hpc_jobs/active/hpc_contemporary_embedding_retry.sh \
  --modality sequence|text|structure|ppi \
  [--benchmark-dir PATH] [--baseline-root PATH] [--plan-dir PATH] \
  [--state-root PATH] [--artifact-catalog PATH] [--results-root PATH]

The wrapper retries only pending pairs for one modality, merges valid arrays
into the one archive-backed SAN state, copies compact reports home, and always
removes job-owned scratch. It never submits another job.
EOF
}

die() { echo "ERROR: $*" >&2; exit 2; }

MODALITY=""
BENCHMARK_DIR=/SAN/bioinf/bmpfp/benchmarks/contemporary/2025_01_to_2026_02_supervisor
BASELINE_ROOT=/SAN/bioinf/bmpfp/embeddings/contemporary/2025_01_to_2026_02_supervisor
PLAN_DIR=""
STATE_ROOT=""
CLI_ARTIFACT_CATALOG="${ARTIFACT_CATALOG:-}"
CLI_RESULTS_ROOT=""
TEXT_CUTOFF_DATE="2025-03-08"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --modality) MODALITY="$2"; shift 2 ;;
    --benchmark-dir) BENCHMARK_DIR="$2"; shift 2 ;;
    --baseline-root) BASELINE_ROOT="$2"; shift 2 ;;
    --plan-dir) PLAN_DIR="$2"; shift 2 ;;
    --state-root) STATE_ROOT="$2"; shift 2 ;;
    --artifact-catalog) CLI_ARTIFACT_CATALOG="$2"; shift 2 ;;
    --results-root) CLI_RESULTS_ROOT="$2"; shift 2 ;;
    --text-cutoff-date) TEXT_CUTOFF_DATE="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) usage >&2; die "Unknown argument: $1" ;;
  esac
done
case "$MODALITY" in sequence|text|structure|ppi) ;; *) usage >&2; die "--modality is required" ;; esac
PLAN_DIR="${PLAN_DIR:-$BASELINE_ROOT/reuse_plan}"
STATE_ROOT="${STATE_ROOT:-$BASELINE_ROOT/retry_state}"
[[ -f "$STATE_ROOT/contract.json" ]] || die "State is not initialized: $STATE_ROOT"

JOB_TOKEN="${JOB_ID:-manual_$$}"
RUN_TAG="${JOB_TOKEN}_$(date +%Y%m%d_%H%M%S)_${MODALITY}"
WORK="/scratch0/contemporary_embedding_retry_${JOB_TOKEN}"
RESULTS_ROOT="${CLI_RESULTS_ROOT:-${RESULTS_ROOT:-$HOME/contemporary_embedding_retry_results}}"
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
  [[ ! -d "$SCRATCH_RESULT_ROOT" ]] || \
    cp -a "$SCRATCH_RESULT_ROOT/." "$staging/" || copy_status=$?
  [[ ! -f "$WORKFLOW_LOG" ]] || \
    cp -p "$WORKFLOW_LOG" "$staging/logs/workflow.log" || copy_status=$?
  if [[ "$workflow_status" == "0" ]]; then
    [[ -f "$staging/RETRY_COMPLETE.json" ]] || copy_status=1
  else
    rm -f "$staging/RETRY_COMPLETE.json" "$staging/EMBEDDING_GATE_PASSED.json"
    printf '{"complete":false,"workflow_exit_status":%s}\n' "$workflow_status" \
      > "$staging/WORKFLOW_FAILED.json" || copy_status=$?
  fi
  if [[ "$copy_status" == "0" ]]; then mv "$staging" "$destination" || copy_status=$?; fi
  if [[ "$copy_status" == "0" ]]; then
    RESULTS_COPIED=1
    echo "Published retry report: $destination"
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
    if [[ "$WORK" == /scratch0/contemporary_embedding_retry_* && ! -L "$WORK" ]]; then
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

[[ -d "$BENCHMARK_DIR" ]] || die "Missing benchmark: $BENCHMARK_DIR"
[[ -d "$PLAN_DIR" ]] || die "Missing reuse plan: $PLAN_DIR"
[[ ! -e "$WORK" ]] || die "Scratch path exists: $WORK"
mkdir -p "$WORK/tmp" "$RESULTS_ROOT"
WORK_OWNED=1
export TMPDIR="$WORK/tmp" TMP="$WORK/tmp" TEMP="$WORK/tmp"

if [[ -z "$FRAMEWORK_COMMIT" ]]; then
  [[ -d "$SUBMISSION_DIR/.git" ]] || die "Submit from a clean framework checkout"
  [[ -z "$(git_in_dir "$SUBMISSION_DIR" status --porcelain)" ]] || \
    die "Submission checkout is dirty"
  FRAMEWORK_COMMIT="$(git_in_dir "$SUBMISSION_DIR" rev-parse HEAD)"
fi
[[ "$FRAMEWORK_COMMIT" =~ ^[0-9a-fA-F]{40}$ ]] || die "Invalid FRAMEWORK_COMMIT"
[[ "$PFP_COMMIT" =~ ^[0-9a-fA-F]{40}$ ]] || die "Invalid PFP_COMMIT"

echo "Host       : $(hostname)"
echo "Modality   : $MODALITY"
echo "State      : $STATE_ROOT"
echo "Scratch    : $WORK"

git clone --no-checkout "$FRAMEWORK_REPO_URL" "$FRAMEWORK_DIR"
git_in_dir "$FRAMEWORK_DIR" checkout --detach "$FRAMEWORK_COMMIT"
git clone --no-checkout "$PFP_REPO_URL" "$PFP_DIR"
git_in_dir "$PFP_DIR" checkout --detach "$PFP_COMMIT"

cd "$FRAMEWORK_DIR"
source scripts/reproduction_common.sh
export ARTIFACT_CATALOG="$CLI_ARTIFACT_CATALOG"
load_framework_paths "$FRAMEWORK_DIR"
# Retry bookkeeping and the baseline archive live on persistent SAN. Expose
# only the caller-selected directories to the immutable MMFP container.
add_mmfp_singularity_bind "$BENCHMARK_DIR"
add_mmfp_singularity_bind "$BASELINE_ROOT"
add_mmfp_singularity_bind "$PLAN_DIR"
add_mmfp_singularity_bind "$(dirname "$STATE_ROOT")"
if [[ "$MODALITY" == "ppi" ]]; then
  artifact_catalog_bind_parent string_embeddings "${STRING_H5_FILE:-}"
fi
activate_or_create_mmfp_env
PYTHON_BIN="$(command -v python)"

command=(
  bash "$FRAMEWORK_DIR/scripts/embeddings/run_contemporary_embedding_retry.sh"
  --pfp-root "$PFP_DIR"
  --work-dir "$WORKFLOW_WORK_DIR"
  --output-dir "$SCRATCH_RESULT_ROOT"
  --benchmark-dir "$BENCHMARK_DIR"
  --plan-dir "$PLAN_DIR"
  --state-root "$STATE_ROOT"
  --modality "$MODALITY"
  --text-cutoff-date "$TEXT_CUTOFF_DATE"
)
if [[ -n "${ARTIFACT_CATALOG:-}" ]]; then
  command+=(--artifact-catalog "$ARTIFACT_CATALOG")
fi
printf 'Command:'; printf ' %q' "${command[@]}"; printf '\n'
set +e
PYTHON_BIN="$PYTHON_BIN" FRAMEWORK_COMMIT="$FRAMEWORK_COMMIT" \
  "${command[@]}" 2>&1 | tee "$WORKFLOW_LOG"
WORKFLOW_STATUS=${PIPESTATUS[0]}
set -e
[[ "$WORKFLOW_STATUS" == "0" ]] || exit "$WORKFLOW_STATUS"
copy_results 0
echo "Finished: $(date)"
