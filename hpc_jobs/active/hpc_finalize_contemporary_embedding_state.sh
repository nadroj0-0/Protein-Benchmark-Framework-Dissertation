#!/usr/bin/env bash
# Consolidate a contemporary baseline archive plus retry delta into one SAN archive.

#$ -l tmem=8G
#$ -l tscratch=40G
#$ -l scratch0free=50G
#$ -l h_rt=48:0:0
#$ -j y
#$ -N cont_emb_final
#$ -V
#$ -notify

set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1

usage() {
  cat <<'EOF'
Usage: qsub hpc_jobs/active/hpc_finalize_contemporary_embedding_state.sh \
  --state-root PATH --benchmark-dir PATH --obo-file FILE --final-root PATH \
  --confirm-retries-finished --retire-source-embeddings \
  [--results-root PATH] [--config FILE]

The job hydrates the authenticated baseline plus retry delta in job-owned
scratch, validates every accepted array, creates one archive, reads that copied
archive back from SAN and validates it again. Only then does it remove the old
baseline archive and retry-state cache. Compact evidence remains on SAN.
EOF
}

die() { echo "ERROR: $*" >&2; exit 2; }
require_value() { [[ $# -ge 2 && -n "$2" ]] || die "$1 requires a value"; }
git_in_dir() {
  local directory="$1"
  shift
  (cd "$directory" && git "$@")
}

STATE_ROOT=""
BENCHMARK_DIR=""
OBO_FILE=""
FINAL_ROOT=""
CLI_RESULTS_ROOT=""
CLI_CONFIG=""
RETRIES_FINISHED=0
RETIRE_SOURCES=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --state-root) require_value "$@"; STATE_ROOT="$2"; shift 2 ;;
    --benchmark-dir) require_value "$@"; BENCHMARK_DIR="$2"; shift 2 ;;
    --obo-file) require_value "$@"; OBO_FILE="$2"; shift 2 ;;
    --final-root) require_value "$@"; FINAL_ROOT="$2"; shift 2 ;;
    --results-root) require_value "$@"; CLI_RESULTS_ROOT="$2"; shift 2 ;;
    --config) require_value "$@"; CLI_CONFIG="$2"; shift 2 ;;
    --confirm-retries-finished) RETRIES_FINISHED=1; shift ;;
    --retire-source-embeddings) RETIRE_SOURCES=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) usage >&2; die "Unknown argument: $1" ;;
  esac
done

for value in STATE_ROOT BENCHMARK_DIR OBO_FILE FINAL_ROOT; do
  [[ -n "${!value}" ]] || die "$value is required"
  [[ "${!value}" == /* ]] || die "$value must be an absolute path"
done
[[ "$RETRIES_FINISHED" == "1" ]] || die "--confirm-retries-finished is required"
[[ "$RETIRE_SOURCES" == "1" ]] || die "--retire-source-embeddings is required"
[[ -d "$STATE_ROOT" ]] || die "State root is missing: $STATE_ROOT"
[[ -d "$BENCHMARK_DIR" ]] || die "Benchmark directory is missing: $BENCHMARK_DIR"
[[ -f "$OBO_FILE" ]] || die "Ontology is missing: $OBO_FILE"
[[ ! -e "$FINAL_ROOT" ]] || die "Final root already exists: $FINAL_ROOT"
if [[ -n "$CLI_RESULTS_ROOT" ]]; then
  [[ "$CLI_RESULTS_ROOT" == /* ]] || die "--results-root must be absolute"
fi

JOB_TOKEN="${JOB_ID:-manual_$$}"
RUN_TAG="${JOB_TOKEN}_$(date +%Y%m%d_%H%M%S)"
WORK="/scratch0/contemporary_embedding_finalize_${JOB_TOKEN}"
RESULTS_ROOT="${CLI_RESULTS_ROOT:-${RESULTS_ROOT:-$HOME/contemporary_embedding_finalization_results}}"
FINAL_REPORT="$RESULTS_ROOT/$RUN_TAG"
FAILED_REPORT="${FINAL_REPORT}.failed"
FRAMEWORK_REPO_URL="${FRAMEWORK_REPO_URL:-https://github.com/nadroj0-0/Protein-Benchmark-Framework-Dissertation.git}"
PFP_REPO_URL="${PFP_REPO_URL:-https://github.com/psipred/PFP.git}"
FRAMEWORK_COMMIT="${FRAMEWORK_COMMIT:-}"
PFP_COMMIT="${PFP_COMMIT:-1e04fd6d6d3c40458fd41ec1a881ed6e24de768e}"
FRAMEWORK_DIR="$WORK/Protein-Benchmark-Framework-Dissertation"
PFP_DIR="$WORK/PFP"
SCRATCH_REPORT="$WORK/report"
WORKFLOW_LOG="$WORK/workflow.log"
SUBMISSION_DIR="${SGE_O_WORKDIR:-$PWD}"
WORK_OWNED=0

publish_report() {
  local status="$1" destination="$FINAL_REPORT" staging="${FINAL_REPORT}.staging-${JOB_TOKEN}"
  [[ "$status" == "0" ]] || {
    destination="$FAILED_REPORT"
    staging="${FAILED_REPORT}.staging-${JOB_TOKEN}"
  }
  [[ ! -e "$destination" && ! -e "$staging" ]] || return 1
  mkdir -p "$staging/logs" || return 1
  [[ ! -d "$SCRATCH_REPORT" ]] || cp -a "$SCRATCH_REPORT/." "$staging/" || return 1
  [[ ! -f "$WORKFLOW_LOG" ]] || cp -p "$WORKFLOW_LOG" "$staging/logs/workflow.log" || return 1
  if [[ "$status" == "0" ]]; then
    [[ -f "$staging/finalization_report.json" ]] || return 1
  else
    printf '{"complete":false,"workflow_exit_status":%s}\n' "$status" \
      > "$staging/WORKFLOW_FAILED.json" || return 1
  fi
  mv "$staging" "$destination"
  echo "Published compact finalization report: $destination"
}

cleanup() {
  local status=$? publish_status=0
  trap - EXIT
  set +e
  publish_report "$status" || publish_status=$?
  if [[ "$WORK_OWNED" == "1" && "$WORK" == /scratch0/contemporary_embedding_finalize_* && ! -L "$WORK" ]]; then
    cd "$HOME"
    rm -rf -- "$WORK"
  else
    echo "Refusing unsafe scratch cleanup: $WORK" >&2
    [[ "$status" != "0" ]] || status=1
  fi
  if [[ "$status" == "0" && "$publish_status" != "0" ]]; then status="$publish_status"; fi
  exit "$status"
}
trap cleanup EXIT
trap 'echo "Received termination signal"; exit 130' INT TERM

[[ ! -e "$WORK" ]] || die "Scratch path already exists: $WORK"
mkdir -p "$WORK/tmp" "$SCRATCH_REPORT" "$RESULTS_ROOT"
WORK_OWNED=1
export TMPDIR="$WORK/tmp" TMP="$WORK/tmp" TEMP="$WORK/tmp"

if [[ -z "$FRAMEWORK_COMMIT" ]]; then
  [[ -d "$SUBMISSION_DIR/.git" ]] || die "Submit from the framework checkout or pass FRAMEWORK_COMMIT"
  [[ -z "$(git_in_dir "$SUBMISSION_DIR" status --porcelain)" ]] || die "Submission checkout is dirty"
  FRAMEWORK_COMMIT="$(git_in_dir "$SUBMISSION_DIR" rev-parse HEAD)"
fi
[[ "$FRAMEWORK_COMMIT" =~ ^[0-9a-fA-F]{40}$ ]] || die "Invalid FRAMEWORK_COMMIT"
[[ "$PFP_COMMIT" =~ ^[0-9a-fA-F]{40}$ ]] || die "Invalid PFP_COMMIT"

echo "Host          : $(hostname)"
echo "State         : $STATE_ROOT"
echo "Benchmark     : $BENCHMARK_DIR"
echo "Final archive : $FINAL_ROOT/contemporary_embedding_cache.tar.gz"
echo "Scratch       : $WORK"

git clone --no-checkout "$FRAMEWORK_REPO_URL" "$FRAMEWORK_DIR"
git_in_dir "$FRAMEWORK_DIR" checkout --detach "$FRAMEWORK_COMMIT"
git clone --no-checkout "$PFP_REPO_URL" "$PFP_DIR"
git_in_dir "$PFP_DIR" checkout --detach "$PFP_COMMIT"

if [[ -n "$CLI_CONFIG" ]]; then
  [[ -f "$CLI_CONFIG" ]] || die "Config is missing: $CLI_CONFIG"
  cp -p "$CLI_CONFIG" "$WORK/run_config.json"
  EFFECTIVE_CONFIG="$WORK/run_config.json"
else
  EFFECTIVE_CONFIG="$FRAMEWORK_DIR/configs/pfp_benchmark_run.temporal.json"
fi

cd "$FRAMEWORK_DIR"
source scripts/reproduction_common.sh
load_framework_paths "$FRAMEWORK_DIR"
add_mmfp_singularity_bind "$WORK"
add_mmfp_singularity_bind "$STATE_ROOT"
add_mmfp_singularity_bind "$BENCHMARK_DIR"
add_mmfp_singularity_bind "$OBO_FILE"
add_mmfp_singularity_bind "$(dirname "$FINAL_ROOT")"
activate_or_create_mmfp_env
PYTHON_BIN="$(command -v python)"

command=(
  "$PYTHON_BIN" "$FRAMEWORK_DIR/scripts/embeddings/finalize_embedding_state.py"
  --state-root "$STATE_ROOT"
  --benchmark-dir "$BENCHMARK_DIR"
  --obo-file "$OBO_FILE"
  --pfp-root "$PFP_DIR"
  --config "$EFFECTIVE_CONFIG"
  --work-dir "$WORK/finalization_work"
  --final-root "$FINAL_ROOT"
  --report-dir "$SCRATCH_REPORT"
  --confirm-retries-finished
  --retire-source-embeddings
)
printf 'Command:'; printf ' %q' "${command[@]}"; printf '\n'
set +e
"${command[@]}" 2>&1 | tee "$WORKFLOW_LOG"
status=${PIPESTATUS[0]}
set -e
exit "$status"
