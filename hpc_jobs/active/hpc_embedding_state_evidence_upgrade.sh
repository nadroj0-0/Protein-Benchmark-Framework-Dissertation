#!/usr/bin/env bash
# Upgrade an existing SAN embedding state with per-array evidence hashes.

#$ -l tmem=8G
#$ -l tscratch=10G
#$ -l h_rt=24:0:0
#$ -j y
#$ -N emb_evidence
#$ -V

set -euo pipefail

usage() {
  cat <<'EOF'
Usage: qsub hpc_jobs/active/hpc_embedding_state_evidence_upgrade.sh \
  --state-root PATH --confirm-retries-finished [--results-root PATH]

Do not submit while any embedding retry job for this state is still running.
The job never regenerates embeddings and never changes accepted membership.
Job-owned scratch is removed on success, failure, or termination.
EOF
}

die() { echo "ERROR: $*" >&2; exit 2; }
require_value() { [[ $# -ge 2 && -n "$2" ]] || die "$1 requires a value"; }

STATE_ROOT=""
CLI_RESULTS_ROOT=""
RETRIES_FINISHED=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --state-root) require_value "$@"; STATE_ROOT="$2"; shift 2 ;;
    --results-root) require_value "$@"; CLI_RESULTS_ROOT="$2"; shift 2 ;;
    --confirm-retries-finished) RETRIES_FINISHED=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) usage >&2; die "Unknown argument: $1" ;;
  esac
done

[[ -d "$STATE_ROOT" ]] || die "State root does not exist: $STATE_ROOT"
[[ "$RETRIES_FINISHED" == "1" ]] || \
  die "Refusing submission without --confirm-retries-finished"
[[ "$STATE_ROOT" == /* ]] || die "--state-root must be an absolute path"
if [[ -n "$CLI_RESULTS_ROOT" ]]; then
  [[ "$CLI_RESULTS_ROOT" == /* ]] || die "--results-root must be an absolute path"
fi
STATE_ROOT="$(cd "$STATE_ROOT" && pwd)"

JOB_TOKEN="${JOB_ID:-manual_$$}"
RUN_TAG="${JOB_TOKEN}_$(date +%Y%m%d_%H%M%S)"
WORK="/scratch0/embedding_state_evidence_${JOB_TOKEN}"
RESULTS_ROOT="${CLI_RESULTS_ROOT:-${RESULTS_ROOT:-$HOME/embedding_state_evidence_results}}"
[[ "$RESULTS_ROOT" == /* ]] || die "RESULTS_ROOT must be an absolute path"
FINAL_OUTPUT="$RESULTS_ROOT/$RUN_TAG"
FAILED_OUTPUT="${FINAL_OUTPUT}.failed"
FRAMEWORK_REPO_URL="${FRAMEWORK_REPO_URL:-https://github.com/nadroj0-0/Protein-Benchmark-Framework-Dissertation.git}"
FRAMEWORK_COMMIT="${FRAMEWORK_COMMIT:-}"
FRAMEWORK_DIR="$WORK/Protein-Benchmark-Framework-Dissertation"
SCRATCH_OUTPUT="$WORK/output"
WORKFLOW_LOG="$WORK/workflow.log"
SUBMISSION_DIR="${SGE_O_WORKDIR:-$PWD}"
WORK_OWNED=0

publish() {
  local status="$1" destination="$FINAL_OUTPUT" staging="${FINAL_OUTPUT}.staging-${JOB_TOKEN}"
  [[ "$status" == "0" ]] || {
    destination="$FAILED_OUTPUT"
    staging="${FAILED_OUTPUT}.staging-${JOB_TOKEN}"
  }
  [[ ! -e "$destination" && ! -e "$staging" ]] || return 1
  mkdir -p "$staging/logs" || return 1
  [[ ! -d "$SCRATCH_OUTPUT" ]] || cp -a "$SCRATCH_OUTPUT/." "$staging/" || return 1
  [[ ! -f "$WORKFLOW_LOG" ]] || cp -p "$WORKFLOW_LOG" "$staging/logs/" || return 1
  if [[ "$status" == "0" ]]; then
    [[ -f "$staging/EVIDENCE_UPGRADE_COMPLETE.json" ]] || return 1
  else
    rm -f "$staging/EVIDENCE_UPGRADE_COMPLETE.json"
    printf '{"complete":false,"workflow_exit_status":%s}\n' "$status" \
      > "$staging/WORKFLOW_FAILED.json" || return 1
  fi
  mv "$staging" "$destination"
  echo "Published evidence-upgrade report: $destination"
}

cleanup() {
  local status=$? publish_status=0
  trap - EXIT
  set +e
  publish "$status" || publish_status=$?
  if [[ "$WORK_OWNED" == "1" && "$WORK" == /scratch0/embedding_state_evidence_* && ! -L "$WORK" ]]; then
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
mkdir -p "$WORK/tmp" "$RESULTS_ROOT" "$SCRATCH_OUTPUT"
WORK_OWNED=1
export TMPDIR="$WORK/tmp" TMP="$WORK/tmp" TEMP="$WORK/tmp"

if [[ -z "$FRAMEWORK_COMMIT" ]]; then
  [[ -d "$SUBMISSION_DIR/.git" ]] || die "Submit from a clean framework checkout"
  [[ -z "$(git -C "$SUBMISSION_DIR" status --porcelain)" ]] || \
    die "Submission checkout is dirty"
  FRAMEWORK_COMMIT="$(git -C "$SUBMISSION_DIR" rev-parse HEAD)"
fi
[[ "$FRAMEWORK_COMMIT" =~ ^[0-9a-fA-F]{40}$ ]] || die "Invalid FRAMEWORK_COMMIT"

git clone --no-checkout "$FRAMEWORK_REPO_URL" "$FRAMEWORK_DIR"
git -C "$FRAMEWORK_DIR" checkout --detach "$FRAMEWORK_COMMIT"

cd "$FRAMEWORK_DIR"
source scripts/reproduction_common.sh
load_framework_paths "$FRAMEWORK_DIR"
add_mmfp_singularity_bind "$STATE_ROOT"
activate_or_create_mmfp_env
PYTHON_BIN="$(command -v python)"

command=(
  bash "$FRAMEWORK_DIR/scripts/embeddings/upgrade_embedding_state_evidence.sh"
  --state-root "$STATE_ROOT"
  --output-dir "$SCRATCH_OUTPUT"
  --confirm-retries-finished
)
printf 'Command:'; printf ' %q' "${command[@]}"; printf '\n'
set +e
PYTHON_BIN="$PYTHON_BIN" "${command[@]}" 2>&1 | tee "$WORKFLOW_LOG"
status=${PIPESTATUS[0]}
set -e
exit "$status"
