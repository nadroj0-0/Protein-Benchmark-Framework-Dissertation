#!/usr/bin/env bash
# Shared Grid Engine body for one homology embedding modality.

set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
REPO_HINT="$(cd "${HERE}/../.." && pwd)"

die() { echo "ERROR: $*" >&2; exit 2; }
git_in_dir() { local directory="$1"; shift; (cd "$directory" && git "$@"); }

MODALITY="${HOMOLOGY_EMBEDDING_MODALITY:-}"
BENCHMARK_DIR=""
LEDGER_DIR=""
BATCH_ROOT=""
CLI_ARTIFACT_CATALOG="${ARTIFACT_CATALOG:-}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --modality) MODALITY="$2"; shift 2 ;;
    --benchmark-dir) BENCHMARK_DIR="$2"; shift 2 ;;
    --ledger-dir) LEDGER_DIR="$2"; shift 2 ;;
    --batch-root) BATCH_ROOT="$2"; shift 2 ;;
    --artifact-catalog) CLI_ARTIFACT_CATALOG="$2"; shift 2 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

case "$MODALITY" in sequence|text|structure|ppi) ;; *) die "Invalid modality: $MODALITY" ;; esac
[[ -d "$BENCHMARK_DIR" ]] || die "Missing benchmark: $BENCHMARK_DIR"
[[ -d "$LEDGER_DIR" ]] || die "Missing ledger: $LEDGER_DIR"
[[ -n "$BATCH_ROOT" ]] || die "--batch-root is required"

JOB_TOKEN="${JOB_ID:-manual_$$}"
TASK_TOKEN="${SGE_TASK_ID:-1}"
WORK="/scratch0/homology_embedding_${MODALITY}_${JOB_TOKEN}_${TASK_TOKEN}"
FRAMEWORK_REPO_URL="${FRAMEWORK_REPO_URL:-https://github.com/nadroj0-0/Protein-Benchmark-Framework-Dissertation.git}"
FRAMEWORK_COMMIT="${FRAMEWORK_COMMIT:-}"
FRAMEWORK_DIR="$WORK/Protein-Benchmark-Framework-Dissertation"
PFP_REPO_URL="${PFP_REPO_URL:-https://github.com/psipred/PFP.git}"
PFP_COMMIT="${PFP_COMMIT:-1e04fd6d6d3c40458fd41ec1a881ed6e24de768e}"
PFP_DIR="$WORK/PFP"
SCRATCH_OUTPUT="$WORK/result"
WORKFLOW_WORK="$WORK/workflow"
WORKFLOW_LOG="$WORK/workflow.log"
INPUT_STAGE="$WORK/input_stage"
BENCHMARK_STAGE="$INPUT_STAGE/benchmark"
LEDGER_STAGE="$INPUT_STAGE/ledger"
FINAL_OUTPUT="$BATCH_ROOT/deltas/$MODALITY"
FAILED_OUTPUT="$BATCH_ROOT/failed/${MODALITY}_${JOB_TOKEN}_${TASK_TOKEN}"
SUBMISSION_DIR="${SGE_O_WORKDIR:-$PWD}"
WORK_OWNED=0
RESULTS_COPIED=0

copy_result() {
  local status="$1"
  local destination="$FINAL_OUTPUT"
  [[ "$RESULTS_COPIED" == "0" ]] || return 0
  if [[ "$status" != "0" ]]; then
    destination="$FAILED_OUTPUT"
  fi
  local staging="${destination}.staging-${JOB_TOKEN}-${TASK_TOKEN}"
  [[ ! -e "$destination" && ! -e "$staging" ]] || return 1
  mkdir -p "$staging/logs"
  if [[ -d "$SCRATCH_OUTPUT" ]]; then
    cp -a "$SCRATCH_OUTPUT/." "$staging/"
  fi
  if [[ -f "$WORKFLOW_LOG" ]]; then
    cp -p "$WORKFLOW_LOG" "$staging/logs/workflow.log"
  fi
  if [[ "$status" == "0" ]]; then
    [[ -f "$staging/WORKFLOW_COMPLETE.json" ]] || return 1
    [[ -f "$staging/artifacts/generated_${MODALITY}.tar.gz" ]] || return 1
  else
    rm -f "$staging/WORKFLOW_COMPLETE.json"
    printf '{"complete":false,"exit_status":%s,"modality":"%s"}\n' \
      "$status" "$MODALITY" > "$staging/WORKFLOW_FAILED.json"
  fi
  mkdir -p "$(dirname "$destination")"
  mv "$staging" "$destination"
  RESULTS_COPIED=1
  echo "==> Published result directory: $destination"
}

cleanup() {
  local status=$?
  local copy_status=0
  trap - EXIT
  set +e
  copy_result "$status" || copy_status=$?
  if [[ "$WORK_OWNED" == "1" && "$WORK" == /scratch0/homology_embedding_* && ! -L "$WORK" ]]; then
    cd "$HOME"
    rm -rf "$WORK"
  fi
  if [[ "$status" == "0" && "$copy_status" != "0" ]]; then
    status="$copy_status"
  fi
  exit "$status"
}
trap cleanup EXIT
trap 'echo "Received termination signal"; exit 130' INT TERM

[[ ! -e "$FINAL_OUTPUT" ]] || die "Completed modality output already exists: $FINAL_OUTPUT"
[[ ! -e "$WORK" ]] || die "Scratch path exists: $WORK"
mkdir -p "$WORK/tmp" "$BATCH_ROOT/deltas" "$BATCH_ROOT/failed"
WORK_OWNED=1
export TMPDIR="$WORK/tmp" TMP="$WORK/tmp" TEMP="$WORK/tmp"

echo "==> Staging immutable benchmark CSVs and source-resolved ledger into scratch"
mkdir -p "$BENCHMARK_STAGE" "$LEDGER_STAGE"
for name in \
  bp-training.csv bp-validation.csv bp-test.csv \
  cc-training.csv cc-validation.csv cc-test.csv \
  mf-training.csv mf-validation.csv mf-test.csv; do
  [[ -f "$BENCHMARK_DIR/$name" ]] || die "Benchmark is missing $name"
  cp -p "$BENCHMARK_DIR/$name" "$BENCHMARK_STAGE/$name"
done
cp -a "$LEDGER_DIR/." "$LEDGER_STAGE/"
[[ -f "$LEDGER_STAGE/output_manifest.json" ]] || die "Staged ledger lacks output_manifest.json"
[[ -f "$LEDGER_STAGE/RUN_COMPLETE.json" ]] || die "Staged ledger lacks RUN_COMPLETE.json"

if [[ -z "$FRAMEWORK_COMMIT" ]]; then
  [[ -d "$SUBMISSION_DIR/.git" ]] || die "Pass FRAMEWORK_COMMIT outside a framework checkout"
  FRAMEWORK_COMMIT="$(git_in_dir "$SUBMISSION_DIR" rev-parse HEAD)"
fi
[[ "$FRAMEWORK_COMMIT" =~ ^[0-9a-fA-F]{40}$ ]] || die "FRAMEWORK_COMMIT must be a full commit"
[[ "$PFP_COMMIT" =~ ^[0-9a-fA-F]{40}$ ]] || die "PFP_COMMIT must be a full commit"

echo "Host              : $(hostname)"
echo "Job/task          : ${JOB_ID:-manual}/${SGE_TASK_ID:-1}"
echo "Modality          : $MODALITY"
echo "Benchmark         : $BENCHMARK_DIR"
echo "Ledger            : $LEDGER_DIR"
echo "Final output      : $FINAL_OUTPUT"
echo "Framework commit  : $FRAMEWORK_COMMIT"
echo "PFP commit        : $PFP_COMMIT"
echo "Started           : $(date)"

git clone --no-checkout "$FRAMEWORK_REPO_URL" "$FRAMEWORK_DIR"
git_in_dir "$FRAMEWORK_DIR" checkout --detach "$FRAMEWORK_COMMIT"
git clone --no-checkout "$PFP_REPO_URL" "$PFP_DIR"
git_in_dir "$PFP_DIR" checkout --detach "$PFP_COMMIT"

cd "$FRAMEWORK_DIR"
source scripts/reproduction_common.sh
export ARTIFACT_CATALOG="$CLI_ARTIFACT_CATALOG"
load_framework_paths "$FRAMEWORK_DIR"
artifact_catalog_bind_parent string_embeddings "${STRING_H5_FILE:-}"
activate_or_create_mmfp_env
PYTHON_BIN="$(command -v python)"

COMMAND=(
  bash "$FRAMEWORK_DIR/scripts/embeddings/run_homology_embedding_modality.sh"
  --pfp-root "$PFP_DIR"
  --work-dir "$WORKFLOW_WORK"
  --output-dir "$SCRATCH_OUTPUT"
  --benchmark-dir "$BENCHMARK_STAGE"
  --ledger-dir "$LEDGER_STAGE"
  --modality "$MODALITY"
)
if [[ -n "$ARTIFACT_CATALOG" ]]; then
  COMMAND+=(--artifact-catalog "$ARTIFACT_CATALOG")
fi

set +e
PYTHON_BIN="$PYTHON_BIN" PFP_COMMIT="$PFP_COMMIT" FRAMEWORK_COMMIT="$FRAMEWORK_COMMIT" \
  "${COMMAND[@]}" 2>&1 | tee "$WORKFLOW_LOG"
status=${PIPESTATUS[0]}
set -e
[[ "$status" == "0" ]] || exit "$status"

copy_result 0
echo "Finished successfully: $(date)"
