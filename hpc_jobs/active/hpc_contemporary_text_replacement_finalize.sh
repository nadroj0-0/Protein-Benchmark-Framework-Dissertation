#!/usr/bin/env bash
# Replace the legacy contemporary text layer and publish one validated variant.

#$ -l tmem=16G
#$ -l tscratch=40G
#$ -l scratch0free=60G
#$ -l h_rt=48:0:0
#$ -j y
#$ -N cont_text_hydrate
#$ -V
#$ -notify

set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1

die() { echo "ERROR: $*" >&2; exit 2; }
require_value() { [[ $# -ge 2 && -n "$2" ]] || die "$1 requires a value"; }
git_in_dir() { local directory="$1"; shift; (cd "$directory" && git "$@"); }

BASE_ROOT=/SAN/bioinf/bmpfp/embeddings/contemporary/2025_01_to_2026_02_supervisor/finalized_pfp_cache
TEXT_RUN_ROOT=/SAN/bioinf/bmpfp/embeddings/contemporary/2025_01_to_2026_02_supervisor/corrected_text_cutoff_2025_03_08/7113674_20260726_111711_text
BENCHMARK_DIR=/SAN/bioinf/bmpfp/benchmarks/contemporary/2025_01_to_2026_02_supervisor
PLAN_DIR=/SAN/bioinf/bmpfp/embeddings/contemporary/2025_01_to_2026_02_supervisor/reuse_plan
INPUT_ACQUISITION=/SAN/bioinf/bmpfp/embeddings/contemporary/2025_01_to_2026_02_supervisor/reports/input_acquisition.tsv
OBO_FILE=/SAN/bioinf/bmpfp/frozen_inputs/ontology/2025-02-06/go-basic.obo
VARIANT_NAME=text-cutoff-2025-03-08__ppi-paper-faithful
VARIANT_ROOT=/SAN/bioinf/bmpfp/embeddings/contemporary/2025_01_to_2026_02_supervisor/variants/text-cutoff-2025-03-08__ppi-paper-faithful
EXPECTED_CUTOFF=2025-03-08
MINIMUM_SAN_FREE_GB="${MINIMUM_SAN_FREE_GB:-45}"
CLI_RESULTS_ROOT=""

usage() {
  cat <<'EOF'
Usage: qsub hpc_jobs/active/hpc_contemporary_text_replacement_finalize.sh [options]

Options may override --base-root, --text-run-root, --benchmark-dir, --plan-dir,
--input-acquisition, --obo-file, --variant-name, --variant-root,
--expected-cutoff, and --results-root.

The CPU-only job extracts the validated old full cache into scratch, removes
its text directory, installs only the corrected cutoff-bound text archive,
builds a fresh baseline and evidence state, then round-trip validates and
publishes a new final archive. It never edits the old final cache.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --base-root) require_value "$@"; BASE_ROOT="$2"; shift 2 ;;
    --text-run-root) require_value "$@"; TEXT_RUN_ROOT="$2"; shift 2 ;;
    --benchmark-dir) require_value "$@"; BENCHMARK_DIR="$2"; shift 2 ;;
    --plan-dir) require_value "$@"; PLAN_DIR="$2"; shift 2 ;;
    --input-acquisition) require_value "$@"; INPUT_ACQUISITION="$2"; shift 2 ;;
    --obo-file) require_value "$@"; OBO_FILE="$2"; shift 2 ;;
    --variant-name) require_value "$@"; VARIANT_NAME="$2"; shift 2 ;;
    --variant-root) require_value "$@"; VARIANT_ROOT="$2"; shift 2 ;;
    --expected-cutoff) require_value "$@"; EXPECTED_CUTOFF="$2"; shift 2 ;;
    --results-root) require_value "$@"; CLI_RESULTS_ROOT="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) usage >&2; die "Unknown argument: $1" ;;
  esac
done

for path in "$BASE_ROOT" "$TEXT_RUN_ROOT" "$BENCHMARK_DIR" "$PLAN_DIR"; do
  [[ -d "$path" ]] || die "Required directory is missing: $path"
done
for path in "$INPUT_ACQUISITION" "$OBO_FILE"; do
  [[ -f "$path" ]] || die "Required file is missing: $path"
done
[[ "$VARIANT_ROOT" == /SAN/* ]] || die "Variant root must be on SAN"
[[ ! -e "$VARIANT_ROOT" ]] || die "Variant root already exists: $VARIANT_ROOT"
[[ "$EXPECTED_CUTOFF" =~ ^[0-9]{4}-[0-9]{2}-[0-9]{2}$ ]] || \
  die "Invalid cutoff: $EXPECTED_CUTOFF"
[[ "$MINIMUM_SAN_FREE_GB" =~ ^[1-9][0-9]*$ ]] || \
  die "MINIMUM_SAN_FREE_GB must be a positive integer"
SAN_FREE_KB="$(df -Pk "$BASE_ROOT" | awk 'NR == 2 {print $4}')"
[[ "$SAN_FREE_KB" =~ ^[0-9]+$ ]] || die "Could not measure SAN free space"
(( SAN_FREE_KB >= MINIMUM_SAN_FREE_GB * 1024 * 1024 )) || \
  die "SAN free space is below the ${MINIMUM_SAN_FREE_GB}G replacement floor"

JOB_TOKEN="${JOB_ID:-manual_$$}"
RUN_TAG="${JOB_TOKEN}_$(date +%Y%m%d_%H%M%S)"
WORK="/scratch0/contemporary_text_replacement_${JOB_TOKEN}"
REPORTS_ROOT="${CLI_RESULTS_ROOT:-${VARIANT_ROOT}_job_reports}"
FINAL_REPORT="$REPORTS_ROOT/$RUN_TAG"
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
RESULTS_COPIED=0

publish_report() {
  local status="$1" destination="$FINAL_REPORT" staging="${FINAL_REPORT}.staging-${JOB_TOKEN}"
  local copy_status=0
  [[ "$RESULTS_COPIED" == "0" ]] || return 0
  if [[ "$status" != "0" ]]; then
    destination="$FAILED_REPORT"
    staging="${FAILED_REPORT}.staging-${JOB_TOKEN}"
  fi
  [[ ! -e "$destination" && ! -e "$staging" ]] || return 1
  mkdir -p "$staging/logs" || return 1
  [[ ! -d "$SCRATCH_REPORT" ]] || cp -a "$SCRATCH_REPORT/." "$staging/" || copy_status=$?
  [[ ! -f "$WORKFLOW_LOG" ]] || cp -p "$WORKFLOW_LOG" "$staging/logs/workflow.log" || copy_status=$?
  if [[ "$status" == "0" ]]; then
    [[ -f "$staging/VARIANT_COMPLETE.json" ]] || copy_status=1
  else
    printf '{"complete":false,"workflow_exit_status":%s}\n' "$status" \
      > "$staging/WORKFLOW_FAILED.json" || copy_status=$?
  fi
  if [[ "$copy_status" == "0" ]]; then mv "$staging" "$destination" || copy_status=$?; fi
  if [[ "$copy_status" == "0" ]]; then
    RESULTS_COPIED=1
    echo "Published replacement report: $destination"
  elif [[ -d "$staging" && ! -L "$staging" ]]; then
    rm -rf "$staging"
  fi
  return "$copy_status"
}

cleanup() {
  local status=$? publish_status=0
  trap - EXIT
  set +e
  publish_report "$status" || publish_status=$?
  if [[ "$WORK_OWNED" == "1" && "$WORK" == /scratch0/contemporary_text_replacement_* && ! -L "$WORK" ]]; then
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
mkdir -p "$WORK/tmp" "$SCRATCH_REPORT" "$REPORTS_ROOT"
WORK_OWNED=1
export TMPDIR="$WORK/tmp" TMP="$WORK/tmp" TEMP="$WORK/tmp"

if [[ -z "$FRAMEWORK_COMMIT" ]]; then
  [[ -d "$SUBMISSION_DIR/.git" ]] || die "Submit from the framework checkout or pass FRAMEWORK_COMMIT"
  [[ -z "$(git_in_dir "$SUBMISSION_DIR" status --porcelain)" ]] || die "Submission checkout is dirty"
  FRAMEWORK_COMMIT="$(git_in_dir "$SUBMISSION_DIR" rev-parse HEAD)"
fi
[[ "$FRAMEWORK_COMMIT" =~ ^[0-9a-fA-F]{40}$ ]] || die "Invalid FRAMEWORK_COMMIT"
[[ "$PFP_COMMIT" =~ ^[0-9a-fA-F]{40}$ ]] || die "Invalid PFP_COMMIT"

echo "Base final cache : $BASE_ROOT"
echo "Corrected text   : $TEXT_RUN_ROOT"
echo "Variant          : $VARIANT_NAME"
echo "Variant root     : $VARIANT_ROOT"
echo "Scratch          : $WORK"

git clone --no-checkout "$FRAMEWORK_REPO_URL" "$FRAMEWORK_DIR"
git_in_dir "$FRAMEWORK_DIR" checkout --detach "$FRAMEWORK_COMMIT"
git clone --no-checkout "$PFP_REPO_URL" "$PFP_DIR"
git_in_dir "$PFP_DIR" checkout --detach "$PFP_COMMIT"

cd "$FRAMEWORK_DIR"
source scripts/reproduction_common.sh
load_framework_paths "$FRAMEWORK_DIR"
for path in \
  "$WORK" "$BASE_ROOT" "$TEXT_RUN_ROOT" "$BENCHMARK_DIR" "$PLAN_DIR" \
  "$(dirname "$INPUT_ACQUISITION")" "$(dirname "$OBO_FILE")" \
  "$(dirname "$VARIANT_ROOT")"; do
  add_mmfp_singularity_bind "$path"
done
activate_or_create_mmfp_env
PYTHON_BIN="$(command -v python)"
POLICY="$FRAMEWORK_DIR/configs/contemporary_embedding_resume.json"
CONFIG="$FRAMEWORK_DIR/configs/pfp_benchmark_run.temporal.json"
BASELINE_ROOT="$VARIANT_ROOT/source_baseline"
STATE_ROOT="$VARIANT_ROOT/retry_state"
FINAL_ROOT="$VARIANT_ROOT/finalized_pfp_cache"

echo "==> [1/3] Compose corrected full baseline with hard text replacement"
"$PYTHON_BIN" "$FRAMEWORK_DIR/scripts/embeddings/compose_contemporary_text_replacement.py" \
  --base-final-root "$BASE_ROOT" \
  --replacement-run-root "$TEXT_RUN_ROOT" \
  --plan-dir "$PLAN_DIR" \
  --policy "$POLICY" \
  --config "$CONFIG" \
  --input-acquisition "$INPUT_ACQUISITION" \
  --expected-cutoff "$EXPECTED_CUTOFF" \
  --variant-name "$VARIANT_NAME" \
  --work-dir "$WORK/composition_work" \
  --output-root "$BASELINE_ROOT" \
  --report "$SCRATCH_REPORT/composition.json"

echo "==> [2/3] Initialise fresh evidence bound to the corrected baseline"
PYTHON_BIN="$PYTHON_BIN" FRAMEWORK_COMMIT="$FRAMEWORK_COMMIT" \
  bash "$FRAMEWORK_DIR/scripts/embeddings/initialize_contemporary_embedding_state.sh" \
    --benchmark-dir "$BENCHMARK_DIR" \
    --plan-dir "$PLAN_DIR" \
    --baseline-root "$BASELINE_ROOT" \
    --state-root "$STATE_ROOT" \
    --pfp-root "$PFP_DIR" \
    --output-dir "$SCRATCH_REPORT/initialization" \
    --text-cutoff-date "$EXPECTED_CUTOFF"

echo "==> [3/3] Exhaustively validate, round-trip, and finalize the variant"
"$PYTHON_BIN" "$FRAMEWORK_DIR/scripts/embeddings/finalize_embedding_state.py" \
  --state-root "$STATE_ROOT" \
  --benchmark-dir "$BENCHMARK_DIR" \
  --obo-file "$OBO_FILE" \
  --pfp-root "$PFP_DIR" \
  --config "$CONFIG" \
  --work-dir "$WORK/finalization_work" \
  --final-root "$FINAL_ROOT" \
  --report-dir "$SCRATCH_REPORT/finalization" \
  --confirm-retries-finished \
  --retire-source-embeddings

"$PYTHON_BIN" - "$VARIANT_ROOT" "$VARIANT_NAME" "$EXPECTED_CUTOFF" \
  "$FRAMEWORK_COMMIT" "$PFP_COMMIT" <<'PY'
import json
import os
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

root = Path(sys.argv[1])
composition = json.loads((root / "source_baseline/COMPOSITION_COMPLETE.json").read_text())
final = json.loads((root / "finalized_pfp_cache/FINAL_CACHE_COMPLETE.json").read_text())
payload = {
    "schema_version": 1,
    "complete": True,
    "completed_at": datetime.now(timezone.utc).isoformat(),
    "variant": sys.argv[2],
    "text_cutoff": sys.argv[3],
    "ppi_policy": "paper-faithful",
    "old_text_carried_forward": False,
    "framework_commit": sys.argv[4],
    "pfp_commit": sys.argv[5],
    "composition_archive_sha256": composition["combined_archive_sha256"],
    "final_archive_sha256": final["archive_sha256"],
    "accepted_counts": final["accepted_counts"],
    "final_root": str((root / "finalized_pfp_cache").resolve()),
}
descriptor, temporary = tempfile.mkstemp(prefix=".VARIANT_COMPLETE.", dir=str(root))
with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
    json.dump(payload, handle, indent=2, sort_keys=True)
    handle.write("\n")
    handle.flush()
    os.fsync(handle.fileno())
os.replace(temporary, root / "VARIANT_COMPLETE.json")
PY

cp -p "$VARIANT_ROOT/VARIANT_COMPLETE.json" "$SCRATCH_REPORT/"
echo "Validated corrected variant: $FINAL_ROOT"
