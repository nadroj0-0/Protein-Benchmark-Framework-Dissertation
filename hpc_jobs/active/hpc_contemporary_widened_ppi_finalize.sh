#!/usr/bin/env bash
# Add the validated widened PPI delta to the corrected contemporary cache.

#$ -l tmem=16G
#$ -l tscratch=40G
#$ -l scratch0free=60G
#$ -l h_rt=48:0:0
#$ -j y
#$ -N cont_ppi_hydrate
#$ -V
#$ -notify

set -Eeuo pipefail
export PYTHONDONTWRITEBYTECODE=1

die() { echo "ERROR: $*" >&2; exit 2; }
require_value() { [[ $# -ge 2 && -n "$2" ]] || die "$1 requires a value"; }
git_in_dir() { local directory="$1"; shift; (cd "$directory" && git "$@"); }

BASE_ROOT=/SAN/bioinf/bmpfp/embeddings/contemporary/2025_01_to_2026_02_supervisor/variants/text-cutoff-2025-03-08__ppi-paper-faithful/finalized_pfp_cache
DELTA_ROOT=/SAN/bioinf/bmpfp/embeddings/contemporary/2025_01_to_2026_02_supervisor/ppi_deltas/validated-wider-unambiguous/7119396_20260729_025638/delta
BENCHMARK_DIR=/SAN/bioinf/bmpfp/benchmarks/contemporary/2025_01_to_2026_02_supervisor
PLAN_DIR=/SAN/bioinf/bmpfp/embeddings/contemporary/2025_01_to_2026_02_supervisor/reuse_plan
INPUT_ACQUISITION=/SAN/bioinf/bmpfp/embeddings/contemporary/2025_01_to_2026_02_supervisor/reports/input_acquisition.tsv
OBO_FILE=/SAN/bioinf/bmpfp/frozen_inputs/ontology/2025-02-06/go-basic.obo
VARIANT_NAME=text-cutoff-2025-03-08__ppi-widened-unambiguous
VARIANT_ROOT=/SAN/bioinf/bmpfp/embeddings/contemporary/2025_01_to_2026_02_supervisor/variants/text-cutoff-2025-03-08__ppi-widened-unambiguous
EXPECTED_CUTOFF=2025-03-08
EXPECTED_TARGET_COUNT=156421
EXPECTED_BASE_PPI_COUNT=100334
EXPECTED_DELTA_COUNT=26330
EXPECTED_FINAL_PPI_COUNT=126664
MINIMUM_SAN_FREE_GB="${MINIMUM_SAN_FREE_GB:-40}"
CLI_RESULTS_ROOT=""

usage() {
  cat <<'EOF'
Usage: qsub hpc_jobs/active/hpc_contemporary_widened_ppi_finalize.sh [options]

Options may override --base-root, --delta-root, --benchmark-dir, --plan-dir,
--input-acquisition, --obo-file, --variant-name, --variant-root,
--expected-cutoff, and --results-root.

The CPU-only job extracts the corrected paper-faithful cache into job-owned
scratch, adds exactly the validated non-overlapping PPI delta, rebuilds all
evidence, and publishes a separately named variant only after exhaustive and
round-trip validation. The paper-faithful cache is opened read-only and is
never modified.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --base-root) require_value "$@"; BASE_ROOT="$2"; shift 2 ;;
    --delta-root) require_value "$@"; DELTA_ROOT="$2"; shift 2 ;;
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

for path in "$BASE_ROOT" "$DELTA_ROOT" "$BENCHMARK_DIR" "$PLAN_DIR"; do
  [[ -d "$path" ]] || die "Required directory is missing: $path"
done
for path in "$INPUT_ACQUISITION" "$OBO_FILE"; do
  [[ -f "$path" ]] || die "Required file is missing: $path"
done
[[ "$BASE_ROOT" == /SAN/* ]] || die "Base root must be on SAN"
[[ "$DELTA_ROOT" == /SAN/* ]] || die "Delta root must be on SAN"
[[ "$VARIANT_ROOT" == /SAN/* ]] || die "Variant root must be on SAN"
[[ "$VARIANT_ROOT" != "$BASE_ROOT" ]] || die "Variant root must differ from base root"
[[ ! -e "$VARIANT_ROOT" ]] || die "Variant root already exists: $VARIANT_ROOT"
[[ ! -e "${VARIANT_ROOT}.staging" ]] || die "Stale publication staging exists"
[[ "$EXPECTED_CUTOFF" =~ ^[0-9]{4}-[0-9]{2}-[0-9]{2}$ ]] || \
  die "Invalid cutoff: $EXPECTED_CUTOFF"
[[ "$MINIMUM_SAN_FREE_GB" =~ ^[1-9][0-9]*$ ]] || \
  die "MINIMUM_SAN_FREE_GB must be a positive integer"
SAN_FREE_KB="$(df -Pk "$BASE_ROOT" | awk 'NR == 2 {print $4}')"
[[ "$SAN_FREE_KB" =~ ^[0-9]+$ ]] || die "Could not measure SAN free space"
(( SAN_FREE_KB >= MINIMUM_SAN_FREE_GB * 1024 * 1024 )) || \
  die "SAN free space is below the ${MINIMUM_SAN_FREE_GB}G publication floor"

JOB_TOKEN="${JOB_ID:-manual_$$}"
RUN_TAG="${JOB_TOKEN}_$(date +%Y%m%d_%H%M%S)"
WORK="/scratch0/contemporary_widened_ppi_${JOB_TOKEN}"
SCRATCH_VARIANT="$WORK/variant"
SCRATCH_REPORT="$WORK/report"
PUBLICATION_STAGING="${VARIANT_ROOT}.staging-${JOB_TOKEN}"
REPORTS_ROOT="${CLI_RESULTS_ROOT:-${VARIANT_ROOT}_job_reports}"
FINAL_REPORT="$REPORTS_ROOT/$RUN_TAG"
FAILED_REPORT="${FINAL_REPORT}.failed"
FRAMEWORK_REPO_URL="${FRAMEWORK_REPO_URL:-https://github.com/nadroj0-0/Protein-Benchmark-Framework-Dissertation.git}"
PFP_REPO_URL="${PFP_REPO_URL:-https://github.com/psipred/PFP.git}"
FRAMEWORK_COMMIT="${FRAMEWORK_COMMIT:-}"
PFP_COMMIT="${PFP_COMMIT:-1e04fd6d6d3c40458fd41ec1a881ed6e24de768e}"
FRAMEWORK_DIR="$WORK/Protein-Benchmark-Framework-Dissertation"
PFP_DIR="$WORK/PFP"
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
    echo "Published widened-PPI report: $destination"
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
  if [[ -d "$PUBLICATION_STAGING" && ! -L "$PUBLICATION_STAGING" ]]; then
    rm -rf -- "$PUBLICATION_STAGING"
  fi
  if [[ "$WORK_OWNED" == "1" && "$WORK" == /scratch0/contemporary_widened_ppi_* && ! -L "$WORK" ]]; then
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
trap 'echo "Received termination signal"; exit 130' INT TERM USR1 USR2

[[ ! -e "$WORK" ]] || die "Scratch path already exists: $WORK"
[[ ! -e "$PUBLICATION_STAGING" ]] || die "Publication staging already exists"
mkdir -p "$WORK/tmp" "$SCRATCH_REPORT" "$REPORTS_ROOT"
WORK_OWNED=1
export TMPDIR="$WORK/tmp" TMP="$WORK/tmp" TEMP="$WORK/tmp"

if [[ -z "$FRAMEWORK_COMMIT" ]]; then
  [[ -d "$SUBMISSION_DIR/.git" ]] || die "Submit from the framework checkout or pass FRAMEWORK_COMMIT"
  FRAMEWORK_COMMIT="$(git_in_dir "$SUBMISSION_DIR" rev-parse HEAD)"
fi
[[ "$FRAMEWORK_COMMIT" =~ ^[0-9a-fA-F]{40}$ ]] || die "Invalid FRAMEWORK_COMMIT"
[[ "$PFP_COMMIT" =~ ^[0-9a-fA-F]{40}$ ]] || die "Invalid PFP_COMMIT"

BASE_ARCHIVE="$BASE_ROOT/contemporary_embedding_cache.tar.gz"
BASE_MARKER="$BASE_ROOT/FINAL_CACHE_COMPLETE.json"
DELTA_ARCHIVE="$DELTA_ROOT/ppi_delta.tar.gz"
DELTA_MARKER="$DELTA_ROOT/DELTA_COMPLETE.json"
for path in "$BASE_ARCHIVE" "$BASE_MARKER" "$DELTA_ARCHIVE" "$DELTA_MARKER"; do
  [[ -f "$path" ]] || die "Required publication input is missing: $path"
done
BASE_ARCHIVE_SHA_BEFORE="$(sha256sum "$BASE_ARCHIVE" | awk '{print $1}')"
BASE_MARKER_SHA_BEFORE="$(sha256sum "$BASE_MARKER" | awk '{print $1}')"

echo "Base final cache : $BASE_ROOT"
echo "PPI delta        : $DELTA_ROOT"
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
  "$WORK" "$BASE_ROOT" "$DELTA_ROOT" "$BENCHMARK_DIR" "$PLAN_DIR" \
  "$(dirname "$INPUT_ACQUISITION")" "$(dirname "$OBO_FILE")" \
  "$(dirname "$VARIANT_ROOT")"; do
  add_mmfp_singularity_bind "$path"
done
activate_or_create_mmfp_env
PYTHON_BIN="$(command -v python)"
POLICY="$FRAMEWORK_DIR/configs/contemporary_embedding_resume.json"
CONFIG="$FRAMEWORK_DIR/configs/pfp_benchmark_run.temporal.json"
BASELINE_ROOT="$SCRATCH_VARIANT/source_baseline"
STATE_ROOT="$SCRATCH_VARIANT/retry_state"
FINAL_ROOT="$SCRATCH_VARIANT/finalized_pfp_cache"

{
  echo "==> [1/5] Compose the corrected cache with the non-overlapping PPI delta"
  "$PYTHON_BIN" "$FRAMEWORK_DIR/scripts/embeddings/compose_contemporary_ppi_delta.py" \
    --base-final-root "$BASE_ROOT" \
    --delta-root "$DELTA_ROOT" \
    --plan-dir "$PLAN_DIR" \
    --policy "$POLICY" \
    --config "$CONFIG" \
    --input-acquisition "$INPUT_ACQUISITION" \
    --variant-name "$VARIANT_NAME" \
    --work-dir "$WORK/composition_work" \
    --output-root "$BASELINE_ROOT" \
    --expected-target-count "$EXPECTED_TARGET_COUNT" \
    --expected-base-ppi-count "$EXPECTED_BASE_PPI_COUNT" \
    --expected-delta-count "$EXPECTED_DELTA_COUNT" \
    --expected-final-ppi-count "$EXPECTED_FINAL_PPI_COUNT" \
    --report "$SCRATCH_REPORT/composition.json"

  echo "==> [2/5] Initialise fresh evidence bound to the widened baseline"
  PYTHON_BIN="$PYTHON_BIN" FRAMEWORK_COMMIT="$FRAMEWORK_COMMIT" \
    bash "$FRAMEWORK_DIR/scripts/embeddings/initialize_contemporary_embedding_state.sh" \
      --benchmark-dir "$BENCHMARK_DIR" \
      --plan-dir "$PLAN_DIR" \
      --baseline-root "$BASELINE_ROOT" \
      --state-root "$STATE_ROOT" \
      --pfp-root "$PFP_DIR" \
      --output-dir "$SCRATCH_REPORT/initialization" \
      --text-cutoff-date "$EXPECTED_CUTOFF"

  echo "==> [3/5] Exhaustively validate and finalize the new cache"
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

  echo "==> [4/5] Assert the paper-faithful source was not modified"
  [[ "$(sha256sum "$BASE_ARCHIVE" | awk '{print $1}')" == "$BASE_ARCHIVE_SHA_BEFORE" ]] || \
    die "Paper-faithful archive changed during hydration"
  [[ "$(sha256sum "$BASE_MARKER" | awk '{print $1}')" == "$BASE_MARKER_SHA_BEFORE" ]] || \
    die "Paper-faithful completion marker changed during hydration"

  "$PYTHON_BIN" - "$SCRATCH_VARIANT" "$VARIANT_ROOT" "$VARIANT_NAME" \
    "$EXPECTED_CUTOFF" "$FRAMEWORK_COMMIT" "$PFP_COMMIT" \
    "$BASE_ARCHIVE_SHA_BEFORE" "$EXPECTED_FINAL_PPI_COUNT" <<'PY'
import json
import os
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

root = Path(sys.argv[1])
publication_root = Path(sys.argv[2])
composition = json.loads((root / "source_baseline/COMPOSITION_COMPLETE.json").read_text())
final = json.loads((root / "finalized_pfp_cache/FINAL_CACHE_COMPLETE.json").read_text())
expected_ppi = int(sys.argv[8])
if final["accepted_counts"]["ppi"] != expected_ppi:
    raise SystemExit(
        f"Final PPI count differs: {final['accepted_counts']['ppi']} != {expected_ppi}"
    )
payload = {
    "schema_version": 1,
    "complete": True,
    "validated": True,
    "completed_at": datetime.now(timezone.utc).isoformat(),
    "variant": sys.argv[3],
    "text_cutoff": sys.argv[4],
    "ppi_policy": "widened-unambiguous",
    "ppi_delta_semantics": "additive-only; no existing PPI arrays replaced",
    "old_text_carried_forward": False,
    "paper_faithful_source_modified": False,
    "framework_commit": sys.argv[5],
    "pfp_commit": sys.argv[6],
    "base_archive_sha256": sys.argv[7],
    "delta_count": composition["delta_count"],
    "replacement_count": composition["replacement_count"],
    "composition_archive_sha256": composition["combined_archive_sha256"],
    "final_archive_sha256": final["archive_sha256"],
    "accepted_counts": final["accepted_counts"],
    "final_root": str(publication_root / "finalized_pfp_cache"),
}
descriptor, temporary = tempfile.mkstemp(prefix=".VARIANT_COMPLETE.", dir=str(root))
with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
    json.dump(payload, handle, indent=2, sort_keys=True)
    handle.write("\n")
    handle.flush()
    os.fsync(handle.fileno())
os.replace(temporary, root / "VARIANT_COMPLETE.json")
PY

  cp -p "$SCRATCH_VARIANT/VARIANT_COMPLETE.json" "$SCRATCH_REPORT/"

  echo "==> [5/5] Atomically publish the validated variant to SAN"
  [[ ! -e "$PUBLICATION_STAGING" && ! -e "$VARIANT_ROOT" ]] || \
    die "Publication destination appeared during validation"
  mkdir -p "$PUBLICATION_STAGING"
  cp -a "$SCRATCH_VARIANT/." "$PUBLICATION_STAGING/"
  [[ -f "$PUBLICATION_STAGING/VARIANT_COMPLETE.json" ]] || \
    die "Staged publication lacks VARIANT_COMPLETE.json"
  [[ -f "$PUBLICATION_STAGING/finalized_pfp_cache/FINAL_CACHE_COMPLETE.json" ]] || \
    die "Staged publication lacks FINAL_CACHE_COMPLETE.json"
  mv "$PUBLICATION_STAGING" "$VARIANT_ROOT"
} 2>&1 | tee "$WORKFLOW_LOG"

echo "Validated widened-PPI variant: $VARIANT_ROOT/finalized_pfp_cache"
