#!/usr/bin/env bash
# Publish only newly available PPI arrays under the validated wider alias policy.

#$ -l tmem=30G
#$ -l tscratch=15G
#$ -l scratch0free=30G
#$ -l h_rt=24:0:0
#$ -j y
#$ -N cont_ppi_wide
#$ -V
#$ -notify

set -Eeuo pipefail
export PYTHONDONTWRITEBYTECODE=1

die() { echo "ERROR: $*" >&2; exit 2; }
require_value() { [[ $# -ge 2 && -n "$2" ]] || die "$1 requires a value"; }
git_in_dir() { local directory="$1"; shift; (cd "$directory" && git "$@"); }

AUDIT_ROOT=/SAN/bioinf/bmpfp/diagnostics/string_alias_policy_coverage/7114128_20260726T192824Z
BASE_PAIR_STATUS=/SAN/bioinf/bmpfp/embeddings/contemporary/2025_01_to_2026_02_supervisor/variants/text-cutoff-2025-03-08__ppi-paper-faithful/finalized_pfp_cache/evidence/pair_status.tsv
STRING_H5=/SAN/bioinf/bmpfp/frozen_inputs/string/v12.0/protein.network.embeddings.v12.0.h5
OUTPUT_PARENT=/SAN/bioinf/bmpfp/embeddings/contemporary/2025_01_to_2026_02_supervisor/ppi_deltas/validated-wider-unambiguous
EXPECTED_POLICY_DETAILS_SHA256=5c6e8f180e8a39d056ce76f0f80ce2896dcc14f88e5740c6567f2ea2c9d0c96c
EXPECTED_AUDIT_SUMMARY_SHA256=5d39795bd294ed4d9644b8968962ddd17bf94fa8202fc7e8d763b369d4739d37
EXPECTED_BASE_PAIR_STATUS_SHA256=cdbff3243209d4bae9f62b9df2527aa29bb49682f60e868dd78459f08ad0408f
EXPECTED_STRING_H5_SHA256=a3a5875df30ec4f0568b9f9d6ecc06565659c59befb221d018e819f3ce5add72
CLI_RESULTS_ROOT=""

usage() {
  cat <<'EOF'
Usage: qsub hpc_jobs/active/hpc_contemporary_unambiguous_ppi_delta.sh [options]

Options may override --audit-root, --base-pair-status, --string-h5,
--output-parent, and --results-root.

The CPU-only job consumes the completed alias-policy audit, excludes every PPI
pair already accepted in the paper-faithful cache, extracts exactly the 26,330
new unambiguous STRING v12 vectors, and publishes one compact round-trip
validated PPI-only delta. It never edits or hydrates an existing cache.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --audit-root) require_value "$@"; AUDIT_ROOT="$2"; shift 2 ;;
    --base-pair-status) require_value "$@"; BASE_PAIR_STATUS="$2"; shift 2 ;;
    --string-h5) require_value "$@"; STRING_H5="$2"; shift 2 ;;
    --output-parent) require_value "$@"; OUTPUT_PARENT="$2"; shift 2 ;;
    --results-root) require_value "$@"; CLI_RESULTS_ROOT="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) usage >&2; die "Unknown argument: $1" ;;
  esac
done

POLICY_DETAILS="$AUDIT_ROOT/protein_policy_details.tsv.gz"
AUDIT_SUMMARY="$AUDIT_ROOT/summary.json"
for path in "$POLICY_DETAILS" "$AUDIT_SUMMARY" "$BASE_PAIR_STATUS" "$STRING_H5"; do
  [[ -f "$path" ]] || die "Required input is missing: $path"
done
[[ "$OUTPUT_PARENT" == /SAN/* ]] || die "Output parent must be on SAN"

JOB_TOKEN="${JOB_ID:-manual_$$}"
RUN_TAG="${JOB_TOKEN}_$(date +%Y%m%d_%H%M%S)"
RUN_ROOT="$OUTPUT_PARENT/$RUN_TAG"
DELTA_ROOT="$RUN_ROOT/delta"
REPORT_ROOT="${CLI_RESULTS_ROOT:-$RUN_ROOT/job_report}"
FAILED_REPORT="${REPORT_ROOT}.failed"
WORK="/scratch0/contemporary_unambiguous_ppi_${JOB_TOKEN}"
FRAMEWORK_REPO_URL="${FRAMEWORK_REPO_URL:-https://github.com/nadroj0-0/Protein-Benchmark-Framework-Dissertation.git}"
FRAMEWORK_COMMIT="${FRAMEWORK_COMMIT:-}"
FRAMEWORK_DIR="$WORK/Protein-Benchmark-Framework-Dissertation"
WORKFLOW_LOG="$WORK/workflow.log"
SUBMISSION_DIR="${SGE_O_WORKDIR:-$PWD}"
WORK_OWNED=0
RESULTS_COPIED=0

publish_report() {
  local status="$1" destination="$REPORT_ROOT" staging="${REPORT_ROOT}.staging-${JOB_TOKEN}"
  local copy_status=0
  [[ "$RESULTS_COPIED" == "0" ]] || return 0
  if [[ "$status" != "0" ]]; then
    destination="$FAILED_REPORT"
    staging="${FAILED_REPORT}.staging-${JOB_TOKEN}"
  fi
  [[ ! -e "$destination" && ! -e "$staging" ]] || return 1
  mkdir -p "$staging" || return 1
  [[ ! -f "$WORKFLOW_LOG" ]] || cp -p "$WORKFLOW_LOG" "$staging/workflow.log" || copy_status=$?
  [[ ! -f "$DELTA_ROOT/DELTA_COMPLETE.json" ]] || \
    cp -p "$DELTA_ROOT/DELTA_COMPLETE.json" "$staging/" || copy_status=$?
  if [[ "$status" == "0" && ! -f "$staging/DELTA_COMPLETE.json" ]]; then copy_status=1; fi
  if [[ "$status" != "0" ]]; then
    printf '{"complete":false,"workflow_exit_status":%s}\n' "$status" \
      > "$staging/WORKFLOW_FAILED.json" || copy_status=$?
  fi
  if [[ "$copy_status" == "0" ]]; then mv "$staging" "$destination" || copy_status=$?; fi
  if [[ "$copy_status" == "0" ]]; then
    RESULTS_COPIED=1
    echo "Published PPI delta job report: $destination"
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
  if [[ "$WORK_OWNED" == "1" && "$WORK" == /scratch0/contemporary_unambiguous_ppi_* && ! -L "$WORK" ]]; then
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
[[ ! -e "$RUN_ROOT" ]] || die "Run root already exists: $RUN_ROOT"
mkdir -p "$WORK" "$OUTPUT_PARENT"
WORK_OWNED=1

if [[ -z "$FRAMEWORK_COMMIT" ]]; then
  [[ -d "$SUBMISSION_DIR/.git" ]] || die "Submit from the framework checkout or pass FRAMEWORK_COMMIT"
  FRAMEWORK_COMMIT="$(git_in_dir "$SUBMISSION_DIR" rev-parse HEAD)"
fi
[[ "$FRAMEWORK_COMMIT" =~ ^[0-9a-fA-F]{40}$ ]] || die "Invalid FRAMEWORK_COMMIT"

git clone --no-checkout "$FRAMEWORK_REPO_URL" "$FRAMEWORK_DIR"
git_in_dir "$FRAMEWORK_DIR" checkout --detach "$FRAMEWORK_COMMIT"

cd "$FRAMEWORK_DIR"
source scripts/reproduction_common.sh
load_framework_paths "$FRAMEWORK_DIR"
for path in "$WORK" "$AUDIT_ROOT" "$(dirname "$BASE_PAIR_STATUS")" \
  "$(dirname "$STRING_H5")" "$OUTPUT_PARENT"; do
  add_mmfp_singularity_bind "$path"
done

# This CPU-only extractor needs only NumPy and h5py. Avoid importing the full
# GPU/PFP stack during startup; Grid Engine terminated the previous attempt in
# that unrelated heavyweight validation before extraction began.
[[ -x "$CONDA_EXE" ]] || die "Missing conda executable: $CONDA_EXE"
[[ -x "$MMFP_ENV_DIR/bin/python" ]] || die "Missing mmfp Python: $MMFP_ENV_DIR/bin/python"
eval "$("$CONDA_EXE" shell.bash hook)"
conda activate "$MMFP_ENV_DIR"
PYTHON_BIN="$MMFP_ENV_DIR/bin/python"
"$PYTHON_BIN" - <<'PY'
import h5py
import numpy

print(f"Validated lightweight extraction runtime: numpy={numpy.__version__} h5py={h5py.__version__}")
PY

{
  echo "Host             : $(hostname)"
  echo "Framework commit : $FRAMEWORK_COMMIT"
  echo "Audit root       : $AUDIT_ROOT"
  echo "Base pair status : $BASE_PAIR_STATUS"
  echo "STRING H5        : $STRING_H5"
  echo "Delta root       : $DELTA_ROOT"
  "$PYTHON_BIN" scripts/embeddings/generate_unambiguous_ppi_delta.py \
    --policy-details "$POLICY_DETAILS" \
    --audit-summary "$AUDIT_SUMMARY" \
    --base-pair-status "$BASE_PAIR_STATUS" \
    --string-h5 "$STRING_H5" \
    --work-dir "$WORK/extraction" \
    --output-root "$DELTA_ROOT" \
    --expected-policy-details-sha256 "$EXPECTED_POLICY_DETAILS_SHA256" \
    --expected-audit-summary-sha256 "$EXPECTED_AUDIT_SUMMARY_SHA256" \
    --expected-base-pair-status-sha256 "$EXPECTED_BASE_PAIR_STATUS_SHA256" \
    --expected-string-h5-sha256 "$EXPECTED_STRING_H5_SHA256" \
    --expected-target-count 156421 \
    --expected-base-count 100334 \
    --expected-direct-count 126114 \
    --expected-delta-count 26330 \
    --expected-final-count 126664
} 2>&1 | tee "$WORKFLOW_LOG"

"$PYTHON_BIN" - "$RUN_ROOT" "$DELTA_ROOT" "$FRAMEWORK_COMMIT" <<'PY'
import json
import os
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

run_root = Path(sys.argv[1])
delta_root = Path(sys.argv[2])
delta = json.loads((delta_root / "DELTA_COMPLETE.json").read_text())
payload = {
    "schema_name": "contemporary-unambiguous-ppi-delta-run",
    "schema_version": 1,
    "complete": True,
    "completed_at": datetime.now(timezone.utc).isoformat(),
    "framework_commit": sys.argv[3],
    "delta_root": str(delta_root.resolve()),
    "delta_count": delta["delta_count"],
    "final_union_count": delta["final_union_count"],
    "archive_sha256": delta["archive_sha256"],
}
descriptor, temporary = tempfile.mkstemp(prefix=".RUN_COMPLETE.", dir=run_root)
with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
    json.dump(payload, handle, indent=2, sort_keys=True)
    handle.write("\n")
    handle.flush()
    os.fsync(handle.fileno())
os.replace(temporary, run_root / "RUN_COMPLETE.json")
PY

echo "Validated PPI-only delta: $DELTA_ROOT"
