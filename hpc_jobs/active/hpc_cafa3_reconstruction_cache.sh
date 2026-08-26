#!/usr/bin/env bash
# Prepare a CAFA3-reconstruction cache from the authenticated regenerated cache.

#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -l tmem=32G
#$ -l scratch0free=80G
#$ -l tscratch=80G
#$ -l h_rt=24:0:0
#$ -l gpu=true
#$ -pe gpu 1
#$ -N c3_recache
#$ -V
#$ -notify

set -Eeuo pipefail
export PYTHONDONTWRITEBYTECODE=1

SOURCE_ARCHIVE_SHA256=c6cdcafd00b0cb871a50beb8cd649ce5c13d5882fcde9c3663a19a86905b9e87
SOURCE_TARGETS_SHA256=34874b6c3f70c6c6ece20019c57a932ad6027ea1bcde113ad07ff86348cfcec6
PFP_COMMIT=1e04fd6d6d3c40458fd41ec1a881ed6e24de768e

usage() {
  cat <<'EOF'
Usage: qsub hpc_jobs/active/hpc_cafa3_reconstruction_cache.sh \
  --benchmark-id ID --benchmark-dir DIR --benchmark-checksums FILE \
  --source-cache-archive FILE --source-targets FILE \
  --results-root DIR --run-tag NAME

The source cache must be the authenticated regenerated CAFA3 hydrated archive.
ProtT5 arrays are reused only after exact ID+sequence or exact-sequence proof.
Only genuinely absent ProtT5 arrays are generated. Auxiliary modalities are
never transferred across protein IDs.
EOF
}

die() { echo "ERROR: $*" >&2; exit 2; }
require_value() { [[ $# -ge 2 && -n "$2" ]] || die "$1 requires a value"; }
git_in_dir() { local directory="$1"; shift; (cd "$directory" && git "$@"); }
sha256_file() {
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$1" | awk '{print $1}'
  else
    shasum -a 256 "$1" | awk '{print $1}'
  fi
}

BENCHMARK_ID=""
BENCHMARK_DIR=""
BENCHMARK_CHECKSUMS=""
SOURCE_CACHE_ARCHIVE=""
SOURCE_TARGETS=""
RESULTS_ROOT=""
RUN_TAG=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --benchmark-id) require_value "$@"; BENCHMARK_ID="$2"; shift 2 ;;
    --benchmark-dir) require_value "$@"; BENCHMARK_DIR="$2"; shift 2 ;;
    --benchmark-checksums) require_value "$@"; BENCHMARK_CHECKSUMS="$2"; shift 2 ;;
    --source-cache-archive) require_value "$@"; SOURCE_CACHE_ARCHIVE="$2"; shift 2 ;;
    --source-targets) require_value "$@"; SOURCE_TARGETS="$2"; shift 2 ;;
    --results-root) require_value "$@"; RESULTS_ROOT="$2"; shift 2 ;;
    --run-tag) require_value "$@"; RUN_TAG="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) usage >&2; die "Unknown argument: $1" ;;
  esac
done

for value in BENCHMARK_ID BENCHMARK_DIR BENCHMARK_CHECKSUMS SOURCE_CACHE_ARCHIVE SOURCE_TARGETS RESULTS_ROOT RUN_TAG; do
  [[ -n "${!value}" ]] || die "$value is required"
done
[[ "$RUN_TAG" =~ ^[A-Za-z0-9._-]+$ && "$RUN_TAG" =~ [A-Za-z0-9] ]] || die "Unsafe run tag"
[[ -d "$BENCHMARK_DIR" && ! -L "$BENCHMARK_DIR" ]] || die "Benchmark directory is missing or unsafe"
[[ -f "$BENCHMARK_CHECKSUMS" && ! -L "$BENCHMARK_CHECKSUMS" ]] || die "Benchmark checksum ledger is missing or unsafe"
[[ -f "$SOURCE_CACHE_ARCHIVE" && ! -L "$SOURCE_CACHE_ARCHIVE" ]] || die "Source cache archive is missing or unsafe"
[[ -f "$SOURCE_TARGETS" && ! -L "$SOURCE_TARGETS" ]] || die "Source target manifest is missing or unsafe"
[[ "$RESULTS_ROOT" == /* && "$RESULTS_ROOT" != / ]] || die "Results root must be a safe absolute path"
[[ "$(sha256_file "$SOURCE_CACHE_ARCHIVE")" == "$SOURCE_ARCHIVE_SHA256" ]] || die "Source cache archive checksum mismatch"
[[ "$(sha256_file "$SOURCE_TARGETS")" == "$SOURCE_TARGETS_SHA256" ]] || die "Source target manifest checksum mismatch"

JOB_TOKEN="${JOB_ID:-manual_$$}"
WORK="/scratch0/cafa3_reconstruction_cache_${JOB_TOKEN}"
FRAMEWORK_DIR="$WORK/Protein-Benchmark-Framework-Dissertation"
PFP_DIR="$WORK/PFP"
INPUTS="$WORK/inputs"
CACHE_ROOT="$WORK/cache"
ROUNDTRIP_ROOT="$WORK/roundtrip"
SCRATCH_OUTPUT="$WORK/output"
FINAL_OUTPUT="$RESULTS_ROOT/$RUN_TAG"
FAILED_OUTPUT="${FINAL_OUTPUT}.failed"
WORKFLOW_LOG="$WORK/workflow.log"
SUBMISSION_DIR="${SGE_O_WORKDIR:-$PWD}"
FRAMEWORK_REPO_URL="${FRAMEWORK_REPO_URL:-https://github.com/nadroj0-0/Protein-Benchmark-Framework-Dissertation.git}"
PFP_REPO_URL="${PFP_REPO_URL:-https://github.com/psipred/PFP.git}"
FRAMEWORK_COMMIT="${FRAMEWORK_COMMIT:-}"
WORK_OWNED=0
PUBLISHED=0

publish() {
  local status="$1" destination="$FINAL_OUTPUT" staging="${FINAL_OUTPUT}.staging-${JOB_TOKEN}" copy_status=0
  [[ "$PUBLISHED" == 0 ]] || return 0
  if [[ "$status" != 0 ]]; then
    destination="$FAILED_OUTPUT"
    staging="${FAILED_OUTPUT}.staging-${JOB_TOKEN}"
  fi
  [[ ! -e "$destination" && ! -e "$staging" ]] || return 1
  mkdir -p "$staging/logs" || return 1
  [[ ! -d "$SCRATCH_OUTPUT" ]] || cp -a "$SCRATCH_OUTPUT/." "$staging/" || copy_status=$?
  [[ ! -f "$WORKFLOW_LOG" ]] || cp -p "$WORKFLOW_LOG" "$staging/logs/workflow.log" || copy_status=$?
  if [[ "$status" == 0 ]]; then
    [[ -f "$staging/CACHE_PREPARATION_COMPLETE.json" ]] || copy_status=1
    [[ -f "$staging/artifacts/cafa3_reconstruction_embedding_cache.tar.gz" ]] || copy_status=1
  else
    rm -f "$staging/CACHE_PREPARATION_COMPLETE.json"
    printf '{"complete":false,"workflow_exit_status":%s}\n' "$status" > "$staging/WORKFLOW_FAILED.json" || copy_status=$?
  fi
  if [[ "$copy_status" == 0 ]]; then mv "$staging" "$destination" || copy_status=$?; fi
  if [[ "$copy_status" == 0 ]]; then
    PUBLISHED=1
    echo "Published CAFA3 reconstruction cache: $destination"
  elif [[ -d "$staging" && ! -L "$staging" ]]; then
    rm -rf -- "$staging"
  fi
  return "$copy_status"
}

cleanup() {
  local status=$? publish_status=0
  trap - EXIT
  set +e
  publish "$status" || publish_status=$?
  if [[ "$WORK_OWNED" == 1 && "$WORK" == /scratch0/cafa3_reconstruction_cache_* && ! -L "$WORK" ]]; then
    cd "$HOME"
    rm -rf -- "$WORK"
  else
    echo "Refusing unsafe scratch cleanup: $WORK" >&2
    [[ "$status" != 0 ]] || status=1
  fi
  if [[ "$status" == 0 && "$publish_status" != 0 ]]; then status=$publish_status; fi
  exit "$status"
}
trap cleanup EXIT
trap 'exit 130' INT TERM

[[ ! -e "$WORK" ]] || die "Scratch path already exists: $WORK"
[[ ! -e "$FINAL_OUTPUT" && ! -e "$FAILED_OUTPUT" ]] || die "Run tag already exists in results root"
mkdir -p "$WORK/tmp" "$INPUTS/benchmark" "$SCRATCH_OUTPUT/reports" "$SCRATCH_OUTPUT/artifacts" "$RESULTS_ROOT"
WORK_OWNED=1
export TMPDIR="$WORK/tmp" TMP="$WORK/tmp" TEMP="$WORK/tmp"

if [[ -z "$FRAMEWORK_COMMIT" ]]; then
  [[ -d "$SUBMISSION_DIR/.git" ]] || die "Submit from a clean framework checkout"
  [[ -z "$(git_in_dir "$SUBMISSION_DIR" status --porcelain)" ]] || die "Submission checkout is dirty"
  FRAMEWORK_COMMIT="$(git_in_dir "$SUBMISSION_DIR" rev-parse HEAD)"
fi
[[ "$FRAMEWORK_COMMIT" =~ ^[0-9a-f]{40}$ ]] || die "Invalid framework commit"

for name in bp-training.csv bp-validation.csv bp-test.csv cc-training.csv cc-validation.csv cc-test.csv mf-training.csv mf-validation.csv mf-test.csv; do
  cp -p "$BENCHMARK_DIR/$name" "$INPUTS/benchmark/$name"
done
cp -p "$BENCHMARK_CHECKSUMS" "$INPUTS/output_checksums.sha256"
cp -p "$SOURCE_CACHE_ARCHIVE" "$INPUTS/source_cache.tar.gz"
cp -p "$SOURCE_TARGETS" "$INPUTS/source_targets.tsv"

git clone --no-checkout "$FRAMEWORK_REPO_URL" "$FRAMEWORK_DIR"
git_in_dir "$FRAMEWORK_DIR" checkout --detach "$FRAMEWORK_COMMIT"
git clone --no-checkout "$PFP_REPO_URL" "$PFP_DIR"
git_in_dir "$PFP_DIR" checkout --detach "$PFP_COMMIT"
[[ "$(git_in_dir "$PFP_DIR" rev-parse HEAD)" == "$PFP_COMMIT" ]] || die "PFP checkout is not pinned"
[[ -z "$(git_in_dir "$PFP_DIR" status --porcelain)" ]] || die "PFP checkout is dirty"

cd "$FRAMEWORK_DIR"
source scripts/reproduction_common.sh
load_framework_paths "$FRAMEWORK_DIR"
add_mmfp_singularity_bind "$WORK"
add_mmfp_singularity_bind "$RESULTS_ROOT"
activate_or_create_mmfp_env
PYTHON_BIN="$(command -v python)"
CONFIG="$FRAMEWORK_DIR/configs/pfp_benchmark_run.cafa3_reconstructed.json"
HELPER="$FRAMEWORK_DIR/scripts/embeddings/prepare_cafa3_reconstruction_sequence_cache.py"

echo "==> [1/7] Safely extract the authenticated regenerated CAFA3 cache" | tee -a "$WORKFLOW_LOG"
"$PYTHON_BIN" scripts/embeddings/manage_embedding_archive.py extract \
  --archive "$INPUTS/source_cache.tar.gz" \
  --output-cache-root "$CACHE_ROOT" \
  --config "$CONFIG" \
  --report "$SCRATCH_OUTPUT/reports/source_archive_extraction.json" >> "$WORKFLOW_LOG" 2>&1

echo "==> [2/7] Bind exact ID/sequence reuse and identify genuine ProtT5 misses" | tee -a "$WORKFLOW_LOG"
"$PYTHON_BIN" "$HELPER" prepare \
  --benchmark-dir "$INPUTS/benchmark" \
  --benchmark-checksums "$INPUTS/output_checksums.sha256" \
  --source-targets "$INPUTS/source_targets.tsv" \
  --expected-source-targets-sha256 "$SOURCE_TARGETS_SHA256" \
  --cache-root "$CACHE_ROOT" \
  --missing-fasta "$PFP_DIR/data/proteins.fasta" \
  --target-manifest "$SCRATCH_OUTPUT/reports/targets.tsv" \
  --actions-tsv "$SCRATCH_OUTPUT/reports/sequence_actions.tsv" \
  --report "$SCRATCH_OUTPUT/reports/cache_prepare.json" >> "$WORKFLOW_LOG" 2>&1

missing_count="$("$PYTHON_BIN" -c 'import json,sys; print(json.load(open(sys.argv[1]))["actions"].get("generate",0))' "$SCRATCH_OUTPUT/reports/cache_prepare.json")"
if [[ "$missing_count" -gt 0 ]]; then
  echo "==> [3/7] Generate $missing_count genuinely missing ProtT5 arrays" | tee -a "$WORKFLOW_LOG"
  "$PYTHON_BIN" -c 'import torch; assert torch.cuda.is_available(); print(torch.cuda.get_device_name(0))' >> "$WORKFLOW_LOG" 2>&1
  rm -rf -- "$PFP_DIR/data/embedding_cache"
  ln -s "$CACHE_ROOT" "$PFP_DIR/data/embedding_cache"
  (
    cd "$PFP_DIR"
    DEVICE=cuda bash "$FRAMEWORK_DIR/scripts/embeddings/generate_embeddings_sequence.sh"
  ) >> "$WORKFLOW_LOG" 2>&1
else
  echo "==> [3/7] No new ProtT5 arrays are required" | tee -a "$WORKFLOW_LOG"
fi

echo "==> [4/7] Validate every target array and auxiliary-coverage boundary" | tee -a "$WORKFLOW_LOG"
"$PYTHON_BIN" "$HELPER" finalize \
  --benchmark-dir "$INPUTS/benchmark" \
  --benchmark-checksums "$INPUTS/output_checksums.sha256" \
  --cache-root "$CACHE_ROOT" \
  --actions-tsv "$SCRATCH_OUTPUT/reports/sequence_actions.tsv" \
  --report "$SCRATCH_OUTPUT/reports/cache_validation.json" >> "$WORKFLOW_LOG" 2>&1

echo "==> [5/7] Create the deterministic benchmark cache archive" | tee -a "$WORKFLOW_LOG"
"$PYTHON_BIN" scripts/embeddings/manage_embedding_archive.py create \
  --cache-root "$CACHE_ROOT" \
  --archive "$SCRATCH_OUTPUT/artifacts/cafa3_reconstruction_embedding_cache.tar.gz" \
  --config "$CONFIG" \
  --report "$SCRATCH_OUTPUT/reports/archive_creation.json" >> "$WORKFLOW_LOG" 2>&1

echo "==> [6/7] Round-trip and revalidate the exact published archive" | tee -a "$WORKFLOW_LOG"
"$PYTHON_BIN" scripts/embeddings/manage_embedding_archive.py extract \
  --archive "$SCRATCH_OUTPUT/artifacts/cafa3_reconstruction_embedding_cache.tar.gz" \
  --output-cache-root "$ROUNDTRIP_ROOT" \
  --config "$CONFIG" \
  --report "$SCRATCH_OUTPUT/reports/archive_roundtrip.json" >> "$WORKFLOW_LOG" 2>&1
"$PYTHON_BIN" "$HELPER" finalize \
  --benchmark-dir "$INPUTS/benchmark" \
  --benchmark-checksums "$INPUTS/output_checksums.sha256" \
  --cache-root "$ROUNDTRIP_ROOT" \
  --actions-tsv "$SCRATCH_OUTPUT/reports/sequence_actions.tsv" \
  --report "$SCRATCH_OUTPUT/reports/roundtrip_validation.json" >> "$WORKFLOW_LOG" 2>&1

echo "==> [7/7] Bind completion to all accepted reports and archive bytes" | tee -a "$WORKFLOW_LOG"
"$PYTHON_BIN" - \
  "$BENCHMARK_ID" "$FRAMEWORK_COMMIT" "$PFP_COMMIT" \
  "$SCRATCH_OUTPUT/reports/cache_prepare.json" \
  "$SCRATCH_OUTPUT/reports/cache_validation.json" \
  "$SCRATCH_OUTPUT/reports/archive_creation.json" \
  "$SCRATCH_OUTPUT/reports/archive_roundtrip.json" \
  "$SCRATCH_OUTPUT/CACHE_PREPARATION_COMPLETE.json" <<'PY'
import hashlib
import json
import pathlib
import sys
from datetime import datetime, timezone

benchmark_id, framework_commit, pfp_commit, prepare_path, validation_path, creation_path, roundtrip_path, output_path = sys.argv[1:]
load = lambda value: json.loads(pathlib.Path(value).read_text(encoding="utf-8"))
creation = load(creation_path)
roundtrip = load(roundtrip_path)
for key in ("archive_sha256", "member_count", "members_by_directory", "member_content_sha256"):
    if creation[key] != roundtrip[key]:
        raise SystemExit(f"Archive round-trip mismatch: {key}")
payload = {
    "schema_version": 1,
    "complete": True,
    "completed_at": datetime.now(timezone.utc).isoformat(),
    "benchmark_id": benchmark_id,
    "framework_commit": framework_commit,
    "pfp_commit": pfp_commit,
    "source_cache_archive_sha256": "c6cdcafd00b0cb871a50beb8cd649ce5c13d5882fcde9c3663a19a86905b9e87",
    "source_targets_sha256": "34874b6c3f70c6c6ece20019c57a932ad6027ea1bcde113ad07ff86348cfcec6",
    "cache_prepare": load(prepare_path),
    "cache_validation": load(validation_path),
    "archive": creation,
    "roundtrip_validated": True,
}
pathlib.Path(output_path).write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
PY

"$PYTHON_BIN" scripts/model_execution/manage_output_manifest.py write --root "$SCRATCH_OUTPUT" >/dev/null
"$PYTHON_BIN" scripts/model_execution/manage_output_manifest.py verify --root "$SCRATCH_OUTPUT" >/dev/null
rm -rf -- "$CACHE_ROOT" "$ROUNDTRIP_ROOT"
publish 0
echo "Finished: $(date)"
