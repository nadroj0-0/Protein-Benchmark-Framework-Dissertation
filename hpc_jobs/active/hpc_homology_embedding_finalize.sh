#!/usr/bin/env bash
#$ -l tmem=16G
#$ -l tscratch=40G
#$ -l scratch0free=60G
#$ -l h_rt=24:0:0
#$ -pe smp 2
#$ -j y
#$ -N hom30_emb_final
#$ -V

set -Eeuo pipefail

die() { echo "ERROR: $*" >&2; exit 2; }
git_in_dir() { local directory="$1"; shift; (cd "$directory" && git "$@"); }

LEDGER_DIR=""
BATCH_ROOT=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --ledger-dir) LEDGER_DIR="$2"; shift 2 ;;
    --batch-root) BATCH_ROOT="$2"; shift 2 ;;
    *) die "Unknown argument: $1" ;;
  esac
done
[[ -d "$LEDGER_DIR" ]] || die "Missing ledger: $LEDGER_DIR"
[[ -n "$BATCH_ROOT" ]] || die "--batch-root is required"

JOB_TOKEN="${JOB_ID:-manual_$$}"
WORK="/scratch0/homology_embedding_finalize_${JOB_TOKEN}"
FRAMEWORK_REPO_URL="${FRAMEWORK_REPO_URL:-https://github.com/nadroj0-0/Protein-Benchmark-Framework-Dissertation.git}"
FRAMEWORK_COMMIT="${FRAMEWORK_COMMIT:-}"
FRAMEWORK_DIR="$WORK/Protein-Benchmark-Framework-Dissertation"
SCRATCH_RESULT="$WORK/result"
STAGED_LEDGER="$WORK/ledger"
STAGED_DELTAS="$WORK/deltas"
STAGED_SOURCES="$WORK/sources"
FINAL_RESULT="$BATCH_ROOT/final"
FAILED_RESULT="$BATCH_ROOT/failed/finalizer_${JOB_TOKEN}"
SUBMISSION_DIR="${SGE_O_WORKDIR:-$PWD}"
WORK_OWNED=0
RESULTS_COPIED=0

wait_for_file() {
  local path="$1" attempts="${2:-30}" delay="${3:-10}"
  local attempt
  for ((attempt=1; attempt<=attempts; attempt++)); do
    [[ -f "$path" ]] && return 0
    echo "Waiting for SAN file ($attempt/$attempts): $path"
    sleep "$delay"
  done
  return 1
}

stage_ledger_with_retry() {
  local attempts="${1:-12}" delay="${2:-10}"
  local attempt
  for ((attempt=1; attempt<=attempts; attempt++)); do
    rm -rf -- "$STAGED_LEDGER"
    mkdir -p "$STAGED_LEDGER"
    if cp -a "$LEDGER_DIR/." "$STAGED_LEDGER/" \
      && [[ -r "$STAGED_LEDGER/output_manifest.json" ]] \
      && [[ -r "$STAGED_LEDGER/RUN_COMPLETE.json" ]] \
      && [[ -r "$STAGED_LEDGER/resolved_embedding_pairs.tsv.gz" ]]; then
      echo "Staged source-resolved ledger locally on attempt $attempt"
      return 0
    fi
    echo "SAN ledger staging failed ($attempt/$attempts); retrying" >&2
    sleep "$delay"
  done
  return 1
}

stage_file_with_retry() {
  local source="$1" destination="$2" attempts="${3:-12}" delay="${4:-10}"
  local attempt temporary="${destination}.partial"
  mkdir -p "$(dirname "$destination")"
  for ((attempt=1; attempt<=attempts; attempt++)); do
    rm -f -- "$temporary" "$destination"
    if cp -p "$source" "$temporary" && mv "$temporary" "$destination"; then
      echo "Staged SAN file locally on attempt $attempt: $source"
      return 0
    fi
    echo "SAN file staging failed ($attempt/$attempts): $source" >&2
    sleep "$delay"
  done
  rm -f -- "$temporary"
  return 1
}

copy_tree_with_retry() {
  local source="$1" destination="$2" attempts="${3:-6}" delay="${4:-10}"
  local attempt
  for ((attempt=1; attempt<=attempts; attempt++)); do
    rm -rf -- "$destination"
    mkdir -p "$destination"
    if cp -a "$source/." "$destination/"; then
      echo "Published result payload on attempt $attempt"
      return 0
    fi
    echo "Result publication failed ($attempt/$attempts); retrying" >&2
    sleep "$delay"
  done
  return 1
}

copy_result() {
  local status="$1"
  local destination="$FINAL_RESULT"
  [[ "$RESULTS_COPIED" == "0" ]] || return 0
  [[ "$status" == "0" ]] || destination="$FAILED_RESULT"
  local staging="${destination}.staging-${JOB_TOKEN}"
  [[ ! -e "$destination" && ! -e "$staging" ]] || return 1
  if [[ -d "$SCRATCH_RESULT" ]]; then
    copy_tree_with_retry "$SCRATCH_RESULT" "$staging" || return 1
  else
    mkdir -p "$staging"
  fi
  if [[ "$status" == "0" ]]; then
    [[ -f "$staging/WORKFLOW_COMPLETE.json" ]] || return 1
    [[ -f "$staging/homology_30_embedding_cache.tar.gz" ]] || return 1
  else
    rm -f "$staging/WORKFLOW_COMPLETE.json"
    printf '{"complete":false,"exit_status":%s}\n' "$status" \
      > "$staging/WORKFLOW_FAILED.json"
  fi
  mkdir -p "$(dirname "$destination")"
  mv "$staging" "$destination"
  RESULTS_COPIED=1
}

cleanup() {
  local status=$?
  local copy_status=0
  trap - EXIT
  set +e
  copy_result "$status" || copy_status=$?
  if [[ "$WORK_OWNED" == "1" && "$WORK" == /scratch0/homology_embedding_finalize_* && ! -L "$WORK" ]]; then
    cd "$HOME"
    rm -rf "$WORK"
  fi
  if [[ "$status" == "0" && "$copy_status" != "0" ]]; then status="$copy_status"; fi
  exit "$status"
}
trap cleanup EXIT
trap 'echo "Received termination signal"; exit 130' INT TERM

[[ ! -e "$FINAL_RESULT" ]] || die "Final result already exists: $FINAL_RESULT"
[[ ! -e "$WORK" ]] || die "Scratch path exists: $WORK"
mkdir -p "$WORK/tmp" "$SCRATCH_RESULT" "$BATCH_ROOT/failed"
WORK_OWNED=1
export TMPDIR="$WORK/tmp" TMP="$WORK/tmp" TEMP="$WORK/tmp"

if [[ -z "$FRAMEWORK_COMMIT" ]]; then
  [[ -d "$SUBMISSION_DIR/.git" ]] || die "Pass FRAMEWORK_COMMIT outside a framework checkout"
  FRAMEWORK_COMMIT="$(git_in_dir "$SUBMISSION_DIR" rev-parse HEAD)"
fi
[[ "$FRAMEWORK_COMMIT" =~ ^[0-9a-fA-F]{40}$ ]] || die "FRAMEWORK_COMMIT must be a full commit"

git clone --no-checkout "$FRAMEWORK_REPO_URL" "$FRAMEWORK_DIR"
git_in_dir "$FRAMEWORK_DIR" checkout --detach "$FRAMEWORK_COMMIT"
cd "$FRAMEWORK_DIR"

source scripts/reproduction_common.sh
load_framework_paths "$FRAMEWORK_DIR"
activate_or_create_mmfp_env
python_bin="$(command -v python)"
wait_for_file "$LEDGER_DIR/output_manifest.json" || \
  die "Ledger output manifest did not become visible: $LEDGER_DIR/output_manifest.json"
stage_ledger_with_retry || die "Could not stage a complete readable ledger from SAN"
ledger_sha="$($python_bin - "$STAGED_LEDGER/output_manifest.json" <<'PY'
import hashlib
import sys
print(hashlib.sha256(open(sys.argv[1], "rb").read()).hexdigest())
PY
)"

generated_args=()
for modality in sequence text structure ppi; do
  delta="$BATCH_ROOT/deltas/$modality"
  source_marker="$delta/WORKFLOW_COMPLETE.json"
  source_archive="$delta/artifacts/generated_${modality}.tar.gz"
  marker="$STAGED_DELTAS/$modality/WORKFLOW_COMPLETE.json"
  archive="$STAGED_DELTAS/$modality/generated_${modality}.tar.gz"
  stage_file_with_retry "$source_marker" "$marker" || \
    die "Could not stage $modality completion marker"
  stage_file_with_retry "$source_archive" "$archive" || \
    die "Could not stage $modality archive"
  "$python_bin" - "$marker" "$archive" "$ledger_sha" "$modality" <<'PY'
import hashlib
import json
import sys
marker = json.load(open(sys.argv[1], encoding="utf-8"))
observed = hashlib.sha256(open(sys.argv[2], "rb").read()).hexdigest()
if marker.get("complete") is not True or marker.get("modality") != sys.argv[4]:
    raise SystemExit("Invalid modality completion marker")
if marker.get("ledger_output_manifest_sha256") != sys.argv[3]:
    raise SystemExit("Delta is bound to a different pair ledger")
if marker.get("archive_sha256") != observed:
    raise SystemExit("Delta archive hash does not match its completion marker")
PY
  generated_args+=(--generated-archive "$modality=$archive")
done

source_table="$WORK/source_archives.tsv"
"$python_bin" - "$STAGED_LEDGER/summary.json" "$source_table" <<'PY'
import json
import sys
from pathlib import Path

summary = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
lines = []
for source in summary["sources"]:
    archive = source["archive"]
    digest = source["archive_sha256"]
    if "\t" in archive or "\n" in archive:
        raise SystemExit("Unsafe source archive path in ledger")
    lines.append(f"{archive}\t{digest}")
Path(sys.argv[2]).write_text("\n".join(lines) + "\n", encoding="utf-8")
PY

source_args=()
source_index=0
while IFS=$'\t' read -r source_archive expected_sha; do
  [[ -n "$source_archive" && -n "$expected_sha" ]] || die "Invalid source archive row"
  source_index=$((source_index + 1))
  staged_archive="$STAGED_SOURCES/source_${source_index}.tar.gz"
  stage_file_with_retry "$source_archive" "$staged_archive" || \
    die "Could not stage reusable source archive: $source_archive"
  observed_sha="$($python_bin - "$staged_archive" <<'PY'
import hashlib
import sys

digest = hashlib.sha256()
with open(sys.argv[1], "rb") as handle:
    for block in iter(lambda: handle.read(1024 * 1024), b""):
        digest.update(block)
print(digest.hexdigest())
PY
)"
  [[ "$observed_sha" == "$expected_sha" ]] || \
    die "Staged reusable source archive hash mismatch: $source_archive"
  source_args+=(--source-archive-override "$source_archive=$staged_archive")
done < "$source_table"
[[ "$source_index" -gt 0 ]] || die "Ledger defines no reusable source archives"

mkdir -p "$SCRATCH_RESULT/reports"
"$python_bin" scripts/embeddings/assemble_pair_resolved_embedding_cache.py \
  --ledger-dir "$STAGED_LEDGER" \
  "${generated_args[@]}" \
  "${source_args[@]}" \
  --policy configs/homology_embedding_generation.json \
  --output-archive "$SCRATCH_RESULT/homology_30_embedding_cache.tar.gz" \
  --report-dir "$SCRATCH_RESULT/reports/assembly"

cp -p "$SCRATCH_RESULT/reports/assembly/RUN_COMPLETE.json" \
  "$SCRATCH_RESULT/WORKFLOW_COMPLETE.json"
copy_result 0
echo "Final homology 30% embedding cache: $FINAL_RESULT/homology_30_embedding_cache.tar.gz"
