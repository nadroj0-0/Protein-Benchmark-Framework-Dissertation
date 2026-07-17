#!/usr/bin/env bash
# Generate the contemporary benchmark's missing embeddings, combine exact
# planner-approved reuse, validate the cache, and publish archive artifacts.

set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
FRAMEWORK_ROOT="$(cd "${HERE}/../.." && pwd)"
# shellcheck source=../reproduction_common.sh
source "$FRAMEWORK_ROOT/scripts/reproduction_common.sh"
PYTHON_BIN="${PYTHON_BIN:-python}"
MMFP_BASE_URL="${MMFP_BASE_URL:-https://zenodo.org/records/19498341/files}"
PREFLIGHT_PER_SPLIT="${PREFLIGHT_PER_SPLIT:-2}"
DISK_POLL_SECONDS="${DISK_POLL_SECONDS:-120}"

TARGET_BENCHMARK_DIR=""
REUSE_PLAN_DIR=""
PFP_ROOT=""
WORK_DIR=""
OUTPUT_DIR=""
TEXT_CUTOFF_DATE=""

usage() {
  cat <<'EOF'
Usage: run_contemporary_embedding_generation.sh \
  --target-benchmark-dir PATH \
  --reuse-plan-dir PATH \
  --pfp-root PATH \
  --work-dir PATH \
  --output-dir PATH \
  --text-cutoff-date YYYY-MM-DD

The PFP checkout must be a disposable, pinned scratch clone. The workflow
never edits tracked upstream PFP source; its one UniProt-only PPI compatibility
change is made in a separate copied script and recorded in the reports. Runtime
data and caches are written only inside the disposable clone.
EOF
}

die() {
  echo "ERROR: $*" >&2
  exit 2
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --target-benchmark-dir) TARGET_BENCHMARK_DIR="$2"; shift 2 ;;
    --reuse-plan-dir) REUSE_PLAN_DIR="$2"; shift 2 ;;
    --pfp-root) PFP_ROOT="$2"; shift 2 ;;
    --work-dir) WORK_DIR="$2"; shift 2 ;;
    --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
    --text-cutoff-date) TEXT_CUTOFF_DATE="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

[[ -d "$TARGET_BENCHMARK_DIR" ]] || die "Missing target benchmark: $TARGET_BENCHMARK_DIR"
[[ -d "$REUSE_PLAN_DIR" ]] || die "Missing reuse plan: $REUSE_PLAN_DIR"
[[ -d "$PFP_ROOT/.git" ]] || die "PFP root is not a Git checkout: $PFP_ROOT"
[[ -n "$WORK_DIR" ]] || die "--work-dir is required"
[[ -n "$OUTPUT_DIR" ]] || die "--output-dir is required"
[[ "$TEXT_CUTOFF_DATE" =~ ^[0-9]{4}-[0-9]{2}-[0-9]{2}$ ]] || \
  die "--text-cutoff-date must use YYYY-MM-DD"
[[ ! -e "$WORK_DIR" ]] || die "Work directory already exists: $WORK_DIR"
[[ ! -e "$OUTPUT_DIR" ]] || die "Output directory already exists: $OUTPUT_DIR"
command -v "$PYTHON_BIN" >/dev/null 2>&1 || die "Python not found: $PYTHON_BIN"

mkdir -p "$WORK_DIR" "$OUTPUT_DIR/logs" "$OUTPUT_DIR/reports" "$OUTPUT_DIR/archives"
WORK_DIR="$(cd "$WORK_DIR" && pwd)"
OUTPUT_DIR="$(cd "$OUTPUT_DIR" && pwd)"
PFP_ROOT="$(cd "$PFP_ROOT" && pwd)"

TARGET_STAGE="$WORK_DIR/target_benchmark"
PLAN_STAGE="$WORK_DIR/reuse_plan"
ARCHIVE_STAGE="$WORK_DIR/published_archives"
PUBLISHED_ROOT="$WORK_DIR/published_cache_root"
RUNTIME_COMPAT="$WORK_DIR/runtime_compat"
FINAL_PACKAGE="$WORK_DIR/final_package"
FINAL_CACHE="$FINAL_PACKAGE/data/embedding_cache"
ACQUISITION_LOG="$OUTPUT_DIR/reports/input_acquisition.tsv"
MODALITY_STATUS="$OUTPUT_DIR/reports/modality_status.tsv"

mkdir -p "$TARGET_STAGE" "$PLAN_STAGE" "$ARCHIVE_STAGE" "$PUBLISHED_ROOT" "$RUNTIME_COMPAT"
printf 'role\tname\tsource\tstaged_path\tsha256\n' > "$ACQUISITION_LOG"
printf 'phase\tmodality\texit_status\n' > "$MODALITY_STATUS"

sha256_file() {
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$1" | awk '{print $1}'
  else
    shasum -a 256 "$1" | awk '{print $1}'
  fi
}

file_count() {
  local directory="$1"
  [[ -d "$directory" ]] || { printf '0\n'; return; }
  find "$directory" -maxdepth 1 -type f -name '*.npy' -print | wc -l | tr -d ' '
}

download_file() {
  local url="$1"
  local destination="$2"
  local partial="${destination}.part"
  echo "Downloading: $url"
  if command -v wget >/dev/null 2>&1; then
    wget --tries=5 --timeout=60 -c "$url" -O "$partial"
  elif command -v curl >/dev/null 2>&1; then
    curl --fail --location --retry 5 --continue-at - --output "$partial" "$url"
  else
    die "Neither wget nor curl is available"
  fi
  mv "$partial" "$destination"
}

archive_sha256() {
  case "$1" in
    mmfp_embeddings_prott5.tar.gz)
      printf '%s\n' '30dd88fc4acbe3bc267bd8d5ae05e4d967fa7c169a6f063f12d2395fb0ffb00f' ;;
    mmfp_embeddings_struct_ppi.tar.gz)
      printf '%s\n' '6d7a243f2c5e2149c162698b4e6a5e297731a4fea835d57b9049bb31f4af32de' ;;
    mmfp_embeddings_text_temporal.tar.gz)
      printf '%s\n' 'df1bf558fab1c018286a5b389665245917ba951f2c6cd8558546d0b1a3b47e36' ;;
    *) return 1 ;;
  esac
}

stop_disk_monitor() {
  if [[ -n "${DISK_MONITOR_PID:-}" ]]; then
    kill "$DISK_MONITOR_PID" 2>/dev/null || true
    wait "$DISK_MONITOR_PID" 2>/dev/null || true
    DISK_MONITOR_PID=""
  fi
}

disk_monitor() {
  while true; do
    {
      printf '%s\t' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
      du -sk "$WORK_DIR" 2>/dev/null | awk '{printf "work_kib=%s", $1}'
      printf '\t'
      df -Pk "$WORK_DIR" | awk 'NR==2 {printf "filesystem_available_kib=%s", $4}'
      printf '\n'
    } >> "$OUTPUT_DIR/reports/scratch_usage.tsv"
    sleep "$DISK_POLL_SECONDS"
  done
}

trap stop_disk_monitor EXIT
trap 'stop_disk_monitor; exit 130' INT TERM
printf 'timestamp\twork_kib\tfilesystem_available_kib\n' > "$OUTPUT_DIR/reports/scratch_usage.tsv"
disk_monitor &
DISK_MONITOR_PID=$!

echo "==> [1/10] Stage and bind the reuse plan to the target CSVs"
for name in \
  bp-training.csv bp-validation.csv bp-test.csv \
  cc-training.csv cc-validation.csv cc-test.csv \
  mf-training.csv mf-validation.csv mf-test.csv; do
  [[ -f "$TARGET_BENCHMARK_DIR/$name" ]] || die "Target benchmark is missing $name"
  cp -p "$TARGET_BENCHMARK_DIR/$name" "$TARGET_STAGE/$name"
  digest="$(sha256_file "$TARGET_STAGE/$name")"
  printf 'target-csv\t%s\t%s\t%s\t%s\n' \
    "$name" "$TARGET_BENCHMARK_DIR/$name" "$TARGET_STAGE/$name" "$digest" >> "$ACQUISITION_LOG"
done
cp -a "$REUSE_PLAN_DIR/." "$PLAN_STAGE/"
[[ -f "$PLAN_STAGE/RUN_COMPLETE.json" ]] || die "Reuse plan lacks RUN_COMPLETE.json"
[[ -f "$PLAN_STAGE/output_manifest.json" ]] || die "Reuse plan lacks output_manifest.json"
mkdir -p "$OUTPUT_DIR/reports/reuse_plan"
for name in RUN_COMPLETE.json output_manifest.json run_manifest.json summary.json; do
  [[ -f "$PLAN_STAGE/$name" ]] || die "Reuse plan lacks $name"
  cp -p "$PLAN_STAGE/$name" "$OUTPUT_DIR/reports/reuse_plan/$name"
  printf 'reuse-plan-manifest\t%s\t%s\t%s\t%s\n' \
    "$name" "$REUSE_PLAN_DIR/$name" "$PLAN_STAGE/$name" \
    "$(sha256_file "$PLAN_STAGE/$name")" >> "$ACQUISITION_LOG"
done

echo "==> [2/10] Prepare deterministic preflight protein views"
"$PYTHON_BIN" "$HERE/prepare_regeneration_workspace.py" \
  --plan-dir "$PLAN_STAGE" \
  --target-benchmark-dir "$TARGET_STAGE" \
  --data-dir "$PFP_ROOT/data" \
  --limit-per-split "$PREFLIGHT_PER_SPLIT" \
  > "$OUTPUT_DIR/reports/preflight_workspace.json"

echo "==> [3/10] Download the same external dependencies used by the CAFA3 workflow"
cd "$PFP_ROOT"
PFP_ROOT="$PFP_ROOT" \
  bash "$HERE/generate_embeddings_dependencies.sh" \
  > "$OUTPUT_DIR/logs/dependencies.log" 2>&1
[[ -f external/dependency_env.sh ]] || die "Dependency script did not write dependency_env.sh"
# shellcheck disable=SC1091
source external/dependency_env.sh

echo "==> [4/10] Authenticate and extract Zijian's published embedding archives"
for name in \
  mmfp_embeddings_prott5.tar.gz \
  mmfp_embeddings_struct_ppi.tar.gz \
  mmfp_embeddings_text_temporal.tar.gz; do
  destination="$ARCHIVE_STAGE/$name"
  download_file "$MMFP_BASE_URL/$name?download=1" "$destination"
  observed="$(sha256_file "$destination")"
  wanted="$(archive_sha256 "$name")"
  [[ "$observed" == "$wanted" ]] || die "Published archive checksum mismatch: $name"
  tar -tzf "$destination" >/dev/null
  tar -xzf "$destination" -C "$PUBLISHED_ROOT"
  printf 'published-embedding-archive\t%s\t%s\t%s\t%s\n' \
    "$name" "$MMFP_BASE_URL/$name" "$destination" "$observed" >> "$ACQUISITION_LOG"
done
PUBLISHED_CACHE="$PUBLISHED_ROOT/data/embedding_cache"
[[ "$(file_count "$PUBLISHED_CACHE/prott5")" == "69811" ]] || die "Published ProtT5 count mismatch"
[[ "$(file_count "$PUBLISHED_CACHE/exp_text_embeddings_temporal")" == "69517" ]] || die "Published text count mismatch"
[[ "$(file_count "$PUBLISHED_CACHE/IF1")" == "67948" ]] || die "Published IF1 count mismatch"
[[ "$(file_count "$PUBLISHED_CACHE/ppi")" == "58294" ]] || die "Published PPI count mismatch"

echo "==> [5/10] Prepare the isolated IF1 runtime and compatibility copies"
IF1_NUMPY_OVERLAY="$WORK_DIR/if1_numpy_1_26_4"
install_mmfp_if1_numpy_overlay "$PYTHON_BIN" "$IF1_NUMPY_OVERLAY"
validate_mmfp_if1_env "$PYTHON_BIN" "$IF1_NUMPY_OVERLAY" \
  > "$OUTPUT_DIR/reports/if1_environment.json"

"$PYTHON_BIN" "$HERE/build_pfp_ppi_compat_copy.py" \
  --source "$PFP_ROOT/scripts/extract_ppi_embeddings.py" \
  --output "$RUNTIME_COMPAT/extract_ppi_embeddings.py" \
  --report "$OUTPUT_DIR/reports/pfp_ppi_compatibility.json"
"$PYTHON_BIN" "$HERE/build_pfp_if1_compat_copy.py" \
  --source "$PFP_ROOT/scripts/extract_esm_if1_embeddings.py" \
  --output "$RUNTIME_COMPAT/extract_esm_if1_embeddings.py" \
  --report "$OUTPUT_DIR/reports/pfp_if1_compatibility.json"

export CAFA_ASSESSMENT_DIR STRING_H5_FILE STRING_ALIAS_FILE
export PPI_EXTRACT_SCRIPT="$RUNTIME_COMPAT/extract_ppi_embeddings.py"
export IF1_EXTRACT_SCRIPT="$RUNTIME_COMPAT/extract_esm_if1_embeddings.py"
export IF1_PYTHON_BIN="$PYTHON_BIN"
export IF1_PYTHONPATH="$IF1_NUMPY_OVERLAY"
export TEXT_CUTOFF_DATE
export TEXT_REPORT_DIR="$PFP_ROOT/results/embedding_reports/text"
export HF_HOME="$WORK_DIR/model_cache/huggingface"
export TORCH_HOME="$WORK_DIR/model_cache/torch"
mkdir -p "$HF_HOME" "$TORCH_HOME"

run_parallel_modalities() {
  local phase="$1"
  local log_dir="$OUTPUT_DIR/logs/$phase"
  local rc=0
  mkdir -p "$log_dir"

  IFS=',' read -ra GPUS <<< "${CUDA_VISIBLE_DEVICES:-0,1,2}"
  gpu() { echo "${GPUS[$1]:-${GPUS[0]}}"; }
  echo "==> [$phase] Parallel embeddings on GPUs: ${GPUS[*]} (PPI on CPU)"

  CUDA_VISIBLE_DEVICES="$(gpu 0)" DEVICE=cuda \
    bash "$HERE/generate_embeddings_sequence.sh" > "$log_dir/sequence.log" 2>&1 &
  local pid_sequence=$!
  CUDA_VISIBLE_DEVICES="$(gpu 1)" \
    bash "$HERE/generate_embeddings_text_temporal_cls.sh" > "$log_dir/text.log" 2>&1 &
  local pid_text=$!
  CUDA_VISIBLE_DEVICES="$(gpu 2)" DEVICE=cuda \
    bash "$HERE/generate_embeddings_structure.sh" > "$log_dir/structure.log" 2>&1 &
  local pid_structure=$!
  CUDA_VISIBLE_DEVICES="" \
    bash "$HERE/generate_embeddings_ppi.sh" > "$log_dir/ppi.log" 2>&1 &
  local pid_ppi=$!

  for specification in \
    "sequence:$pid_sequence" "text:$pid_text" \
    "structure:$pid_structure" "ppi:$pid_ppi"; do
    local name="${specification%%:*}"
    local pid="${specification##*:}"
    local status=0
    wait "$pid" || status=$?
    printf '%s\t%s\t%s\n' "$phase" "$name" "$status" >> "$MODALITY_STATUS"
    if [[ "$status" == "0" ]]; then
      echo "==> [$phase/$name] OK"
    else
      echo "==> [$phase/$name] FAILED with status $status"
      rc=1
    fi
  done
  return "$rc"
}

echo "==> [6/10] Run the bounded preflight in parallel"
run_parallel_modalities preflight || die "A preflight modality failed; inspect logs/preflight"
preflight_count="$($PYTHON_BIN -c 'import json; print(json.load(open("data/regeneration_workspace_manifest.json"))["protein_count"])')"
[[ "$(file_count data/embedding_cache/prott5)" == "$preflight_count" ]] || \
  die "Preflight did not create every ProtT5 embedding"
[[ "$(file_count data/embedding_cache/exp_text_embeddings_temporal)" -gt 0 ]] || \
  die "Preflight produced no temporal text embeddings"
preflight_pdb="$(find data/alphafold_structures -maxdepth 1 -type f -name '*.pdb' -print 2>/dev/null | wc -l | tr -d ' ')"
if [[ "$preflight_pdb" -gt 0 && "$(file_count data/embedding_cache/IF1)" == "0" ]]; then
  die "Preflight downloaded AlphaFold structures but IF1 saved zero embeddings"
fi

# Keep downloaded models and the valid sequence/structure/PPI preflight files,
# but start temporal text extraction cleanly for the full population. PFP marks
# failed API attempts as processed in its checkpoint, so carrying that tiny
# preflight checkpoint forward could prevent a transiently failed protein from
# being retried during the real run.
rm -rf \
  data/embedding_cache/exp_text_embeddings \
  data/embedding_cache/exp_text_embeddings_temporal \
  data/embedding_cache/uniprot_text

echo "==> [7/10] Expand the workspace to all regeneration proteins"
"$PYTHON_BIN" "$HERE/prepare_regeneration_workspace.py" \
  --plan-dir "$PLAN_STAGE" \
  --target-benchmark-dir "$TARGET_STAGE" \
  --data-dir "$PFP_ROOT/data" \
  > "$OUTPUT_DIR/reports/full_workspace.json"

echo "==> [8/10] Generate all four modalities in parallel"
generation_status=0
run_parallel_modalities full || generation_status=$?

echo "==> [9/10] Archive generated outputs before strict assembly"
GENERATED_CACHE="$PFP_ROOT/data/embedding_cache"
for specification in \
  'prott5:prott5' \
  'text:exp_text_embeddings_temporal' \
  'structure:IF1' \
  'ppi:ppi'; do
  modality="${specification%%:*}"
  directory="${specification##*:}"
  if [[ -d "$GENERATED_CACHE/$directory" ]]; then
    tar -czf "$OUTPUT_DIR/archives/generated_${modality}.tar.gz" \
      -C "$PFP_ROOT" "data/embedding_cache/$directory"
  fi
done
if [[ -d data/embedding_cache/uniprot_text ]]; then
  tar -czf "$OUTPUT_DIR/archives/generated_text_provenance.tar.gz" \
    -C "$PFP_ROOT" data/embedding_cache/uniprot_text
fi
[[ "$generation_status" == "0" ]] || die "A full modality failed; partial archives were retained"

echo "==> [10/10] Assemble, validate, and package the complete cache"
"$PYTHON_BIN" "$HERE/assemble_contemporary_embedding_cache.py" \
  --plan-dir "$PLAN_STAGE" \
  --published-cache "$PUBLISHED_CACHE" \
  --generated-cache "$GENERATED_CACHE" \
  --output-cache "$FINAL_CACHE" \
  --report-dir "$OUTPUT_DIR/reports/assembly"

cp -p "$PFP_ROOT/data/regeneration_workspace_manifest.json" \
  "$OUTPUT_DIR/reports/regeneration_workspace_manifest.json"
if [[ -f "$PFP_ROOT/data/alphafold_coverage_results.txt" ]]; then
  cp -p "$PFP_ROOT/data/alphafold_coverage_results.txt" \
    "$OUTPUT_DIR/reports/alphafold_coverage_results.txt"
fi
if [[ -d "$PFP_ROOT/results/embedding_reports" ]]; then
  cp -a "$PFP_ROOT/results/embedding_reports" "$OUTPUT_DIR/reports/"
fi

tar -czf "$OUTPUT_DIR/archives/contemporary_embedding_cache.tar.gz" \
  -C "$FINAL_PACKAGE" data/embedding_cache

{
  printf 'sha256\tsize_bytes\tpath\n'
  for archive in "$OUTPUT_DIR"/archives/*.tar.gz; do
    printf '%s\t%s\t%s\n' \
      "$(sha256_file "$archive")" \
      "$(stat -c '%s' "$archive" 2>/dev/null || stat -f '%z' "$archive")" \
      "$(basename "$archive")"
  done
} > "$OUTPUT_DIR/reports/archive_manifest.tsv"

"$PYTHON_BIN" - "$OUTPUT_DIR" "$TEXT_CUTOFF_DATE" <<'PY'
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

output = Path(sys.argv[1])
payload = {
    "complete": True,
    "schema_version": 1,
    "completed_at": datetime.now(timezone.utc).isoformat(),
    "text_cutoff_date": sys.argv[2],
    "pfp_commit": os.environ.get("PFP_COMMIT", "unknown"),
    "framework_commit": os.environ.get("FRAMEWORK_COMMIT", "unknown"),
    "final_cache_archive": "archives/contemporary_embedding_cache.tar.gz",
    "assembly_summary": "reports/assembly/assembly_summary.json",
}
(output / "WORKFLOW_COMPLETE.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
PY

stop_disk_monitor
echo "==> Contemporary embedding cache complete: $OUTPUT_DIR"
