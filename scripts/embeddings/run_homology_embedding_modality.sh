#!/usr/bin/env bash
# Generate one source-resolved homology embedding delta in a disposable PFP clone.

set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
FRAMEWORK_ROOT="$(cd "${HERE}/../.." && pwd)"
# shellcheck source=../reproduction_common.sh
source "$FRAMEWORK_ROOT/scripts/reproduction_common.sh"

PYTHON_BIN="${PYTHON_BIN:-python}"
PREFLIGHT_PER_SPLIT="${PREFLIGHT_PER_SPLIT:-1}"
PFP_ROOT=""
WORK_DIR=""
OUTPUT_DIR=""
BENCHMARK_DIR=""
LEDGER_DIR=""
MODALITY=""
TEXT_CUTOFF_DATE=""

usage() {
  cat <<'EOF'
Usage: run_homology_embedding_modality.sh \
  --pfp-root PATH --work-dir PATH --output-dir PATH \
  --benchmark-dir PATH --ledger-dir PATH \
  --modality sequence|text|structure|ppi [--text-cutoff-date YYYY-MM-DD] \
  [--artifact-catalog PATH]
EOF
}

die() { echo "ERROR: $*" >&2; exit 2; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --pfp-root) PFP_ROOT="$2"; shift 2 ;;
    --work-dir) WORK_DIR="$2"; shift 2 ;;
    --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
    --benchmark-dir) BENCHMARK_DIR="$2"; shift 2 ;;
    --ledger-dir) LEDGER_DIR="$2"; shift 2 ;;
    --modality) MODALITY="$2"; shift 2 ;;
    --text-cutoff-date) TEXT_CUTOFF_DATE="$2"; shift 2 ;;
    --artifact-catalog) ARTIFACT_CATALOG="$2"; export ARTIFACT_CATALOG; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) usage >&2; die "Unknown argument: $1" ;;
  esac
done

[[ -d "$PFP_ROOT/.git" ]] || die "PFP root is not a Git checkout: $PFP_ROOT"
[[ -d "$BENCHMARK_DIR" ]] || die "Missing benchmark: $BENCHMARK_DIR"
[[ -d "$LEDGER_DIR" ]] || die "Missing source-resolved ledger: $LEDGER_DIR"
[[ -n "$WORK_DIR" && ! -e "$WORK_DIR" ]] || die "Work directory is missing or exists"
[[ -n "$OUTPUT_DIR" && ! -e "$OUTPUT_DIR" ]] || die "Output directory is missing or exists"
case "$MODALITY" in sequence|text|structure|ppi) ;; *) die "Invalid modality: $MODALITY" ;; esac
if [[ -n "$TEXT_CUTOFF_DATE" ]]; then
  [[ "$MODALITY" == "text" ]] || die "--text-cutoff-date is valid only for text"
  [[ "$TEXT_CUTOFF_DATE" =~ ^[0-9]{4}-[0-9]{2}-[0-9]{2}$ ]] || \
    die "--text-cutoff-date must be YYYY-MM-DD"
fi
[[ "$PREFLIGHT_PER_SPLIT" =~ ^[1-9][0-9]*$ ]] || die "PREFLIGHT_PER_SPLIT must be positive"
command -v "$PYTHON_BIN" >/dev/null 2>&1 || die "Python not found: $PYTHON_BIN"

artifact_catalog_configure "$FRAMEWORK_ROOT" "${ARTIFACT_CATALOG:-}"
mkdir -p "$WORK_DIR" "$OUTPUT_DIR/logs" "$OUTPUT_DIR/reports" "$OUTPUT_DIR/artifacts"
WORK_DIR="$(cd "$WORK_DIR" && pwd)"
OUTPUT_DIR="$(cd "$OUTPUT_DIR" && pwd)"
PFP_ROOT="$(cd "$PFP_ROOT" && pwd)"
BENCHMARK_DIR="$(cd "$BENCHMARK_DIR" && pwd)"
LEDGER_DIR="$(cd "$LEDGER_DIR" && pwd)"
RUNTIME_COMPAT="$WORK_DIR/runtime_compat"
mkdir -p "$RUNTIME_COMPAT"

echo "==> [1/7] Stage the modality's external dependencies"
cd "$PFP_ROOT"
if [[ "$MODALITY" != "sequence" ]]; then
  PFP_ROOT="$PFP_ROOT" EMBEDDING_DEPENDENCY_PROFILE="$MODALITY" \
    bash "$HERE/generate_embeddings_dependencies.sh" \
    > "$OUTPUT_DIR/logs/dependencies.log" 2>&1
  # shellcheck disable=SC1091
  source external/dependency_env.sh
  export CAFA_ASSESSMENT_DIR STRING_H5_FILE STRING_ALIAS_FILE
else
  mkdir -p data
  printf 'Sequence generation has no static external dependency bundle.\n' \
    > "$OUTPUT_DIR/logs/dependencies.log"
fi

echo "==> [2/7] Create runtime-only compatibility copies"
if [[ "$MODALITY" == "ppi" ]]; then
  "$PYTHON_BIN" "$HERE/build_pfp_ppi_compat_copy.py" \
    --source "$PFP_ROOT/scripts/extract_ppi_embeddings.py" \
    --output "$RUNTIME_COMPAT/extract_ppi_embeddings.py" \
    --report "$OUTPUT_DIR/reports/pfp_ppi_compatibility.json"
  export PPI_EXTRACT_SCRIPT="$RUNTIME_COMPAT/extract_ppi_embeddings.py"
fi
if [[ "$MODALITY" == "structure" ]]; then
  IF1_NUMPY_OVERLAY="$WORK_DIR/if1_numpy_1_26_4"
  install_mmfp_if1_numpy_overlay "$PYTHON_BIN" "$IF1_NUMPY_OVERLAY"
  validate_mmfp_if1_env "$PYTHON_BIN" "$IF1_NUMPY_OVERLAY" \
    > "$OUTPUT_DIR/reports/if1_environment.json"
  "$PYTHON_BIN" "$HERE/build_pfp_if1_compat_copy.py" \
    --source "$PFP_ROOT/scripts/extract_esm_if1_embeddings.py" \
    --output "$RUNTIME_COMPAT/extract_esm_if1_embeddings.py" \
    --report "$OUTPUT_DIR/reports/pfp_if1_compatibility.json"
  export IF1_EXTRACT_SCRIPT="$RUNTIME_COMPAT/extract_esm_if1_embeddings.py"
  export IF1_PYTHON_BIN="$PYTHON_BIN"
  export IF1_PYTHONPATH="$IF1_NUMPY_OVERLAY"
  export ALPHAFOLD_ACQUISITION_MODE=framework-bounded
  export ALPHAFOLD_PERSISTENT_CACHE_DIR="$WORK_DIR/alphafold_source_cache"
  export ALPHAFOLD_API_WORKERS="${ALPHAFOLD_API_WORKERS:-8}"
  export ALPHAFOLD_DOWNLOAD_WORKERS="${ALPHAFOLD_DOWNLOAD_WORKERS:-8}"
fi
export HF_HOME="${HF_HOME:-$WORK_DIR/model_cache/huggingface}"
export TORCH_HOME="${TORCH_HOME:-$WORK_DIR/model_cache/torch}"
mkdir -p "$HF_HOME" "$TORCH_HOME"

run_text() {
  local phase="$1"
  local raw="data/embedding_cache/exp_text_embeddings"
  local final="data/embedding_cache/exp_text_embeddings_temporal"
  rm -rf "$final"
  if [[ -n "$TEXT_CUTOFF_DATE" ]]; then
    CAFA_ASSESSMENT_DIR="$CAFA_ASSESSMENT_DIR" \
      TEXT_CUTOFF_DATE="$TEXT_CUTOFF_DATE" \
      TEXT_REPORT_DIR="$OUTPUT_DIR/reports/${phase}_text" \
      PYTHON_BIN="$PYTHON_BIN" \
      bash "$HERE/generate_embeddings_text_temporal_cls.sh"
    return
  fi
  "$PYTHON_BIN" scripts/extract_uniprot_text.py extract-current
  "$PYTHON_BIN" scripts/embed_uniprot_descriptions.py --data-dir data &
  local embed_pid=$!
  local embed_status=0
  "$PYTHON_BIN" "$HERE/reduce_text_embeddings_to_cls.py" \
    --directory "$raw" --watch-pid "$embed_pid" \
    --report "$OUTPUT_DIR/reports/${phase}_cls_reduction.json" &
  local reducer_pid=$!
  wait "$embed_pid" || embed_status=$?
  wait "$reducer_pid"
  [[ "$embed_status" == "0" ]] || return "$embed_status"
  mv "$raw" "$final"
}

run_modality() {
  local phase="$1"
  export ALPHAFOLD_PREFETCH_REPORT="$OUTPUT_DIR/reports/${phase}_alphafold_prefetch.json"
  case "$MODALITY" in
    sequence)
      DEVICE=cuda bash "$HERE/generate_embeddings_sequence.sh" ;;
    text)
      run_text "$phase" ;;
    structure)
      DEVICE=cuda bash "$HERE/generate_embeddings_structure.sh" ;;
    ppi)
      CUDA_VISIBLE_DEVICES="" bash "$HERE/generate_embeddings_ppi.sh" ;;
  esac
}

echo "==> [3/7] Build and run a bounded preflight"
"$PYTHON_BIN" "$HERE/prepare_homology_embedding_workspace.py" \
  --ledger-dir "$LEDGER_DIR" --target-benchmark-dir "$BENCHMARK_DIR" \
  --data-dir "$PFP_ROOT/data" --modality "$MODALITY" \
  --limit-per-split "$PREFLIGHT_PER_SPLIT" \
  --report "$OUTPUT_DIR/reports/preflight_workspace.json"
run_modality preflight > "$OUTPUT_DIR/logs/preflight.log" 2>&1

if [[ "$MODALITY" == "text" ]]; then
  rm -rf \
    data/embedding_cache/exp_text_embeddings \
    data/embedding_cache/exp_text_embeddings_temporal \
    data/embedding_cache/uniprot_text
fi

echo "==> [4/7] Expand to every ledger-selected $MODALITY pair"
"$PYTHON_BIN" "$HERE/prepare_homology_embedding_workspace.py" \
  --ledger-dir "$LEDGER_DIR" --target-benchmark-dir "$BENCHMARK_DIR" \
  --data-dir "$PFP_ROOT/data" --modality "$MODALITY" \
  --report "$OUTPUT_DIR/reports/full_workspace.json"

echo "==> [5/7] Generate the complete $MODALITY delta"
run_modality full > "$OUTPUT_DIR/logs/full.log" 2>&1

echo "==> [6/7] Validate and archive the generated delta"
archive="$OUTPUT_DIR/artifacts/generated_${MODALITY}.tar.gz"
assembly="$OUTPUT_DIR/artifacts/generated_${MODALITY}_assembly.tsv.gz"
"$PYTHON_BIN" "$HERE/build_embedding_baseline_archive.py" \
  --generated-cache-root "$PFP_ROOT/data/embedding_cache" \
  --data-dir "$PFP_ROOT/data" --policy "$FRAMEWORK_ROOT/configs/homology_embedding_generation.json" \
  --only-modality "$MODALITY" --archive "$archive" \
  --assembly-report "$assembly" --report "$OUTPUT_DIR/reports/delta_archive.json"

"$PYTHON_BIN" - "$OUTPUT_DIR/reports/full_workspace.json" \
  "$OUTPUT_DIR/reports/delta_archive.json" "$MODALITY" <<'PY'
import json
import sys
workspace = json.load(open(sys.argv[1], encoding="utf-8"))
archive = json.load(open(sys.argv[2], encoding="utf-8"))
modality = sys.argv[3]
requested = int(workspace["protein_count"])
available = int(archive["available_pairs"])
if modality == "sequence" and available != requested:
    raise SystemExit(f"ProtT5 delta is incomplete: {available} / {requested}")
if modality != "sequence" and available == 0:
    raise SystemExit(f"{modality} generated no usable arrays")
PY

if [[ "$MODALITY" == "text" && -d data/embedding_cache/uniprot_text ]]; then
  tar -czf "$OUTPUT_DIR/artifacts/generated_text_provenance.tar.gz" \
    -C "$PFP_ROOT" data/embedding_cache/uniprot_text
fi
if [[ -f data/alphafold_coverage_results.txt ]]; then
  cp -p data/alphafold_coverage_results.txt "$OUTPUT_DIR/reports/"
fi

echo "==> [7/7] Publish the completion marker"
"$PYTHON_BIN" - "$OUTPUT_DIR" "$MODALITY" "$LEDGER_DIR" "$BENCHMARK_DIR" \
  "$TEXT_CUTOFF_DATE" <<'PY'
import hashlib
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

root = Path(sys.argv[1])
modality = sys.argv[2]
archive = root / "artifacts" / f"generated_{modality}.tar.gz"
assembly = root / "artifacts" / f"generated_{modality}_assembly.tsv.gz"
def sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()
payload = {
    "schema_version": 1,
    "complete": True,
    "analysis_kind": "homology_embedding_modality_delta",
    "completed_at": datetime.now(timezone.utc).isoformat(),
    "modality": modality,
    "ledger_dir": str(Path(sys.argv[3]).resolve()),
    "ledger_output_manifest_sha256": sha(Path(sys.argv[3]) / "output_manifest.json"),
    "benchmark_dir": str(Path(sys.argv[4]).resolve()),
    "text_cutoff_date": sys.argv[5] or None,
    "framework_commit": os.environ.get("FRAMEWORK_COMMIT", "unknown"),
    "pfp_commit": os.environ.get("PFP_COMMIT", "unknown"),
    "archive": str(archive.relative_to(root)),
    "archive_sha256": sha(archive),
    "assembly_report": str(assembly.relative_to(root)),
    "assembly_report_sha256": sha(assembly),
}
(root / "WORKFLOW_COMPLETE.json").write_text(
    json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
)
PY

echo "==> Homology $MODALITY delta complete: $OUTPUT_DIR"
