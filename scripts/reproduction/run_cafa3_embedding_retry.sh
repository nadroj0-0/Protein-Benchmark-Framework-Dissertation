#!/usr/bin/env bash
# Retry one CAFA3 embedding modality against the cumulative persistent state.

set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
FRAMEWORK_ROOT="$(cd "${HERE}/../.." && pwd)"
# shellcheck source=../reproduction_common.sh
source "$FRAMEWORK_ROOT/scripts/reproduction_common.sh"

PYTHON_BIN="${PYTHON_BIN:-python}"
CAFA_ASSESSMENT_COMMIT="${CAFA_ASSESSMENT_COMMIT:-d72f0a5abb66d3224bd808e2015b55f1c9d18340}"
CONTROL_COUNT="${CONTROL_COUNT:-20}"
EQUIVALENCE_MINIMUM="${EQUIVALENCE_MINIMUM:-5}"
PFP_ROOT=""
WORK_DIR=""
OUTPUT_DIR=""
EMBEDDING_STATE_ROOT=""
MODALITY=""
TEXT_CUTOFF_DATE="2016-02-17"
EMBEDDING_POLICY="$FRAMEWORK_ROOT/configs/cafa3_embedding_resume.json"

usage() {
  cat <<'EOF'
Usage: run_cafa3_embedding_retry.sh \
  --pfp-root PATH \
  --work-dir PATH \
  --output-dir PATH \
  --embedding-state-root PATH \
  --modality sequence|text|structure|ppi \
  [--text-cutoff-date YYYY-MM-DD] \
  [--embedding-policy PATH]

Only missing pairs for the selected modality are generated. Twenty accepted
control proteins are regenerated in the same subset and must be numerically
equivalent before any new arrays are merged into persistent state.
EOF
}

die() {
  echo "ERROR: $*" >&2
  exit 2
}

git_in_dir() {
  local directory="$1"
  shift
  (cd "$directory" && git "$@")
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --pfp-root) PFP_ROOT="$2"; shift 2 ;;
    --work-dir) WORK_DIR="$2"; shift 2 ;;
    --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
    --embedding-state-root) EMBEDDING_STATE_ROOT="$2"; shift 2 ;;
    --modality) MODALITY="$2"; shift 2 ;;
    --text-cutoff-date) TEXT_CUTOFF_DATE="$2"; shift 2 ;;
    --embedding-policy) EMBEDDING_POLICY="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) usage >&2; die "Unknown argument: $1" ;;
  esac
done

[[ -d "$PFP_ROOT/.git" ]] || die "PFP root is not a Git checkout: $PFP_ROOT"
[[ -n "$WORK_DIR" ]] || die "--work-dir is required"
[[ -n "$OUTPUT_DIR" ]] || die "--output-dir is required"
[[ -n "$EMBEDDING_STATE_ROOT" ]] || die "--embedding-state-root is required"
case "$MODALITY" in sequence|text|structure|ppi) ;; *) die "Invalid --modality: $MODALITY" ;; esac
[[ "$CONTROL_COUNT" =~ ^[1-9][0-9]*$ ]] || die "CONTROL_COUNT must be positive"
[[ "$EQUIVALENCE_MINIMUM" =~ ^[1-9][0-9]*$ ]] || die "EQUIVALENCE_MINIMUM must be positive"
[[ "$EQUIVALENCE_MINIMUM" -le "$CONTROL_COUNT" ]] || \
  die "EQUIVALENCE_MINIMUM cannot exceed CONTROL_COUNT"
[[ -f "$EMBEDDING_POLICY" ]] || die "Missing embedding policy: $EMBEDDING_POLICY"
EMBEDDING_POLICY="$(cd "$(dirname "$EMBEDDING_POLICY")" && pwd)/$(basename "$EMBEDDING_POLICY")"
[[ ! -e "$WORK_DIR" ]] || die "Work directory already exists: $WORK_DIR"
[[ ! -e "$OUTPUT_DIR" ]] || die "Output directory already exists: $OUTPUT_DIR"

PFP_ROOT="$(cd "$PFP_ROOT" && pwd)"
EMBEDDING_STATE_ROOT="$(cd "$EMBEDDING_STATE_ROOT" && pwd)"
mkdir -p "$WORK_DIR" "$OUTPUT_DIR/logs" "$OUTPUT_DIR/reports/embedding_state"
WORK_DIR="$(cd "$WORK_DIR" && pwd)"
OUTPUT_DIR="$(cd "$OUTPUT_DIR" && pwd)"
case "$EMBEDDING_STATE_ROOT/" in
  "$WORK_DIR/"*|"$PFP_ROOT/"*) die "Persistent state cannot live in job scratch" ;;
esac

RUNTIME_COMPAT="$WORK_DIR/runtime_compat"
REQUESTED="$WORK_DIR/requested_pairs.tsv"
CONTROLS="$WORK_DIR/control_pairs.tsv"
MODALITY_STATUS="$OUTPUT_DIR/reports/modality_status.tsv"
mkdir -p "$RUNTIME_COMPAT"
printf 'phase\tmodality\texit_status\n' > "$MODALITY_STATUS"

echo "==> [1/7] Validate the author-supplied environment"
validate_mmfp_env "$PYTHON_BIN" > "$OUTPUT_DIR/reports/environment_validation.txt"

echo "==> [2/7] Stage canonical CAFA3 and embedding dependencies"
cd "$PFP_ROOT"
PFP_ROOT="$PFP_ROOT" CAFA_ASSESSMENT_COMMIT="$CAFA_ASSESSMENT_COMMIT" \
  bash "$FRAMEWORK_ROOT/scripts/embeddings/generate_embeddings_dependencies.sh" \
  > "$OUTPUT_DIR/logs/dependencies.log" 2>&1
source external/dependency_env.sh
export CAFA_ASSESSMENT_DIR STRING_H5_FILE STRING_ALIAS_FILE CAFA3_RAW_DIR
bash "$FRAMEWORK_ROOT/scripts/embeddings/generate_embeddings_prepare_data.sh" \
  > "$OUTPUT_DIR/logs/prepare_data.log" 2>&1
"$PYTHON_BIN" "$FRAMEWORK_ROOT/scripts/embeddings/generate_embeddings_fasta.py" \
  --data-dir data > "$OUTPUT_DIR/logs/proteins_fasta.log" 2>&1

echo "==> [3/7] Recreate compatibility paths without editing PFP"
IF1_NUMPY_OVERLAY="$WORK_DIR/if1_numpy_1_26_4"
install_mmfp_if1_numpy_overlay "$PYTHON_BIN" "$IF1_NUMPY_OVERLAY"
validate_mmfp_if1_env "$PYTHON_BIN" "$IF1_NUMPY_OVERLAY" \
  > "$OUTPUT_DIR/reports/if1_environment.json"
"$PYTHON_BIN" "$FRAMEWORK_ROOT/scripts/embeddings/build_pfp_ppi_compat_copy.py" \
  --source "$PFP_ROOT/scripts/extract_ppi_embeddings.py" \
  --output "$RUNTIME_COMPAT/extract_ppi_embeddings.py" \
  --report "$OUTPUT_DIR/reports/pfp_ppi_compatibility.json"
"$PYTHON_BIN" "$FRAMEWORK_ROOT/scripts/embeddings/build_pfp_if1_compat_copy.py" \
  --source "$PFP_ROOT/scripts/extract_esm_if1_embeddings.py" \
  --output "$RUNTIME_COMPAT/extract_esm_if1_embeddings.py" \
  --report "$OUTPUT_DIR/reports/pfp_if1_compatibility.json"

pfp_commit="$(git_in_dir "$PFP_ROOT" rev-parse HEAD)"
framework_commit="$(git_in_dir "$FRAMEWORK_ROOT" rev-parse HEAD)"
"$PYTHON_BIN" "$FRAMEWORK_ROOT/scripts/embeddings/manage_resumable_embedding_state.py" \
  initialize \
  --state-root "$EMBEDDING_STATE_ROOT" \
  --benchmark-id cafa3-zijian-canonical \
  --benchmark-dir "$CAFA3_RAW_DIR" \
  --data-dir "$PFP_ROOT/data" \
  --policy "$EMBEDDING_POLICY" \
  --pfp-commit "$pfp_commit" \
  --framework-commit "$framework_commit" \
  --environment-report "$OUTPUT_DIR/reports/environment_validation.txt" \
  --source-file "go-ontology=$PFP_ROOT/data/go.obo" \
  --source-file "string-alias=$STRING_ALIAS_FILE" \
  --source-file "string-embeddings=$STRING_H5_FILE" \
  --source-file "pfp-prott5=$PFP_ROOT/scripts/extract_prott5_embeddings.py" \
  --source-file "pfp-text-extract=$PFP_ROOT/scripts/extract_uniprot_text.py" \
  --source-file "pfp-text-embed=$PFP_ROOT/scripts/embed_uniprot_descriptions.py" \
  --source-file "pfp-if1=$PFP_ROOT/scripts/extract_esm_if1_embeddings.py" \
  --source-file "pfp-ppi=$PFP_ROOT/scripts/extract_ppi_embeddings.py" \
  --source-file "framework-if1-compat=$RUNTIME_COMPAT/extract_esm_if1_embeddings.py" \
  --source-file "framework-ppi-compat=$RUNTIME_COMPAT/extract_ppi_embeddings.py" \
  --runtime-value "text_cutoff_date=$TEXT_CUTOFF_DATE" \
  --runtime-value "alphafold_acquisition=framework-bounded" \
  --runtime-value "alphafold_api_workers=${ALPHAFOLD_API_WORKERS:-8}" \
  --runtime-value "alphafold_download_workers=${ALPHAFOLD_DOWNLOAD_WORKERS:-8}" \
  > "$OUTPUT_DIR/reports/embedding_state_initialization.json"

echo "==> [4/7] Select only missing pairs and accepted equivalence controls"
"$PYTHON_BIN" "$FRAMEWORK_ROOT/scripts/embeddings/manage_resumable_embedding_state.py" \
  pending --state-root "$EMBEDDING_STATE_ROOT" --modality "$MODALITY" \
  --output "$REQUESTED" > "$OUTPUT_DIR/reports/pending_selection.json"
requested_count="$(($(wc -l < "$REQUESTED") - 1))"
if [[ "$requested_count" == "0" ]]; then
  "$PYTHON_BIN" "$FRAMEWORK_ROOT/scripts/embeddings/manage_resumable_embedding_state.py" \
    summary --state-root "$EMBEDDING_STATE_ROOT" \
    --report-dir "$OUTPUT_DIR/reports/embedding_state" \
    > "$OUTPUT_DIR/reports/embedding_state_summary.json"
  printf '{"complete":true,"no_work":true,"modality":"%s"}\n' "$MODALITY" \
    > "$OUTPUT_DIR/RETRY_COMPLETE.json"
  echo "No $MODALITY pairs need retrying"
  exit 0
fi
"$PYTHON_BIN" "$FRAMEWORK_ROOT/scripts/embeddings/manage_resumable_embedding_state.py" \
  controls --state-root "$EMBEDDING_STATE_ROOT" --modality "$MODALITY" \
  --count "$CONTROL_COUNT" --output "$CONTROLS" \
  > "$OUTPUT_DIR/reports/control_selection.json"
"$PYTHON_BIN" "$FRAMEWORK_ROOT/scripts/embeddings/prepare_embedding_retry_workspace.py" \
  --data-dir "$PFP_ROOT/data" \
  --requested-pairs "$REQUESTED" \
  --control-pairs "$CONTROLS" \
  --modality "$MODALITY" \
  --report "$OUTPUT_DIR/reports/retry_workspace.json"

export PPI_EXTRACT_SCRIPT="$RUNTIME_COMPAT/extract_ppi_embeddings.py"
export IF1_EXTRACT_SCRIPT="$RUNTIME_COMPAT/extract_esm_if1_embeddings.py"
export IF1_PYTHON_BIN="$PYTHON_BIN"
export IF1_PYTHONPATH="$IF1_NUMPY_OVERLAY"
export TEXT_CUTOFF_DATE
export TEXT_REPORT_DIR="$PFP_ROOT/results/embedding_reports/text"
export HF_HOME="$WORK_DIR/model_cache/huggingface"
export TORCH_HOME="$WORK_DIR/model_cache/torch"
export ALPHAFOLD_ACQUISITION_MODE=framework-bounded
export ALPHAFOLD_PERSISTENT_CACHE_DIR="$EMBEDDING_STATE_ROOT/source_cache/alphafold_structures"
export ALPHAFOLD_API_WORKERS="${ALPHAFOLD_API_WORKERS:-8}"
export ALPHAFOLD_DOWNLOAD_WORKERS="${ALPHAFOLD_DOWNLOAD_WORKERS:-8}"
export ALPHAFOLD_PREFETCH_REPORT="$OUTPUT_DIR/reports/alphafold_prefetch_retry.json"
mkdir -p "$HF_HOME" "$TORCH_HOME"

echo "==> [5/7] Run only the selected modality"
generation_status=0
case "$MODALITY" in
  sequence)
    DEVICE=cuda bash "$FRAMEWORK_ROOT/scripts/embeddings/generate_embeddings_sequence.sh" \
      > "$OUTPUT_DIR/logs/sequence.log" 2>&1 || generation_status=$? ;;
  text)
    bash "$FRAMEWORK_ROOT/scripts/embeddings/generate_embeddings_text_temporal_cls.sh" \
      > "$OUTPUT_DIR/logs/text.log" 2>&1 || generation_status=$? ;;
  structure)
    DEVICE=cuda bash "$FRAMEWORK_ROOT/scripts/embeddings/generate_embeddings_structure.sh" \
      > "$OUTPUT_DIR/logs/structure.log" 2>&1 || generation_status=$? ;;
  ppi)
    CUDA_VISIBLE_DEVICES="" bash "$FRAMEWORK_ROOT/scripts/embeddings/generate_embeddings_ppi.sh" \
      > "$OUTPUT_DIR/logs/ppi.log" 2>&1 || generation_status=$? ;;
esac
printf 'retry\t%s\t%s\n' "$MODALITY" "$generation_status" >> "$MODALITY_STATUS"

echo "==> [6/7] Enforce subset-equivalence before accepting retry outputs"
"$PYTHON_BIN" "$FRAMEWORK_ROOT/scripts/embeddings/verify_embedding_subset_equivalence.py" \
  --state-root "$EMBEDDING_STATE_ROOT" \
  --generated-cache-root "$PFP_ROOT/data/embedding_cache" \
  --control-pairs "$CONTROLS" \
  --modality "$MODALITY" \
  --minimum-compared "$EQUIVALENCE_MINIMUM" \
  --report "$OUTPUT_DIR/reports/subset_equivalence.json"

echo "==> [7/7] Atomically merge validated successes and retain every failure"
attempt_id="${JOB_ID:-local}_$(date -u +%Y%m%dT%H%M%SZ)_${MODALITY}"
merge_command=(
  "$PYTHON_BIN" "$FRAMEWORK_ROOT/scripts/embeddings/manage_resumable_embedding_state.py"
  merge
  --state-root "$EMBEDDING_STATE_ROOT"
  --generated-cache-root "$PFP_ROOT/data/embedding_cache"
  --attempt-id "$attempt_id"
  --requested-pairs "$REQUESTED"
  --allowed-extra-pairs "$CONTROLS"
  --modality-status "$MODALITY_STATUS"
  --report-dir "$OUTPUT_DIR/reports/embedding_state"
)
if [[ -f "$PFP_ROOT/data/alphafold_coverage_results.txt" ]]; then
  merge_command+=(--alphafold-report "$PFP_ROOT/data/alphafold_coverage_results.txt")
fi
if [[ -f "$OUTPUT_DIR/reports/alphafold_prefetch_retry.json" ]]; then
  merge_command+=(--alphafold-prefetch-report "$OUTPUT_DIR/reports/alphafold_prefetch_retry.json")
fi
"${merge_command[@]}" > "$OUTPUT_DIR/reports/embedding_state_merge.json"
printf '%s\n' "$generation_status" > "$OUTPUT_DIR/reports/generator_exit_status.txt"

if [[ -f "$EMBEDDING_STATE_ROOT/EMBEDDING_GATE_PASSED.json" ]]; then
  cp -p "$EMBEDDING_STATE_ROOT/EMBEDDING_GATE_PASSED.json" "$OUTPUT_DIR/"
else
  cp -p "$EMBEDDING_STATE_ROOT/GENERATION_INCOMPLETE.json" "$OUTPUT_DIR/"
fi
printf '{"complete":true,"no_work":false,"modality":"%s","generator_exit_status":%s}\n' \
  "$MODALITY" "$generation_status" > "$OUTPUT_DIR/RETRY_COMPLETE.json"
echo "Retry finished. Valid arrays were preserved even if the generator returned non-zero."
