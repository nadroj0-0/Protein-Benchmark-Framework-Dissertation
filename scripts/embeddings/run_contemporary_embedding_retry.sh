#!/usr/bin/env bash
# Retry one contemporary protein/modality subset into the archive-backed state.

set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
FRAMEWORK_ROOT="$(cd "${HERE}/../.." && pwd)"
# shellcheck source=../reproduction_common.sh
source "$FRAMEWORK_ROOT/scripts/reproduction_common.sh"

PYTHON_BIN="${PYTHON_BIN:-python}"
CONTROL_COUNT="${CONTROL_COUNT:-20}"
EQUIVALENCE_MINIMUM="${EQUIVALENCE_MINIMUM:-5}"
PFP_ROOT=""
WORK_DIR=""
OUTPUT_DIR=""
BENCHMARK_DIR=""
PLAN_DIR=""
STATE_ROOT=""
MODALITY=""
TEXT_CUTOFF_DATE="2025-03-08"

usage() {
  cat <<'EOF'
Usage: run_contemporary_embedding_retry.sh \
  --pfp-root PATH --work-dir PATH --output-dir PATH \
  --benchmark-dir PATH --plan-dir PATH --state-root PATH \
  --modality sequence|text|structure|ppi \
  [--text-cutoff-date YYYY-MM-DD]

Only currently missing pairs for one modality are generated. Accepted control
arrays are materialized from the immutable baseline archive or retry delta into
scratch and must reproduce before new arrays are merged.
EOF
}

die() { echo "ERROR: $*" >&2; exit 2; }
git_in_dir() { local directory="$1"; shift; (cd "$directory" && git "$@"); }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --pfp-root) PFP_ROOT="$2"; shift 2 ;;
    --work-dir) WORK_DIR="$2"; shift 2 ;;
    --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
    --benchmark-dir) BENCHMARK_DIR="$2"; shift 2 ;;
    --plan-dir) PLAN_DIR="$2"; shift 2 ;;
    --state-root) STATE_ROOT="$2"; shift 2 ;;
    --modality) MODALITY="$2"; shift 2 ;;
    --text-cutoff-date) TEXT_CUTOFF_DATE="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) usage >&2; die "Unknown argument: $1" ;;
  esac
done

[[ -d "$PFP_ROOT/.git" ]] || die "PFP root is not a Git checkout: $PFP_ROOT"
[[ -n "$WORK_DIR" ]] || die "--work-dir is required"
[[ -n "$OUTPUT_DIR" ]] || die "--output-dir is required"
[[ -d "$BENCHMARK_DIR" ]] || die "Missing benchmark: $BENCHMARK_DIR"
[[ -d "$PLAN_DIR" ]] || die "Missing reuse plan: $PLAN_DIR"
[[ -f "$STATE_ROOT/contract.json" ]] || die "State is not initialized: $STATE_ROOT"
case "$MODALITY" in sequence|text|structure|ppi) ;; *) die "Invalid modality: $MODALITY" ;; esac
[[ "$CONTROL_COUNT" =~ ^[1-9][0-9]*$ ]] || die "CONTROL_COUNT must be positive"
[[ "$EQUIVALENCE_MINIMUM" =~ ^[1-9][0-9]*$ ]] || die "EQUIVALENCE_MINIMUM must be positive"
[[ "$EQUIVALENCE_MINIMUM" -le "$CONTROL_COUNT" ]] || \
  die "EQUIVALENCE_MINIMUM cannot exceed CONTROL_COUNT"
[[ "$TEXT_CUTOFF_DATE" =~ ^[0-9]{4}-[0-9]{2}-[0-9]{2}$ ]] || \
  die "Invalid text cutoff date: $TEXT_CUTOFF_DATE"
[[ ! -e "$WORK_DIR" ]] || die "Work directory exists: $WORK_DIR"
[[ ! -e "$OUTPUT_DIR" ]] || die "Output directory exists: $OUTPUT_DIR"

PFP_ROOT="$(cd "$PFP_ROOT" && pwd)"
BENCHMARK_DIR="$(cd "$BENCHMARK_DIR" && pwd)"
PLAN_DIR="$(cd "$PLAN_DIR" && pwd)"
STATE_ROOT="$(cd "$STATE_ROOT" && pwd)"
mkdir -p "$WORK_DIR" "$OUTPUT_DIR/logs" "$OUTPUT_DIR/reports/embedding_state"
WORK_DIR="$(cd "$WORK_DIR" && pwd)"
OUTPUT_DIR="$(cd "$OUTPUT_DIR" && pwd)"
case "$STATE_ROOT/" in "$WORK_DIR/"*|"$PFP_ROOT/"*) die "State cannot live in scratch" ;; esac

RUNTIME_COMPAT="$WORK_DIR/runtime_compat"
REQUESTED="$WORK_DIR/requested_pairs.tsv"
CONTROLS="$WORK_DIR/control_pairs.tsv"
REFERENCE_CONTROLS="$WORK_DIR/reference_controls"
MODALITY_STATUS="$OUTPUT_DIR/reports/modality_status.tsv"
mkdir -p "$RUNTIME_COMPAT" "$REFERENCE_CONTROLS"
printf 'phase\tmodality\texit_status\n' > "$MODALITY_STATUS"

pfp_commit="$(git_in_dir "$PFP_ROOT" rev-parse HEAD)"
framework_commit="${FRAMEWORK_COMMIT:-$(git_in_dir "$FRAMEWORK_ROOT" rev-parse HEAD)}"
"$PYTHON_BIN" - "$STATE_ROOT/contract.json" "$pfp_commit" "$framework_commit" \
  "$TEXT_CUTOFF_DATE" <<'PY'
import json
import sys
contract = json.load(open(sys.argv[1]))
expected = {
    "pfp_commit": sys.argv[2],
    "framework_commit": sys.argv[3],
}
for key, observed in expected.items():
    if contract[key] != observed:
        raise SystemExit(f"State contract {key} mismatch: {contract[key]} != {observed}")
cutoff = contract.get("runtime", {}).get("text_cutoff_date")
if cutoff != sys.argv[4]:
    raise SystemExit(f"State contract text cutoff mismatch: {cutoff} != {sys.argv[4]}")
PY

echo "==> [1/9] Validate the author-supplied environment"
validate_mmfp_env "$PYTHON_BIN" > "$OUTPUT_DIR/reports/environment_validation.txt"
"$PYTHON_BIN" - "$STATE_ROOT/contract.json" \
  "$OUTPUT_DIR/reports/environment_validation.txt" <<'PY'
import hashlib
import json
import sys
contract = json.load(open(sys.argv[1]))
observed = hashlib.sha256(open(sys.argv[2], "rb").read()).hexdigest()
expected = contract.get("environment", {}).get("sha256")
if observed != expected:
    raise SystemExit(f"State contract environment mismatch: {expected} != {observed}")
PY

echo "==> [2/9] Stage embedding dependencies"
cd "$PFP_ROOT"
PFP_ROOT="$PFP_ROOT" bash "$HERE/generate_embeddings_dependencies.sh" \
  > "$OUTPUT_DIR/logs/dependencies.log" 2>&1
source external/dependency_env.sh
export CAFA_ASSESSMENT_DIR STRING_H5_FILE STRING_ALIAS_FILE CAFA3_RAW_DIR

echo "==> [3/9] Recreate runtime-only PFP compatibility copies"
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

echo "==> [4/9] Select pending pairs and accepted controls"
"$PYTHON_BIN" "$HERE/manage_resumable_embedding_state.py" pending \
  --state-root "$STATE_ROOT" --modality "$MODALITY" --output "$REQUESTED" \
  > "$OUTPUT_DIR/reports/pending_selection.json"
requested_count="$(($(wc -l < "$REQUESTED") - 1))"
if [[ "$requested_count" == "0" ]]; then
  "$PYTHON_BIN" "$HERE/manage_resumable_embedding_state.py" summary \
    --state-root "$STATE_ROOT" --report-dir "$OUTPUT_DIR/reports/embedding_state" \
    > "$OUTPUT_DIR/reports/embedding_state_summary.json"
  printf '{"complete":true,"no_work":true,"modality":"%s"}\n' "$MODALITY" \
    > "$OUTPUT_DIR/RETRY_COMPLETE.json"
  echo "No $MODALITY pairs need retrying"
  exit 0
fi
"$PYTHON_BIN" "$HERE/manage_resumable_embedding_state.py" controls \
  --state-root "$STATE_ROOT" --modality "$MODALITY" --count "$CONTROL_COUNT" \
  --output "$CONTROLS" > "$OUTPUT_DIR/reports/control_selection.json"

echo "==> [5/9] Build the exact contemporary retry workspace"
"$PYTHON_BIN" "$HERE/prepare_contemporary_retry_workspace.py" \
  --plan-dir "$PLAN_DIR" --target-benchmark-dir "$BENCHMARK_DIR" \
  --data-dir "$PFP_ROOT/data" --requested-pairs "$REQUESTED" \
  --control-pairs "$CONTROLS" --modality "$MODALITY" \
  --report "$OUTPUT_DIR/reports/retry_workspace.json"
"$PYTHON_BIN" "$HERE/manage_resumable_embedding_state.py" materialize \
  --state-root "$STATE_ROOT" --pairs "$CONTROLS" \
  --output-cache-root "$REFERENCE_CONTROLS" \
  --report "$OUTPUT_DIR/reports/control_materialization.json"

export PPI_EXTRACT_SCRIPT="$RUNTIME_COMPAT/extract_ppi_embeddings.py"
export IF1_EXTRACT_SCRIPT="$RUNTIME_COMPAT/extract_esm_if1_embeddings.py"
export IF1_PYTHON_BIN="$PYTHON_BIN"
export IF1_PYTHONPATH="$IF1_NUMPY_OVERLAY"
export TEXT_CUTOFF_DATE
export TEXT_REPORT_DIR="$PFP_ROOT/results/embedding_reports/text"
export HF_HOME="$WORK_DIR/model_cache/huggingface"
export TORCH_HOME="$WORK_DIR/model_cache/torch"
export ALPHAFOLD_ACQUISITION_MODE=framework-bounded
export ALPHAFOLD_PERSISTENT_CACHE_DIR="$STATE_ROOT/source_cache/alphafold_structures"
export ALPHAFOLD_API_WORKERS="${ALPHAFOLD_API_WORKERS:-8}"
export ALPHAFOLD_DOWNLOAD_WORKERS="${ALPHAFOLD_DOWNLOAD_WORKERS:-8}"
export ALPHAFOLD_PREFETCH_REPORT="$OUTPUT_DIR/reports/alphafold_prefetch_retry.json"
mkdir -p "$HF_HOME" "$TORCH_HOME"

echo "==> [6/9] Generate only missing $MODALITY pairs"
generation_status=0
case "$MODALITY" in
  sequence)
    DEVICE=cuda bash "$HERE/generate_embeddings_sequence.sh" \
      > "$OUTPUT_DIR/logs/sequence.log" 2>&1 || generation_status=$? ;;
  text)
    bash "$HERE/generate_embeddings_text_temporal_cls.sh" \
      > "$OUTPUT_DIR/logs/text.log" 2>&1 || generation_status=$? ;;
  structure)
    DEVICE=cuda bash "$HERE/generate_embeddings_structure.sh" \
      > "$OUTPUT_DIR/logs/structure.log" 2>&1 || generation_status=$? ;;
  ppi)
    CUDA_VISIBLE_DEVICES="" bash "$HERE/generate_embeddings_ppi.sh" \
      > "$OUTPUT_DIR/logs/ppi.log" 2>&1 || generation_status=$? ;;
esac
printf 'retry\t%s\t%s\n' "$MODALITY" "$generation_status" >> "$MODALITY_STATUS"

echo "==> [7/9] Prove subset generation matches accepted controls"
"$PYTHON_BIN" "$HERE/verify_embedding_subset_equivalence.py" \
  --state-root "$STATE_ROOT" --reference-cache-root "$REFERENCE_CONTROLS" \
  --generated-cache-root "$PFP_ROOT/data/embedding_cache" \
  --control-pairs "$CONTROLS" --modality "$MODALITY" \
  --minimum-compared "$EQUIVALENCE_MINIMUM" \
  --report "$OUTPUT_DIR/reports/subset_equivalence.json"

echo "==> [8/9] Atomically merge valid retry outputs"
attempt_id="${JOB_ID:-local}_$(date -u +%Y%m%dT%H%M%SZ)_${MODALITY}"
merge_command=(
  "$PYTHON_BIN" "$HERE/manage_resumable_embedding_state.py" merge
  --state-root "$STATE_ROOT"
  --generated-cache-root "$PFP_ROOT/data/embedding_cache"
  --attempt-id "$attempt_id"
  --requested-pairs "$REQUESTED"
  --allowed-extra-pairs "$CONTROLS"
  --modality-status "$MODALITY_STATUS"
  --report-dir "$OUTPUT_DIR/reports/embedding_state"
)
[[ ! -f "$PFP_ROOT/data/alphafold_coverage_results.txt" ]] || \
  merge_command+=(--alphafold-report "$PFP_ROOT/data/alphafold_coverage_results.txt")
[[ ! -f "$ALPHAFOLD_PREFETCH_REPORT" ]] || \
  merge_command+=(--alphafold-prefetch-report "$ALPHAFOLD_PREFETCH_REPORT")
"${merge_command[@]}" > "$OUTPUT_DIR/reports/embedding_state_merge.json"

echo "==> [9/9] Publish compact retry status"
printf '%s\n' "$generation_status" > "$OUTPUT_DIR/reports/generator_exit_status.txt"
if [[ -f "$STATE_ROOT/EMBEDDING_GATE_PASSED.json" ]]; then
  cp -p "$STATE_ROOT/EMBEDDING_GATE_PASSED.json" "$OUTPUT_DIR/"
else
  cp -p "$STATE_ROOT/GENERATION_INCOMPLETE.json" "$OUTPUT_DIR/"
fi
printf '{"complete":true,"no_work":false,"modality":"%s","generator_exit_status":%s}\n' \
  "$MODALITY" "$generation_status" > "$OUTPUT_DIR/RETRY_COMPLETE.json"
echo "Retry complete. Valid arrays were retained; missing pairs remain pending."
