#!/usr/bin/env bash
# Diagnose text/structure numerical repeatability without merging any arrays.

set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
FRAMEWORK_ROOT="$(cd "${HERE}/../.." && pwd)"
# shellcheck source=../reproduction_common.sh
source "$FRAMEWORK_ROOT/scripts/reproduction_common.sh"

PYTHON_BIN="${PYTHON_BIN:-python}"
CONTROL_COUNT="${CONTROL_COUNT:-20}"
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
Usage: run_contemporary_embedding_reproducibility.sh \
  --pfp-root PATH --work-dir PATH --output-dir PATH \
  --benchmark-dir PATH --plan-dir PATH --state-root PATH \
  --modality text|structure [--control-count 20] \
  [--text-cutoff-date YYYY-MM-DD] [--artifact-catalog PATH]

The workflow selects accepted controls, prepares one immutable input view, runs
the same controls twice on the allocated device, and compares both repeats with
each other and the accepted baseline. It never merges into the embedding state.
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
    --control-count) CONTROL_COUNT="$2"; shift 2 ;;
    --text-cutoff-date) TEXT_CUTOFF_DATE="$2"; shift 2 ;;
    --artifact-catalog) ARTIFACT_CATALOG="$2"; export ARTIFACT_CATALOG; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) usage >&2; die "Unknown argument: $1" ;;
  esac
done
artifact_catalog_configure "$FRAMEWORK_ROOT" "${ARTIFACT_CATALOG:-}"

[[ -d "$PFP_ROOT/.git" ]] || die "PFP root is not a Git checkout: $PFP_ROOT"
[[ -n "$WORK_DIR" ]] || die "--work-dir is required"
[[ -n "$OUTPUT_DIR" ]] || die "--output-dir is required"
[[ -d "$BENCHMARK_DIR" ]] || die "Missing benchmark: $BENCHMARK_DIR"
[[ -d "$PLAN_DIR" ]] || die "Missing reuse plan: $PLAN_DIR"
[[ -f "$STATE_ROOT/contract.json" ]] || die "State is not initialized: $STATE_ROOT"
case "$MODALITY" in text|structure) ;; *) die "--modality must be text or structure" ;; esac
[[ "$CONTROL_COUNT" =~ ^[1-9][0-9]*$ ]] || die "--control-count must be positive"
[[ "$TEXT_CUTOFF_DATE" =~ ^[0-9]{4}-[0-9]{2}-[0-9]{2}$ ]] || \
  die "Invalid text cutoff date: $TEXT_CUTOFF_DATE"
[[ ! -e "$WORK_DIR" ]] || die "Work directory exists: $WORK_DIR"
[[ ! -e "$OUTPUT_DIR" ]] || die "Output directory exists: $OUTPUT_DIR"

PFP_ROOT="$(cd "$PFP_ROOT" && pwd)"
BENCHMARK_DIR="$(cd "$BENCHMARK_DIR" && pwd)"
PLAN_DIR="$(cd "$PLAN_DIR" && pwd)"
STATE_ROOT="$(cd "$STATE_ROOT" && pwd)"
mkdir -p "$WORK_DIR" "$OUTPUT_DIR/logs" "$OUTPUT_DIR/reports" \
  "$OUTPUT_DIR/inputs" "$OUTPUT_DIR/generated"
WORK_DIR="$(cd "$WORK_DIR" && pwd)"
OUTPUT_DIR="$(cd "$OUTPUT_DIR" && pwd)"
case "$STATE_ROOT/" in "$WORK_DIR/"*|"$PFP_ROOT/"*) die "State cannot live in scratch" ;; esac

RUNTIME_COMPAT="$WORK_DIR/runtime_compat"
CONTROLS="$OUTPUT_DIR/controls.tsv"
EMPTY_CONTROLS="$WORK_DIR/empty_controls.tsv"
REFERENCE_ROOT="$OUTPUT_DIR/reference"
REPEAT_ONE_ROOT="$OUTPUT_DIR/generated/repeat_1"
REPEAT_TWO_ROOT="$OUTPUT_DIR/generated/repeat_2"
mkdir -p "$RUNTIME_COMPAT" "$REFERENCE_ROOT" "$REPEAT_ONE_ROOT" "$REPEAT_TWO_ROOT"
printf 'protein_id\tmodality\tsequence_sha256\n' > "$EMPTY_CONTROLS"

pfp_commit="$(git_in_dir "$PFP_ROOT" rev-parse HEAD)"
framework_commit="${FRAMEWORK_COMMIT:-$(git_in_dir "$FRAMEWORK_ROOT" rev-parse HEAD)}"
"$PYTHON_BIN" - "$STATE_ROOT/contract.json" "$pfp_commit" "$TEXT_CUTOFF_DATE" \
  "pfp-text-extract=$PFP_ROOT/scripts/extract_uniprot_text.py" \
  "pfp-text-embed=$PFP_ROOT/scripts/embed_uniprot_descriptions.py" \
  "pfp-if1=$PFP_ROOT/scripts/extract_esm_if1_embeddings.py" \
  "framework-if1-compat=$HERE/build_pfp_if1_compat_copy.py" <<'PY'
import hashlib
import json
import sys

contract = json.load(open(sys.argv[1]))
if contract["pfp_commit"] != sys.argv[2]:
    raise SystemExit(
        f"State PFP commit mismatch: {contract['pfp_commit']} != {sys.argv[2]}"
    )
cutoff = contract.get("runtime", {}).get("text_cutoff_date")
if cutoff != sys.argv[3]:
    raise SystemExit(f"State text cutoff mismatch: {cutoff} != {sys.argv[3]}")
sources = {entry["label"]: entry["sha256"] for entry in contract["source_files"]}
for specification in sys.argv[4:]:
    label, path = specification.split("=", 1)
    expected = sources.get(label)
    if expected is None:
        raise SystemExit(f"State contract has no source hash for {label}")
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    observed = digest.hexdigest()
    if observed != expected:
        raise SystemExit(f"State source mismatch for {label}: {expected} != {observed}")
PY

echo "==> [1/8] Validate the frozen MMFP environment"
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
    raise SystemExit(f"State environment mismatch: {expected} != {observed}")
PY

echo "==> [2/8] Stage only the selected modality dependencies"
cd "$PFP_ROOT"
PFP_ROOT="$PFP_ROOT" EMBEDDING_DEPENDENCY_PROFILE="$MODALITY" \
  bash "$HERE/generate_embeddings_dependencies.sh" \
  > "$OUTPUT_DIR/logs/dependencies.log" 2>&1
source external/dependency_env.sh
export CAFA_ASSESSMENT_DIR

IF1_EXTRACT_SCRIPT=""
IF1_NUMPY_OVERLAY=""
if [[ "$MODALITY" == "structure" ]]; then
  echo "==> Preparing the validated IF1 runtime compatibility copy"
  IF1_NUMPY_OVERLAY="$WORK_DIR/if1_numpy_1_26_4"
  install_mmfp_if1_numpy_overlay "$PYTHON_BIN" "$IF1_NUMPY_OVERLAY"
  validate_mmfp_if1_env "$PYTHON_BIN" "$IF1_NUMPY_OVERLAY" \
    > "$OUTPUT_DIR/reports/if1_environment.json"
  IF1_EXTRACT_SCRIPT="$RUNTIME_COMPAT/extract_esm_if1_embeddings.py"
  "$PYTHON_BIN" "$HERE/build_pfp_if1_compat_copy.py" \
    --source "$PFP_ROOT/scripts/extract_esm_if1_embeddings.py" \
    --output "$IF1_EXTRACT_SCRIPT" \
    --report "$OUTPUT_DIR/reports/pfp_if1_compatibility.json"
fi

echo "==> [3/8] Select and materialize accepted controls"
"$PYTHON_BIN" "$HERE/manage_resumable_embedding_state.py" controls \
  --state-root "$STATE_ROOT" --modality "$MODALITY" --count "$CONTROL_COUNT" \
  --output "$CONTROLS" > "$OUTPUT_DIR/reports/control_selection.json"
selected_count="$(($(wc -l < "$CONTROLS") - 1))"
[[ "$selected_count" == "$CONTROL_COUNT" ]] || \
  die "Selected $selected_count controls, expected $CONTROL_COUNT"
"$PYTHON_BIN" "$HERE/manage_resumable_embedding_state.py" materialize \
  --state-root "$STATE_ROOT" --pairs "$CONTROLS" \
  --output-cache-root "$REFERENCE_ROOT" \
  --report "$OUTPUT_DIR/reports/control_materialization.json"

echo "==> [4/8] Build one shared 20-protein input workspace"
"$PYTHON_BIN" "$HERE/prepare_contemporary_retry_workspace.py" \
  --plan-dir "$PLAN_DIR" --target-benchmark-dir "$BENCHMARK_DIR" \
  --data-dir "$PFP_ROOT/data" --requested-pairs "$CONTROLS" \
  --control-pairs "$EMPTY_CONTROLS" --modality "$MODALITY" \
  --report "$OUTPUT_DIR/reports/reproducibility_workspace.json"

export HF_HOME="${HF_HOME:-$STATE_ROOT/source_cache/model_weights/huggingface}"
export TORCH_HOME="${TORCH_HOME:-$STATE_ROOT/source_cache/model_weights/torch}"
mkdir -p "$HF_HOME" "$TORCH_HOME"

runtime_sources=(
  --source-file "pfp-text-extract=$PFP_ROOT/scripts/extract_uniprot_text.py"
  --source-file "pfp-text-embed=$PFP_ROOT/scripts/embed_uniprot_descriptions.py"
  --source-file "pfp-if1=$PFP_ROOT/scripts/extract_esm_if1_embeddings.py"
  --source-file "workflow=$HERE/run_contemporary_embedding_reproducibility.sh"
  --source-file "analysis=$HERE/analyze_embedding_reproducibility.py"
  --source-file "runtime-recorder=$HERE/record_embedding_runtime.py"
  --source-file "text-recipe=$HERE/run_pfp_temporal_text.py"
  --source-file "text-cls-reducer=$HERE/reduce_text_embeddings_to_cls.py"
  --source-file "alphafold-prefetch=$HERE/prefetch_alphafold_structures.py"
)
if [[ -n "$IF1_EXTRACT_SCRIPT" ]]; then
  runtime_sources+=(--source-file "framework-if1-compat=$IF1_EXTRACT_SCRIPT")
fi
"$PYTHON_BIN" "$HERE/record_embedding_runtime.py" \
  --output "$OUTPUT_DIR/reports/runtime_hardware.json" \
  "${runtime_sources[@]}" > "$OUTPUT_DIR/logs/runtime_hardware.log"

input_arguments=()
if [[ "$MODALITY" == "text" ]]; then
  echo "==> [5/8] Freeze one text input, then encode it twice"
  "$PYTHON_BIN" "$HERE/run_pfp_temporal_text.py" \
    --pfp-root "$PFP_ROOT" --cafa-assessment-dir "$CAFA_ASSESSMENT_DIR" \
    --cutoff-date "$TEXT_CUTOFF_DATE" --workers "${TEXT_HISTORY_WORKERS:-5}" \
    > "$OUTPUT_DIR/logs/text_input_preparation.log" 2>&1
  TEXT_INPUT="$PFP_ROOT/data/embedding_cache/uniprot_text/protein_descriptions.tsv"
  [[ -s "$TEXT_INPUT" ]] || die "Text input was not produced: $TEXT_INPUT"
  cp -p "$TEXT_INPUT" "$OUTPUT_DIR/inputs/protein_descriptions.tsv"
  input_arguments+=(--input-file "$TEXT_INPUT")

  run_text_repeat() {
    local repeat_name="$1"
    local repeat_root="$2"
    local current_cache="$PFP_ROOT/data/embedding_cache/exp_text_embeddings"
    local temporal_cache="$PFP_ROOT/data/embedding_cache/exp_text_embeddings_temporal"
    rm -rf "$current_cache" "$temporal_cache"
    "$PYTHON_BIN" "$PFP_ROOT/scripts/embed_uniprot_descriptions.py" \
      --data-dir "$PFP_ROOT/data" > "$OUTPUT_DIR/logs/${repeat_name}_text.log" 2>&1
    "$PYTHON_BIN" "$HERE/reduce_text_embeddings_to_cls.py" \
      --directory "$current_cache" \
      --report "$OUTPUT_DIR/reports/${repeat_name}_text_cls_reduction.json" \
      >> "$OUTPUT_DIR/logs/${repeat_name}_text.log" 2>&1
    mkdir -p "$repeat_root"
    mv "$current_cache" "$repeat_root/exp_text_embeddings_temporal"
  }
  run_text_repeat repeat_1 "$REPEAT_ONE_ROOT"
  [[ "$(sha256sum "$TEXT_INPUT" | awk '{print $1}')" == \
      "$(sha256sum "$OUTPUT_DIR/inputs/protein_descriptions.tsv" | awk '{print $1}')" ]] || \
    die "Text input changed during repeat 1"
  run_text_repeat repeat_2 "$REPEAT_TWO_ROOT"
elif [[ "$MODALITY" == "structure" ]]; then
  echo "==> [5/8] Materialize one authenticated PDB view, then encode it twice"
  PDB_CACHE="$STATE_ROOT/source_cache/alphafold_structures"
  PDB_WORKSPACE="$PFP_ROOT/data/alphafold_structures"
  "$PYTHON_BIN" "$HERE/prefetch_alphafold_structures.py" \
    --pfp-root "$PFP_ROOT" --cafa-assessment-dir "$CAFA_ASSESSMENT_DIR" \
    --data-dir "$PFP_ROOT/data" --persistent-cache-dir "$PDB_CACHE" \
    --workspace-pdb-dir "$PDB_WORKSPACE" \
    --coverage-report "$PFP_ROOT/data/alphafold_coverage_results.txt" \
    --report "$OUTPUT_DIR/reports/alphafold_control_inputs.json" \
    --api-workers "${ALPHAFOLD_API_WORKERS:-8}" \
    --download-workers "${ALPHAFOLD_DOWNLOAD_WORKERS:-8}" \
    > "$OUTPUT_DIR/logs/structure_input_preparation.log" 2>&1
  while IFS=$'\t' read -r protein_id _; do
    [[ "$protein_id" != "protein_id" ]] || continue
    pdb="$PDB_WORKSPACE/${protein_id}.pdb"
    [[ -s "$pdb" ]] || die "Missing control PDB: $pdb"
    input_arguments+=(--input-file "$pdb")
  done < "$CONTROLS"

  run_structure_repeat() {
    local repeat_name="$1"
    local repeat_root="$2"
    mkdir -p "$repeat_root/IF1"
    SINGULARITYENV_PYTHONPATH="$IF1_NUMPY_OVERLAY" \
      MMFP_PYTHONPATH="$IF1_NUMPY_OVERLAY" \
      "$PYTHON_BIN" "$IF1_EXTRACT_SCRIPT" \
        --pdb_dir "$PDB_WORKSPACE" --output_dir "$repeat_root/IF1" \
        --pooling mean --device cuda \
        > "$OUTPUT_DIR/logs/${repeat_name}_structure.log" 2>&1
  }
  run_structure_repeat repeat_1 "$REPEAT_ONE_ROOT"
  run_structure_repeat repeat_2 "$REPEAT_TWO_ROOT"
fi

echo "==> [6/8] Confirm the frozen input did not change"
if [[ "$MODALITY" == "text" ]]; then
  [[ "$(sha256sum "$TEXT_INPUT" | awk '{print $1}')" == \
      "$(sha256sum "$OUTPUT_DIR/inputs/protein_descriptions.tsv" | awk '{print $1}')" ]] || \
    die "Text input changed during repeat 2"
fi

echo "==> [7/8] Compare repeat-to-repeat and both repeats to baseline"
"$PYTHON_BIN" "$HERE/analyze_embedding_reproducibility.py" \
  --contract "$STATE_ROOT/contract.json" --controls "$CONTROLS" \
  --modality "$MODALITY" --baseline-root "$REFERENCE_ROOT" \
  --repeat-one-root "$REPEAT_ONE_ROOT" --repeat-two-root "$REPEAT_TWO_ROOT" \
  --minimum-compared "$CONTROL_COUNT" --output-dir "$OUTPUT_DIR/reports" \
  "${input_arguments[@]}" > "$OUTPUT_DIR/logs/analysis.log"

echo "==> [8/8] Publish a non-merging completion marker"
"$PYTHON_BIN" - "$OUTPUT_DIR/DIAGNOSTIC_COMPLETE.json" "$MODALITY" \
  "$CONTROL_COUNT" "$pfp_commit" "$framework_commit" <<'PY'
import json
import sys

path, modality, controls, pfp_commit, framework_commit = sys.argv[1:]
payload = {
    "complete": True,
    "schema_version": 1,
    "diagnostic_only": True,
    "accepted_embedding_state_modified": False,
    "source_cache_writes_allowed": True,
    "modality": modality,
    "control_count": int(controls),
    "pfp_commit": pfp_commit,
    "framework_commit": framework_commit,
}
open(path, "w", encoding="utf-8").write(
    json.dumps(payload, indent=2, sort_keys=True) + "\n"
)
PY

echo "Diagnostic complete: $OUTPUT_DIR"
