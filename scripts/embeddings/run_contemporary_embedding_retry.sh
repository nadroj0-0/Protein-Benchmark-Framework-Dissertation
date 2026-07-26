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
STRICT_FRAMEWORK_COMMIT=0
REFRESH_ALL_TEXT=0
GENERATE_ALL_TEXT_ONLY=0

usage() {
  cat <<'EOF'
Usage: run_contemporary_embedding_retry.sh \
  --pfp-root PATH --work-dir PATH --output-dir PATH \
  --benchmark-dir PATH --plan-dir PATH --state-root PATH \
  --modality sequence|text|structure|ppi \
  [--text-cutoff-date YYYY-MM-DD] [--artifact-catalog PATH] \
  [--strict-framework-commit] [--refresh-all-text|--generate-all-text-only]

Only currently missing pairs for one modality are generated. Accepted control
arrays are materialized from the immutable baseline archive or retry delta into
scratch and must reproduce before new arrays are merged.

`--refresh-all-text` is a separate replacement transaction. It regenerates
text for every target, hydrates accepted non-text arrays into a fresh cache,
installs only the new text layer, and writes a new validated baseline archive.
It never merges into or edits the existing state.

`--generate-all-text-only` regenerates and validates text for every target,
publishes a text-only archive with cutoff provenance, and stops. It performs
no hydration and never merges into or edits the existing state.
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
    --artifact-catalog) ARTIFACT_CATALOG="$2"; export ARTIFACT_CATALOG; shift 2 ;;
    --strict-framework-commit) STRICT_FRAMEWORK_COMMIT=1; shift ;;
    --refresh-all-text) REFRESH_ALL_TEXT=1; shift ;;
    --generate-all-text-only) GENERATE_ALL_TEXT_ONLY=1; shift ;;
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
case "$MODALITY" in sequence|text|structure|ppi) ;; *) die "Invalid modality: $MODALITY" ;; esac
if [[ "$REFRESH_ALL_TEXT" == "1" && "$GENERATE_ALL_TEXT_ONLY" == "1" ]]; then
  die "--refresh-all-text and --generate-all-text-only are mutually exclusive"
fi
if [[ "$REFRESH_ALL_TEXT" == "1" || "$GENERATE_ALL_TEXT_ONLY" == "1" ]]; then
  [[ "$MODALITY" == "text" ]] || die "Full text modes require --modality text"
fi
FULL_TEXT_SELECTION=$((REFRESH_ALL_TEXT || GENERATE_ALL_TEXT_ONLY))
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
export PFP_COMMIT="$pfp_commit"
"$PYTHON_BIN" - "$STATE_ROOT/contract.json" "$pfp_commit" "$framework_commit" \
  "$TEXT_CUTOFF_DATE" "$STRICT_FRAMEWORK_COMMIT" \
  "pfp-prott5=$PFP_ROOT/scripts/extract_prott5_embeddings.py" \
  "pfp-text-extract=$PFP_ROOT/scripts/extract_uniprot_text.py" \
  "pfp-text-embed=$PFP_ROOT/scripts/embed_uniprot_descriptions.py" \
  "pfp-if1=$PFP_ROOT/scripts/extract_esm_if1_embeddings.py" \
  "pfp-ppi=$PFP_ROOT/scripts/extract_ppi_embeddings.py" \
  "framework-if1-compat=$HERE/build_pfp_if1_compat_copy.py" \
  "framework-ppi-compat=$HERE/build_pfp_ppi_compat_copy.py" <<'PY'
import hashlib
import json
import sys
contract = json.load(open(sys.argv[1]))
observed_pfp_commit = sys.argv[2]
observed_framework_commit = sys.argv[3]
strict_framework_commit = sys.argv[5] == "1"
if contract["pfp_commit"] != observed_pfp_commit:
    raise SystemExit(
        f"State contract pfp_commit mismatch: {contract['pfp_commit']} != "
        f"{observed_pfp_commit}"
    )
if contract["framework_commit"] != observed_framework_commit:
    message = (
        "State framework commit differs: "
        f"initialized={contract['framework_commit']} "
        f"retry={observed_framework_commit}"
    )
    if strict_framework_commit:
        raise SystemExit(message)
    print(
        f"WARNING: {message}; continuing because strict framework revision "
        "matching is disabled. Critical source hashes remain enforced.",
        file=sys.stderr,
    )
cutoff = contract.get("runtime", {}).get("text_cutoff_date")
if cutoff != sys.argv[4]:
    raise SystemExit(f"State contract text cutoff mismatch: {cutoff} != {sys.argv[4]}")
sources = {entry["label"]: entry["sha256"] for entry in contract["source_files"]}
for specification in sys.argv[6:]:
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
        raise SystemExit(
            f"State contract source mismatch for {label}: {expected} != {observed}"
        )
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
PFP_ROOT="$PFP_ROOT" EMBEDDING_DEPENDENCY_PROFILE="$MODALITY" \
  bash "$HERE/generate_embeddings_dependencies.sh" \
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

echo "==> [4/9] Select requested pairs and controls"
if [[ "$FULL_TEXT_SELECTION" == "1" ]]; then
  "$PYTHON_BIN" "$HERE/manage_resumable_embedding_state.py" all-pairs \
    --state-root "$STATE_ROOT" --modality text --output "$REQUESTED" \
    > "$OUTPUT_DIR/reports/full_text_selection.json"
  printf 'protein_id\tmodality\n' > "$CONTROLS"
else
  "$PYTHON_BIN" "$HERE/manage_resumable_embedding_state.py" pending \
    --state-root "$STATE_ROOT" --modality "$MODALITY" --output "$REQUESTED" \
    > "$OUTPUT_DIR/reports/pending_selection.json"
fi
requested_count="$(($(wc -l < "$REQUESTED") - 1))"
if [[ "$requested_count" == "0" && "$FULL_TEXT_SELECTION" == "0" ]]; then
  "$PYTHON_BIN" "$HERE/manage_resumable_embedding_state.py" summary \
    --state-root "$STATE_ROOT" --report-dir "$OUTPUT_DIR/reports/embedding_state" \
    > "$OUTPUT_DIR/reports/embedding_state_summary.json"
  printf '{"complete":true,"no_work":true,"modality":"%s"}\n' "$MODALITY" \
    > "$OUTPUT_DIR/RETRY_COMPLETE.json"
  echo "No $MODALITY pairs need retrying"
  exit 0
fi
if [[ "$FULL_TEXT_SELECTION" == "0" ]]; then
  "$PYTHON_BIN" "$HERE/manage_resumable_embedding_state.py" controls \
    --state-root "$STATE_ROOT" --modality "$MODALITY" --count "$CONTROL_COUNT" \
    --output "$CONTROLS" > "$OUTPUT_DIR/reports/control_selection.json"
fi

echo "==> [5/9] Build the exact contemporary retry workspace"
"$PYTHON_BIN" "$HERE/prepare_contemporary_retry_workspace.py" \
  --plan-dir "$PLAN_DIR" --target-benchmark-dir "$BENCHMARK_DIR" \
  --data-dir "$PFP_ROOT/data" --requested-pairs "$REQUESTED" \
  --control-pairs "$CONTROLS" --modality "$MODALITY" \
  --report "$OUTPUT_DIR/reports/retry_workspace.json"
if [[ "$FULL_TEXT_SELECTION" == "0" ]]; then
  "$PYTHON_BIN" "$HERE/manage_resumable_embedding_state.py" materialize \
    --state-root "$STATE_ROOT" --pairs "$CONTROLS" \
    --output-cache-root "$REFERENCE_CONTROLS" \
    --report "$OUTPUT_DIR/reports/control_materialization.json"
fi

export PPI_EXTRACT_SCRIPT="$RUNTIME_COMPAT/extract_ppi_embeddings.py"
export IF1_EXTRACT_SCRIPT="$RUNTIME_COMPAT/extract_esm_if1_embeddings.py"
export IF1_PYTHON_BIN="$PYTHON_BIN"
export IF1_PYTHONPATH="$IF1_NUMPY_OVERLAY"
export TEXT_CUTOFF_DATE
export TEXT_REPORT_DIR="$PFP_ROOT/results/embedding_reports/text"
export HF_HOME="${HF_HOME:-$STATE_ROOT/source_cache/model_weights/huggingface}"
export TORCH_HOME="${TORCH_HOME:-$STATE_ROOT/source_cache/model_weights/torch}"
export ALPHAFOLD_ACQUISITION_MODE=framework-bounded
export ALPHAFOLD_PERSISTENT_CACHE_DIR="$STATE_ROOT/source_cache/alphafold_structures"
export ALPHAFOLD_API_WORKERS="${ALPHAFOLD_API_WORKERS:-8}"
export ALPHAFOLD_DOWNLOAD_WORKERS="${ALPHAFOLD_DOWNLOAD_WORKERS:-8}"
export ALPHAFOLD_PREFETCH_REPORT="$OUTPUT_DIR/reports/alphafold_prefetch_retry.json"
mkdir -p "$HF_HOME" "$TORCH_HOME"

echo "==> [6/9] Generate selected $MODALITY pairs"
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
if [[ "$FULL_TEXT_SELECTION" == "1" ]]; then
  [[ "$generation_status" == "0" ]] || \
    die "Full text generator failed with status $generation_status"
  printf '{"skipped":true,"reason":"old text used an incorrect effective cutoff"}\n' \
    > "$OUTPUT_DIR/reports/subset_equivalence.json"
else
  "$PYTHON_BIN" "$HERE/verify_embedding_subset_equivalence.py" \
    --state-root "$STATE_ROOT" --reference-cache-root "$REFERENCE_CONTROLS" \
    --generated-cache-root "$PFP_ROOT/data/embedding_cache" \
    --control-pairs "$CONTROLS" --modality "$MODALITY" \
    --minimum-compared "$EQUIVALENCE_MINIMUM" \
    --report "$OUTPUT_DIR/reports/subset_equivalence.json"
fi

if [[ "$GENERATE_ALL_TEXT_ONLY" == "1" ]]; then
  echo "==> [8/9] Validate and archive corrected text without hydration"
  replacement_policy="$WORK_DIR/replacement_policy.json"
  text_source="$PFP_ROOT/data/embedding_cache/exp_text_embeddings_temporal"
  [[ -d "$text_source" ]] || die "Corrected temporal text cache is missing"
  "$PYTHON_BIN" - "$STATE_ROOT/contract.json" "$replacement_policy" <<'PY'
import json
import sys
contract = json.load(open(sys.argv[1], encoding="utf-8"))
with open(sys.argv[2], "w", encoding="utf-8") as handle:
    json.dump(contract["policy"], handle, indent=2, sort_keys=True)
    handle.write("\n")
PY
  mkdir -p "$OUTPUT_DIR/artifacts" "$OUTPUT_DIR/reports/text_provenance"
  archive="$OUTPUT_DIR/artifacts/contemporary_text_embeddings_cutoff_${TEXT_CUTOFF_DATE}.tar.gz"
  assembly="$OUTPUT_DIR/artifacts/contemporary_text_embeddings_cutoff_${TEXT_CUTOFF_DATE}_assembly.tsv.gz"
  text_archive_report="$OUTPUT_DIR/reports/text_only_archive.json"
  "$PYTHON_BIN" "$HERE/build_embedding_baseline_archive.py" \
    --generated-cache-root "$PFP_ROOT/data/embedding_cache" \
    --data-dir "$PFP_ROOT/data" --policy "$replacement_policy" \
    --only-modality text --archive "$archive" \
    --assembly-report "$assembly" --report "$text_archive_report"
  "$PYTHON_BIN" - "$text_archive_report" <<'PY'
import json
import sys
report = json.load(open(sys.argv[1], encoding="utf-8"))
if report.get("modalities") != ["text"]:
    raise SystemExit(f"Unexpected archived modalities: {report.get('modalities')}")
if report.get("available_pairs", 0) <= 0:
    raise SystemExit("Corrected text archive contains no embeddings")
PY
  text_run_report="$PFP_ROOT/data/embedding_cache/uniprot_text/temporal_recipe/framework_temporal_text_run.json"
  cp -p "$text_run_report" "$OUTPUT_DIR/reports/text_provenance/"
  "$PYTHON_BIN" - "$text_run_report" "$OUTPUT_DIR/reports/text_provenance" <<'PY'
import json
import shutil
import sys
from pathlib import Path
report = json.load(open(sys.argv[1], encoding="utf-8"))
destination = Path(sys.argv[2])
for key in ("historical_state_contract", "selected_unisave_versions"):
    source = report.get(key)
    if source:
        path = Path(source)
        if path.is_file():
            shutil.copy2(path, destination / path.name)
PY
  [[ ! -f "$TEXT_REPORT_DIR/cls_reduction.json" ]] || \
    cp -p "$TEXT_REPORT_DIR/cls_reduction.json" "$OUTPUT_DIR/reports/text_provenance/"
elif [[ "$REFRESH_ALL_TEXT" == "1" ]]; then
  echo "==> [8/9] Build a fresh cache with corrected text"
  replacement_cache="$WORK_DIR/corrected_combined_cache"
  replacement_policy="$WORK_DIR/replacement_policy.json"
  text_source="$PFP_ROOT/data/embedding_cache/exp_text_embeddings_temporal"
  [[ -d "$text_source" ]] || die "Corrected temporal text cache is missing"
  "$PYTHON_BIN" "$HERE/manage_resumable_embedding_state.py" hydrate \
    --state-root "$STATE_ROOT" --output-cache-root "$replacement_cache" \
    --exclude-modality text --preserve-evidence \
    --report "$OUTPUT_DIR/reports/non_text_hydration.json"
  [[ ! -e "$replacement_cache/exp_text_embeddings_temporal" ]] || \
    die "Old text entered the replacement cache"
  cp -a "$text_source" "$replacement_cache/exp_text_embeddings_temporal"
  "$PYTHON_BIN" - "$STATE_ROOT/contract.json" "$replacement_policy" <<'PY'
import json
import sys
contract = json.load(open(sys.argv[1], encoding="utf-8"))
with open(sys.argv[2], "w", encoding="utf-8") as handle:
    json.dump(contract["policy"], handle, indent=2, sort_keys=True)
    handle.write("\n")
PY
  mkdir -p "$OUTPUT_DIR/artifacts"
  archive="$OUTPUT_DIR/artifacts/contemporary_embeddings_text_cutoff_${TEXT_CUTOFF_DATE}.tar.gz"
  assembly="$OUTPUT_DIR/artifacts/contemporary_embeddings_text_cutoff_${TEXT_CUTOFF_DATE}_assembly.tsv.gz"
  "$PYTHON_BIN" "$HERE/build_embedding_baseline_archive.py" \
    --generated-cache-root "$replacement_cache" --data-dir "$PFP_ROOT/data" \
    --policy "$replacement_policy" --archive "$archive" \
    --assembly-report "$assembly" \
    --report "$OUTPUT_DIR/reports/replacement_baseline.json"
  cp -p \
    "$PFP_ROOT/data/embedding_cache/uniprot_text/temporal_recipe/framework_temporal_text_run.json" \
    "$OUTPUT_DIR/reports/"
else
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
fi

echo "==> [9/9] Publish compact retry status"
printf '%s\n' "$generation_status" > "$OUTPUT_DIR/reports/generator_exit_status.txt"
if [[ "$FULL_TEXT_SELECTION" == "1" ]]; then
  "$PYTHON_BIN" - "$OUTPUT_DIR" "$STATE_ROOT/contract.json" "$TEXT_CUTOFF_DATE" \
    "$framework_commit" "$pfp_commit" "$GENERATE_ALL_TEXT_ONLY" <<'PY'
import hashlib
import json
import sys
from pathlib import Path
root = Path(sys.argv[1])
contract = Path(sys.argv[2])
cutoff = sys.argv[3]
framework_commit = sys.argv[4]
pfp_commit = sys.argv[5]
text_only = sys.argv[6] == "1"
prefix = "contemporary_text_embeddings" if text_only else "contemporary_embeddings_text"
archive_relative = Path("artifacts") / f"{prefix}_cutoff_{cutoff}.tar.gz"
assembly_relative = Path("artifacts") / f"{prefix}_cutoff_{cutoff}_assembly.tsv.gz"
archive = root / archive_relative
assembly = root / assembly_relative
def sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()
run_report = json.loads(
    (root / "reports" / "text_provenance" / "framework_temporal_text_run.json").read_text(
        encoding="utf-8"
    )
) if text_only else json.loads(
    (root / "reports" / "framework_temporal_text_run.json").read_text(encoding="utf-8")
)
if run_report.get("requested_cutoff") != cutoff or run_report.get("effective_cutoff") != cutoff:
    raise SystemExit("Published text cutoff provenance does not match the requested cutoff")
payload = {
    "complete": True,
    "mode": "full-text-generation-only" if text_only else "full-text-replacement",
    "requested_cutoff": cutoff,
    "effective_cutoff": run_report["effective_cutoff"],
    "source_state_contract_sha256": sha(contract),
    "archive": str(archive_relative),
    "archive_sha256": sha(archive),
    "assembly_report": str(assembly_relative),
    "assembly_report_sha256": sha(assembly),
    "framework_commit": framework_commit,
    "pfp_commit": pfp_commit,
    "old_text_carried_forward": False,
    "hydration_performed": not text_only,
    "state_modified": False,
}
if text_only:
    archive_report_relative = Path("reports") / "text_only_archive.json"
    archive_report_path = root / archive_report_relative
    archive_report = json.loads(archive_report_path.read_text(encoding="utf-8"))
    payload.update({
        "text_archive_report": str(archive_report_relative),
        "text_archive_report_sha256": sha(archive_report_path),
        "target_count": archive_report["target_count"],
        "text_available": archive_report["available_pairs"],
        "text_missing": archive_report["missing_pairs"],
    })
(root / ("TEXT_GENERATION_COMPLETE.json" if text_only else "TEXT_REFRESH_COMPLETE.json")).write_text(
    json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
)
PY
elif [[ -f "$STATE_ROOT/EMBEDDING_GATE_PASSED.json" ]]; then
  cp -p "$STATE_ROOT/EMBEDDING_GATE_PASSED.json" "$OUTPUT_DIR/"
else
  cp -p "$STATE_ROOT/GENERATION_INCOMPLETE.json" "$OUTPUT_DIR/"
fi
printf '{"complete":true,"no_work":false,"modality":"%s","generator_exit_status":%s,"refresh_all_text":%s,"generate_all_text_only":%s}\n' \
  "$MODALITY" "$generation_status" "$REFRESH_ALL_TEXT" "$GENERATE_ALL_TEXT_ONLY" \
  > "$OUTPUT_DIR/RETRY_COMPLETE.json"
if [[ "$GENERATE_ALL_TEXT_ONLY" == "1" ]]; then
  echo "Full corrected text-only archive complete. No hydration or state merge was performed."
elif [[ "$REFRESH_ALL_TEXT" == "1" ]]; then
  echo "Full corrected text archive complete. Existing state was not modified."
else
  echo "Retry complete. Valid arrays were retained; missing pairs remain pending."
fi
