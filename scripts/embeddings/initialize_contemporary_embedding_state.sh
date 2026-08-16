#!/usr/bin/env bash
# Initialize an archive-backed, provenance-bound contemporary retry state.

set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
FRAMEWORK_ROOT="$(cd "${HERE}/../.." && pwd)"
# shellcheck source=../reproduction_common.sh
source "$FRAMEWORK_ROOT/scripts/reproduction_common.sh"

PYTHON_BIN="${PYTHON_BIN:-python3}"
BENCHMARK_DIR=""
PLAN_DIR=""
BASELINE_ROOT=""
STATE_ROOT=""
PFP_ROOT=""
OUTPUT_DIR=""
TEXT_CUTOFF_DATE="2025-03-08"
POLICY="$FRAMEWORK_ROOT/configs/contemporary_embedding_resume.json"

usage() {
  cat <<'EOF'
Usage: initialize_contemporary_embedding_state.sh \
  --benchmark-dir PATH \
  --plan-dir PATH \
  --baseline-root PATH \
  --state-root PATH \
  --pfp-root PATH \
  --output-dir PATH \
  [--text-cutoff-date YYYY-MM-DD] \
  [--policy PATH]

The baseline root must contain archive/contemporary_embedding_cache.tar.gz and
reports/assembly/embedding_assembly.tsv.gz. The command indexes the archive
without extracting its hundreds of thousands of arrays onto persistent storage.
EOF
}

die() { echo "ERROR: $*" >&2; exit 2; }
git_in_dir() { local directory="$1"; shift; (cd "$directory" && git "$@"); }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --benchmark-dir) BENCHMARK_DIR="$2"; shift 2 ;;
    --plan-dir) PLAN_DIR="$2"; shift 2 ;;
    --baseline-root) BASELINE_ROOT="$2"; shift 2 ;;
    --state-root) STATE_ROOT="$2"; shift 2 ;;
    --pfp-root) PFP_ROOT="$2"; shift 2 ;;
    --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
    --text-cutoff-date) TEXT_CUTOFF_DATE="$2"; shift 2 ;;
    --policy) POLICY="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) usage >&2; die "Unknown argument: $1" ;;
  esac
done

[[ -d "$BENCHMARK_DIR" ]] || die "Missing benchmark directory: $BENCHMARK_DIR"
[[ -d "$PLAN_DIR" ]] || die "Missing reuse plan: $PLAN_DIR"
[[ -d "$BASELINE_ROOT" ]] || die "Missing baseline root: $BASELINE_ROOT"
[[ -n "$STATE_ROOT" ]] || die "--state-root is required"
[[ -d "$PFP_ROOT/.git" ]] || die "PFP root is not a Git checkout: $PFP_ROOT"
[[ -n "$OUTPUT_DIR" ]] || die "--output-dir is required"
[[ "$TEXT_CUTOFF_DATE" =~ ^[0-9]{4}-[0-9]{2}-[0-9]{2}$ ]] || \
  die "Invalid text cutoff date: $TEXT_CUTOFF_DATE"
[[ -f "$POLICY" ]] || die "Missing policy: $POLICY"
command -v "$PYTHON_BIN" >/dev/null 2>&1 || die "Python not found: $PYTHON_BIN"

BASELINE_ARCHIVE="$BASELINE_ROOT/archive/contemporary_embedding_cache.tar.gz"
BASELINE_REPORT="$BASELINE_ROOT/reports/assembly/embedding_assembly.tsv.gz"
BUILD_MANIFEST="$BENCHMARK_DIR/build_manifest.json"
REUSE_TABLE="$PLAN_DIR/reuse_proteins.tsv"
REGENERATE_TABLE="$PLAN_DIR/regenerate_proteins.tsv"
for path in \
  "$BASELINE_ARCHIVE" "$BASELINE_REPORT" "$BUILD_MANIFEST" \
  "$REUSE_TABLE" "$REGENERATE_TABLE" \
  "$PLAN_DIR/run_manifest.json" "$PLAN_DIR/output_manifest.json" \
  "$BASELINE_ROOT/reports/input_acquisition.tsv"; do
  [[ -f "$path" ]] || die "Missing required input: $path"
done

expected_pfp_commit="1e04fd6d6d3c40458fd41ec1a881ed6e24de768e"
pfp_commit="$(git_in_dir "$PFP_ROOT" rev-parse HEAD)"
[[ "$pfp_commit" == "$expected_pfp_commit" ]] || \
  die "PFP commit mismatch: expected $expected_pfp_commit, found $pfp_commit"

CACHE_ROLE="$("$PYTHON_BIN" - "$BASELINE_ROOT" "$BASELINE_ARCHIVE" \
  "$BASELINE_REPORT" "$TEXT_CUTOFF_DATE" "$REUSE_TABLE" \
  "$REGENERATE_TABLE" "$pfp_commit" <<'PY'
import csv
import hashlib
import json
import sys
from pathlib import Path

root, archive, assembly = map(Path, sys.argv[1:4])
cutoff = sys.argv[4]
target_tables = [Path(sys.argv[5]), Path(sys.argv[6])]
expected_pfp_commit = sys.argv[7]
composition_path = root / "COMPOSITION_COMPLETE.json"
text_path = root / "provenance" / "TEXT_GENERATION_COMPLETE.json"
present = tuple(
    path.exists() or path.is_symlink() for path in (composition_path, text_path)
)
if present == (False, False):
    print("composition-base-only")
    raise SystemExit(0)
if (
    present != (True, True)
    or composition_path.is_symlink()
    or text_path.is_symlink()
    or not composition_path.is_file()
    or not text_path.is_file()
):
    raise SystemExit("Corrected baseline has incomplete or unsafe composition evidence")

def load(path):
    with path.open(encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise SystemExit(f"Invalid evidence object: {path}")
    return value

def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()

def require_sha256(value, label):
    if not isinstance(value, str) or len(value) != 64 or any(
        character not in "0123456789abcdef" for character in value
    ):
        raise SystemExit(f"Invalid SHA-256 for {label}")
    return value

targets = {}
for table in target_tables:
    with table.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required_columns = {"protein_id", "sequence_sha256"}
        if reader.fieldnames is None or not required_columns.issubset(reader.fieldnames):
            raise SystemExit(f"Target table lacks required columns: {table}")
        for row in reader:
            protein_id = row["protein_id"]
            sequence_sha256 = row["sequence_sha256"]
            if not protein_id or protein_id in targets:
                raise SystemExit(f"Invalid or repeated target: {protein_id!r}")
            targets[protein_id] = require_sha256(
                sequence_sha256, f"target sequence {protein_id}"
            )
if not targets:
    raise SystemExit("Planner population is empty")
target_manifest = "protein_id\tsequence_sha256\n" + "".join(
    f"{protein_id}\t{targets[protein_id]}\n" for protein_id in sorted(targets)
)
target_manifest_sha256 = hashlib.sha256(target_manifest.encode("utf-8")).hexdigest()

composition = load(composition_path)
text = load(text_path)
required_composition = {
    "complete": True,
    "operation": "replace-finalized-text-layer",
    "expected_cutoff": cutoff,
    "old_text_carried_forward": False,
    "roundtrip_validated": True,
    "combined_archive": "archive/contemporary_embedding_cache.tar.gz",
}
required_text = {
    "complete": True,
    "mode": "full-text-generation-only",
    "requested_cutoff": cutoff,
    "effective_cutoff": cutoff,
    "old_text_carried_forward": False,
    "hydration_performed": False,
    "state_modified": False,
}
for key, expected in required_composition.items():
    observed = composition.get(key)
    if (type(expected) is bool and observed is not expected) or (
        type(expected) is not bool and observed != expected
    ):
        raise SystemExit(f"Composition evidence mismatch for {key}")
for key, expected in required_text.items():
    observed = text.get(key)
    if (type(expected) is bool and observed is not expected) or (
        type(expected) is not bool and observed != expected
    ):
        raise SystemExit(f"Text-generation evidence mismatch for {key}")
if composition.get("combined_archive_sha256") != sha256(archive):
    raise SystemExit("Composition evidence does not authenticate the baseline archive")
if composition.get("assembly_report_sha256") != sha256(assembly):
    raise SystemExit("Composition evidence does not authenticate the assembly report")
composition_replacement_sha = require_sha256(
    composition.get("replacement_archive_sha256"), "composition replacement archive"
)
text_archive_sha = require_sha256(text.get("archive_sha256"), "text archive")
if composition_replacement_sha != text_archive_sha:
    raise SystemExit("Composition and text-generation evidence bind different text archives")
for label, value in (
    ("composition", composition.get("target_count")),
    ("text generation", text.get("target_count")),
):
    if type(value) is not int or value != len(targets):
        raise SystemExit(f"{label} target count differs from the planner population")

source_contract_path = root / "provenance" / "source_state_contract.json"
if not source_contract_path.is_file() or source_contract_path.is_symlink():
    raise SystemExit("Corrected baseline lacks a safe source-state contract")
source_contract_sha = require_sha256(
    text.get("source_state_contract_sha256"), "source-state contract"
)
if source_contract_sha != sha256(source_contract_path):
    raise SystemExit("Text-generation evidence does not authenticate the source-state contract")
source_contract = load(source_contract_path)
recorded_contract_sha = require_sha256(
    source_contract.pop("contract_sha256", None), "source-state contract self-hash"
)
canonical_contract = json.dumps(
    source_contract, sort_keys=True, separators=(",", ":"), ensure_ascii=True
).encode("utf-8")
if recorded_contract_sha != hashlib.sha256(canonical_contract).hexdigest():
    raise SystemExit("Source-state contract self-hash is invalid")
if source_contract.get("pfp_commit") != expected_pfp_commit:
    raise SystemExit("Source-state contract has the wrong PFP commit")
source_targets = source_contract.get("targets", {})
if (
    not isinstance(source_targets, dict)
    or type(source_targets.get("count")) is not int
    or source_targets.get("count") != len(targets)
    or source_targets.get("manifest_sha256") != target_manifest_sha256
):
    raise SystemExit("Source-state contract targets differ from the planner population")
if text.get("pfp_commit") != expected_pfp_commit:
    raise SystemExit("Text-generation evidence has the wrong PFP commit")
source_runtime = source_contract.get("runtime", {})
if source_runtime.get("cache_role") != "composition-base-only":
    raise SystemExit("Text-generation source was not the composition base")
if source_runtime.get("text_generation_cutoff") != cutoff:
    raise SystemExit("Source-state contract has the wrong text-generation cutoff")
print(f"accepted-corrected-{cutoff}")
PY
)"

mkdir -p "$OUTPUT_DIR" "$STATE_ROOT"
OUTPUT_DIR="$(cd "$OUTPUT_DIR" && pwd)"
STATE_ROOT="$(cd "$STATE_ROOT" && pwd)"
BENCHMARK_DIR="$(cd "$BENCHMARK_DIR" && pwd)"
PLAN_DIR="$(cd "$PLAN_DIR" && pwd)"
BASELINE_ROOT="$(cd "$BASELINE_ROOT" && pwd)"
PFP_ROOT="$(cd "$PFP_ROOT" && pwd)"
POLICY="$(cd "$(dirname "$POLICY")" && pwd)/$(basename "$POLICY")"

validate_mmfp_env "$PYTHON_BIN" > "$OUTPUT_DIR/environment_validation.txt"

command=(
  "$PYTHON_BIN" "$HERE/manage_resumable_embedding_state.py" initialize
  --state-root "$STATE_ROOT"
  --benchmark-id contemporary-2025-01-to-2026-02-supervisor
  --benchmark-dir "$BENCHMARK_DIR"
  --target-table "$REUSE_TABLE"
  --target-table "$REGENERATE_TABLE"
  --policy "$POLICY"
  --pfp-commit "$pfp_commit"
  --environment-report "$OUTPUT_DIR/environment_validation.txt"
  --baseline-archive "$BASELINE_ARCHIVE"
  --baseline-assembly-report "$BASELINE_REPORT"
  --source-file "benchmark-build-manifest=$BUILD_MANIFEST"
  --source-file "reuse-plan-manifest=$PLAN_DIR/run_manifest.json"
  --source-file "reuse-plan-output=$PLAN_DIR/output_manifest.json"
  --source-file "initial-acquisition=$BASELINE_ROOT/reports/input_acquisition.tsv"
  --source-file "pfp-prott5=$PFP_ROOT/scripts/extract_prott5_embeddings.py"
  --source-file "pfp-text-extract=$PFP_ROOT/scripts/extract_uniprot_text.py"
  --source-file "pfp-text-embed=$PFP_ROOT/scripts/embed_uniprot_descriptions.py"
  --source-file "pfp-if1=$PFP_ROOT/scripts/extract_esm_if1_embeddings.py"
  --source-file "pfp-ppi=$PFP_ROOT/scripts/extract_ppi_embeddings.py"
  --source-file "framework-if1-compat=$HERE/build_pfp_if1_compat_copy.py"
  --source-file "framework-ppi-compat=$HERE/build_pfp_ppi_compat_copy.py"
  --runtime-value "cache_role=$CACHE_ROLE"
  --runtime-value "text_generation_cutoff=$TEXT_CUTOFF_DATE"
  --runtime-value "temporal_profile=supervisor"
  --runtime-value "t1_endpoint_policy=snapshot-membership"
  --runtime-value "exclude_t1_backfill=false"
  --runtime-value "alphafold_acquisition=framework-bounded"
  --runtime-value "alphafold_api_workers=8"
  --runtime-value "alphafold_download_workers=8"
)
if [[ -f "$BASELINE_ROOT/COMPOSITION_COMPLETE.json" ]]; then
  command+=(
    --source-file "replacement-composition=$BASELINE_ROOT/COMPOSITION_COMPLETE.json"
    --runtime-value "text_cutoff_date=$TEXT_CUTOFF_DATE"
  )
fi
if [[ -f "$BASELINE_ROOT/provenance/TEXT_GENERATION_COMPLETE.json" ]]; then
  command+=(
    --source-file "corrected-text-generation=$BASELINE_ROOT/provenance/TEXT_GENERATION_COMPLETE.json"
  )
fi

printf 'Command:'; printf ' %q' "${command[@]}"; printf '\n'
"${command[@]}" > "$OUTPUT_DIR/initialization_summary.json"
cp -p \
  "$STATE_ROOT/contract.json" \
  "$STATE_ROOT/coverage.json" \
  "$STATE_ROOT/baseline_validation.json" \
  "$STATE_ROOT/needs_retry.tsv" \
  "$OUTPUT_DIR/"

printf '{"complete":true,"state_root":"%s"}\n' "$STATE_ROOT" \
  > "$OUTPUT_DIR/INITIALIZATION_COMPLETE.json"
echo "Contemporary embedding state initialized: $STATE_ROOT"
