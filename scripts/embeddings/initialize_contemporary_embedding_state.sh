#!/usr/bin/env bash
# Initialize an archive-backed, provenance-bound contemporary retry state.

set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
FRAMEWORK_ROOT="$(cd "${HERE}/../.." && pwd)"
# shellcheck source=../reproduction_common.sh
source "$FRAMEWORK_ROOT/scripts/reproduction_common.sh"

PYTHON_BIN="${PYTHON_BIN:-python}"
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
without extracting its hundreds of thousands of arrays onto persistent SAN.
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

mkdir -p "$OUTPUT_DIR" "$STATE_ROOT"
OUTPUT_DIR="$(cd "$OUTPUT_DIR" && pwd)"
STATE_ROOT="$(cd "$STATE_ROOT" && pwd)"
BENCHMARK_DIR="$(cd "$BENCHMARK_DIR" && pwd)"
PLAN_DIR="$(cd "$PLAN_DIR" && pwd)"
BASELINE_ROOT="$(cd "$BASELINE_ROOT" && pwd)"
PFP_ROOT="$(cd "$PFP_ROOT" && pwd)"
POLICY="$(cd "$(dirname "$POLICY")" && pwd)/$(basename "$POLICY")"

validate_mmfp_env "$PYTHON_BIN" > "$OUTPUT_DIR/environment_validation.txt"
pfp_commit="$(git_in_dir "$PFP_ROOT" rev-parse HEAD)"
framework_commit="${FRAMEWORK_COMMIT:-$(git_in_dir "$FRAMEWORK_ROOT" rev-parse HEAD)}"

command=(
  "$PYTHON_BIN" "$HERE/manage_resumable_embedding_state.py" initialize
  --state-root "$STATE_ROOT"
  --benchmark-id contemporary-2025_01-to-2026_02-supervisor
  --benchmark-dir "$BENCHMARK_DIR"
  --target-table "$REUSE_TABLE"
  --target-table "$REGENERATE_TABLE"
  --policy "$POLICY"
  --pfp-commit "$pfp_commit"
  --framework-commit "$framework_commit"
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
  --runtime-value "text_cutoff_date=$TEXT_CUTOFF_DATE"
  --runtime-value "temporal_profile=supervisor"
  --runtime-value "t1_endpoint_policy=snapshot-membership"
  --runtime-value "exclude_t1_backfill=false"
  --runtime-value "alphafold_acquisition=framework-bounded"
  --runtime-value "alphafold_api_workers=8"
  --runtime-value "alphafold_download_workers=8"
)

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
