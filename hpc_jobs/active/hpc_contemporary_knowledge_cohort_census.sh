#!/usr/bin/env bash
# Reconstruct t0/t1 knowledge cohorts for the accepted contemporary benchmark.

#$ -S /bin/bash
#$ -l tmem=32G
#$ -l scratch0free=20G
#$ -l tscratch=20G
#$ -l h_rt=36:0:0
#$ -pe smp 2
#$ -j y
#$ -N ct25_cohort
#$ -V
#$ -notify

set -Eeuo pipefail
export PYTHONDONTWRITEBYTECODE=1

usage() {
  cat <<'EOF'
Usage:
  qsub hpc_jobs/active/hpc_contemporary_knowledge_cohort_census.sh \
    --output-dir /SAN/.../knowledge_cohort_census

This CPU-only job rebuilds direct and propagated t0/t1 annotation states from
the frozen sources, verifies the accepted test labels exactly, and publishes
the NK/LK/PK-style cohort census plus a reusable temporal annotation ledger.
EOF
}

die() { echo "ERROR: $*" >&2; exit 2; }
require_value() { [[ $# -ge 2 && -n "$2" ]] || die "$1 requires a value"; }
git_in_dir() {
  local directory="$1"
  shift
  (cd "$directory" && git "$@")
}

OUTPUT_DIR=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --output-dir) require_value "$@"; OUTPUT_DIR="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

[[ "$OUTPUT_DIR" == /SAN/* ]] || die "--output-dir must be an absolute SAN path"
[[ ! -e "$OUTPUT_DIR" ]] || die "Output directory already exists: $OUTPUT_DIR"

BENCHMARK_DIR="${BENCHMARK_DIR:-/SAN/bioinf/bmpfp/benchmarks/contemporary/2025_01_to_2026_02_supervisor}"
T0_SPROT_ARCHIVE="${T0_SPROT_ARCHIVE:-/SAN/bioinf/bmpfp/frozen_inputs/uniprot/2025_01/uniprot_sprot-only2025_01.tar.gz}"
T0_TREMBL="${T0_TREMBL:-/SAN/bioinf/bmpfp/derived_inputs/uniprot/cafa3_target_taxa/2025_01/uniprot_trembl_cafa3_targets.dat.gz}"
T1_SPROT="${T1_SPROT:-/SAN/bioinf/bmpfp/frozen_inputs/uniprot/2026_02/uniprot_sprot.dat.gz}"
T1_TREMBL="${T1_TREMBL:-/SAN/bioinf/bmpfp/derived_inputs/uniprot/cafa3_target_taxa/2026_02/uniprot_trembl_cafa3_targets.dat.gz}"
GOA_T0="${GOA_T0:-/SAN/bioinf/bmpfp/frozen_inputs/goa/225/goa_uniprot_all.gaf.225.gz}"
GOA_T1="${GOA_T1:-/SAN/bioinf/bmpfp/frozen_inputs/goa/234/goa_uniprot_all.gaf.234.gz}"
BENCHMARK_OBO="${BENCHMARK_OBO:-/SAN/bioinf/bmpfp/frozen_inputs/ontology/2025-02-06/go-basic.obo}"
T0_SOURCE_OBO="${T0_SOURCE_OBO:-/SAN/bioinf/bmpfp/frozen_inputs/ontology/2025-03-16/go-basic.obo}"
T1_SOURCE_OBO="${T1_SOURCE_OBO:-/SAN/bioinf/bmpfp/frozen_inputs/ontology/2026-06-19/go-basic.obo}"

for path in \
  "$BENCHMARK_DIR/build_manifest.json" \
  "$T0_SPROT_ARCHIVE" "$T0_TREMBL" "$T1_SPROT" "$T1_TREMBL" \
  "$GOA_T0" "$GOA_T1" "$BENCHMARK_OBO" "$T0_SOURCE_OBO" "$T1_SOURCE_OBO"; do
  [[ -s "$path" ]] || die "Required source is missing or empty: $path"
done

JOB_TOKEN="${JOB_ID:-manual_$$}"
WORK="/scratch0/ct25_cohort_${JOB_TOKEN}"
FRAMEWORK_DIR="$WORK/Protein-Benchmark-Framework-Dissertation"
INPUT_DIR="$WORK/inputs"
ANALYSIS_OUTPUT="$WORK/analysis"
LOG_FILE="$WORK/cohort_census.log"
PUBLISH_STAGE="${OUTPUT_DIR}.staging-${JOB_TOKEN}"
PUBLISH_LOCK="${OUTPUT_DIR}.publish-lock"
FRAMEWORK_REPO_URL="${FRAMEWORK_REPO_URL:-https://github.com/nadroj0-0/Protein-Benchmark-Framework-Dissertation.git}"
FRAMEWORK_COMMIT="${FRAMEWORK_COMMIT:-}"
LOCK_HELD=0

cleanup() {
  local status=$?
  trap - EXIT
  set +e
  if [[ "$LOCK_HELD" == "1" && -d "$PUBLISH_LOCK" && ! -L "$PUBLISH_LOCK" ]]; then
    rmdir -- "$PUBLISH_LOCK"
  fi
  [[ ! -d "$PUBLISH_STAGE" || -L "$PUBLISH_STAGE" ]] || rm -rf -- "$PUBLISH_STAGE"
  [[ ! -d "$WORK" || -L "$WORK" || "$WORK" != /scratch0/ct25_cohort_* ]] || \
    rm -rf -- "$WORK"
  exit "$status"
}
trap cleanup EXIT
trap 'echo "Received termination signal"; exit 130' INT TERM

[[ ! -e "$WORK" ]] || die "Scratch path already exists: $WORK"
[[ ! -e "$PUBLISH_STAGE" ]] || die "Publication stage already exists"
[[ ! -e "$PUBLISH_LOCK" ]] || die "Publication lock already exists"
mkdir -p "$INPUT_DIR" "$(dirname "$OUTPUT_DIR")"

echo "Host             : $(hostname)"
echo "Job ID           : ${JOB_ID:-manual}"
echo "Analysis         : contemporary temporal knowledge-cohort census"
echo "Benchmark        : $BENCHMARK_DIR"
echo "SAN output       : $OUTPUT_DIR"
echo "Started          : $(date -Is)"

git clone --no-checkout "$FRAMEWORK_REPO_URL" "$FRAMEWORK_DIR"
if [[ -n "$FRAMEWORK_COMMIT" ]]; then
  [[ "$FRAMEWORK_COMMIT" =~ ^[0-9a-fA-F]{40}$ ]] || \
    die "FRAMEWORK_COMMIT must be a full commit when supplied"
  git_in_dir "$FRAMEWORK_DIR" checkout --detach "$FRAMEWORK_COMMIT"
else
  git_in_dir "$FRAMEWORK_DIR" checkout main
  FRAMEWORK_COMMIT="$(git_in_dir "$FRAMEWORK_DIR" rev-parse HEAD)"
fi
echo "Framework commit : $FRAMEWORK_COMMIT"

cd "$FRAMEWORK_DIR"
source scripts/reproduction_common.sh
load_framework_paths "$FRAMEWORK_DIR"
add_mmfp_singularity_bind "$WORK"
add_mmfp_singularity_bind /SAN/bioinf/bmpfp
activate_or_create_mmfp_env
PYTHON_BIN="$(command -v python)"

echo "==> Extracting the frozen t0 Swiss-Prot DAT"
tar -xzf "$T0_SPROT_ARCHIVE" -C "$INPUT_DIR" uniprot_sprot.dat.gz
[[ -s "$INPUT_DIR/uniprot_sprot.dat.gz" ]] || die "t0 Swiss-Prot extraction failed"

echo "==> Reconstructing direct/closure states and cohort census"
"$PYTHON_BIN" scripts/diagnostics/build_contemporary_knowledge_cohort_census.py \
  --accepted-benchmark-dir "$BENCHMARK_DIR" \
  --t0-sprot "$INPUT_DIR/uniprot_sprot.dat.gz" \
  --t0-trembl "$T0_TREMBL" \
  --t1-sprot "$T1_SPROT" \
  --t1-trembl "$T1_TREMBL" \
  --goa-t0 "$GOA_T0" \
  --goa-t1 "$GOA_T1" \
  --benchmark-obo "$BENCHMARK_OBO" \
  --t0-source-obo "$T0_SOURCE_OBO" \
  --t1-source-obo "$T1_SOURCE_OBO" \
  --output-dir "$ANALYSIS_OUTPUT" \
  >"$LOG_FILE" 2>&1

# The census manifest intentionally excludes its RUN_COMPLETE marker to avoid a
# circular hash. The outer publication manifest below binds both nested control
# files and is the authoritative transport-integrity check.

echo "==> Publishing cohort census atomically"
mkdir -p "$PUBLISH_STAGE/logs"
cp -a "$ANALYSIS_OUTPUT" "$PUBLISH_STAGE/analysis"
cp -p "$LOG_FILE" "$PUBLISH_STAGE/logs/"
"$PYTHON_BIN" scripts/model_execution/manage_output_manifest.py write \
  --root "$PUBLISH_STAGE" --include-nested-control-files
MANIFEST_SHA256="$("$PYTHON_BIN" -c 'import hashlib,pathlib,sys; print(hashlib.sha256(pathlib.Path(sys.argv[1]).read_bytes()).hexdigest())' "$PUBLISH_STAGE/output_manifest.json")"
"$PYTHON_BIN" -c 'import json,pathlib,sys; pathlib.Path(sys.argv[1]).write_text(json.dumps({"complete":True,"analysis_kind":"contemporary_knowledge_cohort_census","framework_commit":sys.argv[2],"manifest":"output_manifest.json","manifest_sha256":sys.argv[3]},indent=2)+"\n")' \
  "$PUBLISH_STAGE/WORKFLOW_COMPLETE.json" "$FRAMEWORK_COMMIT" "$MANIFEST_SHA256"
"$PYTHON_BIN" scripts/model_execution/manage_output_manifest.py verify \
  --root "$PUBLISH_STAGE" --include-nested-control-files
mkdir "$PUBLISH_LOCK"
LOCK_HELD=1
[[ ! -e "$OUTPUT_DIR" ]] || die "Output directory appeared during publication"
mv -T "$PUBLISH_STAGE" "$OUTPUT_DIR"
rmdir "$PUBLISH_LOCK"
LOCK_HELD=0
echo "Published cohort census: $OUTPUT_DIR"
