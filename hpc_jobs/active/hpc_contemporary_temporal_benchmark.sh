#!/usr/bin/env bash
# UCL Grid Engine wrapper for the contemporary temporal benchmark builder.

#$ -l tmem=64G
#$ -l tscratch=200G
#$ -l scratch0free=200G
#$ -l h_rt=72:0:0
#$ -j y
#$ -N contemporary_cafa
#$ -V

set -euo pipefail

WORK="/scratch0/contemporary_cafa_${JOB_ID:-manual}"
SOURCE_DB_ROOT="${PROTEIN_DATABASE_ROOT:-$HOME/protein_databases}"
RESULTS_ROOT="${RESULTS_ROOT:-$HOME/contemporary_cafa_benchmark_results}"
FRAMEWORK_REPO_URL="${FRAMEWORK_REPO_URL:-https://github.com/nadroj0-0/Protein-Benchmark-Framework-Dissertation.git}"
FRAMEWORK_DIR="$WORK/Protein-Benchmark-Framework-Dissertation"
STAGED_DB_ROOT="$WORK/protein_databases"
RUN_TAG="${JOB_ID:-manual}_$(date +%Y%m%d_%H%M%S)"
SCRATCH_RUN_ROOT="$WORK/run"
FINAL_RUN_ROOT="$RESULTS_ROOT/$RUN_TAG"
COPIED_RESULTS=0

copy_results() {
    if [[ "$COPIED_RESULTS" == "0" && -d "$SCRATCH_RUN_ROOT" ]]; then
        echo "Copying outputs and reports to $FINAL_RUN_ROOT"
        mkdir -p "$FINAL_RUN_ROOT"
        cp -a "$SCRATCH_RUN_ROOT/." "$FINAL_RUN_ROOT/"
        COPIED_RESULTS=1
    fi
}

cleanup() {
    local status=$?
    set +e
    copy_results
    echo "Cleaning scratch directory: $WORK"
    cd "$HOME"
    rm -rf "$WORK"
    exit "$status"
}
trap cleanup EXIT
trap 'echo "Received kill signal"; exit 130' SIGINT SIGTERM

stage_file() {
    local relative="$1"
    local source="$SOURCE_DB_ROOT/$relative"
    local destination="$STAGED_DB_ROOT/$relative"
    if [[ ! -f "$source" ]]; then
        return 1
    fi
    mkdir -p "$(dirname "$destination")"
    cp -a "$source" "$destination"
}

echo "Host        : $(hostname)"
echo "Job ID      : ${JOB_ID:-manual}"
echo "Profile     : ${PROFILE:-contemporary-cafa3-style}"
echo "Scratch     : $WORK"
echo "Final output: $FINAL_RUN_ROOT"

mkdir -p "$WORK" "$RESULTS_ROOT" "$STAGED_DB_ROOT"

echo "Cloning dissertation framework into scratch"
git clone "$FRAMEWORK_REPO_URL" "$FRAMEWORK_DIR"
cd "$FRAMEWORK_DIR"
source scripts/reproduction_common.sh
load_framework_paths "$FRAMEWORK_DIR"
activate_or_create_mmfp_env

echo "Staging frozen database inputs"
stage_file "uniprot/release_2025_01/uniprot_sprot-only2025_01.tar.gz"
stage_file "uniprot/release_2026_02/uniprot_sprot.dat.gz"
if ! stage_file "uniprot/release_2025_01/uniprot_trembl-only2025_01.tar.gz"; then
    if [[ "${ALLOW_SPROT_ONLY:-0}" != "1" ]]; then
        echo "Missing t0 TrEMBL archive under $SOURCE_DB_ROOT" >&2
        exit 1
    fi
fi
if ! stage_file "uniprot/release_2026_02/uniprot_trembl.dat.gz"; then
    if [[ "${ALLOW_SPROT_ONLY:-0}" != "1" ]]; then
        echo "Missing t1 TrEMBL DAT under $SOURCE_DB_ROOT" >&2
        exit 1
    fi
fi
stage_file "goa/release_2025_01/goa_uniprot_all.gaf.225.gz"
stage_file "goa/release_2026_02/goa_uniprot_all.gaf.234.gz"
stage_file "ontology/release_2025-03-07/go-basic.obo"
stage_file "ontology/release_2026-06-15/go-basic.obo"

export DB_ROOT="$STAGED_DB_ROOT"
export RUN_ROOT="$SCRATCH_RUN_ROOT"
export WORK_DIR="$WORK/extracted"
export PYTHON_BIN=python
export REMOVE_ARCHIVES_AFTER_EXTRACT=1
export CAFA_BUILDER_USE_PIGZ=1

bash scripts/benchmark_generation/run_contemporary_temporal_benchmark.sh
copy_results

echo "Finished successfully: $(date)"
