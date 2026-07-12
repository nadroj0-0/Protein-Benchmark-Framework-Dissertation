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

stage_if_present() {
    local relative="$1"
    if stage_file "$relative"; then
        printf 'staged\t%s\n' "$relative" >> "$SCRATCH_RUN_ROOT/logs/local_staging.tsv"
    else
        printf 'download-required\t%s\n' "$relative" >> "$SCRATCH_RUN_ROOT/logs/local_staging.tsv"
    fi
}

stage_input_list() {
    local variable_name="$1"
    local label="$2"
    local value="${!variable_name:-}"
    local staged=()
    local source
    local index=0
    IFS=':' read -r -a sources <<< "$value"
    for source in "${sources[@]}"; do
        if [[ ! -f "$source" ]]; then
            echo "Explicit $variable_name input does not exist: $source" >&2
            exit 1
        fi
        local destination="$STAGED_DB_ROOT/custom/$label/${index}_$(basename "$source")"
        mkdir -p "$(dirname "$destination")"
        cp -a "$source" "$destination"
        staged+=("$destination")
        printf 'staged-custom\t%s\n' "$source" >> "$SCRATCH_RUN_ROOT/logs/local_staging.tsv"
        index=$((index + 1))
    done
    local joined
    joined="$(IFS=:; echo "${staged[*]}")"
    printf -v "$variable_name" '%s' "$joined"
    export "$variable_name"
}

echo "Host        : $(hostname)"
echo "Job ID      : ${JOB_ID:-manual}"
echo "Profile     : ${PROFILE:-contemporary-cafa3-style}"
echo "Scratch     : $WORK"
echo "Final output: $FINAL_RUN_ROOT"

mkdir -p "$WORK" "$RESULTS_ROOT" "$STAGED_DB_ROOT" "$SCRATCH_RUN_ROOT/logs"
printf 'action\tpath\n' > "$SCRATCH_RUN_ROOT/logs/local_staging.tsv"

echo "Cloning dissertation framework into scratch"
git clone "$FRAMEWORK_REPO_URL" "$FRAMEWORK_DIR"
cd "$FRAMEWORK_DIR"
source scripts/reproduction_common.sh
load_framework_paths "$FRAMEWORK_DIR"
activate_or_create_mmfp_env
export PYTHON_BIN="$(command -v python)"

echo "Staging frozen database inputs"
if [[ -n "${UNIPROT_T0_INPUTS:-}" || -n "${UNIPROT_T1_INPUTS:-}" ]]; then
    if [[ -z "${UNIPROT_T0_INPUTS:-}" || -z "${UNIPROT_T1_INPUTS:-}" ]]; then
        echo "UNIPROT_T0_INPUTS and UNIPROT_T1_INPUTS must be supplied together" >&2
        exit 1
    fi
    stage_input_list UNIPROT_T0_INPUTS t0
    stage_input_list UNIPROT_T1_INPUTS t1
else
    stage_if_present "uniprot/release_2025_01/uniprot_sprot.dat.gz"
    stage_if_present "uniprot/release_2025_01/uniprot_sprot-only2025_01.tar.gz"
    stage_if_present "uniprot/release_2025_01/uniprot_trembl_cafa3_targets.dat.gz"
    stage_if_present "uniprot/release_2026_02/uniprot_sprot.dat.gz"
    stage_if_present "uniprot/release_2026_02/uniprot_trembl_cafa3_targets.dat.gz"
fi
stage_if_present "goa/release_2025_01/goa_uniprot_all.gaf.225.gz"
stage_if_present "goa/release_2026_02/goa_uniprot_all.gaf.234.gz"
stage_if_present "ontology/release_2025-02-06/go-basic.obo"
stage_if_present "ontology/release_2025-03-16/go-basic.obo"
stage_if_present "ontology/release_2026-06-19/go-basic.obo"

# Full TrEMBL sources are 100+ GB. If present locally, stream-filter them in
# place rather than duplicating them into 200 GB scratch.
if [[ -f "$SOURCE_DB_ROOT/uniprot/release_2025_01/uniprot_trembl.dat.gz" ]]; then
    export T0_TREMBL_DAT_SOURCE="$SOURCE_DB_ROOT/uniprot/release_2025_01/uniprot_trembl.dat.gz"
elif [[ -f "$SOURCE_DB_ROOT/uniprot/release_2025_01/knowledgebase2025_01.tar.gz" ]]; then
    export T0_TREMBL_ARCHIVE_SOURCE="$SOURCE_DB_ROOT/uniprot/release_2025_01/knowledgebase2025_01.tar.gz"
fi
if [[ -f "$SOURCE_DB_ROOT/uniprot/release_2026_02/uniprot_trembl.dat.gz" ]]; then
    export T1_TREMBL_DAT_SOURCE="$SOURCE_DB_ROOT/uniprot/release_2026_02/uniprot_trembl.dat.gz"
fi

export DB_ROOT="$STAGED_DB_ROOT"
export RUN_ROOT="$SCRATCH_RUN_ROOT"
export WORK_DIR="$WORK/extracted"
export PYTHON_BIN
export REMOVE_ARCHIVES_AFTER_EXTRACT=1
export CAFA_BUILDER_USE_PIGZ=1

bash scripts/benchmark_generation/run_contemporary_temporal_benchmark.sh
copy_results

echo "Finished successfully: $(date)"
