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

CLI_ARTIFACT_CATALOG="${ARTIFACT_CATALOG:-}"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --artifact-catalog)
            [[ $# -ge 2 ]] || { echo "--artifact-catalog requires a path" >&2; exit 2; }
            CLI_ARTIFACT_CATALOG="$2"
            shift 2
            ;;
        *) echo "Unknown argument: $1" >&2; exit 2 ;;
    esac
done

WORK="/scratch0/contemporary_cafa_${JOB_ID:-manual}"
SOURCE_DB_ROOT="${PROTEIN_DATABASE_ROOT:-}"
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
    [[ -n "$SOURCE_DB_ROOT" ]] || return 1
    local source="$SOURCE_DB_ROOT/$relative"
    local destination="$STAGED_DB_ROOT/$relative"
    if [[ ! -f "$source" ]]; then
        return 1
    fi
    mkdir -p "$(dirname "$destination")"
    cp -a "$source" "$destination"
}

stage_catalog_or_root() {
    local artifact_id="$1"
    local relative="$2"
    local explicit_path="${3:-}"
    local source=""
    if [[ -n "$explicit_path" ]]; then
        source="$(resolve_artifact_path "$artifact_id" "$explicit_path" || true)"
    elif [[ -n "$SOURCE_DB_ROOT" && -f "$SOURCE_DB_ROOT/$relative" ]]; then
        source="$SOURCE_DB_ROOT/$relative"
    else
        source="$(resolve_artifact_path "$artifact_id" "" || true)"
    fi
    if [[ -s "$source" ]]; then
        local destination="$STAGED_DB_ROOT/$relative"
        mkdir -p "$(dirname "$destination")"
        cp -a "$source" "$destination"
        printf 'staged-artifact\t%s\t%s\n' "$artifact_id" "$source" \
            >> "$SCRATCH_RUN_ROOT/logs/local_staging.tsv"
    else
        printf 'download-required\t%s\t%s\n' "$artifact_id" "$relative" \
            >> "$SCRATCH_RUN_ROOT/logs/local_staging.tsv"
    fi
}

resolve_large_source() {
    local variable_name="$1"
    local artifact_id="$2"
    local root_relative="$3"
    local explicit="${!variable_name:-}"
    local source=""
    if [[ -n "$explicit" ]]; then
        source="$(resolve_artifact_path "$artifact_id" "$explicit" || true)"
    elif [[ -n "$SOURCE_DB_ROOT" && -s "$SOURCE_DB_ROOT/$root_relative" ]]; then
        source="$SOURCE_DB_ROOT/$root_relative"
    else
        source="$(resolve_artifact_path "$artifact_id" "" || true)"
    fi
    if [[ -n "$source" ]]; then
        printf -v "$variable_name" '%s' "$source"
        export "$variable_name"
        add_mmfp_singularity_bind "$(dirname "$source")"
        printf 'streamed-artifact\t%s\t%s\n' "$artifact_id" "$source" \
            >> "$SCRATCH_RUN_ROOT/logs/local_staging.tsv"
    else
        unset "$variable_name"
    fi
}

stage_if_present() {
    local relative="$1"
    if stage_file "$relative"; then
        printf 'staged\t%s\n' "$relative" >> "$SCRATCH_RUN_ROOT/logs/local_staging.tsv"
    else
        printf 'download-required\t%s\n' "$relative" >> "$SCRATCH_RUN_ROOT/logs/local_staging.tsv"
    fi
}

stage_explicit_file() {
    local variable_name="$1"
    local relative="$2"
    local source="${!variable_name:-}"
    if [[ -z "$source" ]]; then
        stage_if_present "$relative"
        return
    fi
    if [[ ! -f "$source" ]]; then
        artifact_catalog_warn "optional $variable_name cache does not exist; trying the local-root/download fallback: $source"
        stage_if_present "$relative"
        return
    fi
    local destination="$STAGED_DB_ROOT/$relative"
    mkdir -p "$(dirname "$destination")"
    cp -a "$source" "$destination"
    printf 'staged-custom\t%s\n' "$source" >> "$SCRATCH_RUN_ROOT/logs/local_staging.tsv"
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
echo "t1 endpoint : ${T1_ENDPOINT_POLICY:-snapshot-membership}"
echo "Backfill    : ${T1_BACKFILL_POLICY:-allow}"
echo "Scratch     : $WORK"
echo "Final output: $FINAL_RUN_ROOT"

mkdir -p "$WORK" "$RESULTS_ROOT" "$STAGED_DB_ROOT" "$SCRATCH_RUN_ROOT/logs"
printf 'action\tartifact_or_path\tsource_or_destination\n' > "$SCRATCH_RUN_ROOT/logs/local_staging.tsv"

echo "Cloning dissertation framework into scratch"
git clone "$FRAMEWORK_REPO_URL" "$FRAMEWORK_DIR"
cd "$FRAMEWORK_DIR"
source scripts/reproduction_common.sh
export ARTIFACT_CATALOG="$CLI_ARTIFACT_CATALOG"
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
    if [[ ! -s "$STAGED_DB_ROOT/uniprot/release_2025_01/uniprot_sprot.dat.gz" ]]; then
        stage_catalog_or_root uniprot_sprot_t0 \
            "uniprot/release_2025_01/uniprot_sprot-only2025_01.tar.gz"
    fi
    stage_explicit_file T0_TREMBL_FILTERED_INPUT "uniprot/release_2025_01/uniprot_trembl_cafa3_targets.dat.gz"
    stage_catalog_or_root uniprot_sprot_t1 \
        "uniprot/release_2026_02/uniprot_sprot.dat.gz" \
        "${T1_SPROT_INPUT:-}"
    stage_explicit_file T1_TREMBL_FILTERED_INPUT "uniprot/release_2026_02/uniprot_trembl_cafa3_targets.dat.gz"
fi
stage_catalog_or_root goa_t0 "goa/release_2025_01/goa_uniprot_all.gaf.225.gz" "${GOA_T0_INPUT:-}"
stage_catalog_or_root goa_t1 "goa/release_2026_02/goa_uniprot_all.gaf.234.gz" "${GOA_T1_INPUT:-}"
stage_catalog_or_root go_basic_t0 "ontology/release_2025-02-06/go-basic.obo" "${GO_BASIC_T0_INPUT:-}"
stage_catalog_or_root go_basic_t0_source_resolution "ontology/release_2025-03-16/go-basic.obo" "${GO_BASIC_T0_SOURCE_INPUT:-}"
stage_catalog_or_root go_basic_t1 "ontology/release_2026-06-19/go-basic.obo" "${GO_BASIC_T1_INPUT:-}"

# Full TrEMBL sources are 100+ GB. If present locally, stream-filter them in
# place rather than duplicating them into 200 GB scratch.
if [[ -n "${T0_TREMBL_DAT_SOURCE:-}" && -s "$T0_TREMBL_DAT_SOURCE" ]]; then
    add_mmfp_singularity_bind "$(dirname "$T0_TREMBL_DAT_SOURCE")"
    printf 'streamed-explicit\t%s\t%s\n' T0_TREMBL_DAT_SOURCE "$T0_TREMBL_DAT_SOURCE" \
        >> "$SCRATCH_RUN_ROOT/logs/local_staging.tsv"
elif [[ -n "$SOURCE_DB_ROOT" && -s "$SOURCE_DB_ROOT/uniprot/release_2025_01/uniprot_trembl.dat.gz" ]]; then
    export T0_TREMBL_DAT_SOURCE="$SOURCE_DB_ROOT/uniprot/release_2025_01/uniprot_trembl.dat.gz"
    add_mmfp_singularity_bind "$(dirname "$T0_TREMBL_DAT_SOURCE")"
    printf 'streamed-root\t%s\t%s\n' T0_TREMBL_DAT_SOURCE "$T0_TREMBL_DAT_SOURCE" \
        >> "$SCRATCH_RUN_ROOT/logs/local_staging.tsv"
else
    if [[ -n "${T0_TREMBL_DAT_SOURCE:-}" ]]; then
        artifact_catalog_warn "T0_TREMBL_DAT_SOURCE is missing or empty; trying the full-archive fallback: $T0_TREMBL_DAT_SOURCE"
        unset T0_TREMBL_DAT_SOURCE
    fi
    resolve_large_source T0_TREMBL_ARCHIVE_SOURCE uniprot_knowledgebase_t0 \
        "uniprot/release_2025_01/knowledgebase2025_01.tar.gz"
fi
resolve_large_source T1_TREMBL_DAT_SOURCE uniprot_trembl_t1 \
    "uniprot/release_2026_02/uniprot_trembl.dat.gz"

export DB_ROOT="$STAGED_DB_ROOT"
export RUN_ROOT="$SCRATCH_RUN_ROOT"
export WORK_DIR="$WORK/extracted"
export PYTHON_BIN
export REMOVE_ARCHIVES_AFTER_EXTRACT=1
export CAFA_BUILDER_USE_PIGZ=1
export T1_ENDPOINT_POLICY="${T1_ENDPOINT_POLICY:-snapshot-membership}"
export T1_BACKFILL_POLICY="${T1_BACKFILL_POLICY:-allow}"

bash scripts/benchmark_generation/run_contemporary_temporal_benchmark.sh
copy_results

echo "Finished successfully: $(date)"
