#!/bin/bash

# CAFA3 DeepGOPlus historical validation wrapper for UCL/SGE.
# This validates the released DeepGOPlus/TEMPROT/PFP CSV path. Heavy temporary
# files stay in scratch; reports/logs are copied home.

#$ -l tmem=24G
#$ -l tscratch=30G
#$ -l scratch0free=30G
#$ -l h_rt=12:0:0
#$ -j y
#$ -N cafa3_dgp_val
#$ -V

set -euo pipefail

CLI_ARTIFACT_CATALOG="${ARTIFACT_CATALOG:-}"
CLI_DEEPGOPLUS_ARCHIVE=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --artifact-catalog) CLI_ARTIFACT_CATALOG="$2"; shift 2 ;;
        --deepgoplus-archive) CLI_DEEPGOPLUS_ARCHIVE="$2"; shift 2 ;;
        *) echo "Unknown argument: $1" >&2; exit 2 ;;
    esac
done

WORK=/scratch0/cafa3_deepgoplus_validation_${JOB_ID}
OUTDIR="$HOME/cafa3_deepgoplus_validation_reports"
FRAMEWORK_REPO_URL="${FRAMEWORK_REPO_URL:-https://github.com/nadroj0-0/Protein-Benchmark-Framework-Dissertation.git}"
FRAMEWORK_DIR="$WORK/Protein-Benchmark-Framework-Dissertation"
RUN_TAG="${JOB_ID:-manual}_$(date +%Y%m%d_%H%M%S)"

cleanup() {
    local status=$?
    echo
    echo "Cleaning scratch directory: $WORK"
    cd ~/
    rm -rf "$WORK"
    exit "$status"
}
trap cleanup EXIT
trap 'echo "Received kill signal"; exit 130' SIGINT SIGTERM

hostname
echo "Job ID      : ${JOB_ID:-manual}"
echo "Working dir : $WORK"
echo "Output dir  : $OUTDIR/$RUN_TAG"
echo "Started at  : $(date)"
echo

mkdir -p "$WORK" "$OUTDIR"
cd "$WORK"

echo "Cloning dissertation framework into local scratch..."
git clone "$FRAMEWORK_REPO_URL" "$FRAMEWORK_DIR"

cd "$FRAMEWORK_DIR"

source scripts/reproduction_common.sh
export ARTIFACT_CATALOG="$CLI_ARTIFACT_CATALOG"
load_framework_paths "$FRAMEWORK_DIR"
activate_or_create_mmfp_env

echo
echo "Running CAFA3 DeepGOPlus historical validation workflow"
echo "Command:"
echo "bash scripts/validation/run_cafa3_deepgoplus_validation.sh"
echo

export SCRATCH_BASE="$WORK"
export TIMESTAMP="$RUN_TAG"
export REPORT_COPY_DIR="$OUTDIR/$RUN_TAG"
export KEEP_SCRATCH=0
export PYTHON_BIN=python

command=(bash scripts/validation/run_cafa3_deepgoplus_validation.sh)
[[ -z "${ARTIFACT_CATALOG:-}" ]] || command+=(--artifact-catalog "$ARTIFACT_CATALOG")
[[ -z "$CLI_DEEPGOPLUS_ARCHIVE" ]] || command+=(--deepgoplus-archive "$CLI_DEEPGOPLUS_ARCHIVE")
"${command[@]}"
STATUS=$?

echo
echo "Finished at: $(date)"
echo "Reports copied to: $OUTDIR/$RUN_TAG"
exit "$STATUS"
