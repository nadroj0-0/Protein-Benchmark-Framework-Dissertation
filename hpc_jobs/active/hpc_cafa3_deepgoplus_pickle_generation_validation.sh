#!/bin/bash

# CAFA3 DeepGOPlus cafa3_data.py-style pickle-generation validation wrapper.
# Heavy temporary files stay in scratch; reports/logs are copied home.

#$ -l tmem=24G
#$ -l tscratch=30G
#$ -l scratch0free=30G
#$ -l h_rt=12:0:0
#$ -j y
#$ -N cafa3_dgp_pkl
#$ -V

set -euo pipefail

WORK=/scratch0/cafa3_deepgoplus_pickle_generation_${JOB_ID}
OUTDIR="$HOME/cafa3_deepgoplus_pickle_generation_reports"
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
load_framework_paths "$FRAMEWORK_DIR"
activate_or_create_mmfp_env

echo
echo "Running CAFA3 DeepGOPlus pickle-generation validation workflow"
echo "Command:"
echo "bash scripts/validation/run_cafa3_deepgoplus_pickle_generation_validation.sh"
echo

export SCRATCH_BASE="$WORK"
export TIMESTAMP="$RUN_TAG"
export REPORT_COPY_DIR="$OUTDIR/$RUN_TAG"
export KEEP_SCRATCH=0
export PYTHON_BIN=python

bash scripts/validation/run_cafa3_deepgoplus_pickle_generation_validation.sh
STATUS=$?

echo
echo "Finished at: $(date)"
echo "Reports copied to: $OUTDIR/$RUN_TAG"
exit "$STATUS"
