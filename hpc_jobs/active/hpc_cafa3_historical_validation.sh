#!/bin/bash

# CAFA3 historical validation wrapper for UCL/SGE.
# Heavy downloads stay in scratch; generated artefacts, reports and logs are copied home.

#$ -l tmem=48G
#$ -l tscratch=200G
#$ -l scratch0free=200G
#$ -l h_rt=72:0:0
#$ -j y
#$ -N cafa3_hist_val
#$ -V

set -euo pipefail

HISTORICAL_TEST_SOURCE="${HISTORICAL_TEST_SOURCE:-official-groundtruth}"
if [ -z "${HISTORICAL_BENCHMARK_ONTOLOGY:-}" ]; then
    if [ "$HISTORICAL_TEST_SOURCE" = "official-groundtruth" ]; then
        HISTORICAL_BENCHMARK_ONTOLOGY="deepgoplus-packaged"
    else
        HISTORICAL_BENCHMARK_ONTOLOGY="february-go-basic"
    fi
fi

WORK=/scratch0/cafa3_historical_validation_${JOB_ID}
OUTDIR="$HOME/cafa3_historical_validation_reports"
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
echo "Training    : ${HISTORICAL_TRAINING_SNAPSHOT:-september-2016}"
echo "Targets     : ${TARGET_UNIVERSE_POLICY:-official-cafa3-targets}"
echo "Test source : ${HISTORICAL_TEST_SOURCE}"
echo "t1 endpoint : ${HISTORICAL_T1_ENDPOINT_POLICY:-assigned-date-proxy}"
echo "Backfill    : ${HISTORICAL_BACKFILL_POLICY:-exclude-pre-t0}"
echo "Ontology    : ${HISTORICAL_BENCHMARK_ONTOLOGY:-february-go-basic}"
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
export PYTHON_BIN="$(command -v python)"

echo
echo "Running CAFA3 historical validation workflow"
echo "Command:"
echo "bash scripts/validation/run_cafa3_historical_validation.sh"
echo

export SCRATCH_BASE="$WORK"
export TIMESTAMP="$RUN_TAG"
export REPORT_COPY_DIR="$OUTDIR/$RUN_TAG"
export KEEP_SCRATCH=0
export PYTHON_BIN
export DECOMPRESS_GOA=1
export USE_PIGZ=1
export HISTORICAL_TRAINING_SNAPSHOT="${HISTORICAL_TRAINING_SNAPSHOT:-september-2016}"
export TARGET_UNIVERSE_POLICY="${TARGET_UNIVERSE_POLICY:-official-cafa3-targets}"
export HISTORICAL_TEST_SOURCE
export HISTORICAL_T1_ENDPOINT_POLICY="${HISTORICAL_T1_ENDPOINT_POLICY:-assigned-date-proxy}"
export HISTORICAL_BACKFILL_POLICY="${HISTORICAL_BACKFILL_POLICY:-exclude-pre-t0}"
export HISTORICAL_BENCHMARK_ONTOLOGY

bash scripts/validation/run_cafa3_historical_validation.sh
STATUS=$?

echo
echo "Finished at: $(date)"
echo "Generated artefacts and reports copied to: $OUTDIR/$RUN_TAG"
exit "$STATUS"
