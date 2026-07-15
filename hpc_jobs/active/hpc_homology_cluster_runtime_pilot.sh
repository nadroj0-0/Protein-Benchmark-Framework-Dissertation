#!/usr/bin/env bash
# Direct qsub entrypoint for the recommended 30% homology pilot.

#$ -l tmem=64G
#$ -l tscratch=300G
#$ -l scratch0free=300G
#$ -l h_rt=96:0:0
#$ -pe smp 8
#$ -t 1
#$ -j y
#$ -N homology_2026_pilot

set -euo pipefail

FRAMEWORK_ROOT="${FRAMEWORK_SOURCE_ROOT:-${SGE_O_WORKDIR:-$PWD}}"
export HOMOLOGY_RUNTIME_KIND=pilot
exec bash "$FRAMEWORK_ROOT/scripts/benchmark_generation/run_homology_cluster_runtime_hpc.sh"
