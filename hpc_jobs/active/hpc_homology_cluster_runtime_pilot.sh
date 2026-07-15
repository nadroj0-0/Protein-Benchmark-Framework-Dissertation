#!/usr/bin/env bash
# Direct qsub entrypoint for the recommended 30% homology pilot.

# UCL Grid Engine charges consumable tmem/tscratch per SMP slot:
# 8 slots x 8G = 64G memory; 8 slots x 38G = 304G scratch.
#$ -l tmem=8G
#$ -l tscratch=38G
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
