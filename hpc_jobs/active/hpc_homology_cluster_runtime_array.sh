#!/usr/bin/env bash
# Direct qsub entrypoint for all six homology identity thresholds.

# UCL Grid Engine charges consumable tmem/tscratch per SMP slot:
# 2 slots x 32G = 64G memory; 2 slots x 150G = 300G scratch per array task.
# Each identity remains an independent array element; at most two run at once.
#$ -l tmem=32G
#$ -l tscratch=150G
#$ -l scratch0free=300G
#$ -l h_rt=96:0:0
#$ -pe smp 2
#$ -t 1-6
#$ -tc 2
#$ -j y
#$ -N homology_2026_array

set -euo pipefail

if [[ "${1:-}" == "--artifact-catalog" ]]; then
  [[ $# -ge 2 ]] || { echo "--artifact-catalog requires a path" >&2; exit 2; }
  export ARTIFACT_CATALOG="$2"
  shift 2
fi
[[ $# -eq 0 ]] || { echo "Unknown array argument: $1" >&2; exit 2; }

FRAMEWORK_ROOT="${FRAMEWORK_SOURCE_ROOT:-${SGE_O_WORKDIR:-$PWD}}"
export HOMOLOGY_RUNTIME_KIND=array
exec bash "$FRAMEWORK_ROOT/scripts/benchmark_generation/run_homology_cluster_runtime_hpc.sh"
