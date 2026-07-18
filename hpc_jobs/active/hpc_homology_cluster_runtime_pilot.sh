#!/usr/bin/env bash
# Direct qsub entrypoint for the recommended 30% homology pilot.

# UCL Grid Engine charges consumable tmem/tscratch per SMP slot:
# 4 slots x 16G = 64G memory; 4 slots x 75G = 300G scratch.
# Four MMseqs2 threads improve clustering throughput without restoring the old 8-slot wait.
# The earlier single-threaded streaming/indexing stages will not scale linearly with this change.
#$ -l tmem=16G
#$ -l tscratch=75G
#$ -l scratch0free=300G
#$ -l h_rt=96:0:0
#$ -pe smp 4
#$ -t 1
#$ -j y
#$ -N homology_2026_pilot

set -euo pipefail

if [[ "${1:-}" == "--artifact-catalog" ]]; then
  [[ $# -ge 2 ]] || { echo "--artifact-catalog requires a path" >&2; exit 2; }
  export ARTIFACT_CATALOG="$2"
  shift 2
fi
[[ $# -eq 0 ]] || { echo "Unknown pilot argument: $1" >&2; exit 2; }

FRAMEWORK_ROOT="${FRAMEWORK_SOURCE_ROOT:-${SGE_O_WORKDIR:-$PWD}}"
export HOMOLOGY_RUNTIME_KIND=pilot
exec bash "$FRAMEWORK_ROOT/scripts/benchmark_generation/run_homology_cluster_runtime_hpc.sh"
