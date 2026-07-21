#!/usr/bin/env bash
# Direct qsub entrypoint for all six homology identity thresholds.

# UCL Grid Engine charges consumable tmem/tscratch per SMP slot:
# 8 slots x 15G = 120G memory; 8 slots x 38G = 304G scratch per array task.
# Each identity is independent and all six may run concurrently when resources permit.
# Interim 30% pilot accounting exceeded the old 72G memory reservation and approached
# the old 96-hour limit, so production carries deliberate memory and walltime headroom.
# Final pilot accounting must still be reviewed before submission.
#$ -l tmem=15G
#$ -l tscratch=38G
#$ -l scratch0free=300G
#$ -l h_rt=168:0:0
#$ -pe smp 8
#$ -t 1-6
#$ -tc 6
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
