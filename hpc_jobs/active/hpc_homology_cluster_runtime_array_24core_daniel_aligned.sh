#!/usr/bin/env bash
# Fresh six-threshold homology array using Daniel-aligned MMseqs2 defaults.
#
# Scientific profile:
# - UniRef90 FASTA is passed directly to createdb; it is not flattened.
# - createdb input shuffle is left at the pinned MMseqs2 binary's default.
# - cluster uses -s 7.5 and the MMseqs2 default E-value.
# - cluster reassignment and -a 1 backtraces are not requested.
# - Daniel's derived align/convertalis and cluster-FASTA inspection exports are
#   omitted from the six-task production batch. They duplicate the frozen input
#   at a scale the 500G SAN cannot safely hold and are not used by the benchmark.
# - The validated representative/member assignments are checkpointed and then
#   published to the threshold-specific SAN cache as soon as clustering finishes.

# UCL Grid Engine charges tmem/tscratch per SMP slot:
# 24 slots x 7G = 168G memory; 24 slots x 18G = 432G scratch per task.
# All six identities are independent and may run concurrently when resources permit.
#$ -l tmem=7G
#$ -l tscratch=18G
#$ -l scratch0free=450G
#$ -l h_rt=168:0:0
#$ -pe smp 24
#$ -t 1-6
#$ -tc 6
#$ -j y
#$ -N homology_24c_daniel

set -euo pipefail

if [[ "${1:-}" == "--artifact-catalog" ]]; then
  [[ $# -ge 2 ]] || { echo "--artifact-catalog requires a path" >&2; exit 2; }
  export ARTIFACT_CATALOG="$2"
  shift 2
fi
[[ $# -eq 0 ]] || { echo "Unknown array argument: $1" >&2; exit 2; }

FRAMEWORK_ROOT="${FRAMEWORK_SOURCE_ROOT:-${SGE_O_WORKDIR:-$PWD}}"
export HOMOLOGY_RUNTIME_KIND=array
export MMSEQS_PROFILE=daniel-aligned-defaults
export UNIREF90_FASTA="${UNIREF90_FASTA:-/SAN/bioinf/bmpfp/frozen_inputs/uniref90/2026_02/uniref90.fasta.gz}"
export GO_OBO="${GO_OBO:-/SAN/bioinf/bmpfp/frozen_inputs/ontology/2026-06-19/go-basic.obo}"
export HOMOLOGY_COMMON_PREPROCESSING_CACHE="${HOMOLOGY_COMMON_PREPROCESSING_CACHE:-/SAN/bioinf/bmpfp/derived_inputs/homology/2026_02/goa_234/sprot-and-trembl/common_preprocessing}"
export HOMOLOGY_CLUSTER_CACHE_ROOT="${HOMOLOGY_CLUSTER_CACHE_ROOT:-/SAN/bioinf/bmpfp/derived_inputs/homology/2026_02/mmseqs_cluster_assignments}"
export RESULTS_ROOT="${RESULTS_ROOT:-/SAN/bioinf/bmpfp/benchmarks/homology/2026_02/daniel_aligned_24core}"
export REQUIRE_HOMOLOGY_COMMON_CACHE=1
export MINIMUM_CLUSTER_CACHE_FREE_GB="${MINIMUM_CLUSTER_CACHE_FREE_GB:-40}"
export MINIMUM_SCRATCH_GB="${MINIMUM_SCRATCH_GB:-400}"
export SCRATCH_SAFETY_MULTIPLIER="${SCRATCH_SAFETY_MULTIPLIER:-1.5}"
export MMSEQS_WORK_MULTIPLIER="${MMSEQS_WORK_MULTIPLIER:-2.0}"
export PUBLICATION_SAFETY_MULTIPLIER="${PUBLICATION_SAFETY_MULTIPLIER:-1.5}"

exec bash "$FRAMEWORK_ROOT/scripts/benchmark_generation/run_homology_cluster_runtime_hpc.sh"
