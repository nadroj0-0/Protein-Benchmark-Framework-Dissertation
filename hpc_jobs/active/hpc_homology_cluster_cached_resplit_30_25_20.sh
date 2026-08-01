#!/usr/bin/env bash
# Rebuild the 30%, 25% and 20% homology benchmarks from validated MMseqs caches.
# REQUIRE_HOMOLOGY_CLUSTER_CACHE makes any cache miss fatal; this job cannot recluster.

#$ -l tmem=10G
#$ -l tscratch=15G
#$ -l scratch0free=50G
#$ -l h_rt=36:0:0
#$ -pe smp 4
#$ -t 1-3
#$ -tc 3
#$ -j y
#$ -N hom_resplit

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
export UNIREF_LEVEL=50
export MMSEQS_SENSITIVITY=4
export SPLIT_POLICY=cluster-count-random
export UNIREF50_FASTA="${UNIREF50_FASTA:-/SAN/bioinf/bmpfp/frozen_inputs/uniref50/2026_02/uniref50.fasta.gz}"
export GO_OBO="${GO_OBO:-/SAN/bioinf/bmpfp/frozen_inputs/ontology/2026-06-19/go-basic.obo}"
export HOMOLOGY_COMMON_PREPROCESSING_CACHE="${HOMOLOGY_COMMON_PREPROCESSING_CACHE:-/SAN/bioinf/bmpfp/derived_inputs/homology/2026_02/goa_234/sprot-and-trembl/uniref50/common_preprocessing}"
export HOMOLOGY_CLUSTER_CACHE_ROOT="${HOMOLOGY_CLUSTER_CACHE_ROOT:-/SAN/bioinf/bmpfp/derived_inputs/homology/2026_02/mmseqs_cluster_assignments/uniref50_sensitivity_4}"
export RESULTS_ROOT="${RESULTS_ROOT:-/SAN/bioinf/bmpfp/benchmarks/homology/2026_02/uniref50_sensitivity_4_daniel_aligned_cached_random_resplit}"
export REQUIRE_HOMOLOGY_COMMON_CACHE=1
export REQUIRE_HOMOLOGY_CLUSTER_CACHE=1
export MINIMUM_CLUSTER_CACHE_FREE_GB="${MINIMUM_CLUSTER_CACHE_FREE_GB:-20}"
export MINIMUM_SCRATCH_GB="${MINIMUM_SCRATCH_GB:-50}"
export SCRATCH_SAFETY_MULTIPLIER="${SCRATCH_SAFETY_MULTIPLIER:-1}"
export MMSEQS_WORK_MULTIPLIER="${MMSEQS_WORK_MULTIPLIER:-1}"
export PUBLICATION_SAFETY_MULTIPLIER="${PUBLICATION_SAFETY_MULTIPLIER:-1}"

exec bash "$FRAMEWORK_ROOT/scripts/benchmark_generation/run_homology_cluster_runtime_hpc.sh"
