#!/usr/bin/env bash
# Six-threshold UniRef50 scaffold run. The legacy UniRef90/7.5 entrypoints remain unchanged.

# UCL Grid Engine charges tmem/tscratch per SMP slot:
# 24 slots x 7G = 168G memory; 24 slots x 18G = 432G scratch per task.
# The array is intentionally not submittable until the separately validated UniRef50
# FASTA and UniRef50 common-preprocessing cache exist at the configured paths.
#$ -l tmem=7G
#$ -l tscratch=18G
#$ -l scratch0free=450G
#$ -l h_rt=168:0:0
#$ -pe smp 24
#$ -t 1-6
#$ -tc 6
#$ -j y
#$ -N homology_24c_u50

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
export MMSEQS_SENSITIVITY="${MMSEQS_SENSITIVITY:-4}"
export UNIREF50_FASTA="${UNIREF50_FASTA:-/SAN/bioinf/bmpfp/frozen_inputs/uniref50/2026_02/uniref50.fasta.gz}"
export GO_OBO="${GO_OBO:-/SAN/bioinf/bmpfp/frozen_inputs/ontology/2026-06-19/go-basic.obo}"
export HOMOLOGY_COMMON_PREPROCESSING_CACHE="${HOMOLOGY_COMMON_PREPROCESSING_CACHE:-/SAN/bioinf/bmpfp/derived_inputs/homology/2026_02/goa_234/sprot-and-trembl/uniref50/common_preprocessing}"
export HOMOLOGY_CLUSTER_CACHE_ROOT="${HOMOLOGY_CLUSTER_CACHE_ROOT:-/SAN/bioinf/bmpfp/derived_inputs/homology/2026_02/mmseqs_cluster_assignments/uniref50_sensitivity_4}"
export RESULTS_ROOT="${RESULTS_ROOT:-/SAN/bioinf/bmpfp/benchmarks/homology/2026_02/uniref50_sensitivity_4_daniel_aligned_24core}"
export REQUIRE_HOMOLOGY_COMMON_CACHE=1
export MINIMUM_CLUSTER_CACHE_FREE_GB="${MINIMUM_CLUSTER_CACHE_FREE_GB:-40}"
export MINIMUM_SCRATCH_GB="${MINIMUM_SCRATCH_GB:-300}"
export SCRATCH_SAFETY_MULTIPLIER="${SCRATCH_SAFETY_MULTIPLIER:-1.5}"
export MMSEQS_WORK_MULTIPLIER="${MMSEQS_WORK_MULTIPLIER:-2.0}"
export PUBLICATION_SAFETY_MULTIPLIER="${PUBLICATION_SAFETY_MULTIPLIER:-1.5}"

exec bash "$FRAMEWORK_ROOT/scripts/benchmark_generation/run_homology_cluster_runtime_hpc.sh"
