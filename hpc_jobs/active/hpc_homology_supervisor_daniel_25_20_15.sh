#!/usr/bin/env bash
# Build production CSVs from Daniel Buchan's validated 25%, 20% and 15% assignments.

#$ -l tmem=10G
#$ -l tscratch=25G
#$ -l scratch0free=100G
#$ -l h_rt=48:0:0
#$ -pe smp 4
#$ -t 2-4
#$ -tc 3
#$ -j y
#$ -N hom_d_csv

set -euo pipefail

FRAMEWORK_ROOT="${FRAMEWORK_SOURCE_ROOT:-${SGE_O_WORKDIR:-$PWD}}"
TASK_ID="${SGE_TASK_ID:-}"
case "$TASK_ID" in
    2) IDENTITY=25 ;;
    3) IDENTITY=20 ;;
    4) IDENTITY=15 ;;
    *) echo "SGE_TASK_ID must be 2, 3 or 4" >&2; exit 2 ;;
esac
SUPERVISOR_ROOT="/SAN/bioinf/bmpfp/derived_inputs/homology/2026_02/supervisor_daniel_buchan/mmseqs_cluster_assignments/uniref50_sensitivity_4/identity_${IDENTITY}/raw"

export HOMOLOGY_RUNTIME_KIND=array
export SPLIT_POLICY=cluster-count-random
export MMSEQS_PROFILE=daniel-aligned-defaults
export UNIREF_LEVEL=50
export MMSEQS_SENSITIVITY=4
export UNIREF50_FASTA=/SAN/bioinf/bmpfp/frozen_inputs/uniref50/2026_02/uniref50.fasta.gz
export GO_OBO=/SAN/bioinf/bmpfp/frozen_inputs/ontology/2026-06-19/go-basic.obo
export HOMOLOGY_COMMON_PREPROCESSING_CACHE=/SAN/bioinf/bmpfp/derived_inputs/homology/2026_02/goa_234/sprot-and-trembl/uniref50/common_preprocessing
export REQUIRE_HOMOLOGY_COMMON_CACHE=1
export EXTERNAL_CLUSTER_ASSIGNMENTS="$SUPERVISOR_ROOT/cluster_assignments.tsv.gz"
export EXTERNAL_CLUSTER_PROVENANCE="$SUPERVISOR_ROOT/SOURCE_PROVENANCE.json"
export RESULTS_ROOT=/SAN/bioinf/bmpfp/benchmarks/homology/2026_02/supervisor_daniel_buchan/uniref50_sensitivity_4
export MINIMUM_SCRATCH_GB=60
export SCRATCH_SAFETY_MULTIPLIER=1
export MMSEQS_WORK_MULTIPLIER=1
export PUBLICATION_SAFETY_MULTIPLIER=1
unset HOMOLOGY_CLUSTER_CACHE_ROOT
unset REQUIRE_HOMOLOGY_CLUSTER_CACHE

exec bash "$FRAMEWORK_ROOT/scripts/benchmark_generation/run_homology_cluster_runtime_hpc.sh"
