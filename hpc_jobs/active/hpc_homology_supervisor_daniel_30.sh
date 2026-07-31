#!/usr/bin/env bash
# Validate Daniel Buchan's external UniRef50 30% assignment and continue the benchmark.

#$ -l tmem=10G
#$ -l tscratch=25G
#$ -l scratch0free=100G
#$ -l h_rt=48:0:0
#$ -pe smp 4
#$ -t 1
#$ -j y
#$ -N homology_daniel_30

set -euo pipefail

FRAMEWORK_ROOT="${FRAMEWORK_SOURCE_ROOT:-${SGE_O_WORKDIR:-$PWD}}"
SUPERVISOR_ROOT="/SAN/bioinf/bmpfp/derived_inputs/homology/2026_02/supervisor_daniel_buchan/mmseqs_cluster_assignments/uniref50_sensitivity_4/identity_30/raw"

export HOMOLOGY_RUNTIME_KIND=pilot
export MMSEQS_PROFILE=daniel-aligned-defaults
export UNIREF_LEVEL=50
export MMSEQS_SENSITIVITY=4
export UNIREF50_FASTA=/SAN/bioinf/bmpfp/frozen_inputs/uniref50/2026_02/uniref50.fasta.gz
export GO_OBO=/SAN/bioinf/bmpfp/frozen_inputs/ontology/2026-06-19/go-basic.obo
export HOMOLOGY_COMMON_PREPROCESSING_CACHE=/SAN/bioinf/bmpfp/derived_inputs/homology/2026_02/goa_234/sprot-and-trembl/uniref50/common_preprocessing
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
