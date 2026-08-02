#!/usr/bin/env bash
#$ -l tmem=24G
#$ -l tscratch=80G
#$ -l scratch0free=100G
#$ -l h_rt=48:0:0
#$ -pe smp 2
#$ -j y
#$ -N hom30_emb_ppi
#$ -V

set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
exec bash "$HERE/../lib/run_homology_embedding_modality_job.sh" \
  --modality ppi "$@"
