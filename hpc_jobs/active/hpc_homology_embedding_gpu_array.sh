#!/usr/bin/env bash
#$ -l tmem=24G
#$ -l tscratch=100G
#$ -l scratch0free=220G
#$ -l h_rt=96:0:0
#$ -l gpu=true
#$ -pe gpu 1
#$ -t 1-3
#$ -j y
#$ -N hom30_emb_gpu
#$ -V

set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
case "${SGE_TASK_ID:-}" in
  1) modality=sequence ;;
  2) modality=text ;;
  3) modality=structure ;;
  *) echo "Invalid array task: ${SGE_TASK_ID:-unset}" >&2; exit 2 ;;
esac
exec bash "$HERE/../lib/run_homology_embedding_modality_job.sh" \
  --modality "$modality" "$@"
