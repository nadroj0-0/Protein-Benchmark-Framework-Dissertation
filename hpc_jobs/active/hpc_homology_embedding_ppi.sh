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
FRAMEWORK_JOB_ROOT="${FRAMEWORK_JOB_ROOT:-$HOME/Protein-Benchmark-Framework-Dissertation}"
COMMON_JOB="$FRAMEWORK_JOB_ROOT/hpc_jobs/lib/run_homology_embedding_modality_job.sh"
[[ -f "$COMMON_JOB" ]] || {
  echo "Missing shared homology embedding job body: $COMMON_JOB" >&2
  exit 2
}
exec bash "$COMMON_JOB" \
  --modality ppi "$@"
