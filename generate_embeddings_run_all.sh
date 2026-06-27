#!/bin/bash
# generate_embeddings_run_all.sh  — SUB-ORCHESTRATOR (not a standalone entry point)
# Called by reproduce_embeddings_retrain_eval.sh AFTER it has cloned PFP, built the
# env, and cd'd into the repo root. Assumes CWD is the PFP repo root and the mmfp
# env is active. Downloads external deps, regenerates splits, then generates all 4
# modality embeddings in serial. (Modalities are separate scripts so they can be
# parallelised later.)
#
# KNOWN GAP: the sequence (ProtT5) step needs data/proteins.fasta, which no current
# step produces. That stage will fail until a FASTA-generation step is added.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"

# --- External database dependencies (writes external/dependency_env.sh) ---
echo "==> [0/5] External dependencies"
bash "${HERE}/generate_embeddings_dependencies.sh"
source external/dependency_env.sh

echo "==> [1/5] Data preparation (raw CSVs -> splits/labels/sequences)"
bash "${HERE}/generate_embeddings_prepare_data.sh"

echo "==> [2/5] Sequence (ProtT5) embeddings"
bash "${HERE}/generate_embeddings_sequence.sh"

echo "==> [3/5] Text embeddings"
bash "${HERE}/generate_embeddings_text.sh"

echo "==> [4/5] Structure embeddings"
bash "${HERE}/generate_embeddings_structure.sh"

echo "==> [5/5] PPI embeddings"
bash "${HERE}/generate_embeddings_ppi.sh"

echo "==> Embedding generation complete. data/embedding_cache/ populated."