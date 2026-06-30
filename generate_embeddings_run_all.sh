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
echo "==> [0/7] External dependencies"
bash "${HERE}/generate_embeddings_dependencies.sh"
source external/dependency_env.sh

echo "==> [1/7] Data preparation (raw CSVs -> splits/labels/sequences)"
bash "${HERE}/generate_embeddings_prepare_data.sh"

echo "[2/7] Verifying generated dataset splits..."
python "${HERE}/verify_splits.py" \
    --data-dir data \
    --strict

echo "==> [3/7] Sequence (ProtT5) embeddings"
bash "${HERE}/generate_embeddings_sequence.sh"

echo "==> [4/7] Text embeddings"
bash "${HERE}/generate_embeddings_text.sh"

echo "==> [5/7] Structure embeddings"
bash "${HERE}/generate_embeddings_structure.sh"

echo "==> [6/7] PPI embeddings"
bash "${HERE}/generate_embeddings_ppi.sh"

echo "[7/7] Verifying generated embeddings..."
python "${HERE}/verify_embeddings.py" \
    --data-dir data \
    --config "${HERE}/configs/cafa3.json" \
    --strict

echo "==> Embedding generation complete. data/embedding_cache/ populated."