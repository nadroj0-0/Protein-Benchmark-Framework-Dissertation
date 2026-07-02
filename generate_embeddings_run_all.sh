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
echo "==> [0/8] External dependencies"
bash "${HERE}/generate_embeddings_dependencies.sh"
source external/dependency_env.sh

echo "==> [1/8] Data preparation (raw CSVs -> splits/labels/sequences)"
bash "${HERE}/generate_embeddings_prepare_data.sh"

echo "[2/8] Verifying generated dataset splits..."
python "${HERE}/verify_splits.py" \
    --data-dir data \
    --strict

echo "==> [3/8] Building proteins.fasta from split sequences"
python "${HERE}/generate_embeddings_fasta.py" --data-dir data --config "${HERE}/configs/cafa3.json"

echo "==> [4/8] Sequence (ProtT5) embeddings"
bash "${HERE}/generate_embeddings_sequence.sh"

echo "==> [5/8] Text embeddings"
bash "${HERE}/generate_embeddings_text.sh"

echo "==> [6/8] Structure embeddings"
bash "${HERE}/generate_embeddings_structure.sh"

echo "==> [7/8] PPI embeddings"
bash "${HERE}/generate_embeddings_ppi.sh"

echo "[8/8] Verifying generated embeddings..."
python "${HERE}/verify_embeddings.py" \
    --data-dir data \
    --config "${HERE}/configs/cafa3.json" \
    --strict

echo "==> Embedding generation complete. data/embedding_cache/ populated."