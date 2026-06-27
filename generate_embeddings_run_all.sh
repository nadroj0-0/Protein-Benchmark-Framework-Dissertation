#!/bin/bash
# generate_embeddings_run_all.sh  — SINGLE ENTRY POINT
# From-scratch embedding generation. Clones PFP, builds env, downloads external
# deps, regenerates splits, then generates all 4 modality embeddings in serial.
# (Modalities are split into separate scripts so they can be parallelised later.)
#
# KNOWN GAP: the sequence (ProtT5) step needs data/proteins.fasta, which no
# current step produces. That stage will fail until a FASTA-generation step is added.
set -euo pipefail
LOGFILE="generate_embeddings_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOGFILE")
exec 2>&1

HERE="$(cd "$(dirname "$0")" && pwd)"

# --- Shared setup (clone + env + deps; NO Zenodo tarballs — we generate those) ---
git clone https://github.com/psipred/PFP.git
cd PFP
REPO="$(pwd)"

micromamba create -y -n mmfp python=3.11
eval "$(micromamba shell hook --shell bash)"
micromamba activate mmfp

pip install -r requirements.txt
pip install requests h5py fair-esm

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
echo "    Next: retrain + eval on these (reproduce_retrain_eval.sh logic)."
