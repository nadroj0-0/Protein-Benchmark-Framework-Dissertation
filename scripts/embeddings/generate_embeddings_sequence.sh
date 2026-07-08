#!/bin/bash
# generate_embeddings_sequence.sh
# ProtT5-XL sequence embeddings (1024-D).
# DEVICE overridable (default cuda; set DEVICE=cpu for local testing).
# NOTE: reads data/proteins.fasta — produced by the FASTA-generation step.
set -euo pipefail

DEVICE="${DEVICE:-cuda}"
python scripts/extract_prott5_embeddings.py \
    --fasta_file data/proteins.fasta \
    --output_dir data/embedding_cache/prott5 \
    --batch_size 8 \
    --device "${DEVICE}"
