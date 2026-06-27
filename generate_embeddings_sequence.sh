#!/bin/bash
# generate_embeddings_sequence.sh
# ProtT5-XL sequence embeddings (1024-D).  README "4. ProtT5 Embeddings".
# NOTE: reads data/proteins.fasta, which is NOT produced by any current step
#       (sequences are in *_sequences.json). This will fail until a FASTA-
#       generation step exists. Left faithful to the README on purpose.
set -euo pipefail

python scripts/extract_prott5_embeddings.py \
    --fasta_file data/proteins.fasta \
    --output_dir data/embedding_cache/prott5 \
    --batch_size 8
