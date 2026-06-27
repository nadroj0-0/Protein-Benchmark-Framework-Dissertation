#!/bin/bash
# generate_embeddings_structure.sh
# Structure embeddings (512-D ESM-IF1).  README "3. Structure Embeddings".
# Step 1 downloads AlphaFold PDBs; Step 2 runs ESM-IF1 over them.
set -euo pipefail

export CAFA_ASSESSMENT_DIR="${CAFA_ASSESSMENT_DIR:-external/CAFA_assessment_tool}"

python scripts/check_alphafold_coverage.py \
  --cafa-assessment-dir "${CAFA_ASSESSMENT_DIR}"

python scripts/extract_esm_if1_embeddings.py \
    --pdb_dir data/alphafold_structures \
    --output_dir data/embedding_cache/IF1 \
    --pooling mean \
    --device cuda
