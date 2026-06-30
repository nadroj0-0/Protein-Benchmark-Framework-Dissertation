#!/bin/bash
# generate_embeddings_structure.sh
# Structure embeddings (512-D ESM-IF1).  README "3. Structure Embeddings".
# Step 1 downloads AlphaFold PDBs; Step 2 runs ESM-IF1 over them.
# Explicit paths so check_alphafold_coverage.py does NOT use its ../data default
# (which assumes CWD=scripts/ and breaks from repo root).
# DEVICE overridable for local testing.
set -euo pipefail

export CAFA_ASSESSMENT_DIR="${CAFA_ASSESSMENT_DIR:-external/CAFA_assessment_tool}"
DEVICE="${DEVICE:-cuda}"

python scripts/check_alphafold_coverage.py \
  --cafa-assessment-dir "${CAFA_ASSESSMENT_DIR}" \
  --data-dir data \
  --pdb-output-dir data/alphafold_structures \
  --output-file data/alphafold_coverage_results.txt

python scripts/extract_esm_if1_embeddings.py \
    --pdb_dir data/alphafold_structures \
    --output_dir data/embedding_cache/IF1 \
    --pooling mean \
    --device "${DEVICE}"