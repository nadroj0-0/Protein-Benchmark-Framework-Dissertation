#!/bin/bash
# generate_embeddings_structure.sh
# Structure embeddings (512-D ESM-IF1).  README "3. Structure Embeddings".
# Step 1 downloads AlphaFold PDBs; Step 2 runs ESM-IF1 over them.
# Explicit paths so check_alphafold_coverage.py does NOT use its ../data default
# (which assumes CWD=scripts/ and breaks from repo root).
# DEVICE overridable for local testing.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
if [ -f "${HERE}/configs/paths.local.sh" ]; then
  # Machine-specific paths are intentionally not committed.
  # shellcheck disable=SC1091
  source "${HERE}/configs/paths.local.sh"
fi

export CAFA_ASSESSMENT_DIR="${CAFA_ASSESSMENT_DIR:-external/CAFA_assessment_tool}"
if [ ! -d "${CAFA_ASSESSMENT_DIR}" ]; then
  echo "Missing CAFA assessment tool directory: ${CAFA_ASSESSMENT_DIR}" >&2
  echo "Set CAFA_ASSESSMENT_DIR in configs/paths.local.sh or the environment." >&2
  exit 1
fi
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
