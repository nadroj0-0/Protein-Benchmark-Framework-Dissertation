#!/bin/bash
# generate_embeddings_ppi.sh
# PPI embeddings (512-D STRING).  README "1. PPI Embeddings".
set -euo pipefail

export STRING_H5_FILE="${STRING_H5_FILE:-external/string/protein.network.embeddings.v12.0.h5}"
export STRING_ALIAS_FILE="${STRING_ALIAS_FILE:-external/string/protein.aliases.v12.0.txt}"
export CAFA_ASSESSMENT_DIR="${CAFA_ASSESSMENT_DIR:-external/CAFA_assessment_tool}"

python scripts/extract_ppi_embeddings.py \
  --string-h5 "${STRING_H5_FILE}" \
  --string-alias "${STRING_ALIAS_FILE}" \
  --cafa-assessment-dir "${CAFA_ASSESSMENT_DIR}"
