#!/bin/bash
# generate_embeddings_ppi.sh
# PPI embeddings (512-D STRING).  README "1. PPI Embeddings".
# Explicit --data-dir/--output-dir so it doesn't depend on the ./data default.
# NOTE: --cafa3-id-mapping is OPTIONAL (FlyBase supplement); script warns and
# continues if absent. Not passed here — investigate only if PPI coverage is low.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
if [ -f "${HERE}/configs/paths.local.sh" ]; then
  # Machine-specific paths are intentionally not committed.
  # shellcheck disable=SC1091
  source "${HERE}/configs/paths.local.sh"
fi

export STRING_H5_FILE="${STRING_H5_FILE:-external/string/protein.network.embeddings.v12.0.h5}"
export STRING_ALIAS_FILE="${STRING_ALIAS_FILE:-external/string/protein.aliases.v12.0.txt}"
export CAFA_ASSESSMENT_DIR="${CAFA_ASSESSMENT_DIR:-external/CAFA_assessment_tool}"
if [ ! -f "${STRING_H5_FILE}" ]; then
  echo "Missing STRING network embeddings file: ${STRING_H5_FILE}" >&2
  echo "Set STRING_H5_FILE in configs/paths.local.sh or the environment." >&2
  exit 1
fi
if [ ! -f "${STRING_ALIAS_FILE}" ]; then
  echo "Missing STRING alias file: ${STRING_ALIAS_FILE}" >&2
  echo "Set STRING_ALIAS_FILE in configs/paths.local.sh or the environment." >&2
  exit 1
fi
if [ ! -d "${CAFA_ASSESSMENT_DIR}" ]; then
  echo "Missing CAFA assessment tool directory: ${CAFA_ASSESSMENT_DIR}" >&2
  echo "Set CAFA_ASSESSMENT_DIR in configs/paths.local.sh or the environment." >&2
  exit 1
fi

python scripts/extract_ppi_embeddings.py \
  --string-h5 "${STRING_H5_FILE}" \
  --string-alias "${STRING_ALIAS_FILE}" \
  --cafa-assessment-dir "${CAFA_ASSESSMENT_DIR}" \
  --data-dir data \
  --output-dir data/embedding_cache/ppi
