#!/bin/bash
# generate_embeddings_text.sh
# Text embeddings (768-D PubMedBERT).  README "2. Text Embeddings", Steps 1-2.
# Uses live UniProt/UniSave APIs + CAFA assessment tool.
# extract_uniprot_text.py computes its own paths from __file__ (no path args).
# embed_uniprot_descriptions.py gets explicit --data-dir.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${HERE}/../.." && pwd)"
if [ -f "${REPO_ROOT}/configs/paths.local.sh" ]; then
  # Machine-specific paths are intentionally not committed.
  # shellcheck disable=SC1091
  source "${REPO_ROOT}/configs/paths.local.sh"
fi

export CAFA_ASSESSMENT_DIR="${CAFA_ASSESSMENT_DIR:-external/CAFA_assessment_tool}"
if [ ! -d "${CAFA_ASSESSMENT_DIR}" ]; then
  echo "Missing CAFA assessment tool directory: ${CAFA_ASSESSMENT_DIR}" >&2
  echo "Set CAFA_ASSESSMENT_DIR in configs/paths.local.sh or the environment." >&2
  exit 1
fi

python scripts/extract_uniprot_text.py extract-current
python scripts/extract_uniprot_text.py extract-historical
python scripts/extract_uniprot_text.py prepare-temporal-text
python scripts/embed_uniprot_descriptions.py --data-dir data
