#!/bin/bash
# generate_embeddings_text.sh
# Text embeddings (768-D PubMedBERT).  README "2. Text Embeddings", Steps 1-2.
# Uses live UniProt/UniSave APIs + CAFA assessment tool.
# extract_uniprot_text.py computes its own paths from __file__ (no path args).
# embed_uniprot_descriptions.py gets explicit --data-dir.
set -euo pipefail

export CAFA_ASSESSMENT_DIR="${CAFA_ASSESSMENT_DIR:-external/CAFA_assessment_tool}"

python scripts/extract_uniprot_text.py extract-current
python scripts/extract_uniprot_text.py extract-historical
python scripts/extract_uniprot_text.py prepare-temporal-text
python scripts/embed_uniprot_descriptions.py --data-dir data
