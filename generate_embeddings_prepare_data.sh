#!/bin/bash
# generate_embeddings_prepare_data.sh
# Regenerate CAFA3 splits/labels/sequences from the raw CSVs (Zenodo 7409660).
# Runs prepare_cafa3_data.py UNMODIFIED — CSVs already normalised by deps script.
# Explicit --output-dir so it writes to repo-root data/ regardless of CWD.
set -euo pipefail

export CAFA3_RAW_DIR="${CAFA3_RAW_DIR:-external/cafa3_raw}"
python scripts/prepare_cafa3_data.py \
  --cafa3-dir "${CAFA3_RAW_DIR}" \
  --output-dir data
