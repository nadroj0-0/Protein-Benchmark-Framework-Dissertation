#!/bin/bash
# generate_embeddings_prepare_data.sh
# Regenerate CAFA3 splits/labels/sequences from the raw CSVs (Zenodo 7409660).
# Runs prepare_cafa3_data.py UNMODIFIED — the CSVs were already normalised
# (MF 'protein' -> 'proteins') by the dependencies script.
set -euo pipefail

export CAFA3_RAW_DIR="${CAFA3_RAW_DIR:-external/cafa3_raw}"
python scripts/prepare_cafa3_data.py --cafa3-dir "${CAFA3_RAW_DIR}"
