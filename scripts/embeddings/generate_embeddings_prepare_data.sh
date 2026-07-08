#!/bin/bash
# generate_embeddings_prepare_data.sh
# Regenerate CAFA3 splits/labels/sequences from the raw CSVs (Zenodo 7409660).
# Runs prepare_cafa3_data.py UNMODIFIED — CSVs already normalised by deps script.
# Explicit --output-dir so it writes to repo-root data/ regardless of CWD.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${HERE}/../.." && pwd)"
if [ -f "${REPO_ROOT}/configs/paths.local.sh" ]; then
  # Machine-specific paths are intentionally not committed.
  # shellcheck disable=SC1091
  source "${REPO_ROOT}/configs/paths.local.sh"
fi

export CAFA3_RAW_DIR="${CAFA3_RAW_DIR:-external/cafa3_raw}"
if [ ! -d "${CAFA3_RAW_DIR}" ]; then
  echo "Missing CAFA3 raw CSV directory: ${CAFA3_RAW_DIR}" >&2
  echo "Set CAFA3_RAW_DIR in configs/paths.local.sh or the environment." >&2
  exit 1
fi

python scripts/prepare_cafa3_data.py \
  --cafa3-dir "${CAFA3_RAW_DIR}" \
  --output-dir data
