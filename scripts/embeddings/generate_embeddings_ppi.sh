#!/bin/bash
# generate_embeddings_ppi.sh
# PPI embeddings (512-D STRING).  README "1. PPI Embeddings".
# Explicit --data-dir/--output-dir so it doesn't depend on the ./data default.
# CAFA3 test targets use internal T... identifiers. When CAFA3_ID_MAPPING is
# supplied, pass it explicitly so the extractor can resolve those targets to
# the UniProt entry names used by STRING. UniProt-native benchmarks omit it.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${HERE}/../.." && pwd)"
if [ -f "${REPO_ROOT}/configs/paths.local.sh" ]; then
  # Machine-specific paths are intentionally not committed.
  # shellcheck disable=SC1091
  source "${REPO_ROOT}/configs/paths.local.sh"
fi

export STRING_H5_FILE="${STRING_H5_FILE:-external/string/protein.network.embeddings.v12.0.h5}"
export STRING_ALIAS_FILE="${STRING_ALIAS_FILE:-external/string/protein.aliases.v12.0.txt}"
export CAFA_ASSESSMENT_DIR="${CAFA_ASSESSMENT_DIR:-external/CAFA_assessment_tool}"
PPI_EXTRACT_SCRIPT="${PPI_EXTRACT_SCRIPT:-scripts/extract_ppi_embeddings.py}"
CAFA3_ID_MAPPING="${CAFA3_ID_MAPPING:-}"
PYTHON_BIN="${PYTHON_BIN:-python}"
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
if [ ! -f "${PPI_EXTRACT_SCRIPT}" ]; then
  echo "Missing PPI extraction script: ${PPI_EXTRACT_SCRIPT}" >&2
  exit 1
fi
if [ -n "${CAFA3_ID_MAPPING}" ] && [ ! -f "${CAFA3_ID_MAPPING}" ]; then
  echo "Missing CAFA3 ID mapping file: ${CAFA3_ID_MAPPING}" >&2
  exit 1
fi

command=(
  "${PYTHON_BIN}" "${PPI_EXTRACT_SCRIPT}"
  --string-h5 "${STRING_H5_FILE}"
  --string-alias "${STRING_ALIAS_FILE}"
  --cafa-assessment-dir "${CAFA_ASSESSMENT_DIR}"
  --data-dir data
  --output-dir data/embedding_cache/ppi
)
if [ -n "${CAFA3_ID_MAPPING}" ]; then
  command+=(--cafa3-id-mapping "${CAFA3_ID_MAPPING}")
fi
printf 'PPI extractor command:'
printf ' %q' "${command[@]}"
printf '\n'
"${command[@]}"
