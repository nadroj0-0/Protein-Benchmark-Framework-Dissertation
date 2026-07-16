#!/bin/bash
# generate_embeddings_structure.sh
# Structure embeddings (512-D ESM-IF1).  README "3. Structure Embeddings".
# Step 1 downloads AlphaFold PDBs; Step 2 runs ESM-IF1 over them.
# Explicit paths so check_alphafold_coverage.py does NOT use its ../data default
# (which assumes CWD=scripts/ and breaks from repo root).
# DEVICE overridable for local testing.
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
DEVICE="${DEVICE:-cuda}"
IF1_EXTRACT_SCRIPT="${IF1_EXTRACT_SCRIPT:-scripts/extract_esm_if1_embeddings.py}"
IF1_PYTHON_BIN="${IF1_PYTHON_BIN:-python}"
IF1_PYTHONPATH="${IF1_PYTHONPATH:-}"
ALPHAFOLD_ACQUISITION_MODE="${ALPHAFOLD_ACQUISITION_MODE:-pfp}"
ALPHAFOLD_PERSISTENT_CACHE_DIR="${ALPHAFOLD_PERSISTENT_CACHE_DIR:-}"
ALPHAFOLD_API_WORKERS="${ALPHAFOLD_API_WORKERS:-8}"
ALPHAFOLD_DOWNLOAD_WORKERS="${ALPHAFOLD_DOWNLOAD_WORKERS:-8}"
ALPHAFOLD_PREFETCH_REPORT="${ALPHAFOLD_PREFETCH_REPORT:-results/embedding_reports/structure/alphafold_prefetch.json}"

[ -f "$IF1_EXTRACT_SCRIPT" ] || {
  echo "Missing IF1 extractor: $IF1_EXTRACT_SCRIPT" >&2
  exit 1
}
command -v "$IF1_PYTHON_BIN" >/dev/null 2>&1 || {
  echo "Missing IF1 Python entrypoint: $IF1_PYTHON_BIN" >&2
  exit 1
}

if [ "$ALPHAFOLD_ACQUISITION_MODE" = "framework-bounded" ]; then
  [ -n "$ALPHAFOLD_PERSISTENT_CACHE_DIR" ] || {
    echo "Set ALPHAFOLD_PERSISTENT_CACHE_DIR for framework-bounded acquisition" >&2
    exit 1
  }
  python "${REPO_ROOT}/scripts/embeddings/prefetch_alphafold_structures.py" \
    --pfp-root "$PWD" \
    --cafa-assessment-dir "${CAFA_ASSESSMENT_DIR}" \
    --data-dir data \
    --persistent-cache-dir "$ALPHAFOLD_PERSISTENT_CACHE_DIR" \
    --workspace-pdb-dir data/alphafold_structures \
    --coverage-report data/alphafold_coverage_results.txt \
    --report "$ALPHAFOLD_PREFETCH_REPORT" \
    --api-workers "$ALPHAFOLD_API_WORKERS" \
    --download-workers "$ALPHAFOLD_DOWNLOAD_WORKERS"
elif [ "$ALPHAFOLD_ACQUISITION_MODE" = "pfp" ]; then
  python scripts/check_alphafold_coverage.py \
    --cafa-assessment-dir "${CAFA_ASSESSMENT_DIR}" \
    --data-dir data \
    --pdb-output-dir data/alphafold_structures \
    --output-file data/alphafold_coverage_results.txt
else
  echo "Unknown ALPHAFOLD_ACQUISITION_MODE: $ALPHAFOLD_ACQUISITION_MODE" >&2
  exit 1
fi

if [ -n "$IF1_PYTHONPATH" ]; then
  SINGULARITYENV_PYTHONPATH="$IF1_PYTHONPATH" \
    MMFP_PYTHONPATH="$IF1_PYTHONPATH" \
    "$IF1_PYTHON_BIN" "$IF1_EXTRACT_SCRIPT" \
      --pdb_dir data/alphafold_structures \
      --output_dir data/embedding_cache/IF1 \
      --pooling mean \
      --device "$DEVICE"
else
  "$IF1_PYTHON_BIN" "$IF1_EXTRACT_SCRIPT" \
      --pdb_dir data/alphafold_structures \
      --output_dir data/embedding_cache/IF1 \
      --pooling mean \
      --device "$DEVICE"
fi
