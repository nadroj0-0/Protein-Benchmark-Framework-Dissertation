#!/usr/bin/env bash
# Run PFP's temporal text recipe and store only the CLS vector it consumes.

set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${HERE}/../.." && pwd)"
CAFA_ASSESSMENT_DIR="${CAFA_ASSESSMENT_DIR:-external/CAFA_assessment_tool}"
TEXT_CUTOFF_DATE="${TEXT_CUTOFF_DATE:?Set TEXT_CUTOFF_DATE to YYYY-MM-DD}"
TEXT_HISTORY_WORKERS="${TEXT_HISTORY_WORKERS:-5}"
TEXT_REPORT_DIR="${TEXT_REPORT_DIR:-results/embedding_reports/text}"
CURRENT_CACHE="data/embedding_cache/exp_text_embeddings"
TEMPORAL_CACHE="data/embedding_cache/exp_text_embeddings_temporal"

[[ -d "$CAFA_ASSESSMENT_DIR" ]] || {
  echo "Missing CAFA assessment tool directory: $CAFA_ASSESSMENT_DIR" >&2
  exit 1
}
[[ -f scripts/extract_uniprot_text.py ]] || {
  echo "Run this wrapper from the PFP repository root" >&2
  exit 1
}

mkdir -p "$TEXT_REPORT_DIR"

# A preflight may already have produced temporal CLS files. Move them back to
# the PFP script's fixed current-cache output so its normal resume check sees
# them during the full run.
if [[ -d "$TEMPORAL_CACHE" ]]; then
  [[ ! -e "$CURRENT_CACHE" ]] || {
    echo "Both current and temporal text caches exist; refusing an ambiguous merge" >&2
    exit 1
  }
  mv "$TEMPORAL_CACHE" "$CURRENT_CACHE"
fi

python "${REPO_ROOT}/scripts/embeddings/run_pfp_temporal_text.py" \
  --pfp-root "$PWD" \
  --cafa-assessment-dir "$CAFA_ASSESSMENT_DIR" \
  --cutoff-date "$TEXT_CUTOFF_DATE" \
  --workers "$TEXT_HISTORY_WORKERS"

python scripts/embed_uniprot_descriptions.py --data-dir data &
EMBED_PID=$!
python "${REPO_ROOT}/scripts/embeddings/reduce_text_embeddings_to_cls.py" \
  --directory "$CURRENT_CACHE" \
  --watch-pid "$EMBED_PID" \
  --report "$TEXT_REPORT_DIR/cls_reduction.json" &
REDUCER_PID=$!

EMBED_STATUS=0
REDUCER_STATUS=0
wait "$EMBED_PID" || EMBED_STATUS=$?
wait "$REDUCER_PID" || REDUCER_STATUS=$?

if [[ -d "$CURRENT_CACHE" ]]; then
  [[ ! -e "$TEMPORAL_CACHE" ]] || {
    echo "Temporal cache reappeared unexpectedly during text generation" >&2
    exit 1
  }
  mv "$CURRENT_CACHE" "$TEMPORAL_CACHE"
fi

# The PFP embedding script requires the mixed recipe at its fixed input path.
# Restore the current-description source afterwards so a preflight checkpoint
# can resume into a clean full run rather than appending to the mixed TSV.
CURRENT_DESCRIPTION="data/embedding_cache/uniprot_text/protein_descriptions.tsv"
CURRENT_BACKUP="data/embedding_cache/uniprot_text/temporal_recipe/protein_descriptions_current_before_mixed.tsv"
if [[ -f "$CURRENT_BACKUP" ]]; then
  cp -p "$CURRENT_BACKUP" "$CURRENT_DESCRIPTION"
fi

[[ "$EMBED_STATUS" == "0" ]] || {
  echo "PFP text embedding process failed with status $EMBED_STATUS" >&2
  exit "$EMBED_STATUS"
}
[[ "$REDUCER_STATUS" == "0" ]] || {
  echo "CLS reduction process failed with status $REDUCER_STATUS" >&2
  exit "$REDUCER_STATUS"
}

echo "==> Temporal CLS text embeddings: $TEMPORAL_CACHE"
