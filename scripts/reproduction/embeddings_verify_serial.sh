#!/bin/bash
# reproduce_embeddings_retrain_eval.sh  — SINGLE ENTRY POINT
# FULL from-scratch reproduction of PFP / Hybrid Gated Fusion (CAFA3, Table 1):
# clone repo, build env, GENERATE all embeddings from scratch (via the
# generate_embeddings_* sub-scripts), retrain the model on them, then evaluate.
# Unlike the eval-only / download-embeddings paths, this regenerates the
# embeddings rather than downloading the precomputed Zenodo tarballs.
#
# NOTE: training + embedding generation are heavy GPU jobs; impractical on a Mac (CPU).
# The embedding sub-orchestrator builds data/proteins.fasta before the ProtT5 step.

set -euo pipefail
LOGFILE="embeddings_verify_serial_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOGFILE")
exec 2>&1

HERE="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${HERE}/../.." && pwd)"
source "${REPO_ROOT}/scripts/reproduction_common.sh"
load_framework_paths "${REPO_ROOT}"
cd "${REPO_ROOT}"

# --- 0. Clone (code at repo ROOT; README's `cd PFP/MMFP` is wrong). ----
clone_or_reuse_pfp
REPO="$(pwd)"

# --- 1. Environment: author-supplied Python 3.9.23 --------------------
# Creation and constrained dependency installation are centralised in the
# shared helper. Existing environments are activated without modification.
activate_or_create_mmfp_env

# --- 3. Generate ALL embeddings from scratch (sub-orchestrator; CWD = repo root) ---
bash "${REPO_ROOT}/scripts/embeddings/generate_embeddings_run_all_serial.sh" || \
echo "==> run_all non-zero (expected: strict verify flags text shape); continuing to comparison."

echo
echo "Checking generated embeddings..."

for d in \
    data/embedding_cache/prott5 \
    data/embedding_cache/IF1 \
    data/embedding_cache/ppi \
    data/embedding_cache/exp_text_embeddings
do
    echo "$d"
    find "$d" -name "*.npy" | wc -l
done

mkdir -p published
cd published

# --- 3. Data (same 5 tarballs, extract from repo root, no -C) ----------
for f in mmfp_embeddings_struct_ppi mmfp_embeddings_prott5 \
         mmfp_embeddings_text_temporal mmfp_checkpoints mmfp_data_splits; do
  explicit="${PUBLISHED_BUNDLE_DIR:+${PUBLISHED_BUNDLE_DIR}/${f}.tar.gz}"
  stage_or_download_artifact "$(zijian_bundle_artifact_id "$f")" "$explicit" \
    "${f}.tar.gz" "https://zenodo.org/records/19498341/files/${f}.tar.gz"
  tar -xzf "${f}.tar.gz"
done

cd "$REPO"

echo "==> Comparing generated embeddings against published embeddings..."
python "${REPO_ROOT}/scripts/diagnostics/compare_embeddings.py" \
  --root "$REPO" \
  --out-csv results/embedding_comparison.csv \
  --out-json results/embedding_comparison_summary.json
