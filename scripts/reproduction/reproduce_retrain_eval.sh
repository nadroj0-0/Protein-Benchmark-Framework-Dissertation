#!/bin/bash
# reproduce_retrain_eval.sh
# Train-then-eval reproduction of the PFP / Hybrid Gated Fusion paper (CAFA3, Table 1).
# Downloads precomputed embeddings, RETRAINS the model from them, then evaluates.
# The authors' published checkpoints are preserved (not overwritten) for comparison.
#
# NOTE: training is a heavy GPU job. On a Mac (CPU) the train step is impractical
# to finish — this script is logically correct but section 5 is meant for a GPU.

set -euo pipefail
LOGFILE="reproduce_retrain_eval_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOGFILE")
exec 2>&1

HERE="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${HERE}/../.." && pwd)"
source "${REPO_ROOT}/scripts/reproduction_common.sh"
load_framework_paths "${REPO_ROOT}"
cd "${REPO_ROOT}"

# --- 0. Clone (code at repo ROOT; README's `cd PFP/MMFP` is wrong). ----
clone_or_reuse_pfp

# --- 1. Environment: micromamba, Python 3.11 --------------------------
#micromamba create -y -n mmfp python=3.11
#eval "$(micromamba shell hook --shell bash)"
#micromamba activate mmfp
#eval "$(/share/apps/miniforge3_mamba/bin/conda shell.bash hook)"
#conda create -y -n mmfp python=3.11
#conda activate mmfp
#
## --- 2. Dependencies (same incomplete-requirements fix as eval-only) --
##pip install -r requirements.txt
##pip install requests h5py fair-esm
#python -m pip install --upgrade pip setuptools wheel
#pip install -r requirements.txt --prefer-binary
#pip install requests fair-esm
#pip install --only-binary=:all: h5py

activate_or_create_mmfp_env

# --- 3. Data (same 5 tarballs, extract from repo root, no -C) ----------
for f in mmfp_embeddings_struct_ppi mmfp_embeddings_prott5 \
         mmfp_embeddings_text_temporal mmfp_checkpoints mmfp_data_splits; do
  wget -c "https://zenodo.org/records/19498341/files/${f}.tar.gz"
  tar -xzf "${f}.tar.gz"
done

# Preserve the authors' published checkpoints BEFORE we train our own,
# so training doesn't overwrite them and we can compare later.
if [ -d results/full_model ]; then
  mv results/full_model results/full_model_published
fi

# --- 4. Verify TRAINING INPUTS exhaustively (embeddings, splits, ontology).
#        Checkpoints are an OUTPUT of training here, so NOT checked at this stage.
echo "==> Verifying training inputs..."
for d in data/embedding_cache/prott5 data/embedding_cache/IF1 \
         data/embedding_cache/ppi data/embedding_cache/exp_text_embeddings_temporal; do
  [ -d "$d" ] || { echo "MISSING dir: $d"; exit 1; }
done
[ -f data/go.obo ] || { echo "MISSING: data/go.obo"; exit 1; }
for a in BPO CCO MFO; do
  for split in train valid test; do
    [ -f "data/${a}_${split}_names.npy" ]      || { echo "MISSING: ${a}_${split}_names.npy"; exit 1; }
    [ -f "data/${a}_${split}_labels.npz" ]     || { echo "MISSING: ${a}_${split}_labels.npz"; exit 1; }
    [ -f "data/${a}_${split}_sequences.json" ] || { echo "MISSING: ${a}_${split}_sequences.json"; exit 1; }
  done
  [ -f "data/${a}_go_terms.json" ] || { echo "MISSING: ${a}_go_terms.json"; exit 1; }
  [ -f "data/${a}_ia.txt" ]        || { echo "MISSING: ${a}_ia.txt"; exit 1; }
done
echo "==> Training inputs present."

# --- 5. TRAIN from the precomputed embeddings.
#        --num-workers 0 matches the README's training command.
#        Writes fresh checkpoints to results/full_model/.
python train.py \
  --seq-model prott5 \
  --fusion-types gated_bilinear \
  --aspects BPO CCO MFO \
  --use-late-fusion \
  --text-embedding-dir data/embedding_cache/exp_text_embeddings_temporal \
  --output-base results/full_model \
  --num-workers 0 \
  --seed 42

# Verify training actually produced checkpoints where eval expects them.
echo "==> Verifying training outputs..."
for a in BPO CCO MFO; do
  ckpt="results/full_model/fusion_comparison/prott5/${a}/gated_bilinear/best_model.pt"
  [ -f "$ckpt" ] || { echo "MISSING trained checkpoint: $ckpt"; exit 1; }
done
echo "==> Training produced all 3 checkpoints."

# --- 6. EVAL the freshly-trained model with CAFA metrics.
#        Reads results/full_model/... (now YOUR checkpoints). Writes
#        results/full_model_eval/reproduction_summary.csv
python scripts/reproduce_full_model.py

echo "==> Done."
echo "    Your retrain summary: results/full_model_eval/reproduction_summary.csv"
echo "    Authors' checkpoints preserved in: results/full_model_published/"
