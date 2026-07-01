#!/bin/bash
# reproduce_eval_only.sh
# Eval-only reproduction of the PFP / Hybrid Gated Fusion paper (CAFA3, Table 1).
# Downloads precomputed embeddings + the authors' pretrained checkpoints, then
# evaluates them with CAFA metrics. NO training.
# Verified working on macOS (CPU): BPO 0.601 / CCO 0.706 / MFO 0.702.

set -euo pipefail
LOGFILE="reproduce_eval_only_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOGFILE")
exec 2>&1

# --- 0. Clone (code is at the repo ROOT; README's `cd PFP/MMFP` is wrong:
#        no MMFP subdir exists, and that path errors on Linux). ----------
git clone https://github.com/psipred/PFP.git
cd PFP

# --- 1. Environment: micromamba, Python 3.11 --------------------------
#micromamba create -y -n mmfp python=3.11
#eval "$(micromamba shell hook --shell bash)"
#micromamba activate mmfp
eval "$(/share/apps/miniforge3_mamba/bin/conda shell.bash hook)"
conda create -y -n mmfp python=3.11
conda activate mmfp

# --- 2. Dependencies. requirements.txt is INCOMPLETE: the eval script
#        imports extract_uniprot_text.py (needs `requests`); structure/PPI
#        paths need h5py and fair-esm. Those three are the only missing ones.
#pip install -r requirements.txt
#pip install requests h5py fair-esm
python -m pip install --upgrade pip setuptools wheel
pip install -r requirements.txt --prefer-binary
pip install requests fair-esm
pip install --only-binary=:all: h5py

# --- 3. Data. README's single bundle (mmfp_cafa3_data.tar.gz) 404s; it was
#        split into 5 tarballs on Zenodo record 19498341. Each carries its own
#        data/ (and results/) prefix, so they extract from the REPO ROOT with
#        NO -C flag.
for f in mmfp_embeddings_struct_ppi mmfp_embeddings_prott5 \
         mmfp_embeddings_text_temporal mmfp_checkpoints mmfp_data_splits; do
  wget -c "https://zenodo.org/records/19498341/files/${f}.tar.gz"
  tar -xzf "${f}.tar.gz"
done

# --- 4. Verify required INPUTS exhaustively (eval-only: checkpoints ARE inputs).
echo "==> Verifying required inputs..."
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
  ckpt="results/full_model/fusion_comparison/prott5/${a}/gated_bilinear/best_model.pt"
  [ -f "$ckpt" ] || { echo "MISSING checkpoint: $ckpt"; exit 1; }
done
echo "==> All required inputs present."

# --- 5. Evaluate. Loads the 3 checkpoints, scores all aspects with CAFA
#        metrics. Stays on the cached text branch (prints "Text embedding
#        cache present ... skipping"), so NO network calls during eval.
#        Writes results/full_model_eval/reproduction_summary.csv
python scripts/reproduce_full_model.py

echo "==> Done. Summary: results/full_model_eval/reproduction_summary.csv"
