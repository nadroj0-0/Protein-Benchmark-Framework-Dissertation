#!/bin/bash
# reproduce_paper.sh
# Eval-only reproduction of the PFP / Hybrid Gated Fusion paper (CAFA3, Table 1).
# Reconstructed from the exact steps debugged on 2026-06-24/25.
# Uses pretrained checkpoints + precomputed embeddings (NO training).

set -euo pipefail


LOGFILE="reproduce_$(date +%Y%m%d_%H%M%S).log"

exec > >(tee -a "$LOGFILE")
exec 2>&1

# ----------------------------------------------------------------------
# 0. Clone the repo (code lives at the ROOT: train.py, scripts/, mmfp/)
#    NOTE: the README's `cd PFP/MMFP` is wrong — there is no MMFP subdir.
#    On Linux (HPC) that path errors; the real working dir is the repo root.
# ----------------------------------------------------------------------
git clone https://github.com/psipred/PFP.git
cd PFP

# ----------------------------------------------------------------------
# 1. Environment: micromamba, Python 3.11
# ----------------------------------------------------------------------
micromamba create -y -n mmfp python=3.11
eval "$(micromamba shell hook --shell bash)"
micromamba activate mmfp

# ----------------------------------------------------------------------
# 2. Dependencies
#    requirements.txt is INCOMPLETE — reproduce_full_model.py imports
#    extract_uniprot_text.py which needs `requests`; the structure/PPI
#    paths pull in h5py and fair-esm. Install all three explicitly.
# ----------------------------------------------------------------------
pip install -r requirements.txt
pip install torch transformers numpy tqdm requests h5py fair-esm

# ----------------------------------------------------------------------
# 3. Data: the README's single bundle (mmfp_cafa3_data.tar.gz) 404s.
#    It was split into 5 tarballs on the same Zenodo record (19498341).
#    They each carry an internal data/ (and results/) prefix, so they
#    extract correctly from the REPO ROOT with no -C flag.
# ----------------------------------------------------------------------
for f in mmfp_embeddings_struct_ppi \
         mmfp_embeddings_prott5 \
         mmfp_embeddings_text_temporal \
         mmfp_checkpoints \
         mmfp_data_splits; do
  wget -c "https://zenodo.org/records/19498341/files/${f}.tar.gz"
  tar -xzf "${f}.tar.gz"
done


# ----------------------------------------------------------------------
# 4. Verify expected layout before running anything
# ----------------------------------------------------------------------
echo "==> Verifying required inputs..."
for d in data/embedding_cache/prott5 \
         data/embedding_cache/IF1 \
         data/embedding_cache/ppi \
         data/embedding_cache/exp_text_embeddings_temporal; do
  [ -d "$d" ] || { echo "MISSING dir: $d"; exit 1; }
done
[ -f data/go.obo ] || { echo "MISSING: data/go.obo"; exit 1; }
for a in BPO CCO MFO; do
  ckpt="results/full_model/fusion_comparison/prott5/${a}/gated_bilinear/best_model.pt"
  [ -f "$ckpt" ] || { echo "MISSING checkpoint: $ckpt"; exit 1; }
done
echo "==> All required inputs present."


# ----------------------------------------------------------------------
# 5. Training from Scratch
# ----------------------------------------------------------------------

# Full model (Table 1)
python train.py \
  --seq-model prott5 \
  --fusion-types gated_bilinear \
  --aspects BPO CCO MFO \
  --use-late-fusion \
  --text-embedding-dir data/embedding_cache/exp_text_embeddings_temporal \
  --output-base results/full_model \
  --seed 42

# ----------------------------------------------------------------------
# 6. Run the eval-only reproduction (loads checkpoints, scores all 3
#    aspects with CAFA metrics; stays on the cached text branch so no
#    network calls). Writes results/full_model_eval/reproduction_summary.csv
# ----------------------------------------------------------------------
python scripts/reproduce_full_model.py

echo "==> Done. Summary: results/full_model_eval/reproduction_summary.csv"
