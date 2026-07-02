#!/bin/bash
# reproduce_embeddings_retrain_eval.sh  — SINGLE ENTRY POINT
# FULL from-scratch reproduction of PFP / Hybrid Gated Fusion (CAFA3, Table 1):
# clone repo, build env, GENERATE all embeddings from scratch (via the
# generate_embeddings_* sub-scripts), retrain the model on them, then evaluate.
# Unlike the eval-only / download-embeddings paths, this regenerates the
# embeddings rather than downloading the precomputed Zenodo tarballs.
#
# NOTE: training + embedding generation are heavy GPU jobs; impractical on a Mac (CPU).
# KNOWN GAP: the ProtT5 step needs data/proteins.fasta, not yet produced by any step.

set -euo pipefail
LOGFILE="embeddings_verify_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOGFILE")
exec 2>&1

HERE="$(cd "$(dirname "$0")" && pwd)"

# --- 0. Clone (code at repo ROOT; README's `cd PFP/MMFP` is wrong). ----
git clone https://github.com/psipred/PFP.git
cd PFP
REPO="$(pwd)"

# --- 1. Environment: micromamba, Python 3.11 --------------------------
#micromamba create -y -n mmfp python=3.11
#eval "$(micromamba shell hook --shell bash)"
#micromamba activate mmfp

#eval "$(/share/apps/miniforge3_mamba/bin/conda shell.bash hook)"
#conda create -y -n mmfp python=3.11
#conda activate mmfp
#
## --- 2. Dependencies (requirements.txt is incomplete: add requests/h5py/fair-esm) --
##pip install -r requirements.txt
##pip install requests h5py fair-esm
#python -m pip install --upgrade pip setuptools wheel
#pip install -r requirements.txt --prefer-binary
#pip install requests fair-esm
#pip install --only-binary=:all: h5py

eval "$(/share/apps/miniforge3_mamba/bin/conda shell.bash hook)"

ENV_DIR="$HOME/.conda/envs/mmfp"

if [ ! -d "$ENV_DIR" ]; then
    echo "==> Creating Conda environment: mmfp"

    conda create -y -n mmfp python=3.11
    conda activate /home/jsydneyd/.conda/envs/mmfp

    python -m pip install --upgrade pip setuptools wheel
    pip install -r requirements.txt --prefer-binary
    pip install requests fair-esm biopython protobuf sentencepiece torch-geometric "biotite==0.41.2"
    pip install torch-scatter torch-sparse -f https://data.pyg.org/whl/torch-2.6.0+cu124.html
    pip install --only-binary=:all: h5py
    pip install --only-binary=:all: tiktoken

    echo "==> Environment created."
else
    echo "==> Using existing Conda environment."
fi

conda activate /home/jsydneyd/.conda/envs/mmfp

# --- 3. Generate ALL embeddings from scratch (sub-orchestrator; CWD = repo root) ---
bash "${HERE}/generate_embeddings_run_all.sh" || \
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
  wget -c "https://zenodo.org/records/19498341/files/${f}.tar.gz"
  tar -xzf "${f}.tar.gz"
done

cd "$REPO"

echo "==> Comparing generated embeddings against published embeddings..."
python "${HERE}/compare_embeddings.py" \
  --root "$REPO" \
  --out-csv results/embedding_comparison.csv \
  --out-json results/embedding_comparison_summary.json