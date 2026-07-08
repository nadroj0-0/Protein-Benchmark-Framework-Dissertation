#!/bin/bash

set -euo pipefail

LOGFILE="reproduce_$(date +%Y%m%d_%H%M%S).log"

exec > >(tee -a "$LOGFILE")
exec 2>&1


git clone https://github.com/psipred/PFP.git
cd PFP/MMFP

# Create environment
micromamba create -y -n mmfp python=3.11

# Activate environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate mmfp

# Install dependencies
pip install -r requirements.txt


## Data Preparation

# Download precomputed data from Zenodo: https://zenodo.org/records/19498341
wget https://zenodo.org/records/19498341/files/mmfp_cafa3_data.tar.gz

# Extract to data directory
tar -xzf mmfp_cafa3_data.tar.gz -C ./data


### Training from Scratch

# Full model (Table 1)
python train.py \
  --seq-model prott5 \
  --fusion-types gated_bilinear \
  --aspects BPO CCO MFO \
  --use-late-fusion \
  --text-embedding-dir data/embedding_cache/exp_text_embeddings_temporal \
  --output-base results/full_model \
  --seed 42

# Evaluation with CAFA metrics
python scripts/reproduce_full_model.py
