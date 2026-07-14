#!/usr/bin/env bash
# Example local path configuration for Protein-Benchmark-Framework.
#
# Copy this file to configs/paths.local.sh and edit it for your machine.
# Do not commit paths.local.sh; it is intentionally machine-specific.

# Local/reference PFP checkout. On Jordan's Mac this may currently resolve via
# ~/PFP, which is a compatibility symlink into the supplementary archive.
export PFP_DIR="${PFP_DIR:-$HOME/PFP}"

# CAFA assessment tool used by the original evaluation/embedding workflow.
export CAFA_ASSESSMENT_DIR="${CAFA_ASSESSMENT_DIR:-$HOME/CAFA_assessment_tool}"

# Canonical nine CAFA3 CSVs used to prepare PFP-compatible splits.
export CAFA3_RAW_DIR="${CAFA3_RAW_DIR:-$HOME/cafa3_csv/cafa3_raw}"

# Contemporary raw database archive, usually on HPC/SAN rather than in git.
export PROTEIN_DATABASES_DIR="${PROTEIN_DATABASES_DIR:-$HOME/protein_databases}"

# STRING inputs for PPI embedding extraction.
export STRING_H5_FILE="${STRING_H5_FILE:-$HOME/protein_databases/string/protein.network.embeddings.v12.0.h5}"
export STRING_ALIAS_FILE="${STRING_ALIAS_FILE:-$HOME/protein_databases/string/protein.aliases.v12.0.txt}"

# Optional overrides for generate_embeddings_dependencies.sh. Leave unset to
# use the current PFP checkout's external/ and data/ directories.
# export PFP_EXTERNAL_DIR="/path/to/PFP/external"
# export PFP_DATA_DIR="/path/to/PFP/data"
# export PFP_CAFA3_RAW_DIR="/path/to/PFP/external/cafa3_raw"
# export PFP_STRING_DIR="/path/to/PFP/external/string"
# export DEPENDENCY_ENV="/path/to/PFP/external/dependency_env.sh"

# PFP/MMFP Python environment. Use a path or environment name appropriate for
# the active machine; cluster scripts may override this.
export MMFP_ENV="${MMFP_ENV:-mmfp}"
export MMFP_ENV_DIR="${MMFP_ENV_DIR:-$HOME/.conda/envs/$MMFP_ENV}"
export MMFP_PYTHON="${MMFP_PYTHON:-3.9.23}"
export CONDA_EXE="${CONDA_EXE:-/share/apps/miniforge3_mamba/bin/conda}"

# Root-level reproduction wrappers clone/reuse PFP here.
export PFP_GIT_URL="${PFP_GIT_URL:-https://github.com/psipred/PFP.git}"
export PFP_CLONE_DIR="${PFP_CLONE_DIR:-PFP}"

# Scratch directory for verify_csv.sh. The script writes downloaded Zenodo CSVs,
# generated split artefacts, and a patched temporary prepare_cafa3_data.py here.
export VERIFY_CSV_WORKDIR="${VERIFY_CSV_WORKDIR:-$HOME/mmfp_csv_verify}"
