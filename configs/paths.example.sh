#!/usr/bin/env bash
# Example local path configuration for Protein-Benchmark-Framework.
#
# Copy this file to configs/paths.local.sh and edit it for your machine.
# Do not commit paths.local.sh; it is intentionally machine-specific.

# Optional catalogue of already-downloaded large artifacts. Every workflow
# still accepts its individual explicit path overrides. Resolution order is:
# explicit path, catalogue entry, then the workflow's original download.
# export ARTIFACT_CATALOG="/path/to/artifact_paths.tsv"

# Local/reference PFP checkout.
export PFP_DIR="${PFP_DIR:-$HOME/PFP}"

# CAFA assessment tool used by the original evaluation/embedding workflow.
# Keep this checkout outside PFP so PFP source authentication remains isolated.
export CAFA_ASSESSMENT_DIR="${CAFA_ASSESSMENT_DIR:-$HOME/CAFA_assessment_tool}"

# Canonical nine CAFA3 CSVs used to prepare PFP-compatible splits.
export CAFA3_RAW_DIR="${CAFA3_RAW_DIR:-$HOME/cafa3_csv/cafa3_raw}"

# Contemporary raw database archive, kept outside Git.
export PROTEIN_DATABASES_DIR="${PROTEIN_DATABASES_DIR:-$HOME/pfp_inputs}"

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

# Upstream PFP location when a workflow needs to clone it explicitly.
export PFP_GIT_URL="${PFP_GIT_URL:-https://github.com/psipred/PFP.git}"
export PFP_CLONE_DIR="${PFP_CLONE_DIR:-PFP}"
