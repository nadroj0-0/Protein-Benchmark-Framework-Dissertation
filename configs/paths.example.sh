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

# PFP/MMFP Python environment. Use a path or environment name appropriate for
# the active machine; cluster scripts may override this.
export MMFP_ENV="${MMFP_ENV:-mmfp}"

# Scratch directory for verify_csv.sh. The script writes downloaded Zenodo CSVs,
# generated split artefacts, and a patched temporary prepare_cafa3_data.py here.
export VERIFY_CSV_WORKDIR="${VERIFY_CSV_WORKDIR:-$HOME/mmfp_csv_verify}"
