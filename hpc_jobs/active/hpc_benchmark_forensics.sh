#!/usr/bin/env bash
# Run one benchmark-forensics analysis on UCL Grid Engine and publish to SAN.

#$ -S /bin/bash
#$ -l tmem=32G
#$ -l scratch0free=12G
#$ -l tscratch=12G
#$ -l h_rt=24:0:0
#$ -pe smp 1
#$ -j y
#$ -V
#$ -notify

set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1

usage() {
  cat <<'EOF'
Usage: qsub hpc_jobs/active/hpc_benchmark_forensics.sh \
  --dataset cafa3-published|contemporary \
  --output-dir /SAN/.../unique-run-directory

The CAFA3 mode uses Zijian's published nine CSVs, the published DeepGOPlus
pre-projection pickles, and the validated paper-faithful inventory of Zijian's
published embeddings. The contemporary mode uses the completed contemporary
CSVs and pickles plus the finalized hydrated-cache pair-status evidence.
EOF
}

die() { echo "ERROR: $*" >&2; exit 2; }
require_value() { [[ $# -ge 2 && -n "$2" ]] || die "$1 requires a value"; }
git_in_dir() {
  local directory="$1"
  shift
  (cd "$directory" && git "$@")
}

DATASET=""
OUTPUT_DIR=""
SUBMISSION_DIR="${SGE_O_WORKDIR:-$PWD}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dataset) require_value "$@"; DATASET="$2"; shift 2 ;;
    --output-dir) require_value "$@"; OUTPUT_DIR="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

[[ "$DATASET" == "cafa3-published" || "$DATASET" == "contemporary" ]] || \
  die "--dataset must be cafa3-published or contemporary"
[[ "$OUTPUT_DIR" == /SAN/* ]] || die "--output-dir must be an absolute SAN path"
[[ ! -e "$OUTPUT_DIR" ]] || die "Output directory already exists: $OUTPUT_DIR"

JOB_TOKEN="${JOB_ID:-manual_$$}"
WORK="/scratch0/benchmark_forensics_${DATASET}_${JOB_TOKEN}"
FRAMEWORK_DIR="$WORK/Protein-Benchmark-Framework-Dissertation"
SOURCE_DIR="$WORK/source_annotations"
CONFIG_FILE="$WORK/benchmark_forensics.json"
FRAMEWORK_REPO_URL="${FRAMEWORK_REPO_URL:-https://github.com/nadroj0-0/Protein-Benchmark-Framework-Dissertation.git}"
FRAMEWORK_COMMIT="${FRAMEWORK_COMMIT:-}"

CAFA3_BENCHMARK_DIR="${CAFA3_BENCHMARK_DIR:-/SAN/bioinf/bmpfp/reference_artifacts/canonical_cafa3}"
CAFA3_OBO_FILE="${CAFA3_OBO_FILE:-/SAN/bioinf/bmpfp/reference_artifacts/canonical_cafa3/go.obo}"
CAFA3_SOURCE_ARCHIVE="${CAFA3_SOURCE_ARCHIVE:-/SAN/bioinf/bmpfp/reference_artifacts/deepgoplus/data-cafa.tar.gz}"
CAFA3_MODALITY_INVENTORY="${CAFA3_MODALITY_INVENTORY:-/SAN/bioinf/bmpfp/diagnostics/benchmark_forensics/inputs/cafa3_published/embedding_inventory.tsv.gz}"

CONTEMPORARY_BENCHMARK_DIR="${CONTEMPORARY_BENCHMARK_DIR:-/SAN/bioinf/bmpfp/benchmarks/contemporary/2025_01_to_2026_02_supervisor}"
CONTEMPORARY_OBO_FILE="${CONTEMPORARY_OBO_FILE:-/SAN/bioinf/bmpfp/frozen_inputs/ontology/2025-02-06/go-basic.obo}"
CONTEMPORARY_MODALITY_STATUS="${CONTEMPORARY_MODALITY_STATUS:-/SAN/bioinf/bmpfp/embeddings/contemporary/2025_01_to_2026_02_supervisor/finalized_pfp_cache/evidence/pair_status.tsv}"

T0_SPROT_ARCHIVE="${T0_SPROT_ARCHIVE:-/SAN/bioinf/bmpfp/frozen_inputs/uniprot/2025_01/uniprot_sprot-only2025_01.tar.gz}"
T0_TARGET_TREMBL="${T0_TARGET_TREMBL:-/SAN/bioinf/bmpfp/derived_inputs/uniprot/cafa3_target_taxa/2025_01/uniprot_trembl_cafa3_targets.dat.gz}"
T1_SPROT="${T1_SPROT:-/SAN/bioinf/bmpfp/frozen_inputs/uniprot/2026_02/uniprot_sprot.dat.gz}"
T1_TARGET_TREMBL="${T1_TARGET_TREMBL:-/SAN/bioinf/bmpfp/derived_inputs/uniprot/cafa3_target_taxa/2026_02/uniprot_trembl_cafa3_targets.dat.gz}"
T0_SPROT_DAT="$WORK/uniprot_sprot_2025_01.dat.gz"

cleanup() {
  local status=$?
  trap - EXIT
  set +e
  if [[ -d "$WORK" && ! -L "$WORK" && "$WORK" == /scratch0/benchmark_forensics_* ]]; then
    rm -rf -- "$WORK"
  fi
  exit "$status"
}
trap cleanup EXIT
trap 'echo "Received termination signal"; exit 130' INT TERM

for path in "$T0_SPROT_ARCHIVE" "$T0_TARGET_TREMBL" "$T1_SPROT" "$T1_TARGET_TREMBL"; do
  [[ -f "$path" ]] || die "Required taxonomy input is missing: $path"
done
mkdir -p "$WORK" "$SOURCE_DIR" "$(dirname "$OUTPUT_DIR")"

if [[ -z "$FRAMEWORK_COMMIT" ]]; then
  [[ -d "$SUBMISSION_DIR/.git" ]] || \
    die "Submit from a clean framework checkout or pass FRAMEWORK_COMMIT"
  [[ -z "$(git_in_dir "$SUBMISSION_DIR" status --porcelain)" ]] || \
    die "Submission checkout has uncommitted changes"
  FRAMEWORK_COMMIT="$(git_in_dir "$SUBMISSION_DIR" rev-parse HEAD)"
fi
[[ "$FRAMEWORK_COMMIT" =~ ^[0-9a-fA-F]{40}$ ]] || \
  die "FRAMEWORK_COMMIT must be a complete 40-character commit"

echo "Host             : $(hostname)"
echo "Job ID           : ${JOB_ID:-manual}"
echo "Dataset          : $DATASET"
echo "Framework commit : $FRAMEWORK_COMMIT"
echo "Scratch          : $WORK"
echo "SAN output       : $OUTPUT_DIR"
echo "Started          : $(date -Is)"

git clone --no-checkout "$FRAMEWORK_REPO_URL" "$FRAMEWORK_DIR"
git_in_dir "$FRAMEWORK_DIR" checkout --detach "$FRAMEWORK_COMMIT"
[[ "$(git_in_dir "$FRAMEWORK_DIR" rev-parse HEAD)" == "$FRAMEWORK_COMMIT" ]] || \
  die "Scratch checkout does not match FRAMEWORK_COMMIT"

cd "$FRAMEWORK_DIR"
source scripts/reproduction_common.sh
load_framework_paths "$FRAMEWORK_DIR"
add_mmfp_singularity_bind "$WORK"
add_mmfp_singularity_bind /SAN/bioinf/bmpfp
activate_or_create_mmfp_env
PYTHON_BIN="$(command -v python)"

echo "==> Extracting the exact 2025 Swiss-Prot taxonomy source"
tar -xOzf "$T0_SPROT_ARCHIVE" uniprot_sprot.dat.gz > "$T0_SPROT_DAT"
gzip -t "$T0_SPROT_DAT"

if [[ "$DATASET" == "cafa3-published" ]]; then
  for path in "$CAFA3_BENCHMARK_DIR" "$CAFA3_OBO_FILE" "$CAFA3_SOURCE_ARCHIVE" "$CAFA3_MODALITY_INVENTORY"; do
    [[ -e "$path" ]] || die "Required CAFA3 input is missing: $path"
  done
  echo "==> Extracting published DeepGOPlus pre-projection pickles"
  tar -xzf "$CAFA3_SOURCE_ARCHIVE" -C "$SOURCE_DIR" \
    train_data_train.pkl train_data_valid.pkl test_data.pkl
  "$PYTHON_BIN" - "$CONFIG_FILE" "$SOURCE_DIR" "$T0_SPROT_DAT" \
    "$T0_TARGET_TREMBL" "$T1_SPROT" "$T1_TARGET_TREMBL" \
    "$CAFA3_BENCHMARK_DIR" "$CAFA3_OBO_FILE" "$CAFA3_MODALITY_INVENTORY" <<'PY'
import json
import pathlib
import sys

(config, source, t0_sprot, t0_trembl, t1_sprot, t1_trembl,
 benchmark, obo, inventory) = sys.argv[1:]
payload = {
    "schema_version": 1,
    "run_name": "cafa3-published-forensics",
    "top_n": 20,
    "datasets": [{
        "id": "cafa3-published",
        "benchmark_dir": benchmark,
        "obo_file": obo,
        "allow_legacy_singular_protein_header": True,
        "allow_all_zero_rows": False,
        "split_overlap_policy": "disallow",
        "source_annotations": {
            "type": "pfp-pickle-directory",
            "path": source,
            "projection_policy": "Published DeepGOPlus/PFP min_count=50 label-universe projection",
        },
        "taxonomy_sources": [
            {
                "type": "uniprot-dat",
                "path": path,
                "name": name,
                "priority": priority,
            }
            for name, priority, path in (
                ("uniprot-2025-01-trembl", 100, t0_trembl),
                ("uniprot-2025-01-swissprot", 200, t0_sprot),
                ("uniprot-2026-02-trembl", 300, t1_trembl),
                ("uniprot-2026-02-swissprot", 400, t1_sprot),
            )
        ],
        "modality_inventory": {
            "type": "long-table",
            "path": inventory,
            "protein_id_column": "protein_id",
            "modality_column": "modality",
            "states": {
                "published_embedding_available": {
                    "column": "scientifically_eligible",
                    "true_values": ["true"],
                }
            },
        },
    }],
}
pathlib.Path(config).write_text(json.dumps(payload, indent=2) + "\n")
PY
else
  for path in "$CONTEMPORARY_BENCHMARK_DIR" "$CONTEMPORARY_OBO_FILE" "$CONTEMPORARY_MODALITY_STATUS"; do
    [[ -e "$path" ]] || die "Required contemporary input is missing: $path"
  done
  "$PYTHON_BIN" - "$CONFIG_FILE" "$T0_SPROT_DAT" "$T0_TARGET_TREMBL" \
    "$T1_SPROT" "$T1_TARGET_TREMBL" "$CONTEMPORARY_BENCHMARK_DIR" \
    "$CONTEMPORARY_OBO_FILE" "$CONTEMPORARY_MODALITY_STATUS" <<'PY'
import json
import pathlib
import sys

(config, t0_sprot, t0_trembl, t1_sprot, t1_trembl,
 benchmark, obo, modality_status) = sys.argv[1:]
payload = {
    "schema_version": 1,
    "run_name": "contemporary-2025-01-to-2026-02-forensics",
    "top_n": 20,
    "datasets": [{
        "id": "contemporary-2025-01-to-2026-02",
        "benchmark_dir": benchmark,
        "obo_file": obo,
        "allow_legacy_singular_protein_header": False,
        "allow_all_zero_rows": False,
        "split_overlap_policy": "disallow",
        "source_annotations": {
            "type": "pfp-pickle-directory",
            "path": benchmark,
            "projection_policy": "Supervisor-profile min_count=50 label-universe projection",
        },
        "taxonomy_sources": [
            {
                "type": "uniprot-dat",
                "path": path,
                "name": name,
                "priority": priority,
            }
            for name, priority, path in (
                ("uniprot-2025-01-trembl", 100, t0_trembl),
                ("uniprot-2025-01-swissprot", 200, t0_sprot),
                ("uniprot-2026-02-trembl", 300, t1_trembl),
                ("uniprot-2026-02-swissprot", 400, t1_sprot),
            )
        ],
        "modality_inventory": {
            "type": "long-table",
            "path": modality_status,
            "protein_id_column": "protein_id",
            "modality_column": "modality",
            "states": {
                "hydrated_model_input_available": {
                    "column": "state",
                    "true_values": ["accepted"],
                }
            },
        },
    }],
}
pathlib.Path(config).write_text(json.dumps(payload, indent=2) + "\n")
PY
fi

echo "==> Running strict benchmark forensics"
PYTHONPATH="$FRAMEWORK_DIR/benchmark_forensics/src" \
  "$PYTHON_BIN" -m pfp_benchmark_forensics.cli \
    --config "$CONFIG_FILE" \
    --output-dir "$OUTPUT_DIR"

[[ -f "$OUTPUT_DIR/RUN_COMPLETE.json" ]] || die "Completion marker was not published"
[[ -f "$OUTPUT_DIR/output_manifest.json" ]] || die "Output manifest was not published"
echo "Completed         : $(date -Is)"
echo "Published results : $OUTPUT_DIR"
