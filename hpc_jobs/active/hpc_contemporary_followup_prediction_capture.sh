#!/usr/bin/env bash
# Capture validation and test predictions from the accepted corrected contemporary full model.

#$ -S /bin/bash
#$ -l tmem=32G
#$ -l scratch0free=220G
#$ -l tscratch=220G
#$ -l h_rt=24:0:0
#$ -q gpu.q@zeus1.local,gpu.q@zeus2.local
#$ -l gpu=true
#$ -pe gpu 1
#$ -j y
#$ -N ct25_capture
#$ -V
#$ -notify

set -Eeuo pipefail
export PYTHONDONTWRITEBYTECODE=1

usage() {
  cat <<'EOF'
Usage: qsub hpc_jobs/active/hpc_contemporary_followup_prediction_capture.sh \
  --output-dir /SAN/.../unique-capture-directory

Runs inference only, without retraining. It captures validation and test arrays
from the accepted full-model checkpoints trained with the corrected 2025-03-08
text cache and paper-faithful PPI policy. The fresh test capture must exactly
match the previously accepted test prediction content before publication.
EOF
}

die() { echo "ERROR: $*" >&2; exit 2; }
require_value() { [[ $# -ge 2 && -n "$2" ]] || die "$1 requires a value"; }
git_in_dir() {
  local directory="$1"
  shift
  (cd "$directory" && git "$@")
}

OUTPUT_DIR=""
SUBMISSION_DIR="${SGE_O_WORKDIR:-$PWD}"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --output-dir) require_value "$@"; OUTPUT_DIR="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

[[ "$OUTPUT_DIR" == /SAN/* ]] || die "--output-dir must be an absolute SAN path"
[[ ! -e "$OUTPUT_DIR" ]] || die "Output directory already exists: $OUTPUT_DIR"

SOURCE_RUN="${SOURCE_RUN:-/SAN/bioinf/bmpfp/model_runs/contemporary/2025_01_to_2026_02_supervisor/variants/text-cutoff-2025-03-08__ppi-paper-faithful/full/7118745_20260728_164527}"
CACHE_ARCHIVE="${CACHE_ARCHIVE:-/SAN/bioinf/bmpfp/embeddings/contemporary/2025_01_to_2026_02_supervisor/variants/text-cutoff-2025-03-08__ppi-paper-faithful/finalized_pfp_cache/contemporary_embedding_cache.tar.gz}"
OBO_FILE="${OBO_FILE:-/SAN/bioinf/bmpfp/frozen_inputs/ontology/2025-02-06/go-basic.obo}"
SOURCE_TEST_MANIFEST="$SOURCE_RUN/evaluation/prediction_artifacts/prediction_artifact_manifest.json"
SOURCE_IA_DIR="$SOURCE_RUN/evaluation/prediction_artifacts"
SOURCE_CONFIG="$SOURCE_RUN/run_config.json"
SOURCE_PREPARATION_REPORT="$SOURCE_RUN/reports/preparation.json"
SOURCE_EMBEDDING_REPORT="$SOURCE_RUN/reports/embedding_cache.json"

for path in \
  "$SOURCE_RUN/WORKFLOW_COMPLETE.json" \
  "$SOURCE_TEST_MANIFEST" \
  "$SOURCE_CONFIG" \
  "$SOURCE_PREPARATION_REPORT" \
  "$SOURCE_EMBEDDING_REPORT" \
  "$CACHE_ARCHIVE" \
  "$OBO_FILE"; do
  [[ -f "$path" ]] || die "Required immutable input is missing: $path"
done
[[ -d "$SOURCE_RUN/models" ]] || die "Checkpoint directory is missing"
[[ -d "$SOURCE_RUN/prepared_data" ]] || die "Prepared-data directory is missing"
for aspect in BPO CCO MFO; do
  [[ -f "$SOURCE_RUN/models/fusion_comparison/prott5/$aspect/gated_bilinear/best_model.pt" ]] || \
    die "Saved full-model checkpoint is missing for $aspect"
  [[ -f "$SOURCE_IA_DIR/${aspect}_ia.txt" ]] || die "Accepted IA file is missing for $aspect"
done

JOB_TOKEN="${JOB_ID:-manual_$$}"
WORK="/scratch0/ct25_followup_capture_${JOB_TOKEN}"
FRAMEWORK_DIR="$WORK/Protein-Benchmark-Framework-Dissertation"
PFP_DIR="$WORK/PFP"
CACHE_ROOT="$WORK/embedding_cache"
DATA_DIR="$WORK/prepared_data"
CHECKPOINT_ROOT="$WORK/models"
PAIR_ROOT="$WORK/capture_pair"
LOG_DIR="$WORK/logs"
PUBLISH_STAGE="${OUTPUT_DIR}.staging-${JOB_TOKEN}"
PUBLISH_LOCK="${OUTPUT_DIR}.publish-lock"
FRAMEWORK_REPO_URL="${FRAMEWORK_REPO_URL:-https://github.com/nadroj0-0/Protein-Benchmark-Framework-Dissertation.git}"
FRAMEWORK_COMMIT="${FRAMEWORK_COMMIT:-}"
PFP_REPO_URL="${PFP_REPO_URL:-https://github.com/psipred/PFP.git}"
PFP_COMMIT="${PFP_COMMIT:-1e04fd6d6d3c40458fd41ec1a881ed6e24de768e}"
LOCK_HELD=0

cleanup() {
  local status=$?
  trap - EXIT
  set +e
  if [[ "$LOCK_HELD" == "1" && -d "$PUBLISH_LOCK" && ! -L "$PUBLISH_LOCK" ]]; then
    rmdir -- "$PUBLISH_LOCK"
  fi
  [[ ! -d "$PUBLISH_STAGE" || -L "$PUBLISH_STAGE" ]] || rm -rf -- "$PUBLISH_STAGE"
  [[ ! -d "$WORK" || -L "$WORK" || "$WORK" != /scratch0/ct25_followup_capture_* ]] || \
    rm -rf -- "$WORK"
  exit "$status"
}
trap cleanup EXIT
trap 'echo "Received termination signal"; exit 130' INT TERM

[[ ! -e "$WORK" ]] || die "Scratch path already exists: $WORK"
[[ ! -e "$PUBLISH_STAGE" ]] || die "Publication stage already exists: $PUBLISH_STAGE"
[[ ! -e "$PUBLISH_LOCK" ]] || die "Publication lock already exists: $PUBLISH_LOCK"
mkdir -p "$WORK" "$PAIR_ROOT" "$LOG_DIR" "$(dirname "$OUTPUT_DIR")"

if [[ -z "$FRAMEWORK_COMMIT" ]]; then
  [[ -d "$SUBMISSION_DIR/.git" ]] || die "Submit from a clean framework checkout"
  [[ -z "$(git_in_dir "$SUBMISSION_DIR" status --porcelain)" ]] || \
    die "Submission checkout has uncommitted changes"
  FRAMEWORK_COMMIT="$(git_in_dir "$SUBMISSION_DIR" rev-parse HEAD)"
fi
[[ "$FRAMEWORK_COMMIT" =~ ^[0-9a-fA-F]{40}$ ]] || \
  die "FRAMEWORK_COMMIT must be a full commit"

echo "Host             : $(hostname)"
echo "Job ID           : ${JOB_ID:-manual}"
echo "Source run       : $SOURCE_RUN"
echo "Embedding policy : text-cutoff-2025-03-08__ppi-paper-faithful"
echo "Model mode       : full"
echo "Framework commit : $FRAMEWORK_COMMIT"
echo "SAN output       : $OUTPUT_DIR"
echo "Started          : $(date -Is)"

git clone --no-checkout "$FRAMEWORK_REPO_URL" "$FRAMEWORK_DIR"
git_in_dir "$FRAMEWORK_DIR" checkout --detach "$FRAMEWORK_COMMIT"
git clone --no-checkout "$PFP_REPO_URL" "$PFP_DIR"
git_in_dir "$PFP_DIR" checkout --detach "$PFP_COMMIT"
[[ "$(git_in_dir "$FRAMEWORK_DIR" rev-parse HEAD)" == "$FRAMEWORK_COMMIT" ]] || \
  die "Framework checkout differs from requested commit"
[[ "$(git_in_dir "$PFP_DIR" rev-parse HEAD)" == "$PFP_COMMIT" ]] || \
  die "PFP checkout differs from requested commit"

cd "$FRAMEWORK_DIR"
source scripts/reproduction_common.sh
load_framework_paths "$FRAMEWORK_DIR"
add_mmfp_singularity_bind "$WORK"
add_mmfp_singularity_bind /SAN/bioinf/bmpfp
activate_or_create_mmfp_env
PYTHON_BIN="$(command -v python)"
"$PYTHON_BIN" -c 'import torch; assert torch.cuda.is_available(); assert torch.cuda.device_count() == 1; print(torch.cuda.get_device_name(0))'

echo "==> Extracting the immutable corrected paper-faithful cache"
"$PYTHON_BIN" scripts/embeddings/manage_embedding_archive.py extract \
  --archive "$CACHE_ARCHIVE" \
  --output-cache-root "$CACHE_ROOT" \
  --config "$SOURCE_CONFIG" \
  --report "$WORK/cache_extraction.json" \
  >"$LOG_DIR/cache_extraction.log" 2>&1

echo "==> Staging accepted prepared data and full-model checkpoints"
cp -a "$SOURCE_RUN/prepared_data" "$DATA_DIR"
cp -a "$SOURCE_RUN/models" "$CHECKPOINT_ROOT"

COMMON_ARGS=(
  --pfp-root "$PFP_DIR"
  --data-dir "$DATA_DIR"
  --cache-root "$CACHE_ROOT"
  --obo-file "$OBO_FILE"
  --checkpoint-root "$CHECKPOINT_ROOT"
  --config "$SOURCE_CONFIG"
  --mode full
  --aspect BPO --aspect CCO --aspect MFO
  --ia-file-dir "$SOURCE_IA_DIR"
  --benchmark-id contemporary-2025_01-to-2026_02-supervisor
  --framework-commit "$FRAMEWORK_COMMIT"
  --pfp-commit "$PFP_COMMIT"
  --preparation-report "$SOURCE_PREPARATION_REPORT"
  --embedding-report "$SOURCE_EMBEDDING_REPORT"
  --num-workers 0
  --seed 42
)

for split in valid test; do
  echo "==> Capturing $split predictions"
  mkdir -p "$PAIR_ROOT/$split"
  "$PYTHON_BIN" scripts/model_execution/evaluate_pfp_checkpoints.py \
    "${COMMON_ARGS[@]}" \
    --evaluation-split "$split" \
    --output-dir "$PAIR_ROOT/$split/evaluation" \
    --prediction-artifact-dir "$PAIR_ROOT/$split/evaluation/prediction_artifacts" \
    >"$LOG_DIR/${split}_capture.log" 2>&1
done

echo "==> Proving the fresh test capture exactly reproduces the accepted test arrays"
"$PYTHON_BIN" - "$SOURCE_TEST_MANIFEST" \
  "$PAIR_ROOT/test/evaluation/prediction_artifacts/prediction_artifact_manifest.json" \
  "$PAIR_ROOT/test_reproduction_check.json" <<'PY'
import json
import pathlib
import sys

old_path, new_path, output_path = map(pathlib.Path, sys.argv[1:])
old = json.loads(old_path.read_text())
new = json.loads(new_path.read_text())
checks = {
    "benchmark_id": (old["benchmark_id"], new["benchmark_id"]),
    "mode": (old["mode"], new["mode"]),
    "seed": (old["seed"], new["seed"]),
    "selected_aspects": (old["selected_aspects"], new["selected_aspects"]),
    "config_sha256": (old["config"]["sha256"], new["config"]["sha256"]),
    "obo_sha256": (old["obo"]["sha256"], new["obo"]["sha256"]),
    "pfp_commit": (old["provenance"]["pfp_commit"], new["provenance"]["pfp_commit"]),
    "benchmark_fingerprint": (
        old["provenance"]["benchmark_fingerprint"],
        new["provenance"]["benchmark_fingerprint"],
    ),
    "source_csv_sha256": (
        old["provenance"]["source_csv_sha256"],
        new["provenance"]["source_csv_sha256"],
    ),
}
aspect_fields = (
    "shape",
    "scores_content_sha256",
    "truth_content_sha256",
    "protein_ids_sha256",
    "go_terms_sha256",
    "checkpoint_sha256",
    "ia_file_sha256",
    "canonical_cafa_metrics",
)
for aspect in old["selected_aspects"]:
    for field in aspect_fields:
        checks[f"{aspect}.{field}"] = (
            old["aspects"][aspect][field],
            new["aspects"][aspect][field],
        )
mismatches = {
    key: {"accepted": left, "fresh": right}
    for key, (left, right) in checks.items()
    if left != right
}
report = {
    "schema_version": 1,
    "status": "passed" if not mismatches else "failed",
    "accepted_test_manifest": str(old_path),
    "fresh_test_manifest": str(new_path),
    "checks": len(checks),
    "mismatches": mismatches,
}
output_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
if mismatches:
    raise SystemExit(f"Fresh test capture differs from accepted arrays: {sorted(mismatches)}")
PY

echo "==> Publishing paired captures atomically"
mkdir -p "$PUBLISH_STAGE"
cp -a "$PAIR_ROOT/valid" "$PAIR_ROOT/test" "$PUBLISH_STAGE/"
cp -p "$PAIR_ROOT/test_reproduction_check.json" "$WORK/cache_extraction.json" \
  "$PUBLISH_STAGE/"
cp -a "$LOG_DIR" "$PUBLISH_STAGE/logs"
"$PYTHON_BIN" scripts/model_execution/manage_output_manifest.py write \
  --root "$PUBLISH_STAGE" --include-nested-control-files
MANIFEST_SHA256="$(
  "$PYTHON_BIN" -c 'import hashlib,pathlib,sys; print(hashlib.sha256(pathlib.Path(sys.argv[1]).read_bytes()).hexdigest())' \
    "$PUBLISH_STAGE/output_manifest.json"
)"
"$PYTHON_BIN" -c 'import json,pathlib,sys; pathlib.Path(sys.argv[1]).write_text(json.dumps({"complete":True,"analysis_kind":"validation_test_prediction_capture","embedding_policy":"text-cutoff-2025-03-08__ppi-paper-faithful","mode":"full","manifest":"output_manifest.json","manifest_sha256":sys.argv[2]},indent=2)+"\n")' \
  "$PUBLISH_STAGE/WORKFLOW_COMPLETE.json" "$MANIFEST_SHA256"
"$PYTHON_BIN" scripts/model_execution/manage_output_manifest.py verify \
  --root "$PUBLISH_STAGE" --include-nested-control-files
mkdir "$PUBLISH_LOCK"
LOCK_HELD=1
[[ ! -e "$OUTPUT_DIR" ]] || die "Output directory appeared during publication"
mv -T "$PUBLISH_STAGE" "$OUTPUT_DIR"
rmdir "$PUBLISH_LOCK"
LOCK_HELD=0
echo "Published paired validation/test captures: $OUTPUT_DIR"
