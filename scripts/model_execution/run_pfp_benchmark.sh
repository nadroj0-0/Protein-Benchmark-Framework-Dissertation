#!/usr/bin/env bash
# Prepare, validate, train and/or evaluate PFP on any conforming nine-CSV benchmark.
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FRAMEWORK_ROOT="$(cd "$HERE/../.." && pwd)"
PYTHON_BIN="${PYTHON_BIN:-$(command -v python3 || command -v python || true)}"
[[ -n "$PYTHON_BIN" ]] || { echo "ERROR: Python is required" >&2; exit 2; }

usage() {
  cat <<'EOF'
Usage:
  bash scripts/model_execution/run_pfp_benchmark.sh \
    --benchmark-id ID --benchmark-dir DIR --obo-file FILE --pfp-root DIR \
    --work-dir DIR --output-dir DIR --config FILE \
    --execution-mode prepare-only|eval-only|train-eval \
    [--embedding-cache-root DIR] [--checkpoint-root DIR] \
    [--modality-mode full|sequence-only] [--aspect BPO|CCO|MFO] \
    [--seed N] [--num-workers N] [--ia-file-dir DIR] \
    [--reference-data-dir DIR] [--reference-source-archive FILE] \
    [--expected-metrics FILE] \
    [--benchmark-evidence FILE] \
    [--embedding-evidence FILE] [--require-embedding-evidence] \
    [--reference-tolerance FLOAT] [--require-reference-match] \
    [--expected-pfp-commit SHA] [--allow-unversioned-pfp] \
    [--allow-dirty-framework]

This entrypoint never downloads inputs and never modifies the PFP checkout.
Every path is explicit so the same workflow can run locally or under an HPC
wrapper. Repeating --aspect selects a subset; the default is all three.
EOF
}

die() { echo "ERROR: $*" >&2; exit 2; }
require_value() { [[ $# -ge 2 && -n "$2" ]] || die "$1 requires a value"; }

BENCHMARK_ID=""
BENCHMARK_DIR=""
OBO_FILE=""
PFP_ROOT=""
WORK_DIR=""
OUTPUT_DIR=""
CONFIG=""
EXECUTION_MODE=""
EMBEDDING_CACHE_ROOT=""
CHECKPOINT_ROOT=""
MODALITY_MODE="full"
SEED=42
NUM_WORKERS=0
IA_FILE_DIR=""
REFERENCE_DATA_DIR=""
REFERENCE_SOURCE_ARCHIVE=""
EXPECTED_METRICS=""
REFERENCE_TOLERANCE=""
REQUIRE_REFERENCE_MATCH=0
REQUIRE_EMBEDDING_EVIDENCE=0
EXPECTED_PFP_COMMIT="1e04fd6d6d3c40458fd41ec1a881ed6e24de768e"
ALLOW_UNVERSIONED_PFP=0
ALLOW_DIRTY_FRAMEWORK=0
ASPECTS=()
BENCHMARK_EVIDENCE=()
EMBEDDING_EVIDENCE=()
ASPECT_COUNT=0
BENCHMARK_EVIDENCE_COUNT=0
EMBEDDING_EVIDENCE_COUNT=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --benchmark-id) require_value "$@"; BENCHMARK_ID="$2"; shift 2 ;;
    --benchmark-dir) require_value "$@"; BENCHMARK_DIR="$2"; shift 2 ;;
    --obo-file) require_value "$@"; OBO_FILE="$2"; shift 2 ;;
    --pfp-root) require_value "$@"; PFP_ROOT="$2"; shift 2 ;;
    --work-dir) require_value "$@"; WORK_DIR="$2"; shift 2 ;;
    --output-dir) require_value "$@"; OUTPUT_DIR="$2"; shift 2 ;;
    --config) require_value "$@"; CONFIG="$2"; shift 2 ;;
    --execution-mode) require_value "$@"; EXECUTION_MODE="$2"; shift 2 ;;
    --embedding-cache-root) require_value "$@"; EMBEDDING_CACHE_ROOT="$2"; shift 2 ;;
    --checkpoint-root) require_value "$@"; CHECKPOINT_ROOT="$2"; shift 2 ;;
    --modality-mode) require_value "$@"; MODALITY_MODE="$2"; shift 2 ;;
    --aspect) require_value "$@"; ASPECTS+=("$2"); ASPECT_COUNT=$((ASPECT_COUNT + 1)); shift 2 ;;
    --seed) require_value "$@"; SEED="$2"; shift 2 ;;
    --num-workers) require_value "$@"; NUM_WORKERS="$2"; shift 2 ;;
    --ia-file-dir) require_value "$@"; IA_FILE_DIR="$2"; shift 2 ;;
    --reference-data-dir) require_value "$@"; REFERENCE_DATA_DIR="$2"; shift 2 ;;
    --reference-source-archive) require_value "$@"; REFERENCE_SOURCE_ARCHIVE="$2"; shift 2 ;;
    --expected-metrics) require_value "$@"; EXPECTED_METRICS="$2"; shift 2 ;;
    --benchmark-evidence) require_value "$@"; BENCHMARK_EVIDENCE+=("$2"); BENCHMARK_EVIDENCE_COUNT=$((BENCHMARK_EVIDENCE_COUNT + 1)); shift 2 ;;
    --embedding-evidence) require_value "$@"; EMBEDDING_EVIDENCE+=("$2"); EMBEDDING_EVIDENCE_COUNT=$((EMBEDDING_EVIDENCE_COUNT + 1)); shift 2 ;;
    --require-embedding-evidence) REQUIRE_EMBEDDING_EVIDENCE=1; shift ;;
    --reference-tolerance) require_value "$@"; REFERENCE_TOLERANCE="$2"; shift 2 ;;
    --require-reference-match) REQUIRE_REFERENCE_MATCH=1; shift ;;
    --expected-pfp-commit) require_value "$@"; EXPECTED_PFP_COMMIT="$2"; shift 2 ;;
    --allow-unversioned-pfp) ALLOW_UNVERSIONED_PFP=1; shift ;;
    --allow-dirty-framework) ALLOW_DIRTY_FRAMEWORK=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

[[ -n "$BENCHMARK_ID" ]] || die "--benchmark-id is required"
[[ -n "$BENCHMARK_DIR" ]] || die "--benchmark-dir is required"
[[ -n "$OBO_FILE" ]] || die "--obo-file is required"
[[ -n "$PFP_ROOT" ]] || die "--pfp-root is required"
[[ -n "$WORK_DIR" ]] || die "--work-dir is required"
[[ -n "$OUTPUT_DIR" ]] || die "--output-dir is required"
[[ -n "$CONFIG" ]] || die "--config is required"
[[ -n "$EXECUTION_MODE" ]] || die "--execution-mode is required"
[[ "$EXECUTION_MODE" =~ ^(prepare-only|eval-only|train-eval)$ ]] || die "Invalid --execution-mode"
[[ "$MODALITY_MODE" =~ ^(full|sequence-only)$ ]] || die "Invalid --modality-mode"
if [[ "$ALLOW_DIRTY_FRAMEWORK" == "1" ]] && \
   [[ "$EXECUTION_MODE" != "prepare-only" || "$ALLOW_UNVERSIONED_PFP" != "1" ]]; then
  die "--allow-dirty-framework is restricted to prepare-only unversioned test fixtures"
fi
[[ "$SEED" =~ ^[0-9]+$ ]] || die "--seed must be a non-negative integer"
[[ "$NUM_WORKERS" =~ ^[0-9]+$ ]] || die "--num-workers must be a non-negative integer"
[[ -d "$BENCHMARK_DIR" ]] || die "Benchmark directory does not exist: $BENCHMARK_DIR"
[[ -f "$OBO_FILE" ]] || die "GO OBO file does not exist: $OBO_FILE"
[[ -f "$PFP_ROOT/train.py" && -f "$PFP_ROOT/scripts/prepare_cafa3_data.py" ]] || \
  die "PFP root is not a compatible checkout: $PFP_ROOT"
if git -C "$PFP_ROOT" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  OBSERVED_PFP_COMMIT="$(git -C "$PFP_ROOT" rev-parse HEAD)"
  [[ "$OBSERVED_PFP_COMMIT" == "$EXPECTED_PFP_COMMIT" ]] || \
    die "PFP commit mismatch: expected $EXPECTED_PFP_COMMIT, found $OBSERVED_PFP_COMMIT"
  [[ -z "$(git -C "$PFP_ROOT" status --porcelain --untracked-files=no)" ]] || \
    die "PFP has tracked modifications; use an immutable checkout"
  while IFS= read -r untracked; do
    case "$untracked" in
      *.py|*.pyc|*.pyo|*.so|*.pth|*.egg-info/*)
        die "PFP contains an untracked executable/importable file: $untracked" ;;
    esac
  done < <(git -C "$PFP_ROOT" ls-files --others --exclude-standard)
elif [[ "$ALLOW_UNVERSIONED_PFP" != "1" ]]; then
  die "PFP root is not a Git checkout; use --allow-unversioned-pfp only for fixtures"
fi
[[ -f "$CONFIG" ]] || die "Run config does not exist: $CONFIG"
if git -C "$FRAMEWORK_ROOT" rev-parse --is-inside-work-tree >/dev/null 2>&1 && \
   [[ -n "$(git -C "$FRAMEWORK_ROOT" status --porcelain)" ]] && \
   [[ "$ALLOW_DIRTY_FRAMEWORK" != "1" ]]; then
  die "Framework checkout has uncommitted changes; commit them before a model run"
fi
if [[ -n "$IA_FILE_DIR" ]]; then
  [[ -d "$IA_FILE_DIR" ]] || die "IA file directory does not exist: $IA_FILE_DIR"
fi
[[ ! -e "$WORK_DIR" ]] || die "Work directory must not already exist: $WORK_DIR"
[[ ! -e "$OUTPUT_DIR" ]] || die "Output directory must not already exist: $OUTPUT_DIR"
if [[ "$EXECUTION_MODE" != "prepare-only" ]]; then
  [[ -d "$EMBEDDING_CACHE_ROOT" ]] || die "--embedding-cache-root is required and must exist"
fi
if [[ "$EXECUTION_MODE" == "eval-only" ]]; then
  [[ -d "$CHECKPOINT_ROOT" ]] || die "--checkpoint-root is required and must exist"
  [[ -d "$REFERENCE_DATA_DIR" ]] || \
    die "--reference-data-dir is required for eval-only checkpoint/term binding"
fi
if [[ -n "$REFERENCE_SOURCE_ARCHIVE" ]]; then
  [[ -f "$REFERENCE_SOURCE_ARCHIVE" ]] || \
    die "--reference-source-archive does not exist: $REFERENCE_SOURCE_ARCHIVE"
  [[ -n "$REFERENCE_DATA_DIR" ]] || \
    die "--reference-source-archive requires --reference-data-dir"
fi
if [[ "$REQUIRE_REFERENCE_MATCH" == "1" ]]; then
  [[ -f "$EXPECTED_METRICS" ]] || \
    die "--require-reference-match also requires --expected-metrics"
  [[ "$EXECUTION_MODE" != "prepare-only" ]] || \
    die "Reference metrics cannot be checked in prepare-only mode"
fi

mkdir -p "$WORK_DIR" "$WORK_DIR/logs" "$WORK_DIR/reports"
DATA_DIR="$WORK_DIR/data"
MODEL_OUTPUT="$WORK_DIR/model_output"
EVALUATION_OUTPUT="$WORK_DIR/evaluation"
OUTPUT_STAGE="${OUTPUT_DIR}.staging-$$"
[[ ! -e "$OUTPUT_STAGE" ]] || die "Output staging path already exists: $OUTPUT_STAGE"
trap 'rm -rf -- "$OUTPUT_STAGE"' EXIT

ASPECT_ARGS=()
if [[ "$ASPECT_COUNT" -eq 0 ]]; then
  ASPECTS=(BPO CCO MFO)
fi
SEEN_ASPECTS=" "
for aspect in "${ASPECTS[@]}"; do
  [[ "$aspect" =~ ^(BPO|CCO|MFO)$ ]] || die "Invalid aspect: $aspect"
  case "$SEEN_ASPECTS" in
    *" $aspect "*) die "Duplicate aspect: $aspect" ;;
  esac
  SEEN_ASPECTS="${SEEN_ASPECTS}${aspect} "
  ASPECT_ARGS+=(--aspect "$aspect")
done

PREPARE_COMMAND=(
  "$PYTHON_BIN" "$HERE/prepare_pfp_benchmark.py"
  --benchmark-dir "$BENCHMARK_DIR"
  --data-dir "$DATA_DIR"
  --obo-file "$OBO_FILE"
  --pfp-root "$PFP_ROOT"
  --config "$CONFIG"
  --report "$WORK_DIR/reports/preparation.json"
  --log-dir "$WORK_DIR/logs"
)
if [[ -n "$REFERENCE_DATA_DIR" ]]; then
  PREPARE_COMMAND+=(--reference-data-dir "$REFERENCE_DATA_DIR")
fi
if [[ -n "$REFERENCE_SOURCE_ARCHIVE" ]]; then
  PREPARE_COMMAND+=(--reference-source-archive "$REFERENCE_SOURCE_ARCHIVE")
fi
if [[ "$BENCHMARK_EVIDENCE_COUNT" -gt 0 ]]; then
  for evidence in "${BENCHMARK_EVIDENCE[@]}"; do PREPARE_COMMAND+=(--benchmark-evidence "$evidence"); done
fi
echo "==> Validate and prepare benchmark"
"${PREPARE_COMMAND[@]}" >"$WORK_DIR/logs/preparation_driver.log" 2>&1

if [[ "$EXECUTION_MODE" != "prepare-only" ]]; then
  echo "==> Exhaustively validate embedding cache"
  EMBEDDING_COMMAND=(
    "$PYTHON_BIN" "$HERE/validate_pfp_embedding_cache.py"
    --data-dir "$DATA_DIR"
    --cache-root "$EMBEDDING_CACHE_ROOT"
    --config "$CONFIG"
    --mode "$MODALITY_MODE"
    --report "$WORK_DIR/reports/embedding_cache.json"
    --issues-tsv "$WORK_DIR/reports/embedding_cache_issues.tsv"
    --preparation-report "$WORK_DIR/reports/preparation.json"
    "${ASPECT_ARGS[@]}"
  )
  if [[ -n "$IA_FILE_DIR" ]]; then EMBEDDING_COMMAND+=(--ia-file-dir "$IA_FILE_DIR"); fi
  if [[ "$EMBEDDING_EVIDENCE_COUNT" -gt 0 ]]; then
    for evidence in "${EMBEDDING_EVIDENCE[@]}"; do EMBEDDING_COMMAND+=(--embedding-evidence "$evidence"); done
  fi
  if [[ "$REQUIRE_EMBEDDING_EVIDENCE" == "1" ]]; then EMBEDDING_COMMAND+=(--require-embedding-evidence); fi
  "${EMBEDDING_COMMAND[@]}" \
    >"$WORK_DIR/logs/embedding_validation.log" 2>&1

  mkdir -p "$DATA_DIR/embedding_cache"
  config_path() {
    "$PYTHON_BIN" -c 'import json, pathlib, sys; c=json.load(open(sys.argv[1])); p=pathlib.Path(c["modalities"][sys.argv[2]]["directory"]); print(p if p.is_absolute() else pathlib.Path(sys.argv[3])/p)' \
      "$CONFIG" "$1" "$EMBEDDING_CACHE_ROOT"
  }
  SEQUENCE_DIR="$(config_path sequence)"
  TEXT_DIR="$(config_path text)"
  STRUCTURE_DIR="$(config_path structure)"
  PPI_DIR="$(config_path ppi)"
  ln -s "$SEQUENCE_DIR" "$DATA_DIR/embedding_cache/prott5"
  ln -s "$STRUCTURE_DIR" "$DATA_DIR/embedding_cache/IF1"
  ln -s "$PPI_DIR" "$DATA_DIR/embedding_cache/ppi"

  if [[ "$EXECUTION_MODE" == "eval-only" ]]; then
    echo "==> Evaluate existing checkpoints"
    EVAL_COMMAND=(
      "$PYTHON_BIN" "$HERE/evaluate_pfp_checkpoints.py"
      --pfp-root "$PFP_ROOT"
      --data-dir "$DATA_DIR"
      --cache-root "$EMBEDDING_CACHE_ROOT"
      --obo-file "$OBO_FILE"
      --checkpoint-root "$CHECKPOINT_ROOT"
      --output-dir "$EVALUATION_OUTPUT"
      --config "$CONFIG"
      --mode "$MODALITY_MODE"
      --num-workers "$NUM_WORKERS"
      --seed "$SEED"
      "${ASPECT_ARGS[@]}"
    )
    if [[ -n "$IA_FILE_DIR" ]]; then EVAL_COMMAND+=(--ia-file-dir "$IA_FILE_DIR"); fi
    "${EVAL_COMMAND[@]}" >"$WORK_DIR/logs/evaluation.log" 2>&1
  else
    echo "==> Train and evaluate fresh PFP models"
    mkdir -p "$MODEL_OUTPUT"
    EMPTY_DIR="$WORK_DIR/empty_modality"
    mkdir -p "$EMPTY_DIR"
    for aspect in "${ASPECTS[@]}"; do
      TRAIN_COMMAND=(
        "$PYTHON_BIN" "$PFP_ROOT/train.py"
        --single
        --seq-model prott5
        --fusion-types gated_bilinear
        --aspects "$aspect"
        --data-dir "$DATA_DIR"
        --obo-file "$OBO_FILE"
        --output-base "$MODEL_OUTPUT"
        --text-embedding-dir "$TEXT_DIR"
        --use-late-fusion
        --late-output-mode hybrid
        --aux-loss-weight 0.8
        --modality-dropout 0.1
        --num-workers "$NUM_WORKERS"
        --seed "$SEED"
      )
      if [[ "$MODALITY_MODE" == "sequence-only" ]]; then
        TRAIN_COMMAND+=(--text-embedding-dir "$EMPTY_DIR" --no-struct --no-ppi)
      fi
      if [[ -n "$IA_FILE_DIR" ]]; then TRAIN_COMMAND+=(--ia-file-dir "$IA_FILE_DIR"); fi
      "${TRAIN_COMMAND[@]}" >"$WORK_DIR/logs/train_${aspect}.log" 2>&1
      RESULT_DIR="$MODEL_OUTPUT/fusion_comparison/prott5/$aspect/gated_bilinear"
      [[ -s "$RESULT_DIR/best_model.pt" && -s "$RESULT_DIR/results.json" ]] || \
        die "PFP did not produce a complete result for $aspect"
    done
    echo "==> Strictly re-evaluate fresh checkpoints"
    EVAL_COMMAND=(
      "$PYTHON_BIN" "$HERE/evaluate_pfp_checkpoints.py"
      --pfp-root "$PFP_ROOT"
      --data-dir "$DATA_DIR"
      --cache-root "$EMBEDDING_CACHE_ROOT"
      --obo-file "$OBO_FILE"
      --checkpoint-root "$MODEL_OUTPUT"
      --output-dir "$EVALUATION_OUTPUT"
      --config "$CONFIG"
      --mode "$MODALITY_MODE"
      --num-workers "$NUM_WORKERS"
      --seed "$SEED"
      "${ASPECT_ARGS[@]}"
    )
    if [[ -n "$IA_FILE_DIR" ]]; then EVAL_COMMAND+=(--ia-file-dir "$IA_FILE_DIR"); fi
    "${EVAL_COMMAND[@]}" >"$WORK_DIR/logs/strict_post_training_evaluation.log" 2>&1
  fi

  echo "==> Revalidate embedding bytes after model execution"
  POST_EMBEDDING_COMMAND=("${EMBEDDING_COMMAND[@]}")
  for ((index = 0; index < ${#POST_EMBEDDING_COMMAND[@]}; index++)); do
    if [[ "${POST_EMBEDDING_COMMAND[$index]}" == "--report" ]]; then
      POST_EMBEDDING_COMMAND[$((index + 1))]="$WORK_DIR/reports/embedding_cache_post.json"
    elif [[ "${POST_EMBEDDING_COMMAND[$index]}" == "--issues-tsv" ]]; then
      POST_EMBEDDING_COMMAND[$((index + 1))]="$WORK_DIR/reports/embedding_cache_post_issues.tsv"
    fi
  done
  "${POST_EMBEDDING_COMMAND[@]}" >"$WORK_DIR/logs/embedding_validation_post.log" 2>&1
  "$PYTHON_BIN" - "$WORK_DIR/reports/embedding_cache.json" \
    "$WORK_DIR/reports/embedding_cache_post.json" <<'PY'
import json
import sys

before = json.load(open(sys.argv[1], encoding="utf-8"))
after = json.load(open(sys.argv[2], encoding="utf-8"))
for modality, summary in before["modalities"].items():
    if summary["valid_content_sha256"] != after["modalities"][modality]["valid_content_sha256"]:
        raise SystemExit(f"Embedding cache changed during model execution: {modality}")
PY
fi

FRAMEWORK_COMMIT="$(git -C "$FRAMEWORK_ROOT" rev-parse HEAD 2>/dev/null || echo unknown)"
if [[ "$FRAMEWORK_COMMIT" != "unknown" ]] && \
   [[ -n "$(git -C "$FRAMEWORK_ROOT" status --porcelain 2>/dev/null)" ]]; then
  FRAMEWORK_COMMIT="${FRAMEWORK_COMMIT}-dirty"
fi
PFP_COMMIT="$(git -C "$PFP_ROOT" rev-parse HEAD 2>/dev/null || echo unversioned-fixture)"
SUMMARY_COMMAND=(
  "$PYTHON_BIN" "$HERE/summarize_pfp_benchmark_run.py"
  --benchmark-id "$BENCHMARK_ID"
  --execution-mode "$EXECUTION_MODE"
  --modality-mode "$MODALITY_MODE"
  --seed "$SEED"
  --framework-commit "$FRAMEWORK_COMMIT"
  --pfp-commit "$PFP_COMMIT"
  --preparation-report "$WORK_DIR/reports/preparation.json"
  --output-json "$WORK_DIR/reports/run_report.json"
  --output-md "$WORK_DIR/reports/run_report.md"
  "${ASPECT_ARGS[@]}"
)
if [[ "$EXECUTION_MODE" != "prepare-only" ]]; then
  SUMMARY_COMMAND+=(--embedding-report "$WORK_DIR/reports/embedding_cache.json")
  SUMMARY_COMMAND+=(--embedding-post-report "$WORK_DIR/reports/embedding_cache_post.json")
fi
if [[ "$EXECUTION_MODE" == "train-eval" ]]; then
  SUMMARY_COMMAND+=(--result-root "$MODEL_OUTPUT" --evaluation-root "$EVALUATION_OUTPUT")
elif [[ "$EXECUTION_MODE" == "eval-only" ]]; then
  SUMMARY_COMMAND+=(--evaluation-root "$EVALUATION_OUTPUT")
fi
if [[ -n "$EXPECTED_METRICS" ]]; then SUMMARY_COMMAND+=(--expected-metrics "$EXPECTED_METRICS"); fi
if [[ -n "$REFERENCE_TOLERANCE" ]]; then SUMMARY_COMMAND+=(--reference-tolerance "$REFERENCE_TOLERANCE"); fi
if [[ "$REQUIRE_REFERENCE_MATCH" == "1" ]]; then SUMMARY_COMMAND+=(--require-reference-match); fi
"${SUMMARY_COMMAND[@]}" >"$WORK_DIR/logs/summary.log" 2>&1

echo "==> Publish completed run"
mkdir -p "$OUTPUT_STAGE/logs" "$OUTPUT_STAGE/reports" "$OUTPUT_STAGE/prepared_data"
cp -p "$WORK_DIR/logs/"* "$OUTPUT_STAGE/logs/"
cp -p "$WORK_DIR/reports/"* "$OUTPUT_STAGE/reports/"
find "$DATA_DIR" -maxdepth 1 -type f -exec cp -p {} "$OUTPUT_STAGE/prepared_data/" \;
cp -p "$CONFIG" "$OUTPUT_STAGE/run_config.json"
if [[ -d "$MODEL_OUTPUT" ]]; then cp -a "$MODEL_OUTPUT" "$OUTPUT_STAGE/models"; fi
if [[ -d "$EVALUATION_OUTPUT" ]]; then cp -a "$EVALUATION_OUTPUT" "$OUTPUT_STAGE/evaluation"; fi
"$PYTHON_BIN" "$HERE/manage_output_manifest.py" write --root "$OUTPUT_STAGE"
MANIFEST_SHA256="$("$PYTHON_BIN" -c 'import hashlib, pathlib, sys; print(hashlib.sha256(pathlib.Path(sys.argv[1]).read_bytes()).hexdigest())' "$OUTPUT_STAGE/output_manifest.json")"
"$PYTHON_BIN" -c 'import json, pathlib, sys; pathlib.Path(sys.argv[1]).write_text(json.dumps({"complete": True, "manifest": "output_manifest.json", "manifest_sha256": sys.argv[2], "run_report": "reports/run_report.json"}, indent=2)+"\n")' \
  "$OUTPUT_STAGE/WORKFLOW_COMPLETE.json" "$MANIFEST_SHA256"
"$PYTHON_BIN" "$HERE/manage_output_manifest.py" verify --root "$OUTPUT_STAGE"
mv "$OUTPUT_STAGE" "$OUTPUT_DIR"
trap - EXIT
echo "Completed PFP benchmark run: $OUTPUT_DIR"
