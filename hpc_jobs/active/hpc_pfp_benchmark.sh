#!/usr/bin/env bash
# Generic one-GPU UCL Grid Engine wrapper for prepared PFP-compatible benchmarks.

#$ -l tmem=32G
#$ -l scratch0free=250G
#$ -l tscratch=250G
#$ -l h_rt=96:0:0
#$ -l gpu=true
#$ -pe gpu 1
#$ -j y
#$ -N pfp_benchmark
#$ -V
#$ -notify

set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1
PYTHON_BIN="${PYTHON_BIN:-$(command -v python3 || command -v python || true)}"

usage() {
  cat <<'EOF'
Usage: qsub hpc_jobs/active/hpc_pfp_benchmark.sh \
  --benchmark-id ID --benchmark-dir DIR \
  (--embedding-cache-root DIR | --embedding-cache-archive FILE) \
  --obo-file FILE --results-root DIR \
  --execution-mode eval-only|train-eval \
  [--checkpoint-root DIR] [--config FILE] \
  [--modality-mode full|sequence-only] [--cache-staging copy|direct] \
  [--aspect BPO|CCO|MFO] [--seed N] [--num-workers N] \
  [--ia-file-dir DIR] [--expected-metrics FILE] \
  [--reference-data-dir DIR] [--reference-source-archive FILE] \
  [--benchmark-evidence FILE] \
  [--embedding-evidence FILE] [--require-embedding-evidence] \
  [--reference-tolerance FLOAT] [--require-reference-match]

The benchmark, ontology and cache are existing prerequisites. This wrapper
does not download or generate embeddings. A cache archive is safely extracted
into job-owned scratch. A cache directory is copied by default, while
--cache-staging direct reads that directory in place. Repeat --aspect to select
multiple aspects.
EOF
}

die() { echo "ERROR: $*" >&2; exit 2; }
require_value() { [[ $# -ge 2 && -n "$2" ]] || die "$1 requires a value"; }

BENCHMARK_ID=""
BENCHMARK_DIR=""
CACHE_ROOT=""
CACHE_ARCHIVE=""
OBO_FILE=""
RESULTS_ROOT=""
EXECUTION_MODE=""
CHECKPOINT_ROOT=""
CLI_CONFIG=""
MODALITY_MODE="full"
CACHE_STAGING="copy"
SEED=42
NUM_WORKERS=0
IA_FILE_DIR=""
EXPECTED_METRICS=""
REFERENCE_DATA_DIR=""
REFERENCE_SOURCE_ARCHIVE=""
REFERENCE_TOLERANCE=""
REQUIRE_REFERENCE_MATCH=0
REQUIRE_EMBEDDING_EVIDENCE=0
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
    --embedding-cache-root) require_value "$@"; CACHE_ROOT="$2"; shift 2 ;;
    --embedding-cache-archive) require_value "$@"; CACHE_ARCHIVE="$2"; shift 2 ;;
    --obo-file) require_value "$@"; OBO_FILE="$2"; shift 2 ;;
    --results-root) require_value "$@"; RESULTS_ROOT="$2"; shift 2 ;;
    --execution-mode) require_value "$@"; EXECUTION_MODE="$2"; shift 2 ;;
    --checkpoint-root) require_value "$@"; CHECKPOINT_ROOT="$2"; shift 2 ;;
    --config) require_value "$@"; CLI_CONFIG="$2"; shift 2 ;;
    --modality-mode) require_value "$@"; MODALITY_MODE="$2"; shift 2 ;;
    --cache-staging) require_value "$@"; CACHE_STAGING="$2"; shift 2 ;;
    --aspect) require_value "$@"; ASPECTS+=("$2"); ASPECT_COUNT=$((ASPECT_COUNT + 1)); shift 2 ;;
    --seed) require_value "$@"; SEED="$2"; shift 2 ;;
    --num-workers) require_value "$@"; NUM_WORKERS="$2"; shift 2 ;;
    --ia-file-dir) require_value "$@"; IA_FILE_DIR="$2"; shift 2 ;;
    --expected-metrics) require_value "$@"; EXPECTED_METRICS="$2"; shift 2 ;;
    --reference-data-dir) require_value "$@"; REFERENCE_DATA_DIR="$2"; shift 2 ;;
    --reference-source-archive) require_value "$@"; REFERENCE_SOURCE_ARCHIVE="$2"; shift 2 ;;
    --benchmark-evidence) require_value "$@"; BENCHMARK_EVIDENCE+=("$2"); BENCHMARK_EVIDENCE_COUNT=$((BENCHMARK_EVIDENCE_COUNT + 1)); shift 2 ;;
    --embedding-evidence) require_value "$@"; EMBEDDING_EVIDENCE+=("$2"); EMBEDDING_EVIDENCE_COUNT=$((EMBEDDING_EVIDENCE_COUNT + 1)); shift 2 ;;
    --require-embedding-evidence) REQUIRE_EMBEDDING_EVIDENCE=1; shift ;;
    --reference-tolerance) require_value "$@"; REFERENCE_TOLERANCE="$2"; shift 2 ;;
    --require-reference-match) REQUIRE_REFERENCE_MATCH=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

SUBMISSION_DIR="${SGE_O_WORKDIR:-$PWD}"
resolve_submission_path() {
  if [[ -z "$1" ]]; then
    printf '\n'
    return
  fi
  case "$1" in
    /*) printf '%s\n' "$1" ;;
    *) printf '%s/%s\n' "$SUBMISSION_DIR" "$1" ;;
  esac
}

BENCHMARK_DIR="$(resolve_submission_path "$BENCHMARK_DIR")"
CACHE_ROOT="$(resolve_submission_path "$CACHE_ROOT")"
CACHE_ARCHIVE="$(resolve_submission_path "$CACHE_ARCHIVE")"
OBO_FILE="$(resolve_submission_path "$OBO_FILE")"
RESULTS_ROOT="$(resolve_submission_path "$RESULTS_ROOT")"
[[ -z "$CHECKPOINT_ROOT" ]] || CHECKPOINT_ROOT="$(resolve_submission_path "$CHECKPOINT_ROOT")"
[[ -z "$CLI_CONFIG" ]] || CLI_CONFIG="$(resolve_submission_path "$CLI_CONFIG")"
[[ -z "$IA_FILE_DIR" ]] || IA_FILE_DIR="$(resolve_submission_path "$IA_FILE_DIR")"
[[ -z "$EXPECTED_METRICS" ]] || EXPECTED_METRICS="$(resolve_submission_path "$EXPECTED_METRICS")"
[[ -z "$REFERENCE_DATA_DIR" ]] || REFERENCE_DATA_DIR="$(resolve_submission_path "$REFERENCE_DATA_DIR")"
[[ -z "$REFERENCE_SOURCE_ARCHIVE" ]] || REFERENCE_SOURCE_ARCHIVE="$(resolve_submission_path "$REFERENCE_SOURCE_ARCHIVE")"
for index in "${!BENCHMARK_EVIDENCE[@]}"; do
  BENCHMARK_EVIDENCE[$index]="$(resolve_submission_path "${BENCHMARK_EVIDENCE[$index]}")"
done
for index in "${!EMBEDDING_EVIDENCE[@]}"; do
  EMBEDDING_EVIDENCE[$index]="$(resolve_submission_path "${EMBEDDING_EVIDENCE[$index]}")"
done

for value in BENCHMARK_ID BENCHMARK_DIR OBO_FILE RESULTS_ROOT EXECUTION_MODE; do
  [[ -n "${!value}" ]] || die "$value is required"
done
if [[ -n "$CACHE_ROOT" && -n "$CACHE_ARCHIVE" ]]; then
  die "Use exactly one of --embedding-cache-root or --embedding-cache-archive"
fi
if [[ -z "$CACHE_ROOT" && -z "$CACHE_ARCHIVE" ]]; then
  die "One of --embedding-cache-root or --embedding-cache-archive is required"
fi
[[ "$EXECUTION_MODE" =~ ^(eval-only|train-eval)$ ]] || die "Invalid --execution-mode"
[[ "$MODALITY_MODE" =~ ^(full|sequence-only)$ ]] || die "Invalid --modality-mode"
[[ "$CACHE_STAGING" =~ ^(copy|direct)$ ]] || die "Invalid --cache-staging"
[[ "$NUM_WORKERS" =~ ^[0-9]+$ && "$NUM_WORKERS" -le 8 ]] || \
  die "--num-workers must be an integer from 0 to 8 for this one-GPU wrapper"
ALLOCATED_SLOTS="${NSLOTS:-1}"
[[ "$ALLOCATED_SLOTS" =~ ^[1-9][0-9]*$ ]] || die "Invalid NSLOTS: $ALLOCATED_SLOTS"
MAX_WORKERS=$((ALLOCATED_SLOTS - 1))
[[ "$NUM_WORKERS" -le "$MAX_WORKERS" ]] || \
  die "--num-workers=$NUM_WORKERS leaves no allocated CPU slot for the main process; max=$MAX_WORKERS"
[[ -d "$BENCHMARK_DIR" ]] || die "Benchmark directory is missing: $BENCHMARK_DIR"
if [[ -n "$CACHE_ROOT" ]]; then
  [[ -d "$CACHE_ROOT" ]] || die "Embedding cache is missing: $CACHE_ROOT"
else
  [[ -f "$CACHE_ARCHIVE" ]] || die "Embedding cache archive is missing: $CACHE_ARCHIVE"
  [[ "$CACHE_STAGING" == "copy" ]] || \
    die "--embedding-cache-archive requires --cache-staging copy"
fi
[[ -f "$OBO_FILE" ]] || die "GO OBO is missing: $OBO_FILE"
if [[ "$EXECUTION_MODE" == "eval-only" ]]; then
  [[ -d "$CHECKPOINT_ROOT" ]] || die "--checkpoint-root is required for eval-only"
  [[ -d "$REFERENCE_DATA_DIR" ]] || die "--reference-data-dir is required for eval-only"
fi
if [[ -n "$REFERENCE_SOURCE_ARCHIVE" ]]; then
  [[ -f "$REFERENCE_SOURCE_ARCHIVE" ]] || \
    die "--reference-source-archive is missing: $REFERENCE_SOURCE_ARCHIVE"
  [[ -n "$REFERENCE_DATA_DIR" ]] || \
    die "--reference-source-archive requires --reference-data-dir"
fi
if [[ "$REQUIRE_REFERENCE_MATCH" == "1" ]]; then
  [[ -f "$EXPECTED_METRICS" ]] || \
    die "--require-reference-match also requires --expected-metrics"
fi
if [[ "$BENCHMARK_EVIDENCE_COUNT" -gt 0 ]]; then
  for evidence in "${BENCHMARK_EVIDENCE[@]}"; do
    [[ -f "$evidence" ]] || die "Benchmark evidence is missing: $evidence"
  done
fi
if [[ "$EMBEDDING_EVIDENCE_COUNT" -gt 0 ]]; then
  for evidence in "${EMBEDDING_EVIDENCE[@]}"; do
    [[ -f "$evidence" ]] || die "Embedding evidence is missing: $evidence"
  done
fi

require_unique_basenames() {
  local label="$1"
  shift
  local seen=" " path name
  for path in "$@"; do
    name="${path##*/}"
    case "$seen" in
      *" $name "*) die "$label contains repeated filename: $name" ;;
    esac
    seen="${seen}${name} "
  done
}
if [[ "$BENCHMARK_EVIDENCE_COUNT" -gt 0 ]]; then
  require_unique_basenames "Benchmark evidence" "${BENCHMARK_EVIDENCE[@]}"
fi
if [[ "$EMBEDDING_EVIDENCE_COUNT" -gt 0 ]]; then
  require_unique_basenames "Embedding evidence" "${EMBEDDING_EVIDENCE[@]}"
fi

JOB_TOKEN="${JOB_ID:-manual}"
if [[ -n "${SGE_TASK_ID:-}" && "${SGE_TASK_ID}" != "undefined" ]]; then
  JOB_TOKEN="${JOB_TOKEN}_${SGE_TASK_ID}"
fi
RUN_TAG="${JOB_TOKEN}_$(date +%Y%m%d_%H%M%S)"
WORK="/scratch0/pfp_benchmark_${JOB_TOKEN}"
FRAMEWORK_DIR="$WORK/Protein-Benchmark-Framework-Dissertation"
PFP_DIR="$WORK/PFP"
SCRATCH_INPUTS="$WORK/inputs"
SCRATCH_OUTPUT="$WORK/output"
WORKFLOW_WORK="$WORK/workflow"
WORKFLOW_LOG="$WORK/workflow.log"
FINAL_OUTPUT="$RESULTS_ROOT/$RUN_TAG"
FAILED_OUTPUT="${FINAL_OUTPUT}.failed"
FRAMEWORK_REPO_URL="${FRAMEWORK_REPO_URL:-https://github.com/nadroj0-0/Protein-Benchmark-Framework-Dissertation.git}"
PFP_REPO_URL="${PFP_REPO_URL:-https://github.com/psipred/PFP.git}"
PFP_COMMIT="${PFP_COMMIT:-1e04fd6d6d3c40458fd41ec1a881ed6e24de768e}"
FRAMEWORK_COMMIT="${FRAMEWORK_COMMIT:-}"
WORK_OWNED=0

publish() {
  local status="$1"
  local destination="$FINAL_OUTPUT"
  local staging="${FINAL_OUTPUT}.staging-${JOB_TOKEN}"
  local copy_status=0
  if [[ "$status" != "0" ]]; then
    destination="$FAILED_OUTPUT"
    staging="${FAILED_OUTPUT}.staging-${JOB_TOKEN}"
  fi
  [[ ! -e "$destination" && ! -e "$staging" ]] || return 1
  if [[ "$status" == "0" && -d "$SCRATCH_OUTPUT" ]]; then
    local payload_kb destination_free_kb
    payload_kb="$(du -sk "$SCRATCH_OUTPUT" | awk '{print $1}')" || return 1
    destination_free_kb="$(df -Pk "$RESULTS_ROOT" | awk 'END {print $4}')" || return 1
    if [[ "$destination_free_kb" -lt $((payload_kb + 1024 * 1024)) ]]; then
      echo "Insufficient destination space for publication plus 1 GiB margin" >&2
      return 1
    fi
  fi
  mkdir -p "$staging/logs" || return 1
  if [[ -d "$SCRATCH_OUTPUT" ]]; then
    cp -a "$SCRATCH_OUTPUT/." "$staging/" || copy_status=$?
  elif [[ -d "$WORKFLOW_WORK/logs" ]]; then
    cp -a "$WORKFLOW_WORK/logs/." "$staging/logs/" || copy_status=$?
    if [[ -d "$WORKFLOW_WORK/reports" ]]; then
      mkdir -p "$staging/reports" || copy_status=$?
      cp -a "$WORKFLOW_WORK/reports/." "$staging/reports/" || copy_status=$?
    fi
  fi
  if [[ -f "$WORKFLOW_LOG" ]]; then
    cp -p "$WORKFLOW_LOG" "$staging/logs/workflow.log" || copy_status=$?
  fi
  if [[ "$status" == "0" ]]; then
    [[ -f "$staging/WORKFLOW_COMPLETE.json" ]] || copy_status=1
    if [[ "$copy_status" == "0" ]]; then
      rm -f "$staging/output_manifest.json" "$staging/WORKFLOW_COMPLETE.json"
      "$PYTHON_BIN" "$FRAMEWORK_DIR/scripts/model_execution/manage_output_manifest.py" \
        write --root "$staging" >/dev/null || copy_status=$?
      if [[ "$copy_status" == "0" ]]; then
        local manifest_sha
        manifest_sha="$("$PYTHON_BIN" -c 'import hashlib,pathlib,sys; print(hashlib.sha256(pathlib.Path(sys.argv[1]).read_bytes()).hexdigest())' "$staging/output_manifest.json")" || copy_status=$?
        if [[ "$copy_status" == "0" ]]; then
          "$PYTHON_BIN" -c 'import json,pathlib,sys; pathlib.Path(sys.argv[1]).write_text(json.dumps({"complete":True,"manifest":"output_manifest.json","manifest_sha256":sys.argv[2],"run_report":"reports/run_report.json"},indent=2)+"\n")' \
            "$staging/WORKFLOW_COMPLETE.json" "$manifest_sha" || copy_status=$?
        fi
      fi
      if [[ "$copy_status" == "0" ]]; then
        "$PYTHON_BIN" "$FRAMEWORK_DIR/scripts/model_execution/manage_output_manifest.py" \
          verify --root "$staging" >/dev/null || copy_status=$?
      fi
    fi
  else
    rm -f "$staging/output_manifest.json" "$staging/WORKFLOW_COMPLETE.json"
    printf '{"complete":false,"workflow_exit_status":%s}\n' "$status" \
      >"$staging/WORKFLOW_FAILED.json" || copy_status=$?
    if [[ "$copy_status" == "0" && -f "$FRAMEWORK_DIR/scripts/model_execution/manage_output_manifest.py" ]]; then
      "$PYTHON_BIN" "$FRAMEWORK_DIR/scripts/model_execution/manage_output_manifest.py" \
        write --root "$staging" >/dev/null || copy_status=$?
      if [[ "$copy_status" == "0" ]]; then
        "$PYTHON_BIN" "$FRAMEWORK_DIR/scripts/model_execution/manage_output_manifest.py" \
          verify --root "$staging" >/dev/null || copy_status=$?
      fi
    fi
  fi
  if [[ "$copy_status" == "0" ]]; then
    mv "$staging" "$destination" || copy_status=$?
  fi
  if [[ "$copy_status" == "0" ]]; then
    echo "Published: $destination"
  elif [[ -d "$staging" && ! -L "$staging" ]]; then
    rm -rf -- "$staging"
  fi
  return "$copy_status"
}

remove_owned_work() {
  if [[ "$WORK_OWNED" == "1" && "$WORK" == /scratch0/pfp_benchmark_* && ! -L "$WORK" ]]; then
    cd "$HOME"
    rm -rf -- "$WORK"
    [[ ! -e "$WORK" ]] || return 1
    return 0
  fi
  echo "Refusing unsafe scratch cleanup: $WORK" >&2
  return 1
}

emergency_cleanup() {
  local signal_status="$1"
  trap - EXIT
  trap '' INT TERM USR2
  set +e
  remove_owned_work || echo "Scratch cleanup failed after signal: $WORK" >&2
  exit "$signal_status"
}

cleanup() {
  local status=$?
  local publish_status=0
  trap - EXIT
  trap 'emergency_cleanup 130' INT TERM
  trap 'emergency_cleanup 152' USR2
  set +e
  publish "$status" || publish_status=$?
  if ! remove_owned_work; then
    echo "Scratch cleanup failed: $WORK" >&2
    [[ "$status" != "0" ]] || status=1
  fi
  trap '' INT TERM USR2
  if [[ "$status" == "0" && "$publish_status" != "0" ]]; then status="$publish_status"; fi
  exit "$status"
}
trap cleanup EXIT
trap 'emergency_cleanup 130' INT TERM
trap 'emergency_cleanup 152' USR2

[[ ! -e "$WORK" ]] || die "Scratch path already exists: $WORK"
mkdir -p "$WORK"
WORK_OWNED=1
mkdir -p "$WORK/tmp" "$SCRATCH_INPUTS" "$RESULTS_ROOT"
export TMPDIR="$WORK/tmp" TMP="$WORK/tmp" TEMP="$WORK/tmp"

if [[ -z "$FRAMEWORK_COMMIT" ]]; then
  [[ -d "$SUBMISSION_DIR/.git" ]] || die "Submit from the framework checkout or pass FRAMEWORK_COMMIT"
  [[ -z "$(git -C "$SUBMISSION_DIR" status --porcelain)" ]] || die "Submission checkout has uncommitted changes"
  FRAMEWORK_COMMIT="$(git -C "$SUBMISSION_DIR" rev-parse HEAD)"
fi
[[ "$FRAMEWORK_COMMIT" =~ ^[0-9a-fA-F]{40}$ ]] || die "FRAMEWORK_COMMIT must be a full commit"

echo "Host            : $(hostname)"
echo "Benchmark       : $BENCHMARK_ID"
echo "Execution       : $EXECUTION_MODE"
echo "Modalities      : $MODALITY_MODE"
echo "Cache staging   : $CACHE_STAGING"
if [[ -n "$CACHE_ARCHIVE" ]]; then echo "Cache archive   : $CACHE_ARCHIVE"; fi
echo "Scratch         : $WORK"
echo "Final output    : $FINAL_OUTPUT"

git clone --no-checkout "$FRAMEWORK_REPO_URL" "$FRAMEWORK_DIR"
git -C "$FRAMEWORK_DIR" checkout --detach "$FRAMEWORK_COMMIT"
git clone --no-checkout "$PFP_REPO_URL" "$PFP_DIR"
git -C "$PFP_DIR" checkout --detach "$PFP_COMMIT"

input_kb=0
add_input_kb() {
  local value
  value="$(du -sk "$1" | awk '{print $1}')"
  input_kb=$((input_kb + value))
}
add_input_kb "$FRAMEWORK_DIR"
add_input_kb "$PFP_DIR"
add_input_kb "$BENCHMARK_DIR"
add_input_kb "$OBO_FILE"
if [[ -n "$CACHE_ARCHIVE" ]]; then
  add_input_kb "$CACHE_ARCHIVE"
elif [[ "$CACHE_STAGING" == "copy" ]]; then
  add_input_kb "$CACHE_ROOT"
fi
if [[ "$EXECUTION_MODE" == "eval-only" ]]; then
  add_input_kb "$CHECKPOINT_ROOT"
  add_input_kb "$REFERENCE_DATA_DIR"
fi
if [[ -n "$REFERENCE_SOURCE_ARCHIVE" ]]; then
  add_input_kb "$REFERENCE_SOURCE_ARCHIVE"
fi
if [[ -n "$IA_FILE_DIR" ]]; then add_input_kb "$IA_FILE_DIR"; fi
if [[ -n "$CLI_CONFIG" ]]; then add_input_kb "$CLI_CONFIG"; fi
if [[ -n "$EXPECTED_METRICS" ]]; then add_input_kb "$EXPECTED_METRICS"; fi
for evidence in "${BENCHMARK_EVIDENCE[@]}"; do add_input_kb "$evidence"; done
for evidence in "${EMBEDDING_EVIDENCE[@]}"; do add_input_kb "$evidence"; done
scratch_free_kb="$(df -Pk /scratch0 | awk 'END {print $4}')"
required_kb=$((input_kb + 80 * 1024 * 1024))
echo "Staged input estimate: ${input_kb} KiB"
echo "Scratch free space  : ${scratch_free_kb} KiB"
echo "Required with margin: ${required_kb} KiB"
[[ "$scratch_free_kb" -ge "$required_kb" ]] || \
  die "Insufficient /scratch0 space for staged inputs plus 80 GiB output margin"

cp -a "$BENCHMARK_DIR" "$SCRATCH_INPUTS/benchmark"
cp -p "$OBO_FILE" "$SCRATCH_INPUTS/go.obo"
if [[ -n "$CACHE_ARCHIVE" ]]; then
  cp -p "$CACHE_ARCHIVE" "$SCRATCH_INPUTS/embedding_cache.tar.gz"
  "$PYTHON_BIN" "$FRAMEWORK_DIR/scripts/embeddings/manage_embedding_archive.py" \
    extract \
    --archive "$SCRATCH_INPUTS/embedding_cache.tar.gz" \
    --output-cache-root "$SCRATCH_INPUTS/embedding_cache" \
    --config "${CLI_CONFIG:-$FRAMEWORK_DIR/configs/pfp_benchmark_run.temporal.json}" \
    --report "$SCRATCH_INPUTS/embedding_archive_extraction.json"
  EFFECTIVE_CACHE="$SCRATCH_INPUTS/embedding_cache"
elif [[ "$CACHE_STAGING" == "copy" ]]; then
  mkdir -p "$SCRATCH_INPUTS/embedding_cache"
  cp -a "$CACHE_ROOT/." "$SCRATCH_INPUTS/embedding_cache/"
  EFFECTIVE_CACHE="$SCRATCH_INPUTS/embedding_cache"
else
  EFFECTIVE_CACHE="$CACHE_ROOT"
fi
if [[ "$EXECUTION_MODE" == "eval-only" ]]; then
  cp -a "$CHECKPOINT_ROOT" "$SCRATCH_INPUTS/checkpoints"
  EFFECTIVE_CHECKPOINTS="$SCRATCH_INPUTS/checkpoints"
  cp -a "$REFERENCE_DATA_DIR" "$SCRATCH_INPUTS/reference_data"
fi
if [[ -n "$REFERENCE_SOURCE_ARCHIVE" ]]; then
  cp -p "$REFERENCE_SOURCE_ARCHIVE" "$SCRATCH_INPUTS/mmfp_data_splits.tar.gz"
fi
if [[ "$BENCHMARK_EVIDENCE_COUNT" -gt 0 ]]; then
  mkdir -p "$SCRATCH_INPUTS/benchmark_evidence"
  for evidence in "${BENCHMARK_EVIDENCE[@]}"; do
    cp -p "$evidence" "$SCRATCH_INPUTS/benchmark_evidence/"
  done
fi
if [[ "$EMBEDDING_EVIDENCE_COUNT" -gt 0 ]]; then
  mkdir -p "$SCRATCH_INPUTS/embedding_evidence"
  for evidence in "${EMBEDDING_EVIDENCE[@]}"; do
    cp -p "$evidence" "$SCRATCH_INPUTS/embedding_evidence/"
  done
fi
if [[ -n "$CLI_CONFIG" ]]; then
  cp -p "$CLI_CONFIG" "$SCRATCH_INPUTS/run_config.json"
  EFFECTIVE_CONFIG="$SCRATCH_INPUTS/run_config.json"
else
  EFFECTIVE_CONFIG="$FRAMEWORK_DIR/configs/pfp_benchmark_run.temporal.json"
fi
if [[ -n "$EXPECTED_METRICS" ]]; then
  cp -p "$EXPECTED_METRICS" "$SCRATCH_INPUTS/expected_metrics.json"
  EFFECTIVE_EXPECTED="$SCRATCH_INPUTS/expected_metrics.json"
fi
if [[ -n "$IA_FILE_DIR" ]]; then
  cp -a "$IA_FILE_DIR" "$SCRATCH_INPUTS/ia"
  EFFECTIVE_IA="$SCRATCH_INPUTS/ia"
fi

cd "$FRAMEWORK_DIR"
source scripts/reproduction_common.sh
load_framework_paths "$FRAMEWORK_DIR"
add_mmfp_singularity_bind "$WORK"
if [[ "$CACHE_STAGING" == "direct" ]]; then add_mmfp_singularity_bind "$CACHE_ROOT"; fi
activate_or_create_mmfp_env
PYTHON_BIN="$(command -v python)"
"$PYTHON_BIN" -c 'import torch,sys; assert torch.cuda.is_available(), "CUDA is unavailable"; assert torch.cuda.device_count() == 1, f"expected one visible GPU, found {torch.cuda.device_count()}"; value=torch.ones(1, device="cuda") * 2; assert value.item() == 2; print(f"CUDA preflight: {torch.cuda.get_device_name(0)}")'

COMMAND=(
  bash "$FRAMEWORK_DIR/scripts/model_execution/run_pfp_benchmark.sh"
  --benchmark-id "$BENCHMARK_ID"
  --benchmark-dir "$SCRATCH_INPUTS/benchmark"
  --embedding-cache-root "$EFFECTIVE_CACHE"
  --obo-file "$SCRATCH_INPUTS/go.obo"
  --pfp-root "$PFP_DIR"
  --work-dir "$WORKFLOW_WORK"
  --output-dir "$SCRATCH_OUTPUT"
  --config "$EFFECTIVE_CONFIG"
  --execution-mode "$EXECUTION_MODE"
  --modality-mode "$MODALITY_MODE"
  --seed "$SEED"
  --num-workers "$NUM_WORKERS"
  --expected-pfp-commit "$PFP_COMMIT"
)
if [[ "$ASPECT_COUNT" -gt 0 ]]; then
  for aspect in "${ASPECTS[@]}"; do COMMAND+=(--aspect "$aspect"); done
fi
if [[ "$EXECUTION_MODE" == "eval-only" ]]; then COMMAND+=(--checkpoint-root "$EFFECTIVE_CHECKPOINTS"); fi
if [[ "$EXECUTION_MODE" == "eval-only" ]]; then COMMAND+=(--reference-data-dir "$SCRATCH_INPUTS/reference_data"); fi
if [[ -n "$REFERENCE_SOURCE_ARCHIVE" ]]; then
  COMMAND+=(--reference-source-archive "$SCRATCH_INPUTS/mmfp_data_splits.tar.gz")
fi
if [[ "$BENCHMARK_EVIDENCE_COUNT" -gt 0 ]]; then
  for evidence in "$SCRATCH_INPUTS/benchmark_evidence/"*; do COMMAND+=(--benchmark-evidence "$evidence"); done
fi
if [[ "$EMBEDDING_EVIDENCE_COUNT" -gt 0 ]]; then
  for evidence in "$SCRATCH_INPUTS/embedding_evidence/"*; do COMMAND+=(--embedding-evidence "$evidence"); done
fi
if [[ "$REQUIRE_EMBEDDING_EVIDENCE" == "1" ]]; then COMMAND+=(--require-embedding-evidence); fi
if [[ -n "${EFFECTIVE_IA:-}" ]]; then COMMAND+=(--ia-file-dir "$EFFECTIVE_IA"); fi
if [[ -n "${EFFECTIVE_EXPECTED:-}" ]]; then COMMAND+=(--expected-metrics "$EFFECTIVE_EXPECTED"); fi
if [[ -n "$REFERENCE_TOLERANCE" ]]; then COMMAND+=(--reference-tolerance "$REFERENCE_TOLERANCE"); fi
if [[ "$REQUIRE_REFERENCE_MATCH" == "1" ]]; then COMMAND+=(--require-reference-match); fi

set +e
"${COMMAND[@]}" 2>&1 | tee "$WORKFLOW_LOG"
PIPE_RESULTS=("${PIPESTATUS[@]}")
STATUS="${PIPE_RESULTS[0]}"
if [[ "$STATUS" == "0" && "${PIPE_RESULTS[1]}" != "0" ]]; then
  STATUS="${PIPE_RESULTS[1]}"
fi
set -e
exit "$STATUS"
