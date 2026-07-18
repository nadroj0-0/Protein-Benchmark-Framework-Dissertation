#!/usr/bin/env bash
# Rebuild the canonical CAFA3 embeddings, train fresh PFP models, and evaluate.

set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
FRAMEWORK_ROOT="$(cd "${HERE}/../.." && pwd)"
# shellcheck source=../reproduction_common.sh
source "$FRAMEWORK_ROOT/scripts/reproduction_common.sh"

PYTHON_BIN="${PYTHON_BIN:-python}"
MMFP_BASE_URL="${MMFP_BASE_URL:-https://zenodo.org/records/19498341/files}"
PREFLIGHT_PER_SPLIT="${PREFLIGHT_PER_SPLIT:-2}"
DISK_POLL_SECONDS="${DISK_POLL_SECONDS:-120}"
CAFA_ASSESSMENT_COMMIT="${CAFA_ASSESSMENT_COMMIT:-d72f0a5abb66d3224bd808e2015b55f1c9d18340}"

PFP_ROOT=""
WORK_DIR=""
OUTPUT_DIR=""
TEXT_CUTOFF_DATE="2016-02-17"
EMBEDDING_STATE_ROOT=""
EMBEDDING_MODE="initial"
EMBEDDING_POLICY="$FRAMEWORK_ROOT/configs/cafa3_embedding_resume.json"
PUBLISHED_EMBEDDING_ARCHIVE_DIR=""

usage() {
  cat <<'EOF'
Usage: run_cafa3_full_from_scratch_reproduction.sh \
  --pfp-root PATH \
  --work-dir PATH \
  --output-dir PATH \
  --embedding-state-root PATH \
  [--embedding-mode initial|resume] \
  [--embedding-policy PATH] \
  [--text-cutoff-date YYYY-MM-DD] \
  [--published-embedding-archive-dir PATH] [--artifact-catalog PATH]

The PFP root must be a disposable pinned clone. Published embeddings are
staged or downloaded only after all four modalities have been regenerated, compared,
then deleted before training. They are never used as training inputs.

The persistent embedding state must be outside job-owned scratch. In initial
mode, validated arrays are consolidated into one archive-backed baseline. If the historical coverage
gate is not met, the workflow exits successfully with GENERATION_INCOMPLETE.json
and does not train. Resume mode hydrates a passing state and continues at the
published-cache comparison stage.
EOF
}

die() {
  echo "ERROR: $*" >&2
  exit 2
}

git_in_dir() {
  local directory="$1"
  shift
  (cd "$directory" && git "$@")
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --pfp-root)
      [[ $# -ge 2 ]] || die "--pfp-root requires a path"
      PFP_ROOT="$2"; shift 2 ;;
    --work-dir)
      [[ $# -ge 2 ]] || die "--work-dir requires a path"
      WORK_DIR="$2"; shift 2 ;;
    --output-dir)
      [[ $# -ge 2 ]] || die "--output-dir requires a path"
      OUTPUT_DIR="$2"; shift 2 ;;
    --embedding-state-root)
      [[ $# -ge 2 ]] || die "--embedding-state-root requires a path"
      EMBEDDING_STATE_ROOT="$2"; shift 2 ;;
    --embedding-mode)
      [[ $# -ge 2 ]] || die "--embedding-mode requires initial or resume"
      EMBEDDING_MODE="$2"; shift 2 ;;
    --embedding-policy)
      [[ $# -ge 2 ]] || die "--embedding-policy requires a path"
      EMBEDDING_POLICY="$2"; shift 2 ;;
    --text-cutoff-date)
      [[ $# -ge 2 ]] || die "--text-cutoff-date requires YYYY-MM-DD"
      TEXT_CUTOFF_DATE="$2"; shift 2 ;;
    --published-embedding-archive-dir)
      [[ $# -ge 2 ]] || die "--published-embedding-archive-dir requires a path"
      PUBLISHED_EMBEDDING_ARCHIVE_DIR="$2"; shift 2 ;;
    --artifact-catalog)
      [[ $# -ge 2 ]] || die "--artifact-catalog requires a path"
      ARTIFACT_CATALOG="$2"; export ARTIFACT_CATALOG; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

[[ -d "$PFP_ROOT/.git" ]] || die "PFP root is not a Git checkout: $PFP_ROOT"
[[ -n "$WORK_DIR" ]] || die "--work-dir is required"
[[ -n "$OUTPUT_DIR" ]] || die "--output-dir is required"
[[ -n "$EMBEDDING_STATE_ROOT" ]] || die "--embedding-state-root is required"
[[ "$EMBEDDING_MODE" == "initial" || "$EMBEDDING_MODE" == "resume" ]] || \
  die "--embedding-mode must be initial or resume"
[[ -f "$EMBEDDING_POLICY" ]] || die "Missing embedding policy: $EMBEDDING_POLICY"
EMBEDDING_POLICY="$(cd "$(dirname "$EMBEDDING_POLICY")" && pwd)/$(basename "$EMBEDDING_POLICY")"
[[ "$TEXT_CUTOFF_DATE" =~ ^[0-9]{4}-[0-9]{2}-[0-9]{2}$ ]] || \
  die "--text-cutoff-date must use YYYY-MM-DD"
[[ "$PREFLIGHT_PER_SPLIT" =~ ^[1-9][0-9]*$ ]] || \
  die "PREFLIGHT_PER_SPLIT must be a positive integer"
[[ ! -e "$WORK_DIR" ]] || die "Work directory already exists: $WORK_DIR"
[[ ! -e "$OUTPUT_DIR" ]] || die "Output directory already exists: $OUTPUT_DIR"
command -v "$PYTHON_BIN" >/dev/null 2>&1 || die "Python not found: $PYTHON_BIN"
artifact_catalog_configure "$FRAMEWORK_ROOT" "${ARTIFACT_CATALOG:-}"

PFP_ROOT="$(cd "$PFP_ROOT" && pwd)"
[[ -z "$(git_in_dir "$PFP_ROOT" status --porcelain --untracked-files=no)" ]] || \
  die "PFP checkout has tracked changes; use a clean disposable clone"
[[ ! -e "$PFP_ROOT/data/embedding_cache" ]] || \
  die "PFP clone already has an embedding cache"
[[ ! -e "$PFP_ROOT/results/full_model" ]] || \
  die "PFP clone already has full_model results"

mkdir -p "$WORK_DIR" "$OUTPUT_DIR/logs" "$OUTPUT_DIR/reports"
WORK_DIR="$(cd "$WORK_DIR" && pwd)"
OUTPUT_DIR="$(cd "$OUTPUT_DIR" && pwd)"
mkdir -p "$EMBEDDING_STATE_ROOT"
EMBEDDING_STATE_ROOT="$(cd "$EMBEDDING_STATE_ROOT" && pwd)"
case "$EMBEDDING_STATE_ROOT/" in
  "$WORK_DIR/"*|"$PFP_ROOT/"*)
    die "Embedding state must be persistent and outside job-owned scratch: $EMBEDDING_STATE_ROOT" ;;
esac

ARCHIVE_STAGE="$WORK_DIR/published_archives"
PUBLISHED_ROOT="$WORK_DIR/published_cache_root"
RUNTIME_COMPAT="$WORK_DIR/runtime_compat"
PREFLIGHT_BACKUP="$WORK_DIR/preflight_full_split_backup"
ACQUISITION_LOG="$OUTPUT_DIR/reports/input_acquisition.tsv"
MODALITY_STATUS="$OUTPUT_DIR/reports/modality_status.tsv"
STATE_REPORT_DIR="$OUTPUT_DIR/reports/embedding_state"
BASELINE_DIR="$EMBEDDING_STATE_ROOT/baseline"
BASELINE_ARCHIVE="$BASELINE_DIR/cafa3_embedding_cache.tar.gz"
BASELINE_ASSEMBLY_REPORT="$BASELINE_DIR/embedding_assembly.tsv.gz"
mkdir -p "$ARCHIVE_STAGE" "$PUBLISHED_ROOT" "$RUNTIME_COMPAT"
printf 'role\tname\tsource\tstaged_path\tsha256\n' > "$ACQUISITION_LOG"
printf 'phase\tmodality\texit_status\n' > "$MODALITY_STATUS"

sha256_file() {
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$1" | awk '{print $1}'
  else
    shasum -a 256 "$1" | awk '{print $1}'
  fi
}

publish_baseline_pair() {
  local source_archive="$1"
  local source_report="$2"
  local staging="$EMBEDDING_STATE_ROOT/.baseline.staging.${JOB_ID:-$$}"
  if [[ -d "$BASELINE_DIR" ]]; then
    [[ -f "$BASELINE_ARCHIVE" && -f "$BASELINE_ASSEMBLY_REPORT" ]] || \
      die "Persistent baseline directory is incomplete: $BASELINE_DIR"
    [[ "$(sha256_file "$source_archive")" == "$(sha256_file "$BASELINE_ARCHIVE")" ]] || \
      die "Existing persistent archive differs from generated source"
    [[ "$(sha256_file "$source_report")" == "$(sha256_file "$BASELINE_ASSEMBLY_REPORT")" ]] || \
      die "Existing persistent assembly report differs from generated source"
    return
  fi
  [[ ! -e "$BASELINE_DIR" ]] || die "Unsafe persistent baseline path: $BASELINE_DIR"
  rm -rf "$staging"
  mkdir -p "$staging"
  if ! cp -p "$source_archive" "$staging/cafa3_embedding_cache.tar.gz" || \
     ! cp -p "$source_report" "$staging/embedding_assembly.tsv.gz"; then
    rm -rf "$staging"
    die "Failed to stage archive-backed baseline"
  fi
  [[ "$(sha256_file "$source_archive")" == \
      "$(sha256_file "$staging/cafa3_embedding_cache.tar.gz")" ]] || {
    rm -rf "$staging"
    die "Persistent baseline archive failed authentication"
  }
  [[ "$(sha256_file "$source_report")" == \
      "$(sha256_file "$staging/embedding_assembly.tsv.gz")" ]] || {
    rm -rf "$staging"
    die "Persistent baseline assembly report failed authentication"
  }
  mv "$staging" "$BASELINE_DIR"
}

file_count() {
  local directory="$1"
  [[ -d "$directory" ]] || { printf '0\n'; return; }
  find "$directory" -maxdepth 1 -type f -name '*.npy' -print | wc -l | tr -d ' '
}

download_file() {
  local url="$1"
  local destination="$2"
  local partial="${destination}.part"
  echo "Downloading: $url"
  if command -v wget >/dev/null 2>&1; then
    wget --tries=5 --timeout=60 -c "$url" -O "$partial"
  elif command -v curl >/dev/null 2>&1; then
    curl --fail --location --retry 5 --continue-at - --output "$partial" "$url"
  else
    die "Neither wget nor curl is available"
  fi
  mv "$partial" "$destination"
}

archive_sha256() {
  case "$1" in
    mmfp_embeddings_prott5.tar.gz)
      printf '%s\n' '30dd88fc4acbe3bc267bd8d5ae05e4d967fa7c169a6f063f12d2395fb0ffb00f' ;;
    mmfp_embeddings_struct_ppi.tar.gz)
      printf '%s\n' '6d7a243f2c5e2149c162698b4e6a5e297731a4fea835d57b9049bb31f4af32de' ;;
    mmfp_embeddings_text_temporal.tar.gz)
      printf '%s\n' 'df1bf558fab1c018286a5b389665245917ba951f2c6cd8558546d0b1a3b47e36' ;;
    *) return 1 ;;
  esac
}

stop_disk_monitor() {
  if [[ -n "${DISK_MONITOR_PID:-}" ]]; then
    kill "$DISK_MONITOR_PID" 2>/dev/null || true
    wait "$DISK_MONITOR_PID" 2>/dev/null || true
    DISK_MONITOR_PID=""
  fi
}

disk_monitor() {
  while true; do
    {
      printf '%s\t' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
      du -sk "$WORK_DIR" "$PFP_ROOT" 2>/dev/null | \
        awk '{sum += $1} END {printf "run_kib=%s", sum}'
      printf '\t'
      df -Pk "$WORK_DIR" | awk 'NR==2 {printf "filesystem_available_kib=%s", $4}'
      printf '\n'
    } >> "$OUTPUT_DIR/reports/scratch_usage.tsv"
    sleep "$DISK_POLL_SECONDS"
  done
}

trap stop_disk_monitor EXIT
trap 'stop_disk_monitor; exit 130' INT TERM
printf 'timestamp\trun_kib\tfilesystem_available_kib\n' \
  > "$OUTPUT_DIR/reports/scratch_usage.tsv"
disk_monitor &
DISK_MONITOR_PID=$!

run_parallel_modalities() {
  local phase="$1"
  local log_dir="$OUTPUT_DIR/logs/$phase"
  local rc=0
  mkdir -p "$log_dir"

  IFS=',' read -ra GPUS <<< "${CUDA_VISIBLE_DEVICES:-0,1,2}"
  gpu() { echo "${GPUS[$1]:-${GPUS[0]}}"; }
  echo "==> [$phase] GPUs: ${GPUS[*]} (PPI runs on CPU)"

  CUDA_VISIBLE_DEVICES="$(gpu 0)" DEVICE=cuda \
    bash "$FRAMEWORK_ROOT/scripts/embeddings/generate_embeddings_sequence.sh" \
    > "$log_dir/sequence.log" 2>&1 &
  local pid_sequence=$!
  CUDA_VISIBLE_DEVICES="$(gpu 1)" \
    bash "$FRAMEWORK_ROOT/scripts/embeddings/generate_embeddings_text_temporal_cls.sh" \
    > "$log_dir/text.log" 2>&1 &
  local pid_text=$!
  CUDA_VISIBLE_DEVICES="$(gpu 2)" DEVICE=cuda \
    ALPHAFOLD_PREFETCH_REPORT="$OUTPUT_DIR/reports/alphafold_prefetch_${phase}.json" \
    bash "$FRAMEWORK_ROOT/scripts/embeddings/generate_embeddings_structure.sh" \
    > "$log_dir/structure.log" 2>&1 &
  local pid_structure=$!
  CUDA_VISIBLE_DEVICES="" \
    bash "$FRAMEWORK_ROOT/scripts/embeddings/generate_embeddings_ppi.sh" \
    > "$log_dir/ppi.log" 2>&1 &
  local pid_ppi=$!

  for specification in \
    "sequence:$pid_sequence" "text:$pid_text" \
    "structure:$pid_structure" "ppi:$pid_ppi"; do
    local name="${specification%%:*}"
    local pid="${specification##*:}"
    local status=0
    wait "$pid" || status=$?
    printf '%s\t%s\t%s\n' "$phase" "$name" "$status" >> "$MODALITY_STATUS"
    if [[ "$status" == "0" ]]; then
      echo "==> [$phase/$name] OK"
    else
      echo "==> [$phase/$name] FAILED with status $status"
      rc=1
    fi
  done
  return "$rc"
}

initialize_embedding_state() {
  local pfp_commit framework_commit
  pfp_commit="$(git_in_dir "$PFP_ROOT" rev-parse HEAD)"
  framework_commit="$(git_in_dir "$FRAMEWORK_ROOT" rev-parse HEAD)"
  local command=(
    "$PYTHON_BIN" "$FRAMEWORK_ROOT/scripts/embeddings/manage_resumable_embedding_state.py"
    initialize \
    --state-root "$EMBEDDING_STATE_ROOT" \
    --benchmark-id cafa3-zijian-canonical \
    --benchmark-dir "$CAFA3_RAW_DIR" \
    --data-dir "$PFP_ROOT/data" \
    --policy "$EMBEDDING_POLICY" \
    --pfp-commit "$pfp_commit" \
    --framework-commit "$framework_commit" \
    --environment-report "$OUTPUT_DIR/reports/environment_validation.txt" \
    --source-file "go-ontology=$PFP_ROOT/data/go.obo" \
    --source-file "string-alias=$STRING_ALIAS_FILE" \
    --source-file "string-embeddings=$STRING_H5_FILE" \
    --source-file "pfp-prott5=$PFP_ROOT/scripts/extract_prott5_embeddings.py" \
    --source-file "pfp-text-extract=$PFP_ROOT/scripts/extract_uniprot_text.py" \
    --source-file "pfp-text-embed=$PFP_ROOT/scripts/embed_uniprot_descriptions.py" \
    --source-file "pfp-if1=$PFP_ROOT/scripts/extract_esm_if1_embeddings.py" \
    --source-file "pfp-ppi=$PFP_ROOT/scripts/extract_ppi_embeddings.py" \
    --source-file "framework-if1-compat=$RUNTIME_COMPAT/extract_esm_if1_embeddings.py" \
    --source-file "framework-ppi-compat=$RUNTIME_COMPAT/extract_ppi_embeddings.py" \
    --runtime-value "text_cutoff_date=$TEXT_CUTOFF_DATE" \
    --runtime-value "alphafold_acquisition=framework-bounded" \
    --runtime-value "alphafold_api_workers=${ALPHAFOLD_API_WORKERS:-8}" \
    --runtime-value "alphafold_download_workers=${ALPHAFOLD_DOWNLOAD_WORKERS:-8}"
  )
  if [[ -f "$BASELINE_ARCHIVE" && -f "$BASELINE_ASSEMBLY_REPORT" ]]; then
    command+=(
      --baseline-archive "$BASELINE_ARCHIVE"
      --baseline-assembly-report "$BASELINE_ASSEMBLY_REPORT"
    )
  fi
  "${command[@]}" > "$OUTPUT_DIR/reports/embedding_state_initialization.json"
}

state_gate_passed() {
  "$PYTHON_BIN" - "$EMBEDDING_STATE_ROOT/coverage.json" <<'PY'
import json
import sys
raise SystemExit(0 if json.load(open(sys.argv[1]))["embedding_gate_passed"] else 1)
PY
}

publish_incomplete_generation() {
  mkdir -p "$STATE_REPORT_DIR"
  "$PYTHON_BIN" "$FRAMEWORK_ROOT/scripts/embeddings/manage_resumable_embedding_state.py" \
    summary --state-root "$EMBEDDING_STATE_ROOT" --report-dir "$STATE_REPORT_DIR" \
    > "$OUTPUT_DIR/reports/embedding_state_summary.json"
  cp -p "$EMBEDDING_STATE_ROOT/GENERATION_INCOMPLETE.json" \
    "$OUTPUT_DIR/GENERATION_INCOMPLETE.json"
  if [[ -f "$PFP_ROOT/data/alphafold_coverage_results.txt" ]]; then
    cp -p "$PFP_ROOT/data/alphafold_coverage_results.txt" "$OUTPUT_DIR/reports/"
  fi
  if [[ -d "$PFP_ROOT/results/embedding_reports" ]]; then
    cp -a "$PFP_ROOT/results/embedding_reports" "$OUTPUT_DIR/reports/"
  fi
  stop_disk_monitor
  echo "==> Embedding generation remains incomplete. Valid arrays are durable at:"
  echo "    $EMBEDDING_STATE_ROOT"
  echo "==> Retry manifest: $EMBEDDING_STATE_ROOT/needs_retry.tsv"
  echo "==> Training and evaluation were intentionally not started."
}

echo "==> [1/13] Validate and record the author-supplied environment"
validate_mmfp_env "$PYTHON_BIN" > "$OUTPUT_DIR/reports/environment_validation.txt"

echo "==> [2/13] Download canonical CAFA3 and embedding-generation dependencies"
cd "$PFP_ROOT"
PFP_ROOT="$PFP_ROOT" \
CAFA_ASSESSMENT_COMMIT="$CAFA_ASSESSMENT_COMMIT" \
  bash "$FRAMEWORK_ROOT/scripts/embeddings/generate_embeddings_dependencies.sh" \
  > "$OUTPUT_DIR/logs/dependencies.log" 2>&1
[[ -f external/dependency_env.sh ]] || die "Dependency environment was not written"
# shellcheck disable=SC1091
source external/dependency_env.sh
export CAFA_ASSESSMENT_DIR STRING_H5_FILE STRING_ALIAS_FILE CAFA3_RAW_DIR

cp -p "$CAFA3_RAW_DIR/.zenodo_md5.txt" "$OUTPUT_DIR/reports/cafa3_zenodo_md5.txt"
for name in \
  bp-training.csv bp-validation.csv bp-test.csv \
  cc-training.csv cc-validation.csv cc-test.csv \
  mf-training.csv mf-validation.csv mf-test.csv; do
  path="$CAFA3_RAW_DIR/$name"
  [[ -f "$path" ]] || die "Missing canonical CAFA3 CSV: $path"
  printf 'cafa3-csv\t%s\t%s/%s\t%s\t%s\n' \
    "$name" "https://zenodo.org/records/7409660/files" "$name" "$path" \
    "$(sha256_file "$path")" >> "$ACQUISITION_LOG"
done
printf 'go-ontology\tgo.obo\t%s\t%s\t%s\n' \
  "$CAFA_ASSESSMENT_DIR/precrec/go_cafa3.obo" "$PFP_ROOT/data/go.obo" \
  "$(sha256_file "$PFP_ROOT/data/go.obo")" >> "$ACQUISITION_LOG"
cafa_commit_full="$(git_in_dir "$CAFA_ASSESSMENT_DIR" rev-parse HEAD)"
printf '%s\n' "$cafa_commit_full" > "$OUTPUT_DIR/reports/cafa_assessment_commit.txt"
printf 'cafa-assessment-code\tbenchmark_folder.py\t%s@%s\t%s\t%s\n' \
  "https://github.com/ashleyzhou972/CAFA_assessment_tool" "$cafa_commit_full" \
  "$CAFA_ASSESSMENT_DIR/benchmark_folder.py" \
  "$(sha256_file "$CAFA_ASSESSMENT_DIR/benchmark_folder.py")" >> "$ACQUISITION_LOG"
for specification in "string-alias:$STRING_ALIAS_FILE" "string-h5:$STRING_H5_FILE"; do
  role="${specification%%:*}"
  path="${specification#*:}"
  [[ -f "$path" ]] || die "Missing STRING input: $path"
  printf '%s\t%s\t%s\t%s\t%s\n' \
    "$role" "$(basename "$path")" "STRING v12.0" "$path" \
    "$(sha256_file "$path")" >> "$ACQUISITION_LOG"
done

echo "==> [3/13] Prepare and validate the exact PFP split contract"
bash "$FRAMEWORK_ROOT/scripts/embeddings/generate_embeddings_prepare_data.sh" \
  > "$OUTPUT_DIR/logs/prepare_data.log" 2>&1
"$PYTHON_BIN" "$FRAMEWORK_ROOT/scripts/verification/verify_splits.py" \
  --data-dir data --strict > "$OUTPUT_DIR/reports/full_split_validation_before_preflight.txt"
"$PYTHON_BIN" "$FRAMEWORK_ROOT/scripts/embeddings/generate_embeddings_fasta.py" \
  --data-dir data > "$OUTPUT_DIR/logs/proteins_fasta.log" 2>&1
full_fasta_count="$(grep -c '^>' "$PFP_ROOT/data/proteins.fasta")"
[[ "$full_fasta_count" == "69811" ]] || \
  die "Canonical CAFA3 FASTA count mismatch: $full_fasta_count != 69811"
printf 'prepared-fasta\tproteins.fasta\tgenerated-from-nine-csvs\t%s\t%s\n' \
  "$PFP_ROOT/data/proteins.fasta" "$(sha256_file "$PFP_ROOT/data/proteins.fasta")" \
  >> "$ACQUISITION_LOG"

echo "==> [4/13] Prepare IF1/PPI runtime compatibility without editing PFP"
IF1_NUMPY_OVERLAY="$WORK_DIR/if1_numpy_1_26_4"
install_mmfp_if1_numpy_overlay "$PYTHON_BIN" "$IF1_NUMPY_OVERLAY"
validate_mmfp_if1_env "$PYTHON_BIN" "$IF1_NUMPY_OVERLAY" \
  > "$OUTPUT_DIR/reports/if1_environment.json"
"$PYTHON_BIN" "$FRAMEWORK_ROOT/scripts/embeddings/build_pfp_ppi_compat_copy.py" \
  --source "$PFP_ROOT/scripts/extract_ppi_embeddings.py" \
  --output "$RUNTIME_COMPAT/extract_ppi_embeddings.py" \
  --report "$OUTPUT_DIR/reports/pfp_ppi_compatibility.json"
"$PYTHON_BIN" "$FRAMEWORK_ROOT/scripts/embeddings/build_pfp_if1_compat_copy.py" \
  --source "$PFP_ROOT/scripts/extract_esm_if1_embeddings.py" \
  --output "$RUNTIME_COMPAT/extract_esm_if1_embeddings.py" \
  --report "$OUTPUT_DIR/reports/pfp_if1_compatibility.json"

export PPI_EXTRACT_SCRIPT="$RUNTIME_COMPAT/extract_ppi_embeddings.py"
export IF1_EXTRACT_SCRIPT="$RUNTIME_COMPAT/extract_esm_if1_embeddings.py"
export IF1_PYTHON_BIN="$PYTHON_BIN"
export IF1_PYTHONPATH="$IF1_NUMPY_OVERLAY"
export TEXT_CUTOFF_DATE
export TEXT_REPORT_DIR="$PFP_ROOT/results/embedding_reports/text"
export HF_HOME="$WORK_DIR/model_cache/huggingface"
export TORCH_HOME="$WORK_DIR/model_cache/torch"
export ALPHAFOLD_ACQUISITION_MODE=framework-bounded
# Structures are acquisition intermediates. Keeping tens of thousands of PDBs
# as loose SAN files would consume the same inode quota this workflow protects.
export ALPHAFOLD_PERSISTENT_CACHE_DIR="$WORK_DIR/source_cache/alphafold_structures"
export ALPHAFOLD_API_WORKERS="${ALPHAFOLD_API_WORKERS:-8}"
export ALPHAFOLD_DOWNLOAD_WORKERS="${ALPHAFOLD_DOWNLOAD_WORKERS:-8}"
mkdir -p "$HF_HOME" "$TORCH_HOME"

if [[ "$EMBEDDING_MODE" == "resume" ]]; then
  [[ -f "$BASELINE_ARCHIVE" && -f "$BASELINE_ASSEMBLY_REPORT" ]] || \
    die "Resume mode requires the archive-backed baseline under $BASELINE_DIR"
  echo "==> Initialize and authenticate the cumulative embedding state"
  initialize_embedding_state
  echo "==> [5-8/13] Resume from the authenticated persistent embedding state"
  "$PYTHON_BIN" "$FRAMEWORK_ROOT/scripts/embeddings/manage_resumable_embedding_state.py" \
    summary --state-root "$EMBEDDING_STATE_ROOT" --report-dir "$STATE_REPORT_DIR" \
    > "$OUTPUT_DIR/reports/embedding_state_summary.json"
  if ! state_gate_passed; then
    publish_incomplete_generation
    exit 0
  fi
  rm -rf \
    data/embedding_cache/prott5 \
    data/embedding_cache/exp_text_embeddings_temporal \
    data/embedding_cache/IF1 \
    data/embedding_cache/ppi
  "$PYTHON_BIN" "$FRAMEWORK_ROOT/scripts/embeddings/manage_resumable_embedding_state.py" \
    hydrate --state-root "$EMBEDDING_STATE_ROOT" \
    --output-cache-root "$PFP_ROOT/data/embedding_cache" \
    --report "$OUTPUT_DIR/reports/embedding_state_hydration.json"
else
  if find "$EMBEDDING_STATE_ROOT" -mindepth 1 -maxdepth 1 -print -quit | grep -q .; then
    die "Initial mode requires an empty embedding state root; use retry/resume for existing state"
  fi

  echo "==> [5/13] Create a reversible bounded preflight view"
  "$PYTHON_BIN" "$FRAMEWORK_ROOT/scripts/embeddings/prepare_cafa3_embedding_preflight.py" \
    create --data-dir "$PFP_ROOT/data" --backup-dir "$PREFLIGHT_BACKUP" \
    --limit-per-split "$PREFLIGHT_PER_SPLIT" \
    > "$OUTPUT_DIR/reports/preflight_workspace.json"
  "$PYTHON_BIN" "$FRAMEWORK_ROOT/scripts/verification/verify_splits.py" \
    --data-dir data --strict > "$OUTPUT_DIR/reports/preflight_split_validation.txt"

  echo "==> [6/13] Run all four modalities on the bounded preflight"
  run_parallel_modalities preflight || die "A preflight modality failed"
  fasta_count="$(grep -c '^>' data/proteins.fasta)"
  [[ "$(file_count data/embedding_cache/prott5)" == "$fasta_count" ]] || \
    die "Preflight did not create every ProtT5 embedding"
  [[ "$(file_count data/embedding_cache/exp_text_embeddings_temporal)" -gt 0 ]] || \
    die "Preflight produced no temporal text embeddings"
  preflight_pdb="$(find data/alphafold_structures -maxdepth 1 -type f -name '*.pdb' -print 2>/dev/null | wc -l | tr -d ' ')"
  if [[ "$preflight_pdb" -gt 0 && "$(file_count data/embedding_cache/IF1)" == "0" ]]; then
    die "Preflight downloaded structures but IF1 saved zero embeddings"
  fi

  # PFP checkpoints failed text API calls as processed. Restart text cleanly for
  # the full population while retaining valid preflight arrays in other modes.
  rm -rf \
    data/embedding_cache/exp_text_embeddings \
    data/embedding_cache/exp_text_embeddings_temporal \
    data/embedding_cache/uniprot_text

  echo "==> [7/13] Restore and authenticate the complete prepared dataset"
  "$PYTHON_BIN" "$FRAMEWORK_ROOT/scripts/embeddings/prepare_cafa3_embedding_preflight.py" \
    restore --data-dir "$PFP_ROOT/data" --backup-dir "$PREFLIGHT_BACKUP" \
    > "$OUTPUT_DIR/reports/full_workspace_restored.json"
  "$PYTHON_BIN" "$FRAMEWORK_ROOT/scripts/verification/verify_splits.py" \
    --data-dir data --strict > "$OUTPUT_DIR/reports/full_split_validation_after_preflight.txt"
  [[ "$(grep -c '^>' data/proteins.fasta)" == "69811" ]] || \
    die "Full FASTA was not restored to all 69,811 proteins"

  echo "==> [8/13] Regenerate all four complete embedding modalities in parallel"
  generation_status=0
  run_parallel_modalities full || generation_status=$?
  baseline_work="$WORK_DIR/initial_embedding_baseline"
  mkdir -p "$baseline_work"
  "$PYTHON_BIN" "$FRAMEWORK_ROOT/scripts/embeddings/build_embedding_baseline_archive.py" \
    --generated-cache-root "$PFP_ROOT/data/embedding_cache" \
    --data-dir "$PFP_ROOT/data" \
    --policy "$EMBEDDING_POLICY" \
    --archive "$baseline_work/cafa3_embedding_cache.tar.gz" \
    --assembly-report "$baseline_work/embedding_assembly.tsv.gz" \
    --report "$OUTPUT_DIR/reports/initial_embedding_baseline.json"
  publish_baseline_pair \
    "$baseline_work/cafa3_embedding_cache.tar.gz" \
    "$baseline_work/embedding_assembly.tsv.gz"

  echo "==> Authenticate the archive baseline and initialize retry state"
  initialize_embedding_state
  attempt_id="${JOB_ID:-local}_$(date -u +%Y%m%dT%H%M%SZ)_initial"
  # This merge records reasons for missing pairs. Arrays already authenticated
  # in the baseline are recognized as accepted and are not copied into SAN.
  merge_command=(
    "$PYTHON_BIN" "$FRAMEWORK_ROOT/scripts/embeddings/manage_resumable_embedding_state.py"
    merge
    --state-root "$EMBEDDING_STATE_ROOT"
    --generated-cache-root "$PFP_ROOT/data/embedding_cache"
    --attempt-id "$attempt_id"
    --modality-status "$MODALITY_STATUS"
    --report-dir "$STATE_REPORT_DIR"
  )
  if [[ -f "$PFP_ROOT/data/alphafold_coverage_results.txt" ]]; then
    merge_command+=(--alphafold-report "$PFP_ROOT/data/alphafold_coverage_results.txt")
  fi
  if [[ -f "$OUTPUT_DIR/reports/alphafold_prefetch_full.json" ]]; then
    merge_command+=(--alphafold-prefetch-report "$OUTPUT_DIR/reports/alphafold_prefetch_full.json")
  fi
  "${merge_command[@]}" > "$OUTPUT_DIR/reports/embedding_state_merge.json"
  printf '%s\n' "$generation_status" > "$OUTPUT_DIR/reports/full_generation_exit_status.txt"

  # Downstream work sees only arrays that passed the independent state validator.
  rm -rf \
    data/embedding_cache/prott5 \
    data/embedding_cache/exp_text_embeddings_temporal \
    data/embedding_cache/IF1 \
    data/embedding_cache/ppi
  "$PYTHON_BIN" "$FRAMEWORK_ROOT/scripts/embeddings/manage_resumable_embedding_state.py" \
    hydrate --state-root "$EMBEDDING_STATE_ROOT" \
    --output-cache-root "$PFP_ROOT/data/embedding_cache" \
    --report "$OUTPUT_DIR/reports/embedding_state_hydration.json"
  if ! state_gate_passed; then
    publish_incomplete_generation
    exit 0
  fi
fi

echo "==> Validate the hydrated cache before published comparison and training"
"$PYTHON_BIN" "$FRAMEWORK_ROOT/scripts/verification/verify_embeddings.py" \
  --data-dir data --config "$FRAMEWORK_ROOT/configs/cafa3.json" \
  --strict --all-arrays --require-min-coverage \
  > "$OUTPUT_DIR/reports/generated_embedding_validation.txt"

echo "==> [9/13] Download and authenticate Zijian's published embedding cache"
for name in \
  mmfp_embeddings_prott5.tar.gz \
  mmfp_embeddings_struct_ppi.tar.gz \
  mmfp_embeddings_text_temporal.tar.gz; do
  destination="$ARCHIVE_STAGE/$name"
  explicit_archive=""
  if [[ -n "$PUBLISHED_EMBEDDING_ARCHIVE_DIR" ]]; then
    explicit_archive="$PUBLISHED_EMBEDDING_ARCHIVE_DIR/$name"
  fi
  archive_source="$(resolve_artifact_path "$(zijian_embedding_artifact_id "$name")" "$explicit_archive" || true)"
  if [[ -n "$archive_source" ]]; then
    echo "Staging published embedding archive: $archive_source"
    cp -p "$archive_source" "$destination"
  else
    download_file "$MMFP_BASE_URL/$name?download=1" "$destination"
  fi
  observed="$(sha256_file "$destination")"
  wanted="$(archive_sha256 "$name")"
  [[ "$observed" == "$wanted" ]] || die "Published archive checksum mismatch: $name"
  tar -tzf "$destination" >/dev/null
  tar -xzf "$destination" -C "$PUBLISHED_ROOT"
  printf 'published-embedding-archive\t%s\t%s\t%s\t%s\n' \
    "$name" "${archive_source:-$MMFP_BASE_URL/$name}" "$destination" "$observed" >> "$ACQUISITION_LOG"
done
PUBLISHED_CACHE="$PUBLISHED_ROOT/data/embedding_cache"
[[ "$(file_count "$PUBLISHED_CACHE/prott5")" == "69811" ]] || die "Published ProtT5 count mismatch"
[[ "$(file_count "$PUBLISHED_CACHE/exp_text_embeddings_temporal")" == "69517" ]] || die "Published text count mismatch"
[[ "$(file_count "$PUBLISHED_CACHE/IF1")" == "67948" ]] || die "Published IF1 count mismatch"
[[ "$(file_count "$PUBLISHED_CACHE/ppi")" == "58294" ]] || die "Published PPI count mismatch"

echo "==> [10/13] Compare regenerated and published arrays without enforcing equality"
"$PYTHON_BIN" "$FRAMEWORK_ROOT/scripts/diagnostics/compare_embeddings.py" \
  --generated-cache-root "$PFP_ROOT/data/embedding_cache" \
  --published-cache-root "$PUBLISHED_CACHE" \
  --out-csv "$OUTPUT_DIR/reports/embedding_comparison.csv" \
  --out-json "$OUTPUT_DIR/reports/embedding_comparison_summary.json" \
  > "$OUTPUT_DIR/logs/embedding_comparison.log" 2>&1
gzip -f "$OUTPUT_DIR/reports/embedding_comparison.csv"

rm -rf "$PUBLISHED_ROOT" "$ARCHIVE_STAGE"
[[ ! -e "$PUBLISHED_ROOT" && ! -e "$ARCHIVE_STAGE" ]] || \
  die "Published embedding cache was not discarded"
printf '{"discarded":true,"reason":"comparison complete; never used for training"}\n' \
  > "$OUTPUT_DIR/reports/published_cache_discarded.json"

echo "==> [11/13] Train fresh gated-bilinear late-fusion models"
rm -rf results/full_model results/full_model_eval
"$PYTHON_BIN" train.py \
  --seq-model prott5 \
  --fusion-types gated_bilinear \
  --aspects BPO CCO MFO \
  --use-late-fusion \
  --text-embedding-dir data/embedding_cache/exp_text_embeddings_temporal \
  --output-base results/full_model \
  --num-workers 0 \
  --seed 42 > "$OUTPUT_DIR/logs/training.log" 2>&1

mkdir -p "$OUTPUT_DIR/reports/training"
for aspect in BPO CCO MFO; do
  result="results/full_model/fusion_comparison/prott5/$aspect/gated_bilinear/results.json"
  checkpoint="results/full_model/fusion_comparison/prott5/$aspect/gated_bilinear/best_model.pt"
  [[ -s "$result" ]] || die "Training did not produce $aspect results"
  [[ -s "$checkpoint" ]] || die "Training did not produce $aspect checkpoint"
  cp -p "$result" "$OUTPUT_DIR/reports/training/${aspect}_results.json"
done

echo "==> [12/13] Evaluate fresh checkpoints against published paper metrics"
evaluation_status=0
set +e
"$PYTHON_BIN" scripts/reproduce_full_model.py \
  > "$OUTPUT_DIR/logs/evaluation.log" 2>&1
evaluation_status=$?
set -e
[[ "$evaluation_status" == "0" || "$evaluation_status" == "1" ]] || \
  die "Evaluation failed as infrastructure (status $evaluation_status)"
EVAL_DIR="results/full_model_eval"
[[ -s "$EVAL_DIR/reproduction_summary.json" ]] || \
  die "Evaluation did not produce its JSON summary"
[[ -s "$EVAL_DIR/reproduction_summary.csv" ]] || \
  die "Evaluation did not produce its CSV summary"
mkdir -p "$OUTPUT_DIR/reports/evaluation"
cp -p "$EVAL_DIR/reproduction_summary.json" "$OUTPUT_DIR/reports/evaluation/"
cp -p "$EVAL_DIR/reproduction_summary.csv" "$OUTPUT_DIR/reports/evaluation/"
printf '%s\n' "$evaluation_status" > "$OUTPUT_DIR/reports/evaluation/exit_status.txt"
for aspect in BPO CCO MFO; do
  result="$EVAL_DIR/eval_only/fusion_comparison/prott5/$aspect/gated_bilinear/results.json"
  [[ -s "$result" ]] || die "Evaluation did not produce $aspect results"
  cp -p "$result" "$OUTPUT_DIR/reports/evaluation/${aspect}_results.json"
done

echo "==> [13/13] Build the complete compact reproduction report"
"$PYTHON_BIN" "$FRAMEWORK_ROOT/scripts/diagnostics/build_cafa3_full_reproduction_report.py" \
  --pfp-root "$PFP_ROOT" \
  --framework-root "$FRAMEWORK_ROOT" \
  --embedding-summary "$OUTPUT_DIR/reports/embedding_comparison_summary.json" \
  --evaluation-summary "$OUTPUT_DIR/reports/evaluation/reproduction_summary.json" \
  --modality-status "$MODALITY_STATUS" \
  --input-acquisition "$ACQUISITION_LOG" \
  --evaluation-exit-status "$evaluation_status" \
  --text-cutoff-date "$TEXT_CUTOFF_DATE" \
  --published-cache-discarded \
  --output-json "$OUTPUT_DIR/cafa3_full_reproduction_report.json" \
  --output-md "$OUTPUT_DIR/cafa3_full_reproduction_report.md" \
  > "$OUTPUT_DIR/logs/report.log" 2>&1

if [[ -d "$PFP_ROOT/results/embedding_reports" ]]; then
  cp -a "$PFP_ROOT/results/embedding_reports" "$OUTPUT_DIR/reports/"
fi
if [[ -f "$PFP_ROOT/data/alphafold_coverage_results.txt" ]]; then
  cp -p "$PFP_ROOT/data/alphafold_coverage_results.txt" "$OUTPUT_DIR/reports/"
fi
if [[ -f "$PREFLIGHT_BACKUP/preflight_backup_manifest.json" ]]; then
  cp -p "$PREFLIGHT_BACKUP/preflight_backup_manifest.json" \
    "$OUTPUT_DIR/reports/preflight_backup_manifest.json"
else
  printf '{"skipped":true,"reason":"embedding state resume mode"}\n' \
    > "$OUTPUT_DIR/reports/preflight_skipped.json"
fi

"$PYTHON_BIN" - "$OUTPUT_DIR" "$PFP_ROOT" "$FRAMEWORK_ROOT" "$TEXT_CUTOFF_DATE" \
  "$EMBEDDING_MODE" "$EMBEDDING_STATE_ROOT" <<'PY'
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

output, pfp, framework = map(Path, sys.argv[1:4])
payload = {
    "complete": True,
    "schema_version": 1,
    "completed_at": datetime.now(timezone.utc).isoformat(),
    "text_cutoff_date": sys.argv[4],
    "embedding_mode": sys.argv[5],
    "embedding_state_root": sys.argv[6],
    "pfp_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=pfp, text=True).strip(),
    "framework_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=framework, text=True).strip(),
    "report_markdown": "cafa3_full_reproduction_report.md",
    "report_json": "cafa3_full_reproduction_report.json",
    "generated_embeddings_persisted": True,
    "published_embeddings_persisted": False,
}
(output / "WORKFLOW_COMPLETE.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
PY

stop_disk_monitor
echo "==> Full CAFA3 from-scratch reproduction complete: $OUTPUT_DIR"
