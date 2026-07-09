#!/bin/bash
# Validate DeepGOPlus cafa3_data.py-style CAFA files -> pickle generation.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${HERE}/../.." && pwd)"

TIMESTAMP="${TIMESTAMP:-$(date +%Y%m%d_%H%M%S)}"
SCRATCH_BASE="${SCRATCH_BASE:-${TMPDIR:-/tmp}}"
RUN_DIR="${RUN_DIR:-${SCRATCH_BASE}/cafa3_deepgoplus_pickle_generation_${TIMESTAMP}}"
GENERATED="${RUN_DIR}/generated"
REFERENCE="${RUN_DIR}/reference"
REPORTS="${RUN_DIR}/reports"
LOGS="${RUN_DIR}/logs"
REPORT_COPY_DIR="${REPORT_COPY_DIR:-${HOME}/cafa3_deepgoplus_pickle_generation_reports/${TIMESTAMP}}"
KEEP_SCRATCH="${KEEP_SCRATCH:-0}"
PYTHON_BIN="${PYTHON_BIN:-python3}"

DEFAULT_DEEPGOPLUS_CAFA_URLS=(
  "https://deepgo.cbrc.kaust.edu.sa/data/data-cafa.tar.gz"
  "https://deepgo.cbrc.kaust.edu.sa/data/deepgoplus-cafa.tar.gz"
)

REQUIRED_RELATIVE_FILES=(
  "go.obo"
  "CAFA3_training_data/uniprot_sprot_exp.fasta"
  "CAFA3_training_data/uniprot_sprot_exp.txt"
  "CAFA3_targets/targets_all.fasta"
  "benchmark20171115/groundtruth/leafonly_all.txt"
  "train_data.pkl"
  "test_data.pkl"
  "terms.pkl"
)

mkdir -p "$GENERATED" "$REFERENCE" "$REPORTS" "$LOGS" "$REPORT_COPY_DIR"
LOG_FILE="${LOGS}/run.log"

copy_back() {
  mkdir -p "$REPORT_COPY_DIR"
  if [ -d "$REPORTS" ]; then
    cp -R "$REPORTS/." "$REPORT_COPY_DIR"/ 2>/dev/null || true
  fi
  if [ -d "$LOGS" ]; then
    mkdir -p "$REPORT_COPY_DIR/logs"
    cp -R "$LOGS/." "$REPORT_COPY_DIR/logs"/ 2>/dev/null || true
  fi
}

cleanup() {
  local status=$?
  echo
  echo "==> Final status: ${status}"
  echo "==> Copying reports/logs to: ${REPORT_COPY_DIR}"
  copy_back
  if [ "${KEEP_SCRATCH}" = "1" ]; then
    echo "==> KEEP_SCRATCH=1; preserving scratch run directory: ${RUN_DIR}"
  else
    echo "==> Removing scratch run directory: ${RUN_DIR}"
    cd "${HOME}"
    rm -rf "$RUN_DIR"
  fi
  exit "$status"
}
trap cleanup EXIT
trap 'echo "Interrupted"; exit 130' INT TERM

exec > >(tee -a "$LOG_FILE")
exec 2>&1

require_cmd() {
  command -v "$1" >/dev/null 2>&1 || {
    echo "Missing required command: $1" >&2
    exit 1
  }
}

download() {
  local url="$1"
  local out="$2"
  mkdir -p "$(dirname "$out")"
  if [ -s "$out" ]; then
    echo "  exists: $out"
  else
    echo "  downloading: $url"
  fi
  wget -c "$url" -O "$out"
}

extract_archive() {
  local archive="$1"
  local dest="$2"
  mkdir -p "$dest"
  case "$archive" in
    *.tar.gz|*.tgz)
      tar -xzf "$archive" -C "$dest"
      ;;
    *.zip)
      require_cmd unzip
      unzip -q -o "$archive" -d "$dest"
      ;;
    *)
      echo "Unsupported archive type: $archive" >&2
      return 1
      ;;
  esac
}

locate_cafa_root() {
  local root="$1"
  local candidate ok rel
  while IFS= read -r candidate; do
    ok=1
    for rel in "${REQUIRED_RELATIVE_FILES[@]}"; do
      [ -f "${candidate}/${rel}" ] || { ok=0; break; }
    done
    if [ "$ok" -eq 1 ]; then
      printf '%s\n' "$candidate"
      return 0
    fi
  done < <(find "$root" -type d | sort)
  return 1
}

sha256_file() {
  local file="$1"
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$file" | awk '{print $1}'
  else
    shasum -a 256 "$file" | awk '{print $1}'
  fi
}

write_manifest() {
  local builder_command="$1"
  local cafa_status="$2"
  local cafa_root="$3"
  local manifest="${REPORTS}/run_manifest.md"
  {
    echo "# CAFA3 DeepGOPlus Pickle Generation Run Manifest"
    echo
    echo "- Run timestamp: ${TIMESTAMP}"
    echo "- Hostname: $(hostname)"
    echo "- User: ${USER:-unknown}"
    echo "- Repository root: ${REPO_ROOT}"
    echo "- Repository commit: $(git -C "$REPO_ROOT" rev-parse HEAD 2>/dev/null || echo unknown)"
    echo "- Scratch run directory: ${RUN_DIR}"
    echo "- Report copy directory: ${REPORT_COPY_DIR}"
    echo "- Source mode: cafa3-files"
    echo "- DeepGOPlus CAFA source: ${cafa_status}"
    echo "- DeepGOPlus CAFA root: ${cafa_root}"
    echo
    echo "## DeepGOPlus CAFA Archive Candidates"
    echo
    for url in "${DEFAULT_DEEPGOPLUS_CAFA_URLS[@]}"; do
      echo "- ${url}"
    done
    echo
    echo "## Required Inputs"
    echo
    for rel in "${REQUIRED_RELATIVE_FILES[@]}"; do
      echo "- ${rel}"
    done
    echo
    echo "## Builder Command"
    echo
    echo '```bash'
    echo "$builder_command"
    echo '```'
    echo
    echo "## Downloaded/Generated File Inventory"
    echo
    echo "|path|bytes|sha256|"
    echo "|---|---:|---|"
    find "$REFERENCE" "$GENERATED" -type f | sort | while IFS= read -r file; do
      bytes="$(wc -c < "$file" | tr -d ' ')"
      sha="$(sha256_file "$file")"
      echo "|${file#${RUN_DIR}/}|${bytes}|${sha}|"
    done
  } > "$manifest"
}

for cmd in wget tar find "$PYTHON_BIN" awk tee wc; do
  require_cmd "$cmd"
done

echo "=============================================================="
echo "CAFA3 DeepGOPlus Pickle Generation Validation"
echo "=============================================================="
echo "Run dir         : $RUN_DIR"
echo "Report copy dir : $REPORT_COPY_DIR"
echo "Scratch cleanup : $([ "$KEEP_SCRATCH" = "1" ] && echo disabled || echo enabled)"
echo

echo "==> [1/4] Locate or download DeepGOPlus CAFA source archive"
CAFA_STATUS="not found"
CAFA_ROOT=""
if [ -n "${DEEPGOPLUS_CAFA_DIR:-}" ] && [ -d "${DEEPGOPLUS_CAFA_DIR}" ]; then
  mkdir -p "${REFERENCE}/deepgoplus_cafa"
  cp -R "${DEEPGOPLUS_CAFA_DIR}/." "${REFERENCE}/deepgoplus_cafa"/
  CAFA_ROOT="$(locate_cafa_root "${REFERENCE}/deepgoplus_cafa" || true)"
  CAFA_STATUS="copied from DEEPGOPLUS_CAFA_DIR=${DEEPGOPLUS_CAFA_DIR}"
else
  if [ -n "${DEEPGOPLUS_CAFA_URL:-}" ]; then
    CANDIDATE_URLS=("$DEEPGOPLUS_CAFA_URL")
  else
    CANDIDATE_URLS=("${DEFAULT_DEEPGOPLUS_CAFA_URLS[@]}")
  fi

  mkdir -p "${REFERENCE}/deepgoplus_cafa"
  for url in "${CANDIDATE_URLS[@]}"; do
    ARCHIVE="${REFERENCE}/$(basename "$url")"
    if download "$url" "$ARCHIVE"; then
      rm -rf "${REFERENCE}/deepgoplus_cafa/extracted"
      mkdir -p "${REFERENCE}/deepgoplus_cafa/extracted"
      extract_archive "$ARCHIVE" "${REFERENCE}/deepgoplus_cafa/extracted"
      CAFA_ROOT="$(locate_cafa_root "${REFERENCE}/deepgoplus_cafa/extracted" || true)"
      if [ -n "$CAFA_ROOT" ]; then
        CAFA_STATUS="downloaded from ${url}"
        break
      fi
      echo "  Archive did not contain the complete required CAFA file/pickle set: ${url}"
    fi
  done
fi
if [ -z "$CAFA_ROOT" ]; then
  echo "Could not locate the complete DeepGOPlus CAFA source/reference file set." >&2
  exit 1
fi
echo "  CAFA source root: ${CAFA_ROOT}"

echo "==> [2/4] Regenerate DeepGOPlus pickles from CAFA files"
BUILDER_PYTHONPATH="${REPO_ROOT}/benchmark_builders/contemporary_cafa/src${PYTHONPATH:+:${PYTHONPATH}}"
BUILDER_CMD=(
  "$PYTHON_BIN" -m cafa_benchmark_builder
  --source-mode cafa3-files
  --go-obo "${CAFA_ROOT}/go.obo"
  --train-sequences-file "${CAFA_ROOT}/CAFA3_training_data/uniprot_sprot_exp.fasta"
  --train-annotations-file "${CAFA_ROOT}/CAFA3_training_data/uniprot_sprot_exp.txt"
  --test-sequences-file "${CAFA_ROOT}/CAFA3_targets/targets_all.fasta"
  --test-annotations-file "${CAFA_ROOT}/benchmark20171115/groundtruth/leafonly_all.txt"
  --output-dir "$GENERATED"
  --min-count 50
)
printf '%q ' "${BUILDER_CMD[@]}" > "${LOGS}/builder_command.txt"
echo >> "${LOGS}/builder_command.txt"
PYTHONPATH="$BUILDER_PYTHONPATH" "${BUILDER_CMD[@]}" 2>&1 | tee "${LOGS}/builder.log"

echo "==> [3/4] Write run manifest"
BUILDER_COMMAND_TEXT="$(cat "${LOGS}/builder_command.txt")"
write_manifest "$BUILDER_COMMAND_TEXT" "$CAFA_STATUS" "$CAFA_ROOT"

echo "==> [4/4] Compare regenerated pickles with DeepGOPlus references"
COMPARE_CMD=(
  "$PYTHON_BIN" "${REPO_ROOT}/scripts/validation/compare_deepgoplus_pickles.py"
  --generated-dir "$GENERATED"
  --reference-pickle-dir "$CAFA_ROOT"
  --reports-dir "$REPORTS"
  --manifest-md "${REPORTS}/run_manifest.md"
)
"${COMPARE_CMD[@]}" 2>&1 | tee "${LOGS}/comparison.log"

echo
echo "==> DeepGOPlus pickle generation validation complete."
echo "==> Reports will be copied to: ${REPORT_COPY_DIR}"
