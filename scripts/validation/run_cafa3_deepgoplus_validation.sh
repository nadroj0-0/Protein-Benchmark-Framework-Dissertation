#!/bin/bash
# Validate the historical CAFA3 DeepGOPlus -> TEMPROT -> PFP CSV path.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${HERE}/../.." && pwd)"

TIMESTAMP="${TIMESTAMP:-$(date +%Y%m%d_%H%M%S)}"
SCRATCH_BASE="${SCRATCH_BASE:-${TMPDIR:-/tmp}}"
RUN_DIR="${RUN_DIR:-${SCRATCH_BASE}/cafa3_deepgoplus_validation_${TIMESTAMP}}"
GENERATED="${RUN_DIR}/generated"
REFERENCE="${RUN_DIR}/reference"
REPORTS="${RUN_DIR}/reports"
LOGS="${RUN_DIR}/logs"
REPORT_COPY_DIR="${REPORT_COPY_DIR:-${HOME}/cafa3_deepgoplus_validation_reports/${TIMESTAMP}}"
KEEP_SCRATCH="${KEEP_SCRATCH:-0}"
PYTHON_BIN="${PYTHON_BIN:-python3}"

CAFA3_REFERENCE_RECORD_URL="https://zenodo.org/records/7409660"
CAFA3_REFERENCE_CSV_BASE_URL="${CAFA3_REFERENCE_RECORD_URL}/files"
DEFAULT_DEEPGOPLUS_PICKLES_URLS=(
  "https://deepgo.cbrc.kaust.edu.sa/data/data-cafa.tar.gz"
  "https://deepgo.cbrc.kaust.edu.sa/data/deepgoplus-cafa.tar.gz"
)

CSV_FILES=(
  bp-training.csv bp-validation.csv bp-test.csv
  cc-training.csv cc-validation.csv cc-test.csv
  mf-training.csv mf-validation.csv mf-test.csv
)
PICKLE_FILES=(
  train_data.pkl test_data.pkl terms.pkl
  train_data_train.pkl train_data_valid.pkl
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

locate_complete_set() {
  local root="$1"
  shift
  local candidate ok file
  while IFS= read -r candidate; do
    ok=1
    for file in "$@"; do
      [ -f "${candidate}/${file}" ] || { ok=0; break; }
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
  local pickle_status="$2"
  local go_obo="$3"
  local manifest="${REPORTS}/run_manifest.md"
  {
    echo "# CAFA3 DeepGOPlus Historical Validation Run Manifest"
    echo
    echo "- Run timestamp: ${TIMESTAMP}"
    echo "- Hostname: $(hostname)"
    echo "- User: ${USER:-unknown}"
    echo "- Repository root: ${REPO_ROOT}"
    echo "- Repository commit: $(git -C "$REPO_ROOT" rev-parse HEAD 2>/dev/null || echo unknown)"
    echo "- Scratch run directory: ${RUN_DIR}"
    echo "- Report copy directory: ${REPORT_COPY_DIR}"
    echo "- Source mode: DeepGOPlus released pickles"
    echo "- DeepGOPlus pickle source: ${pickle_status}"
    echo "- GO OBO used for export: ${go_obo}"
    echo
    echo "## Reference URLs"
    echo
    echo "- Canonical CAFA3 CSV record: ${CAFA3_REFERENCE_RECORD_URL}"
    echo "- DeepGOPlus CAFA archive candidates:"
    for url in "${DEFAULT_DEEPGOPLUS_PICKLES_URLS[@]}"; do
      echo "  - ${url}"
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
echo "CAFA3 DeepGOPlus Historical Validation"
echo "=============================================================="
echo "Run dir         : $RUN_DIR"
echo "Report copy dir : $REPORT_COPY_DIR"
echo "Scratch cleanup : $([ "$KEEP_SCRATCH" = "1" ] && echo disabled || echo enabled)"
echo

echo "==> [1/5] Locate or download DeepGOPlus reference pickles"
PICKLE_STATUS="not found"
REFERENCE_PICKLE_DIR=""
if [ -n "${DEEPGOPLUS_PICKLES_DIR:-}" ] && [ -d "${DEEPGOPLUS_PICKLES_DIR}" ]; then
  mkdir -p "${REFERENCE}/deepgoplus_pickles"
  cp -R "${DEEPGOPLUS_PICKLES_DIR}/." "${REFERENCE}/deepgoplus_pickles"/
  REFERENCE_PICKLE_DIR="$(locate_complete_set "${REFERENCE}/deepgoplus_pickles" "${PICKLE_FILES[@]}")"
  PICKLE_STATUS="copied from DEEPGOPLUS_PICKLES_DIR=${DEEPGOPLUS_PICKLES_DIR}"
else
  if [ -n "${DEEPGOPLUS_PICKLES_URL:-}" ]; then
    CANDIDATE_URLS=("$DEEPGOPLUS_PICKLES_URL")
  else
    CANDIDATE_URLS=("${DEFAULT_DEEPGOPLUS_PICKLES_URLS[@]}")
  fi

  mkdir -p "${REFERENCE}/deepgoplus_pickles"
  for DEEPGOPLUS_PICKLES_URL in "${CANDIDATE_URLS[@]}"; do
    PICKLE_ARCHIVE="${REFERENCE}/$(basename "$DEEPGOPLUS_PICKLES_URL")"
    if download "$DEEPGOPLUS_PICKLES_URL" "$PICKLE_ARCHIVE"; then
      rm -rf "${REFERENCE}/deepgoplus_pickles/extracted"
      mkdir -p "${REFERENCE}/deepgoplus_pickles/extracted"
      extract_archive "$PICKLE_ARCHIVE" "${REFERENCE}/deepgoplus_pickles/extracted"
      REFERENCE_PICKLE_DIR="$(locate_complete_set "${REFERENCE}/deepgoplus_pickles/extracted" "${PICKLE_FILES[@]}" || true)"
      if [ -n "$REFERENCE_PICKLE_DIR" ]; then
        PICKLE_STATUS="downloaded from ${DEEPGOPLUS_PICKLES_URL}"
        break
      fi
      echo "  Archive did not contain the required DeepGOPlus pickle set: ${DEEPGOPLUS_PICKLES_URL}"
    fi
  done
fi
if [ -z "$REFERENCE_PICKLE_DIR" ]; then
  echo "Could not locate required DeepGOPlus pickle files from candidate archive(s)." >&2
  exit 1
fi
echo "  Reference pickle directory: ${REFERENCE_PICKLE_DIR}"

GO_OBO_INPUT="${GO_OBO:-${REFERENCE_PICKLE_DIR}/go.obo}"
if [ ! -f "$GO_OBO_INPUT" ]; then
  echo "Could not locate GO OBO file. Set GO_OBO or use a DeepGOPlus directory containing go.obo." >&2
  exit 1
fi
echo "  GO OBO input: ${GO_OBO_INPUT}"

echo "==> [2/5] Export PFP CSVs from DeepGOPlus pickles"
BUILDER_PYTHONPATH="${REPO_ROOT}/benchmark_builders/contemporary_cafa/src${PYTHONPATH:+:${PYTHONPATH}}"
BUILDER_CMD=(
  "$PYTHON_BIN" -m cafa_benchmark_builder
  --source-mode deepgoplus
  --deepgoplus-dir "$REFERENCE_PICKLE_DIR"
  --go-obo "$GO_OBO_INPUT"
  --output-dir "$GENERATED"
)
printf '%q ' "${BUILDER_CMD[@]}" > "${LOGS}/builder_command.txt"
echo >> "${LOGS}/builder_command.txt"
PYTHONPATH="$BUILDER_PYTHONPATH" "${BUILDER_CMD[@]}" 2>&1 | tee "${LOGS}/builder.log"

echo "==> [3/5] Download canonical CAFA3 reference CSV artefacts"
REFERENCE_CSV_DIR="${REFERENCE}/cafa3_zenodo_7409660"
mkdir -p "$REFERENCE_CSV_DIR"
for csv_file in "${CSV_FILES[@]}"; do
  download "${CAFA3_REFERENCE_CSV_BASE_URL}/${csv_file}?download=1" "${REFERENCE_CSV_DIR}/${csv_file}"
done
REFERENCE_CSV_DIR="$(locate_complete_set "$REFERENCE_CSV_DIR" "${CSV_FILES[@]}")" || {
  echo "Could not locate all 9 reference CSVs under ${REFERENCE}/cafa3_zenodo_7409660" >&2
  exit 1
}
echo "  Reference CSV directory: ${REFERENCE_CSV_DIR}"

echo "==> [4/5] Write run manifest"
BUILDER_COMMAND_TEXT="$(cat "${LOGS}/builder_command.txt")"
write_manifest "$BUILDER_COMMAND_TEXT" "$PICKLE_STATUS" "$GO_OBO_INPUT"

echo "==> [5/5] Compare generated and reference outputs"
COMPARE_CMD=(
  "$PYTHON_BIN" "${REPO_ROOT}/scripts/validation/compare_cafa3_outputs.py"
  --generated-dir "$GENERATED"
  --reference-csv-dir "$REFERENCE_CSV_DIR"
  --reference-pickle-dir "$REFERENCE_PICKLE_DIR"
  --reports-dir "$REPORTS"
  --manifest-md "${REPORTS}/run_manifest.md"
)
"${COMPARE_CMD[@]}" 2>&1 | tee "${LOGS}/comparison.log"

echo
echo "==> DeepGOPlus validation complete."
echo "==> Reports will be copied to: ${REPORT_COPY_DIR}"
