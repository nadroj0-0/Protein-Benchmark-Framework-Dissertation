#!/bin/bash
# Build and validate a CAFA3 historical benchmark recreation in scratch.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${HERE}/../.." && pwd)"

TIMESTAMP="${TIMESTAMP:-$(date +%Y%m%d_%H%M%S)}"
SCRATCH_BASE="${SCRATCH_BASE:-${TMPDIR:-/scratch0/${USER:-$(whoami)}}}"
RUN_DIR="${RUN_DIR:-${SCRATCH_BASE}/cafa3_historical_validation_${TIMESTAMP}}"
RAW="${RUN_DIR}/raw"
GENERATED="${RUN_DIR}/generated"
REFERENCE="${RUN_DIR}/reference"
REPORTS="${RUN_DIR}/reports"
LOGS="${RUN_DIR}/logs"
REPORT_COPY_DIR="${REPORT_COPY_DIR:-${HOME}/cafa3_historical_validation_reports/${TIMESTAMP}}"
KEEP_SCRATCH="${KEEP_SCRATCH:-0}"
PYTHON_BIN="${PYTHON_BIN:-python}"
USE_PIGZ="${USE_PIGZ:-1}"
DECOMPRESS_GOA="${DECOMPRESS_GOA:-0}"
GOA_PROGRESS_INTERVAL="${GOA_PROGRESS_INTERVAL:-1000000}"

GOA_T0_URL="https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/UNIPROT/goa_uniprot_all.gaf.163.gz"
GOA_T1_URL="https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/UNIPROT/goa_uniprot_all.gaf.172.gz"
UNIPROT_T0_URL="https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2017_02/knowledgebase/uniprot_sprot-only2017_02.tar.gz"
UNIPROT_T1_URL="https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2017_11/knowledgebase/uniprot_sprot-only2017_11.tar.gz"
GO_T0_OBO_URL="https://release.geneontology.org/2017-02-01/ontology/go.obo"
GO_T0_BASIC_URL="https://release.geneontology.org/2017-02-01/ontology/go-basic.obo"
GO_T1_OBO_URL="https://release.geneontology.org/2017-11-01/ontology/go.obo"
GO_T1_BASIC_URL="https://release.geneontology.org/2017-11-01/ontology/go-basic.obo"
CAFA3_REFERENCE_RECORD_URL="https://zenodo.org/records/7409660"
CAFA3_REFERENCE_CSV_BASE_URL="${CAFA3_REFERENCE_RECORD_URL}/files"

UNIPROT_T0_RELEASE_DATE="15-Feb-2017"
UNIPROT_T1_RELEASE_DATE="22-Nov-2017"
CAFA3_T0_DATE="2017-02-13"
CAFA3_T1_DATE="2017-11-15"

CSV_FILES=(
  bp-training.csv bp-validation.csv bp-test.csv
  cc-training.csv cc-validation.csv cc-test.csv
  mf-training.csv mf-validation.csv mf-test.csv
)
PICKLE_FILES=(
  train_data.pkl test_data.pkl terms.pkl
  train_data_train.pkl train_data_valid.pkl
)

mkdir -p "$RAW" "$GENERATED" "$REFERENCE" "$REPORTS" "$LOGS" "$REPORT_COPY_DIR"
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

download_optional() {
  local url="$1"
  local out="$2"
  mkdir -p "$(dirname "$out")"
  if [ -s "$out" ]; then
    echo "  optional exists: $out"
    return 0
  fi
  echo "  optional download: $url"
  wget -c "$url" -O "$out" || {
    echo "  optional download failed: $url"
    rm -f "$out"
  }
}

sha256_file() {
  local file="$1"
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$file" | awk '{print $1}'
  else
    shasum -a 256 "$file" | awk '{print $1}'
  fi
}

extract_tar_once() {
  local tarball="$1"
  local dest="$2"
  local marker="${dest}/.extracted"
  mkdir -p "$dest"
  if [ -f "$marker" ]; then
    echo "  extracted already: $tarball"
    return 0
  fi
  echo "  extracting: $tarball -> $dest"
  tar -xzf "$tarball" -C "$dest"
  date > "$marker"
}

find_uniprot_input() {
  local dir="$1"
  local match
  match="$(find "$dir" -type f -name 'uniprot_sprot.dat.gz' | sort | head -1 || true)"
  if [ -z "$match" ]; then
    match="$(find "$dir" -type f -name '*.dat.gz' | sort | head -1 || true)"
  fi
  if [ -z "$match" ]; then
    match="$(find "$dir" -type f \( -name 'uniprot_sprot.fasta.gz' -o -name '*.fasta.gz' -o -name '*.fa.gz' \) | sort | head -1 || true)"
  fi
  if [ -z "$match" ]; then
    echo "Could not locate a Swiss-Prot DAT/FASTA file under $dir" >&2
    exit 1
  fi
  printf '%s\n' "$match"
}

decompress_gaf_if_requested() {
  local gz="$1"
  local out="${gz%.gz}"
  if [ "$DECOMPRESS_GOA" != "1" ]; then
    printf '%s\n' "$gz"
    return 0
  fi
  if [ -s "$out" ]; then
    echo "  decompressed GOA exists: $out" >&2
    printf '%s\n' "$out"
    return 0
  fi
  echo "  decompressing GOA in scratch: $gz -> $out" >&2
  if [ "$USE_PIGZ" != "0" ] && command -v pigz >/dev/null 2>&1; then
    pigz -dc "$gz" > "$out"
  else
    gzip -dc "$gz" > "$out"
  fi
  printf '%s\n' "$out"
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

discover_deepgoplus_url() {
  local roots=("$REPO_ROOT")
  if [ -n "${TEMPROT_DIR:-}" ] && [ -d "${TEMPROT_DIR}" ]; then
    roots+=("$TEMPROT_DIR")
  fi
  grep -RhoE 'https?://[^[:space:]"'"'"'<>)]*(data-cafa|deepgoplus|deepgo)[^[:space:]"'"'"'<>)]*\.(tar\.gz|tgz|zip)' "${roots[@]}" 2>/dev/null | head -1 || true
}

write_manifest() {
  local builder_command="$1"
  local pickle_status="$2"
  local manifest="${REPORTS}/run_manifest.md"
  {
    echo "# CAFA3 Historical Validation Run Manifest"
    echo
    echo "- Run timestamp: ${TIMESTAMP}"
    echo "- Hostname: $(hostname)"
    echo "- User: ${USER:-unknown}"
    echo "- Repository root: ${REPO_ROOT}"
    echo "- Repository commit: $(git -C "$REPO_ROOT" rev-parse HEAD 2>/dev/null || echo unknown)"
    echo "- Scratch run directory: ${RUN_DIR}"
    echo "- Report copy directory: ${REPORT_COPY_DIR}"
    echo "- GOA pigz stream: ${CAFA_BUILDER_USE_PIGZ:-${USE_PIGZ}}"
    echo "- GOA pre-unzip: ${DECOMPRESS_GOA}"
    echo "- GOA progress interval: ${GOA_PROGRESS_INTERVAL}"
    echo "- CAFA3 t0 date: ${CAFA3_T0_DATE}"
    echo "- CAFA3 t1 date: ${CAFA3_T1_DATE}"
    echo "- UniProt t0 release date: ${UNIPROT_T0_RELEASE_DATE}"
    echo "- UniProt t1 release date: ${UNIPROT_T1_RELEASE_DATE}"
    echo "- DeepGOPlus pickle reference status: ${pickle_status}"
    echo
    echo "## Raw Snapshot URLs"
    echo
    echo "- GOA t0: ${GOA_T0_URL}"
    echo "- GOA t1: ${GOA_T1_URL}"
    echo "- UniProt t0: ${UNIPROT_T0_URL}"
    echo "- UniProt t1: ${UNIPROT_T1_URL}"
    echo "- GO t0 go.obo: ${GO_T0_OBO_URL}"
    echo "- GO t0 go-basic.obo: ${GO_T0_BASIC_URL}"
    echo "- GO t1 go.obo: ${GO_T1_OBO_URL}"
    echo "- GO t1 go-basic.obo: ${GO_T1_BASIC_URL}"
    echo "- Canonical CAFA3 reference CSV record: ${CAFA3_REFERENCE_RECORD_URL}"
    echo
    echo "## Builder Command"
    echo
    echo '```bash'
    echo "$builder_command"
    echo '```'
    echo
    echo "## Artefact Locations"
    echo
    echo "- Generated outputs: ${GENERATED}"
    echo "- Reference outputs: ${REFERENCE}"
    echo "- Reports: ${REPORTS}"
    echo
    echo "## Downloaded/Generated File Inventory"
    echo
    echo "|path|bytes|sha256|"
    echo "|---|---:|---|"
    find "$RAW" "$REFERENCE" "$GENERATED" -type f | sort | while IFS= read -r file; do
      bytes="$(wc -c < "$file" | tr -d ' ')"
      sha="$(sha256_file "$file")"
      echo "|${file#${RUN_DIR}/}|${bytes}|${sha}|"
    done
  } > "$manifest"
}

for cmd in wget tar find "$PYTHON_BIN" awk tee wc; do
  require_cmd "$cmd"
done
if [ "$DECOMPRESS_GOA" = "1" ]; then
  require_cmd gzip
fi
export CAFA_BUILDER_USE_PIGZ="${CAFA_BUILDER_USE_PIGZ:-${USE_PIGZ}}"
export CAFA_BUILDER_GOA_PROGRESS_INTERVAL="${GOA_PROGRESS_INTERVAL}"

echo "=============================================================="
echo "CAFA3 Historical Validation"
echo "=============================================================="
echo "Run dir         : $RUN_DIR"
echo "Report copy dir : $REPORT_COPY_DIR"
echo "Scratch cleanup : $([ "$KEEP_SCRATCH" = "1" ] && echo disabled || echo enabled)"
echo "GOA pigz stream : $([ "${CAFA_BUILDER_USE_PIGZ}" = "0" ] && echo disabled || echo enabled-if-available)"
echo "GOA pre-unzip   : $([ "$DECOMPRESS_GOA" = "1" ] && echo enabled || echo disabled)"
echo "GOA progress    : every ${GOA_PROGRESS_INTERVAL} parsed rows"
echo

echo "==> [1/8] Download historical raw snapshots into scratch"
download "$GOA_T0_URL" "${RAW}/goa/goa_uniprot_all.gaf.163.gz"
download "$GOA_T1_URL" "${RAW}/goa/goa_uniprot_all.gaf.172.gz"
download "$UNIPROT_T0_URL" "${RAW}/uniprot/uniprot_sprot-only2017_02.tar.gz"
download "$UNIPROT_T1_URL" "${RAW}/uniprot/uniprot_sprot-only2017_11.tar.gz"
download "$GO_T0_OBO_URL" "${RAW}/go/2017-02-01/go.obo"
download "$GO_T0_BASIC_URL" "${RAW}/go/2017-02-01/go-basic.obo"
download "$GO_T1_OBO_URL" "${RAW}/go/2017-11-01/go.obo"
download "$GO_T1_BASIC_URL" "${RAW}/go/2017-11-01/go-basic.obo"

echo "==> [2/8] Download best-effort metadata"
download_optional "https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2017_02/relnotes.txt" "${RAW}/metadata/uniprot_2017_02_relnotes.txt"
download_optional "https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2017_02/changes.html" "${RAW}/metadata/uniprot_2017_02_changes.html"
download_optional "https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2017_11/relnotes.txt" "${RAW}/metadata/uniprot_2017_11_relnotes.txt"
download_optional "https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2017_11/changes.html" "${RAW}/metadata/uniprot_2017_11_changes.html"
download_optional "https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/UNIPROT/README" "${RAW}/metadata/goa_old_uniprot_README"
download_optional "https://release.geneontology.org/2017-02-01/summary.txt" "${RAW}/metadata/go_2017-02-01_summary.txt"
download_optional "https://release.geneontology.org/2017-11-01/summary.txt" "${RAW}/metadata/go_2017-11-01_summary.txt"

echo "==> [3/8] Extract UniProt tarballs in scratch"
extract_tar_once "${RAW}/uniprot/uniprot_sprot-only2017_02.tar.gz" "${RAW}/uniprot/release_2017_02"
extract_tar_once "${RAW}/uniprot/uniprot_sprot-only2017_11.tar.gz" "${RAW}/uniprot/release_2017_11"
UNIPROT_T0_INPUT="$(find_uniprot_input "${RAW}/uniprot/release_2017_02")"
UNIPROT_T1_INPUT="$(find_uniprot_input "${RAW}/uniprot/release_2017_11")"
echo "  UniProt t0 input: ${UNIPROT_T0_INPUT}"
echo "  UniProt t1 input: ${UNIPROT_T1_INPUT}"

echo "==> [3b/8] Resolve GOA inputs"
GOA_T0_INPUT="$(decompress_gaf_if_requested "${RAW}/goa/goa_uniprot_all.gaf.163.gz")"
GOA_T1_INPUT="$(decompress_gaf_if_requested "${RAW}/goa/goa_uniprot_all.gaf.172.gz")"
echo "  GOA t0 input: ${GOA_T0_INPUT}"
echo "  GOA t1 input: ${GOA_T1_INPUT}"

echo "==> [4/8] Run benchmark builder"
BUILDER_PYTHONPATH="${REPO_ROOT}/benchmark_builders/contemporary_cafa/src${PYTHONPATH:+:${PYTHONPATH}}"
BUILDER_CMD=(
  "$PYTHON_BIN" -m cafa_benchmark_builder
  --uniprot-t0 "$UNIPROT_T0_INPUT"
  --uniprot-t1 "$UNIPROT_T1_INPUT"
  --goa-t0 "$GOA_T0_INPUT"
  --goa-t1 "$GOA_T1_INPUT"
  --go-obo "${RAW}/go/2017-11-01/go-basic.obo"
  --reviewed-only
  --output-dir "$GENERATED"
)
printf '%q ' "${BUILDER_CMD[@]}" > "${LOGS}/builder_command.txt"
echo >> "${LOGS}/builder_command.txt"
PYTHONPATH="$BUILDER_PYTHONPATH" "${BUILDER_CMD[@]}" 2>&1 | tee "${LOGS}/builder.log"

echo "==> [5/8] Download canonical CAFA3 reference CSV artefacts"
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

echo "==> [6/8] Locate or download DeepGOPlus reference pickles"
PICKLE_STATUS="not found; pickle comparison skipped"
REFERENCE_PICKLE_DIR=""
if [ -n "${DEEPGOPLUS_PICKLES_DIR:-}" ] && [ -d "${DEEPGOPLUS_PICKLES_DIR}" ]; then
  mkdir -p "${REFERENCE}/deepgoplus_pickles"
  cp -R "${DEEPGOPLUS_PICKLES_DIR}/." "${REFERENCE}/deepgoplus_pickles"/
  REFERENCE_PICKLE_DIR="$(locate_complete_set "${REFERENCE}/deepgoplus_pickles" "${PICKLE_FILES[@]}" || true)"
  PICKLE_STATUS="copied from DEEPGOPLUS_PICKLES_DIR=${DEEPGOPLUS_PICKLES_DIR}"
else
  DEEPGOPLUS_PICKLES_URL="${DEEPGOPLUS_PICKLES_URL:-$(discover_deepgoplus_url)}"
  if [ -n "$DEEPGOPLUS_PICKLES_URL" ]; then
    case "$DEEPGOPLUS_PICKLES_URL" in
      *.tar.gz) PICKLE_SUFFIX=".tar.gz" ;;
      *.tgz) PICKLE_SUFFIX=".tgz" ;;
      *.zip) PICKLE_SUFFIX=".zip" ;;
      *) PICKLE_SUFFIX=".archive" ;;
    esac
    PICKLE_ARCHIVE="${REFERENCE}/deepgoplus_pickles_reference${PICKLE_SUFFIX}"
    download "$DEEPGOPLUS_PICKLES_URL" "$PICKLE_ARCHIVE"
    mkdir -p "${REFERENCE}/deepgoplus_pickles"
    extract_archive "$PICKLE_ARCHIVE" "${REFERENCE}/deepgoplus_pickles"
    REFERENCE_PICKLE_DIR="$(locate_complete_set "${REFERENCE}/deepgoplus_pickles" "${PICKLE_FILES[@]}" || true)"
    PICKLE_STATUS="downloaded from ${DEEPGOPLUS_PICKLES_URL}"
  fi
fi
if [ -z "$REFERENCE_PICKLE_DIR" ]; then
  echo "  DeepGOPlus reference pickles not found; comparison will continue CSV-only."
else
  echo "  Reference pickle directory: ${REFERENCE_PICKLE_DIR}"
fi

echo "==> [7/8] Write run manifest"
BUILDER_COMMAND_TEXT="$(cat "${LOGS}/builder_command.txt")"
write_manifest "$BUILDER_COMMAND_TEXT" "$PICKLE_STATUS"

echo "==> [8/8] Compare generated and reference outputs"
COMPARE_CMD=(
  "$PYTHON_BIN" "${REPO_ROOT}/scripts/validation/compare_cafa3_outputs.py"
  --generated-dir "$GENERATED"
  --reference-csv-dir "$REFERENCE_CSV_DIR"
  --reports-dir "$REPORTS"
  --manifest-md "${REPORTS}/run_manifest.md"
)
if [ -n "$REFERENCE_PICKLE_DIR" ]; then
  COMPARE_CMD+=(--reference-pickle-dir "$REFERENCE_PICKLE_DIR")
fi
"${COMPARE_CMD[@]}" 2>&1 | tee "${LOGS}/comparison.log"

echo
echo "==> Validation complete."
echo "==> Reports will be copied to: ${REPORT_COPY_DIR}"
