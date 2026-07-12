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
INCLUDE_TREMBL_TARGETS="${INCLUDE_TREMBL_TARGETS:-1}"
PIGZ_THREADS="${PIGZ_THREADS:-1}"
HISTORICAL_TRAINING_SNAPSHOT="${HISTORICAL_TRAINING_SNAPSHOT:-september-2016}"
TARGET_UNIVERSE_POLICY="${TARGET_UNIVERSE_POLICY:-official-cafa3-targets}"
HISTORICAL_TEST_SOURCE="${HISTORICAL_TEST_SOURCE:-official-groundtruth}"
HISTORICAL_T1_ENDPOINT_POLICY="${HISTORICAL_T1_ENDPOINT_POLICY:-assigned-date-proxy}"
HISTORICAL_BACKFILL_POLICY="${HISTORICAL_BACKFILL_POLICY:-exclude-pre-t0}"
HISTORICAL_BENCHMARK_ONTOLOGY="${HISTORICAL_BENCHMARK_ONTOLOGY:-}"

FILTER_DAT="${REPO_ROOT}/scripts/benchmark_generation/filter_uniprot_dat.py"
EXTRACT_MEMBER="${REPO_ROOT}/scripts/benchmark_generation/extract_tar_member.py"
TARGET_TAXA="${REPO_ROOT}/benchmark_builders/contemporary_cafa/src/cafa_benchmark_builder/resources/cafa3_target_taxa.txt"

GOA_T0_URL="https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/UNIPROT/goa_uniprot_all.gaf.163.gz"
GOA_T1_URL="https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/UNIPROT/goa_uniprot_all.gaf.172.gz"
UNIPROT_T0_URL="https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2017_02/knowledgebase/uniprot_sprot-only2017_02.tar.gz"
UNIPROT_T1_URL="https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2017_11/knowledgebase/uniprot_sprot-only2017_11.tar.gz"
UNIPROT_T0_FULL_URL="https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2017_02/knowledgebase/knowledgebase2017_02.tar.gz"
UNIPROT_T1_FULL_URL="https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2017_11/knowledgebase/knowledgebase2017_11.tar.gz"
GO_T0_OBO_URL="https://release.geneontology.org/2017-02-01/ontology/go.obo"
GO_T0_BASIC_URL="https://release.geneontology.org/2017-02-01/ontology/go-basic.obo"
GO_T1_OBO_URL="https://release.geneontology.org/2017-11-01/ontology/go.obo"
GO_T1_BASIC_URL="https://release.geneontology.org/2017-11-01/ontology/go-basic.obo"
CAFA3_REFERENCE_RECORD_URL="https://zenodo.org/records/7409660"
CAFA3_REFERENCE_CSV_BASE_URL="${CAFA3_REFERENCE_RECORD_URL}/files"
DEFAULT_DEEPGOPLUS_PICKLES_URL="https://deepgo.cbrc.kaust.edu.sa/data/data-cafa.tar.gz"
UNIPROT_2016_08_URL="https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2016_08/knowledgebase/uniprot_sprot-only2016_08.tar.gz"
UNIPROT_2016_08_METALINK_URL="https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2016_08/knowledgebase/RELEASE.metalink"
UNIPROT_2016_08_SIZE="1515243865"
UNIPROT_2016_08_MD5="bad05b535790955ddd5e0d1833915f9f"

UNIPROT_T0_RELEASE_DATE="15-Feb-2017"
UNIPROT_T1_RELEASE_DATE="22-Nov-2017"
CAFA3_T0_DATE="2017-02-13"
CAFA3_T1_DATE="2017-11-15"
TRAINING_SNAPSHOT_ID=""
TRAINING_SNAPSHOT_DATE=""

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
  if [ -d "$GENERATED" ]; then
    mkdir -p "$REPORT_COPY_DIR/generated"
    cp -R "$GENERATED/." "$REPORT_COPY_DIR/generated"/ 2>/dev/null || true
  fi
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
  echo "==> Copying generated artefacts, reports and logs to: ${REPORT_COPY_DIR}"
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

md5_file() {
  local file="$1"
  if command -v md5sum >/dev/null 2>&1; then
    md5sum "$file" | awk '{print $1}'
  else
    md5 -q "$file"
  fi
}

stage_override_or_download() {
  local override="$1"
  local url="$2"
  local out="$3"
  mkdir -p "$(dirname "$out")"
  if [ -n "$override" ]; then
    [ -f "$override" ] || { echo "Explicit input does not exist: $override" >&2; exit 1; }
    echo "  staging explicit input: $override -> $out"
    cp -a "$override" "$out"
  else
    download "$url" "$out"
  fi
}

verify_size_and_md5() {
  local file="$1"
  local expected_size="$2"
  local expected_md5="$3"
  local actual_size actual_md5
  actual_size="$(wc -c < "$file" | tr -d ' ')"
  actual_md5="$(md5_file "$file")"
  if [ "$actual_size" != "$expected_size" ] || [ "$actual_md5" != "$expected_md5" ]; then
    echo "Integrity check failed for $file" >&2
    echo "  expected size/md5: $expected_size $expected_md5" >&2
    echo "  actual size/md5  : $actual_size $actual_md5" >&2
    exit 1
  fi
  echo "  verified official size and MD5: $actual_size $actual_md5"
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

compress_stream() {
  if [ "$USE_PIGZ" != "0" ] && command -v pigz >/dev/null 2>&1; then
    pigz -p "$PIGZ_THREADS" -c
  else
    gzip -c
  fi
}

decompress_stream() {
  if [ "$USE_PIGZ" != "0" ] && command -v pigz >/dev/null 2>&1; then
    pigz -p "$PIGZ_THREADS" -dc
  else
    gzip -dc
  fi
}

filter_trembl_archive_url() {
  local label="$1"
  local url="$2"
  local destination="$3"
  local temporary="${destination}.part"
  mkdir -p "$(dirname "$destination")"
  echo "  streaming ${label} and retaining only CAFA3 target taxa"
  rm -f "$temporary"
  if wget --progress=dot:giga -O - "$url" \
    | "$PYTHON_BIN" "$EXTRACT_MEMBER" --suffix uniprot_trembl.dat.gz \
    | decompress_stream \
    | "$PYTHON_BIN" "$FILTER_DAT" --taxa-file "$TARGET_TAXA" \
    | compress_stream > "$temporary"; then
    mv "$temporary" "$destination"
  else
    rm -f "$temporary"
    return 1
  fi
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
  grep -RhoE 'https?://[^[:space:]"'"'"'<>)]*(data-cafa|deepgoplus|deepgo)[^[:space:]"'"'"'<>)]*\.(tar\.gz|tgz|zip)' "${roots[@]}" 2>/dev/null \
    | grep -v '^https\?://example/' \
    | head -1 || true
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
    echo "- TrEMBL target records included: ${INCLUDE_TREMBL_TARGETS}"
    echo "- Historical training snapshot profile: ${HISTORICAL_TRAINING_SNAPSHOT}"
    echo "- Training snapshot identifier: ${TRAINING_SNAPSHOT_ID}"
    echo "- Training snapshot release date: ${TRAINING_SNAPSHOT_DATE}"
    echo "- Target universe policy: ${TARGET_UNIVERSE_POLICY}"
    echo "- Historical test source: ${HISTORICAL_TEST_SOURCE}"
    echo "- Historical t1 endpoint policy: ${HISTORICAL_T1_ENDPOINT_POLICY}"
    echo "- Historical backfill policy: ${HISTORICAL_BACKFILL_POLICY}"
    echo "- Historical benchmark ontology: ${HISTORICAL_BENCHMARK_ONTOLOGY}"
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
    echo "- Closest public pre-freeze Swiss-Prot snapshot: ${UNIPROT_2016_08_URL}"
    echo "- Published UniProt 2016_08 integrity metadata: ${UNIPROT_2016_08_METALINK_URL}"
    echo "- UniProt t0 complete archive used for target-taxon TrEMBL: ${UNIPROT_T0_FULL_URL}"
    echo "- UniProt t1 complete archive used for target-taxon TrEMBL: ${UNIPROT_T1_FULL_URL}"
    echo "- GO t0 go.obo: ${GO_T0_OBO_URL}"
    echo "- GO t0 go-basic.obo: ${GO_T0_BASIC_URL}"
    echo "- GO t1 go.obo: ${GO_T1_OBO_URL}"
    echo "- GO t1 go-basic.obo: ${GO_T1_BASIC_URL}"
    echo "- Canonical CAFA3 reference CSV record: ${CAFA3_REFERENCE_RECORD_URL}"
    echo "- Default DeepGOPlus CAFA pickle archive: ${DEFAULT_DEEPGOPLUS_PICKLES_URL}"
    echo
    echo "## Historical Interpretation Boundary"
    echo
    echo "- The official CAFA3 training package was dated 24-Sep-2016."
    echo "- UniProtKB 2016_08 (07-Sep-2016) is the last public monthly Swiss-Prot release before that package date; it is the closest defensible public sequence snapshot, not a claim of an exact organiser-internal freeze."
    if [ "$HISTORICAL_TEST_SOURCE" = "official-groundtruth" ]; then
      echo "- Released CAFA target IDs and exact FASTA sequences define the test universe directly; UniProt mapping is intentionally bypassed because it is not part of the artifact-reproduction claim."
      echo "- Released leafonly_all.txt labels define the test benchmark exactly as consumed by DeepGOPlus; this mode validates artifact reproduction, not raw t1 snapshot reconstruction."
    else
      echo "- In raw-GOA official-target mode, conservative UniProt mapping is reported and is not used to replace released target sequences."
      echo "- Raw GOA mode reconstructs the temporal test set from the closest public archives and retains the documented 15-Nov-2017 snapshot limitation."
    fi
    echo "- February 2017 remains the target t0 knowledge baseline and a named legacy training-snapshot option."
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
for helper in "$FILTER_DAT" "$EXTRACT_MEMBER" "$TARGET_TAXA"; do
  [ -f "$helper" ] || { echo "Missing required helper: $helper" >&2; exit 1; }
done
if [ "$DECOMPRESS_GOA" = "1" ] || [ "$INCLUDE_TREMBL_TARGETS" = "1" ]; then
  require_cmd gzip
fi
export CAFA_BUILDER_USE_PIGZ="${CAFA_BUILDER_USE_PIGZ:-${USE_PIGZ}}"
export CAFA_BUILDER_GOA_PROGRESS_INTERVAL="${GOA_PROGRESS_INTERVAL}"

case "$HISTORICAL_TRAINING_SNAPSHOT" in
  september-2016|february-2017-legacy) ;;
  *)
    echo "Unknown HISTORICAL_TRAINING_SNAPSHOT: $HISTORICAL_TRAINING_SNAPSHOT" >&2
    echo "Expected september-2016 or february-2017-legacy" >&2
    exit 1
    ;;
esac
case "$TARGET_UNIVERSE_POLICY" in
  reconstructed-all-qualifying|official-cafa3-targets) ;;
  *)
    echo "Unknown TARGET_UNIVERSE_POLICY: $TARGET_UNIVERSE_POLICY" >&2
    exit 1
    ;;
esac
case "$HISTORICAL_TEST_SOURCE" in
  official-groundtruth|raw-goa) ;;
  *)
    echo "Unknown HISTORICAL_TEST_SOURCE: $HISTORICAL_TEST_SOURCE" >&2
    echo "Expected official-groundtruth or raw-goa" >&2
    exit 1
    ;;
esac
if [ -z "$HISTORICAL_BENCHMARK_ONTOLOGY" ]; then
  if [ "$HISTORICAL_TEST_SOURCE" = "official-groundtruth" ]; then
    HISTORICAL_BENCHMARK_ONTOLOGY="deepgoplus-packaged"
  else
    HISTORICAL_BENCHMARK_ONTOLOGY="february-go-basic"
  fi
fi
case "$HISTORICAL_T1_ENDPOINT_POLICY" in
  assigned-date-proxy|snapshot-membership) ;;
  *) echo "Unknown HISTORICAL_T1_ENDPOINT_POLICY: $HISTORICAL_T1_ENDPOINT_POLICY" >&2; exit 1 ;;
esac
case "$HISTORICAL_BACKFILL_POLICY" in
  exclude-pre-t0|allow) ;;
  *) echo "Unknown HISTORICAL_BACKFILL_POLICY: $HISTORICAL_BACKFILL_POLICY" >&2; exit 1 ;;
esac
case "$HISTORICAL_BENCHMARK_ONTOLOGY" in
  february-go-basic|deepgoplus-packaged) ;;
  *) echo "Unknown HISTORICAL_BENCHMARK_ONTOLOGY: $HISTORICAL_BENCHMARK_ONTOLOGY" >&2; exit 1 ;;
esac
if [ "$HISTORICAL_TEST_SOURCE" = "official-groundtruth" ] \
  && [ "$TARGET_UNIVERSE_POLICY" != "official-cafa3-targets" ]; then
  echo "official-groundtruth requires TARGET_UNIVERSE_POLICY=official-cafa3-targets" >&2
  exit 1
fi
if [ "$HISTORICAL_TEST_SOURCE" = "official-groundtruth" ] \
  && [ "$HISTORICAL_TRAINING_SNAPSHOT" != "september-2016" ]; then
  echo "official-groundtruth validation requires HISTORICAL_TRAINING_SNAPSHOT=september-2016" >&2
  exit 1
fi

echo "=============================================================="
echo "CAFA3 Historical Validation"
echo "=============================================================="
echo "Run dir         : $RUN_DIR"
echo "Report copy dir : $REPORT_COPY_DIR"
echo "Scratch cleanup : $([ "$KEEP_SCRATCH" = "1" ] && echo disabled || echo enabled)"
echo "GOA pigz stream : $([ "${CAFA_BUILDER_USE_PIGZ}" = "0" ] && echo disabled || echo enabled-if-available)"
echo "GOA pre-unzip   : $([ "$DECOMPRESS_GOA" = "1" ] && echo enabled || echo disabled)"
echo "GOA progress    : every ${GOA_PROGRESS_INTERVAL} parsed rows"
echo "Training source : ${HISTORICAL_TRAINING_SNAPSHOT}"
echo "Target universe : ${TARGET_UNIVERSE_POLICY}"
echo "Test source     : ${HISTORICAL_TEST_SOURCE}"
echo "t1 endpoint    : ${HISTORICAL_T1_ENDPOINT_POLICY}"
echo "Backfill       : ${HISTORICAL_BACKFILL_POLICY}"
echo "Label ontology : ${HISTORICAL_BENCHMARK_ONTOLOGY}"
echo

echo "==> [1/8] Download historical raw snapshots into scratch"
if [ "$HISTORICAL_TEST_SOURCE" = "raw-goa" ]; then
  download "$GOA_T0_URL" "${RAW}/goa/goa_uniprot_all.gaf.163.gz"
  download "$GOA_T1_URL" "${RAW}/goa/goa_uniprot_all.gaf.172.gz"
  download "$UNIPROT_T0_URL" "${RAW}/uniprot/uniprot_sprot-only2017_02.tar.gz"
  download "$UNIPROT_T1_URL" "${RAW}/uniprot/uniprot_sprot-only2017_11.tar.gz"
elif [ "$HISTORICAL_TRAINING_SNAPSHOT" = "february-2017-legacy" ]; then
  download "$UNIPROT_T0_URL" "${RAW}/uniprot/uniprot_sprot-only2017_02.tar.gz"
else
  echo "  released ground truth selected; skipping GOA and 2017 target snapshots"
fi
case "$HISTORICAL_TRAINING_SNAPSHOT" in
  september-2016)
    TRAINING_SNAPSHOT_ID="UniProtKB-2016_08"
    TRAINING_SNAPSHOT_DATE="07-Sep-2016"
    TRAINING_ARCHIVE="${RAW}/uniprot/uniprot_sprot-only2016_08.tar.gz"
    stage_override_or_download \
      "${HISTORICAL_TRAINING_UNIPROT_ARCHIVE:-}" \
      "$UNIPROT_2016_08_URL" \
      "$TRAINING_ARCHIVE"
    verify_size_and_md5 "$TRAINING_ARCHIVE" "$UNIPROT_2016_08_SIZE" "$UNIPROT_2016_08_MD5"
    download_optional "$UNIPROT_2016_08_METALINK_URL" "${RAW}/metadata/uniprot_2016_08_RELEASE.metalink"
    ;;
  february-2017-legacy)
    TRAINING_SNAPSHOT_ID="UniProtKB-2017_02-legacy"
    TRAINING_SNAPSHOT_DATE="$UNIPROT_T0_RELEASE_DATE"
    TRAINING_ARCHIVE="${RAW}/uniprot/uniprot_sprot-only2017_02.tar.gz"
    ;;
esac
if [ "$HISTORICAL_TEST_SOURCE" = "raw-goa" ]; then
  download "$GO_T0_OBO_URL" "${RAW}/go/2017-02-01/go.obo"
  download "$GO_T0_BASIC_URL" "${RAW}/go/2017-02-01/go-basic.obo"
  download "$GO_T1_OBO_URL" "${RAW}/go/2017-11-01/go.obo"
  download "$GO_T1_BASIC_URL" "${RAW}/go/2017-11-01/go-basic.obo"
fi
DEEPGOPLUS_PICKLES_URL="${DEEPGOPLUS_PICKLES_URL:-$DEFAULT_DEEPGOPLUS_PICKLES_URL}"
OFFICIAL_CAFA3_ARCHIVE="${REFERENCE}/deepgoplus_pickles_reference.tar.gz"
stage_override_or_download \
  "${OFFICIAL_CAFA3_ARCHIVE_INPUT:-}" \
  "$DEEPGOPLUS_PICKLES_URL" \
  "$OFFICIAL_CAFA3_ARCHIVE"

echo "==> [2/8] Download best-effort metadata"
if [ "$HISTORICAL_TEST_SOURCE" = "raw-goa" ]; then
  download_optional "https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2017_02/relnotes.txt" "${RAW}/metadata/uniprot_2017_02_relnotes.txt"
  download_optional "https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2017_02/changes.html" "${RAW}/metadata/uniprot_2017_02_changes.html"
  download_optional "https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2017_11/relnotes.txt" "${RAW}/metadata/uniprot_2017_11_relnotes.txt"
  download_optional "https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2017_11/changes.html" "${RAW}/metadata/uniprot_2017_11_changes.html"
  download_optional "https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/UNIPROT/README" "${RAW}/metadata/goa_old_uniprot_README"
  download_optional "https://release.geneontology.org/2017-02-01/summary.txt" "${RAW}/metadata/go_2017-02-01_summary.txt"
  download_optional "https://release.geneontology.org/2017-11-01/summary.txt" "${RAW}/metadata/go_2017-11-01_summary.txt"
else
  echo "  raw-snapshot metadata not required for released-groundtruth mode"
fi

echo "==> [3/8] Extract UniProt tarballs in scratch"
if [ "$HISTORICAL_TEST_SOURCE" = "raw-goa" ]; then
  extract_tar_once "${RAW}/uniprot/uniprot_sprot-only2017_02.tar.gz" "${RAW}/uniprot/release_2017_02"
  extract_tar_once "${RAW}/uniprot/uniprot_sprot-only2017_11.tar.gz" "${RAW}/uniprot/release_2017_11"
elif [ "$HISTORICAL_TRAINING_SNAPSHOT" = "february-2017-legacy" ]; then
  extract_tar_once "${RAW}/uniprot/uniprot_sprot-only2017_02.tar.gz" "${RAW}/uniprot/release_2017_02"
fi
if [ "$HISTORICAL_TRAINING_SNAPSHOT" = "september-2016" ]; then
  extract_tar_once "$TRAINING_ARCHIVE" "${RAW}/uniprot/release_2016_08"
  TRAINING_UNIPROT_INPUT="$(find_uniprot_input "${RAW}/uniprot/release_2016_08")"
else
  TRAINING_UNIPROT_INPUT="$(find_uniprot_input "${RAW}/uniprot/release_2017_02")"
fi
echo "  Training UniProt input: ${TRAINING_UNIPROT_INPUT}"
if [ "$HISTORICAL_TEST_SOURCE" = "raw-goa" ]; then
  TARGET_T0_SPROT_INPUT="$(find_uniprot_input "${RAW}/uniprot/release_2017_02")"
  TARGET_T1_SPROT_INPUT="$(find_uniprot_input "${RAW}/uniprot/release_2017_11")"
  echo "  Target-map t0 Swiss-Prot input: ${TARGET_T0_SPROT_INPUT}"
  echo "  Target-map t1 Swiss-Prot input: ${TARGET_T1_SPROT_INPUT}"
fi

echo "==> [3a/8] Extract official CAFA3/DeepGOPlus archive"
mkdir -p "${REFERENCE}/deepgoplus_pickles"
extract_archive "$OFFICIAL_CAFA3_ARCHIVE" "${REFERENCE}/deepgoplus_pickles"
OFFICIAL_TARGET_FASTA="$(find "${REFERENCE}/deepgoplus_pickles" -type f -path '*/CAFA3_targets/targets_all.fasta' | sort | head -1 || true)"
OFFICIAL_TARGET_MAPPING_DIR="$(find "${REFERENCE}/deepgoplus_pickles" -type d -path '*/CAFA3_targets/Mapping files' | sort | head -1 || true)"
OFFICIAL_TRAINING_ANNOTATIONS="$(find "${REFERENCE}/deepgoplus_pickles" -type f -path '*/CAFA3_training_data/uniprot_sprot_exp.txt' | sort | head -1 || true)"
OFFICIAL_TEST_ANNOTATIONS="$(find "${REFERENCE}/deepgoplus_pickles" -type f -path '*/benchmark20171115/groundtruth/leafonly_all.txt' | sort | head -1 || true)"
OFFICIAL_GO_OBO="$(find "${REFERENCE}/deepgoplus_pickles" -type f -name 'go.obo' | sort | head -1 || true)"
for required in "$OFFICIAL_TARGET_FASTA" "$OFFICIAL_TRAINING_ANNOTATIONS" "$OFFICIAL_TEST_ANNOTATIONS" "$OFFICIAL_GO_OBO"; do
  [ -f "$required" ] || { echo "Missing official CAFA3 file in $OFFICIAL_CAFA3_ARCHIVE" >&2; exit 1; }
done
[ -d "$OFFICIAL_TARGET_MAPPING_DIR" ] || {
  echo "Missing official CAFA3 mapping directory in $OFFICIAL_CAFA3_ARCHIVE" >&2
  exit 1
}
echo "  Official target FASTA: ${OFFICIAL_TARGET_FASTA}"
echo "  Official target mappings: ${OFFICIAL_TARGET_MAPPING_DIR}"
echo "  Official training labels: ${OFFICIAL_TRAINING_ANNOTATIONS}"
echo "  Official test labels: ${OFFICIAL_TEST_ANNOTATIONS}"
echo "  DeepGOPlus ontology: ${OFFICIAL_GO_OBO}"

case "$HISTORICAL_BENCHMARK_ONTOLOGY" in
  february-go-basic)
    RAW_BENCHMARK_GO_OBO="${RAW}/go/2017-02-01/go-basic.obo"
    ;;
  deepgoplus-packaged)
    RAW_BENCHMARK_GO_OBO="$OFFICIAL_GO_OBO"
    ;;
esac

UNIPROT_T0_TREMBL="${RAW}/uniprot/release_2017_02/uniprot_trembl_cafa3_targets.dat.gz"
UNIPROT_T1_TREMBL="${RAW}/uniprot/release_2017_11/uniprot_trembl_cafa3_targets.dat.gz"
if [ "$HISTORICAL_TEST_SOURCE" = "official-groundtruth" ]; then
  echo "==> [3b/8] Released ground truth selected; target TrEMBL mapping is not required"
elif [ "$INCLUDE_TREMBL_TARGETS" = "1" ]; then
  echo "==> [3b/8] Stream-filter historical TrEMBL target populations"
  filter_trembl_archive_url "UniProt 2017_02 TrEMBL" "$UNIPROT_T0_FULL_URL" "$UNIPROT_T0_TREMBL"
  filter_trembl_archive_url "UniProt 2017_11 TrEMBL" "$UNIPROT_T1_FULL_URL" "$UNIPROT_T1_TREMBL"
else
  echo "==> [3b/8] INCLUDE_TREMBL_TARGETS=0; running Swiss-Prot-only target diagnostic"
fi

echo "==> [3c/8] Resolve GOA inputs"
if [ "$HISTORICAL_TEST_SOURCE" = "raw-goa" ]; then
  GOA_T0_INPUT="$(decompress_gaf_if_requested "${RAW}/goa/goa_uniprot_all.gaf.163.gz")"
  GOA_T1_INPUT="$(decompress_gaf_if_requested "${RAW}/goa/goa_uniprot_all.gaf.172.gz")"
  echo "  GOA t0 input: ${GOA_T0_INPUT}"
  echo "  GOA t1 input: ${GOA_T1_INPUT}"
else
  echo "  bypassed: released training and test annotation files are authoritative"
fi

echo "==> [4/8] Run benchmark builder"
BUILDER_PYTHONPATH="${REPO_ROOT}/benchmark_builders/contemporary_cafa/src${PYTHONPATH:+:${PYTHONPATH}}"
# The public monthly GO products only bracket the exact GOA snapshot dates.
# Preserve unmapped-term counts for the forensic comparison instead of aborting.
BUILDER_CMD=(
  "$PYTHON_BIN" -m cafa_benchmark_builder
  --source-mode snapshots
  --profile cafa3-reconstructed
  --uniprot-t0 "$TRAINING_UNIPROT_INPUT"
  --target-universe-policy "$TARGET_UNIVERSE_POLICY"
  --training-snapshot-id "$TRAINING_SNAPSHOT_ID"
  --training-snapshot-date "$TRAINING_SNAPSHOT_DATE"
  --training-reviewed-only
  --test-eligibility-policy ontology-no-knowledge
  --no-strict-qc
  --output-dir "$GENERATED"
  --report-dir "${REPORTS}/builder"
)
if [ "$HISTORICAL_TEST_SOURCE" = "official-groundtruth" ]; then
  BUILDER_CMD+=(
    --go-obo "$OFFICIAL_GO_OBO"
    --go-obo-t0 "$OFFICIAL_GO_OBO"
    --go-obo-t1 "$OFFICIAL_GO_OBO"
    --training-annotations-file "$OFFICIAL_TRAINING_ANNOTATIONS"
    --test-annotations-file "$OFFICIAL_TEST_ANNOTATIONS"
    --official-target-fasta "$OFFICIAL_TARGET_FASTA"
  )
else
  BUILDER_CMD+=(
    --uniprot-t1 "$TARGET_T1_SPROT_INPUT"
    --target-uniprot-t0 "$TARGET_T0_SPROT_INPUT"
    --target-uniprot-t1 "$TARGET_T1_SPROT_INPUT"
    --goa-t0 "$GOA_T0_INPUT"
    --goa-t1 "$GOA_T1_INPUT"
    --go-obo "$RAW_BENCHMARK_GO_OBO"
    --go-obo-t0 "${RAW}/go/2017-02-01/go-basic.obo"
    --go-obo-t1 "${RAW}/go/2017-11-01/go-basic.obo"
    --t0-cutoff 20170213
    --t1-cutoff 20171115
    --t1-endpoint-policy "$HISTORICAL_T1_ENDPOINT_POLICY"
  )
  if [ "$HISTORICAL_BACKFILL_POLICY" = "exclude-pre-t0" ]; then
    BUILDER_CMD+=(--exclude-t1-backfill)
  else
    BUILDER_CMD+=(--allow-t1-backfill)
  fi
  if [ "$INCLUDE_TREMBL_TARGETS" = "1" ]; then
    BUILDER_CMD+=(
      --target-uniprot-t0 "$UNIPROT_T0_TREMBL"
      --target-uniprot-t1 "$UNIPROT_T1_TREMBL"
      --include-unreviewed-targets
    )
  else
    BUILDER_CMD+=(--target-reviewed-only)
  fi
  if [ "$HISTORICAL_TRAINING_SNAPSHOT" = "september-2016" ]; then
    BUILDER_CMD+=(--training-annotations-file "$OFFICIAL_TRAINING_ANNOTATIONS")
  fi
  if [ "$TARGET_UNIVERSE_POLICY" = "official-cafa3-targets" ]; then
    BUILDER_CMD+=(
      --official-target-fasta "$OFFICIAL_TARGET_FASTA"
      --official-target-mapping-dir "$OFFICIAL_TARGET_MAPPING_DIR"
    )
  fi
fi
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
  REFERENCE_PICKLE_DIR="$(locate_complete_set "${REFERENCE}/deepgoplus_pickles" "${PICKLE_FILES[@]}" || true)"
  PICKLE_STATUS="downloaded/staged from ${DEEPGOPLUS_PICKLES_URL}"
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
if [ "$HISTORICAL_TEST_SOURCE" = "official-groundtruth" ]; then
  [ -n "$REFERENCE_PICKLE_DIR" ] || {
    echo "Official-groundtruth validation requires the released DeepGOPlus pickles" >&2
    exit 1
  }
  "$PYTHON_BIN" "${REPO_ROOT}/scripts/validation/validate_cafa3_official_test_artifacts.py" \
    --generated-dir "$GENERATED" \
    --reference-pickle-dir "$REFERENCE_PICKLE_DIR" \
    --reference-csv-dir "$REFERENCE_CSV_DIR" \
    --report "${REPORTS}/cafa3_official_test_artifact_gate.md" \
    2>&1 | tee "${LOGS}/official_test_artifact_gate.log"
fi
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
echo "==> Generated artefacts, reports and logs will be copied to: ${REPORT_COPY_DIR}"
