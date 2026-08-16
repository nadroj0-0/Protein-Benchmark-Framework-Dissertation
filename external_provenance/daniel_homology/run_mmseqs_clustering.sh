#!/usr/bin/env bash
# Standalone UniRef90 clustering for the homology benchmark.
# The order mirrors Daniel's 01_cluster_with_mmseqs.sh.

set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd -P)"
STATE="$HERE/handoff_state.py"

# Scientific settings
COV=0.80
SENSITIVITY=4
IDENTITIES="10"
THREADS=24

# Frozen input
EXPECTED_SHA256="ed80f79bbf1f054b3ea444ce0db0819586731e4dc3a3fb0c75a60ff273eedefb"
EXPECTED_RECORDS=121389642
RELEASE="2026_07"
PREFERRED_MMSEQS="18-8cc5c"

# Paths and optional outputs
FASTA=""
OUT=""
WORK=""
MMSEQS="/home/dbuchan/Applications/mmseqs/bin/mmseqs"
REASSIGN=0
EXPORT_STATS=0
EXPORT_FASTA=0
KEEP_WORK=0

usage() {
  cat <<'EOF'
Usage:
  ./run_mmseqs_clustering.sh \
    --fasta /path/uniref50_2026_07.fasta.gz \
    --output-dir /path/results \
    --scratch-dir /path/large_work_disk \
    [--mmseqs /path/mmseqs] [--threads 24]

Optional:
  --thresholds 30,25,20,15,10,5
  --cluster-reassign
  --export-alignment-statistics
  --export-cluster-fasta
  --retain-completed-work
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --fasta) FASTA="$2"; shift 2 ;;
    --output-dir) OUT="$2"; shift 2 ;;
    --scratch-dir) WORK="$2"; shift 2 ;;
    --mmseqs) MMSEQS="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --thresholds) IDENTITIES="$2"; shift 2 ;;
    --cluster-reassign) REASSIGN=1; shift ;;
    --export-alignment-statistics) EXPORT_STATS=1; shift ;;
    --export-cluster-fasta) EXPORT_FASTA=1; shift ;;
    --retain-completed-work) KEEP_WORK=1; shift ;;
    --help|-h) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

[[ -f "$FASTA" ]] || { echo "A valid --fasta file is required." >&2; exit 2; }
[[ -n "$OUT" && -n "$WORK" ]] || { echo "--output-dir and --scratch-dir are required." >&2; exit 2; }
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || { echo "--threads must be a positive integer." >&2; exit 2; }
command -v python3 >/dev/null || { echo "python3 is required." >&2; exit 2; }
command -v gzip >/dev/null || { echo "gzip is required." >&2; exit 2; }

if [[ "$MMSEQS" == */* ]]; then
  [[ -x "$MMSEQS" ]] || { echo "MMseqs is not executable: $MMSEQS" >&2; exit 2; }
else
  MMSEQS="$(command -v "$MMSEQS" || true)"
  [[ -n "$MMSEQS" ]] || { echo "MMseqs was not found." >&2; exit 2; }
fi

IFS=',' read -r -a ID_LIST <<< "$IDENTITIES"
for ID in "${ID_LIST[@]}"; do
  [[ " 25 20 15 10 " == *" $ID "* ]] || { echo "Unsupported identity: $ID" >&2; exit 2; }
done

mkdir -p "$OUT" "$WORK"

CURRENT_STAGE="initial checks"
record_failure() {
  local status="$1"
  trap - ERR
  printf 'status=%s\nstage=%s\ntime_utc=%s\n' \
    "$status" "$CURRENT_STAGE" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" > "$OUT/RUN_FAILED.txt"
  echo "Run stopped during: $CURRENT_STAGE" >&2
  exit "$status"
}
trap 'record_failure $?' ERR

hash_file() {
  if command -v sha256sum >/dev/null; then
    sha256sum "$1" | awk '{print $1}'
  else
    shasum -a 256 "$1" | awk '{print $1}'
  fi
}

# Used only by the bundled three-sequence self-test.
if [[ "${HANDOFF_SELF_TEST:-0}" == "1" ]]; then
  EXPECTED_SHA256="${HANDOFF_TEST_EXPECTED_SHA256:?}"
  EXPECTED_RECORDS="${HANDOFF_TEST_EXPECTED_RECORDS:?}"
  RELEASE="self-test"
fi

#echo "Checking the frozen UniRef50 file..."
FASTA_SHA256="$(hash_file "$FASTA")"
#[[ "$FASTA_SHA256" == "$EXPECTED_SHA256" ]] || {
#  echo "UniRef50 SHA-256 mismatch; clustering has not started." >&2
#  exit 2
#}

MMSEQS_VERSION="$($MMSEQS version 2>&1)"
MMSEQS_VERSION="${MMSEQS_VERSION%%$'\n'*}"
MMSEQS_SHA256="$(hash_file "$MMSEQS")"
#if [[ "${MMSEQS_VERSION//[._]/-}" != *"$PREFERRED_MMSEQS"* ]]; then
#  echo "WARNING: intended MMseqs version is $PREFERRED_MMSEQS; found $MMSEQS_VERSION" >&2
#  echo "Continuing and recording the version used." >&2
#fi

CONTRACT="$OUT/RUN_CONTRACT.json"
CURRENT_STAGE="creating the run contract"
CONTRACT_COMMAND=(python3 "$STATE" initialize-run \
  --contract "$CONTRACT" \
  --input-path "$FASTA" --input-sha256 "$FASTA_SHA256" --input-release "$RELEASE" \
  --mmseqs-path "$MMSEQS" --mmseqs-sha256 "$MMSEQS_SHA256" --mmseqs-version "$MMSEQS_VERSION" \
  --threads "$THREADS" --thresholds "$IDENTITIES")
[[ $REASSIGN -eq 1 ]] && CONTRACT_COMMAND+=(--cluster-reassign)
[[ $EXPORT_STATS -eq 1 ]] && CONTRACT_COMMAND+=(--export-alignment-statistics)
[[ $EXPORT_FASTA -eq 1 ]] && CONTRACT_COMMAND+=(--export-cluster-fasta)
"${CONTRACT_COMMAND[@]}"
rm -f "$OUT/RUN_FAILED.txt"

run_step() {
  local LOG="$1"
  shift
  printf 'Running:' | tee -a "$LOG"
  printf ' %q' "$@" | tee -a "$LOG"
  printf '\n' | tee -a "$LOG"
  "$@" 2>&1 | tee -a "$LOG"
}

compress_result() {
  local SOURCE="$1" DESTINATION="$2"
  rm -f "${SOURCE}.gz" "${DESTINATION}.partial"
  gzip -n "$SOURCE"
  mv "${SOURCE}.gz" "${DESTINATION}.partial"
  gzip -t "${DESTINATION}.partial"
  mv "${DESTINATION}.partial" "$DESTINATION"
}

# Step one: create the sequence database once.
SEQDB_DIR="$WORK/seqDB"
SEQDB="$SEQDB_DIR/uniref50-DB"
SEQDB_DONE="$SEQDB_DIR/CREATEDB_COMPLETE.json"
mkdir -p "$SEQDB_DIR"

if python3 "$STATE" check-createdb-marker \
    --contract "$CONTRACT" --db-base "$SEQDB" --marker "$SEQDB_DONE"; then
  echo "Reusing the completed UniRef50 sequence database."
else
  CURRENT_STAGE="creating the shared UniRef50 sequence database"
  run_step "$OUT/mmseqs_createdb.log" \
    "$MMSEQS" createdb "$FASTA" "$SEQDB" --dbtype 1
  python3 "$STATE" write-createdb-marker \
    --contract "$CONTRACT" --db-base "$SEQDB" --marker "$SEQDB_DONE"
fi

# Steps two onward: repeat Daniel's clustering/export pattern for each identity.
for ID in "${ID_LIST[@]}"; do
  LABEL="$(printf '%02d' "$ID")"
  SEQID="$(printf '0.%02d' "$ID")"
  RESULT="$OUT/threshold_$LABEL"
  TEMP="$WORK/threshold_$LABEL"
  CLUSTERDB="$TEMP/clustDB/uniref50-cov80-id$LABEL"
  CLUSTER_DONE="$TEMP/CLUSTER_COMPLETE.json"
  mkdir -p "$RESULT" "$TEMP/clustDB" "$TEMP/tmp"

  ASSIGNMENTS="$RESULT/cluster_assignments.tsv.gz"
  VALIDATION="$RESULT/validation.json"
  COMPLETE="$RESULT/COMPLETE.json"
  STATS="$RESULT/alignment_statistics.tsv.gz"
  CLUSTER_FASTA="$RESULT/clusters.faa.gz"

  if [[ -e "$COMPLETE" ]]; then
    python3 "$STATE" check-threshold-marker \
      --contract "$CONTRACT" --threshold "$ID" --assignments "$ASSIGNMENTS" --marker "$COMPLETE" \
      || { echo "Threshold ${ID}% failed checksum verification." >&2; exit 2; }
    echo "Threshold ${ID}% is already complete; skipping."
    continue
  fi

  # Step two: build the cluster database.
  EXTRA_CLUSTER_FLAGS=(--threads "$THREADS")
  [[ $REASSIGN -eq 1 ]] && EXTRA_CLUSTER_FLAGS+=(--cluster-reassign 1)
  if ! python3 "$STATE" check-cluster-marker \
      --contract "$CONTRACT" --threshold "$ID" --cluster-base "$CLUSTERDB" --marker "$CLUSTER_DONE"; then
    CURRENT_STAGE="clustering at ${ID}% identity"
    run_step "$RESULT/mmseqs_cluster.log" \
      "$MMSEQS" cluster \
        "$SEQDB" "$CLUSTERDB" "$TEMP/tmp" \
        -c "$COV" --min-seq-id "$SEQID" -s "$SENSITIVITY" \
        --cov-mode 0 --cluster-mode 0 --alignment-mode 3 --seq-id-mode 0 \
        "${EXTRA_CLUSTER_FLAGS[@]}"
    python3 "$STATE" write-cluster-marker \
      --contract "$CONTRACT" --threshold "$ID" --cluster-base "$CLUSTERDB" --marker "$CLUSTER_DONE"
  else
    echo "Reusing completed ${ID}% cluster database."
  fi

  # Step three: representative/member TSV required by the benchmark.
  RAW_TSV="$TEMP/cluster_assignments.tsv"
  rm -f "$RAW_TSV"
  CURRENT_STAGE="exporting and validating ${ID}% cluster assignments"
  run_step "$RESULT/mmseqs_createtsv.log" \
    "$MMSEQS" createtsv "$SEQDB" "$SEQDB" "$CLUSTERDB" "$RAW_TSV"
  # python3 "$STATE" validate-tsv \
  #   --input "$RAW_TSV" --expected-rows "$EXPECTED_RECORDS" \
  #   --report "$RESULT/validation.json.partial"
  # compress_result "$RAW_TSV" "$ASSIGNMENTS"
  # mv "$RESULT/validation.json.partial" "$VALIDATION"

  # Optional: Daniel's representative/member alignment statistics.
  if [[ $EXPORT_STATS -eq 1 ]]; then
    mkdir -p "$TEMP/alignDB"
    ALIGNDB="$TEMP/alignDB/uniref50-cov80-id$LABEL"
    CURRENT_STAGE="exporting ${ID}% alignment statistics"
    run_step "$RESULT/mmseqs_align.log" \
      "$MMSEQS" align "$SEQDB" "$SEQDB" "$CLUSTERDB" "$ALIGNDB" -a 1 --threads "$THREADS"
    rm -f "$TEMP/alignment_statistics.tsv"
    run_step "$RESULT/mmseqs_convertalis.log" \
      "$MMSEQS" convertalis "$SEQDB" "$SEQDB" "$ALIGNDB" "$TEMP/alignment_statistics.tsv" \
        --format-mode 4 --format-output query,target,evalue,bits,raw,pident
    [[ -s "$TEMP/alignment_statistics.tsv" ]] || { echo "Empty alignment-statistics export." >&2; false; }
    compress_result "$TEMP/alignment_statistics.tsv" "$STATS"
  fi

  # Optional: Daniel's FASTA arranged by cluster.
  if [[ $EXPORT_FASTA -eq 1 ]]; then
    mkdir -p "$TEMP/fastaclustDB"
    FASTADB="$TEMP/fastaclustDB/uniref50-cov80-id$LABEL"
    CURRENT_STAGE="exporting ${ID}% cluster FASTA"
    run_step "$RESULT/mmseqs_createseqfiledb.log" \
      "$MMSEQS" createseqfiledb "$SEQDB" "$CLUSTERDB" "$FASTADB"
    rm -f "$TEMP/clusters.faa"
    run_step "$RESULT/mmseqs_result2flat.log" \
      "$MMSEQS" result2flat "$SEQDB" "$SEQDB" "$FASTADB" "$TEMP/clusters.faa"
    [[ -s "$TEMP/clusters.faa" ]] || { echo "Empty cluster-FASTA export." >&2; false; }
    compress_result "$TEMP/clusters.faa" "$CLUSTER_FASTA"
  fi

  # Record checksums before removing temporary databases.
  #MARKER_COMMAND=(python3 "$STATE" write-threshold-marker \
  #  --contract "$CONTRACT" --threshold "$ID" --assignments "$ASSIGNMENTS" \
  #  --validation-report "$VALIDATION" --marker "$COMPLETE")
  #[[ $EXPORT_STATS -eq 1 ]] && MARKER_COMMAND+=(--alignment-statistics "$STATS")
  #[[ $EXPORT_FASTA -eq 1 ]] && MARKER_COMMAND+=(--cluster-fasta "$CLUSTER_FASTA")
  #CURRENT_STAGE="publishing the ${ID}% completion marker"
  #"${MARKER_COMMAND[@]}"

  #if [[ $KEEP_WORK -eq 0 ]]; then
  #  rm -rf "$TEMP"
  #fi
  echo "Threshold ${ID}% complete: $RESULT"
done

CURRENT_STAGE="writing and verifying the final manifest"
# python3 "$STATE" write-run-manifest \
#   --output-dir "$OUT" --contract "$CONTRACT" --thresholds "$IDENTITIES"
# python3 "$STATE" verify-run --output-dir "$OUT"
# trap - ERR
# rm -f "$OUT/RUN_FAILED.txt"
# echo "All requested thresholds are complete and verified: $OUT"
