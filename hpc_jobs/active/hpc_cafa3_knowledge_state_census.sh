#!/usr/bin/env bash
# Classify Zijian's published CAFA3 test rows with the official CAFA3 type lists.

#$ -S /bin/bash
#$ -l tmem=8G
#$ -l scratch0free=8G
#$ -l tscratch=8G
#$ -l h_rt=4:0:0
#$ -pe smp 1
#$ -j y
#$ -N c3_cohort
#$ -V
#$ -notify

set -Eeuo pipefail
export PYTHONDONTWRITEBYTECODE=1

usage() {
  cat <<'EOF'
Usage:
  qsub hpc_jobs/active/hpc_cafa3_knowledge_state_census.sh \
    --output-dir /SAN/.../cafa3_knowledge_state_census

This CPU-only job classifies Zijian's exact published CAFA3 test rows using
the organizer-provided type1 (no-knowledge) and type2 (limited-knowledge)
lists. It also reports root-only/non-root truth counts for all nine CSVs.
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
while [[ $# -gt 0 ]]; do
  case "$1" in
    --output-dir) require_value "$@"; OUTPUT_DIR="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

[[ "$OUTPUT_DIR" == /SAN/* ]] || die "--output-dir must be an absolute SAN path"
[[ ! -e "$OUTPUT_DIR" ]] || die "Output directory already exists: $OUTPUT_DIR"

PUBLISHED_CSV_DIR="${PUBLISHED_CSV_DIR:-/SAN/bioinf/bmpfp/reference_artifacts/canonical_cafa3}"
OFFICIAL_CAFA_ARCHIVE="${OFFICIAL_CAFA_ARCHIVE:-/SAN/bioinf/bmpfp/reference_artifacts/deepgoplus/data-cafa.tar.gz}"

for prefix in bp cc mf; do
  for split in training validation test; do
    [[ -s "$PUBLISHED_CSV_DIR/$prefix-$split.csv" ]] || \
      die "Published CAFA3 CSV is missing or empty: $PUBLISHED_CSV_DIR/$prefix-$split.csv"
  done
done
[[ -s "$OFFICIAL_CAFA_ARCHIVE" ]] || \
  die "Official CAFA3 archive is missing or empty: $OFFICIAL_CAFA_ARCHIVE"

JOB_TOKEN="${JOB_ID:-manual_$$}"
WORK="/scratch0/c3_cohort_${JOB_TOKEN}"
FRAMEWORK_DIR="$WORK/Protein-Benchmark-Framework-Dissertation"
ANALYSIS_OUTPUT="$WORK/analysis"
LOG_FILE="$WORK/cafa3_knowledge_state_census.log"
PUBLISH_STAGE="${OUTPUT_DIR}.staging-${JOB_TOKEN}"
PUBLISH_LOCK="${OUTPUT_DIR}.publish-lock"
FRAMEWORK_REPO_URL="${FRAMEWORK_REPO_URL:-https://github.com/nadroj0-0/Protein-Benchmark-Framework-Dissertation.git}"
FRAMEWORK_COMMIT="${FRAMEWORK_COMMIT:-}"
LOCK_HELD=0

cleanup() {
  local status=$?
  trap - EXIT
  set +e
  if [[ "$LOCK_HELD" == "1" && -d "$PUBLISH_LOCK" && ! -L "$PUBLISH_LOCK" ]]; then
    rmdir -- "$PUBLISH_LOCK"
  fi
  [[ ! -d "$PUBLISH_STAGE" || -L "$PUBLISH_STAGE" ]] || rm -rf -- "$PUBLISH_STAGE"
  [[ ! -d "$WORK" || -L "$WORK" || "$WORK" != /scratch0/c3_cohort_* ]] || rm -rf -- "$WORK"
  exit "$status"
}
trap cleanup EXIT
trap 'echo "Received termination signal"; exit 130' INT TERM

[[ ! -e "$WORK" ]] || die "Scratch path already exists: $WORK"
[[ ! -e "$PUBLISH_STAGE" ]] || die "Publication stage already exists"
[[ ! -e "$PUBLISH_LOCK" ]] || die "Publication lock already exists"
mkdir -p "$WORK" "$(dirname "$OUTPUT_DIR")"

echo "Host             : $(hostname)"
echo "Job ID           : ${JOB_ID:-manual}"
echo "Analysis         : official CAFA3 knowledge-state census"
echo "Published CSVs   : $PUBLISHED_CSV_DIR"
echo "Official archive : $OFFICIAL_CAFA_ARCHIVE"
echo "SAN output       : $OUTPUT_DIR"
echo "Started          : $(date -Is)"

git clone --no-checkout "$FRAMEWORK_REPO_URL" "$FRAMEWORK_DIR"
if [[ -n "$FRAMEWORK_COMMIT" ]]; then
  [[ "$FRAMEWORK_COMMIT" =~ ^[0-9a-fA-F]{40}$ ]] || \
    die "FRAMEWORK_COMMIT must be a full commit when supplied"
  git_in_dir "$FRAMEWORK_DIR" checkout --detach "$FRAMEWORK_COMMIT"
else
  git_in_dir "$FRAMEWORK_DIR" checkout main
  FRAMEWORK_COMMIT="$(git_in_dir "$FRAMEWORK_DIR" rev-parse HEAD)"
fi
echo "Framework commit : $FRAMEWORK_COMMIT"

cd "$FRAMEWORK_DIR"
python3 scripts/diagnostics/build_cafa3_knowledge_state_census.py \
  --published-csv-dir "$PUBLISHED_CSV_DIR" \
  --official-cafa-archive "$OFFICIAL_CAFA_ARCHIVE" \
  --output-dir "$ANALYSIS_OUTPUT" \
  >"$LOG_FILE" 2>&1

echo "==> Publishing census atomically"
mkdir -p "$PUBLISH_STAGE/logs"
cp -a "$ANALYSIS_OUTPUT" "$PUBLISH_STAGE/analysis"
cp -p "$LOG_FILE" "$PUBLISH_STAGE/logs/"
python3 scripts/model_execution/manage_output_manifest.py write \
  --root "$PUBLISH_STAGE" --include-nested-control-files
MANIFEST_SHA256="$(python3 -c 'import hashlib,pathlib,sys; print(hashlib.sha256(pathlib.Path(sys.argv[1]).read_bytes()).hexdigest())' "$PUBLISH_STAGE/output_manifest.json")"
python3 -c 'import json,pathlib,sys; pathlib.Path(sys.argv[1]).write_text(json.dumps({"complete":True,"analysis_kind":"cafa3_official_knowledge_state_census","framework_commit":sys.argv[2],"manifest":"output_manifest.json","manifest_sha256":sys.argv[3]},indent=2)+"\n")' \
  "$PUBLISH_STAGE/WORKFLOW_COMPLETE.json" "$FRAMEWORK_COMMIT" "$MANIFEST_SHA256"
python3 scripts/model_execution/manage_output_manifest.py verify \
  --root "$PUBLISH_STAGE" --include-nested-control-files
mkdir "$PUBLISH_LOCK"
LOCK_HELD=1
[[ ! -e "$OUTPUT_DIR" ]] || die "Output directory appeared during publication"
mv -T "$PUBLISH_STAGE" "$OUTPUT_DIR"
rmdir "$PUBLISH_LOCK"
LOCK_HELD=0

echo "Published CAFA3 knowledge-state census: $OUTPUT_DIR"
echo "Finished         : $(date -Is)"
