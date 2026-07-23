#!/usr/bin/env bash
# Hydrate accepted CAFA3 reproduction embeddings and compare them with Zijian's cache.

#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -l tmem=8G
#$ -l tscratch=40G
#$ -l scratch0free=60G
#$ -l h_rt=24:0:0
#$ -N cafa3_hyd_cmp
#$ -V
#$ -notify

set -Eeuo pipefail
export PYTHONDONTWRITEBYTECODE=1

usage() {
  cat <<'EOF'
Usage: qsub hpc_jobs/active/hpc_cafa3_hydrate_compare.sh \
  [--state-root PATH] [--artifact-catalog FILE] \
  [--destination-root PATH]

Hydrates every accepted array from the authenticated CAFA3 reproduction state
into job-local scratch, creates and round-trip validates one compact archive,
then compares that archived cache against Zijian's authenticated published
embeddings. The source state is read-only and is never retired or modified.
EOF
}

die() { echo "ERROR: $*" >&2; exit 2; }
require_value() { [[ $# -ge 2 && -n "$2" ]] || die "$1 requires a value"; }
git_in_dir() { local directory="$1"; shift; (cd "$directory" && git "$@"); }
sha256_file() {
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$1" | awk '{print $1}'
  else
    shasum -a 256 "$1" | awk '{print $1}'
  fi
}

STATE_ROOT=/SAN/bioinf/bmpfp/embedding_states/cafa3_full_reproduction
ARTIFACT_CATALOG_PATH=/SAN/bioinf/bmpfp/manifests/artifact_paths.tsv
DESTINATION_ROOT=/SAN/bioinf/bmpfp/diagnostics/cafa3_embedding_hydration_comparison
while [[ $# -gt 0 ]]; do
  case "$1" in
    --state-root) require_value "$@"; STATE_ROOT="$2"; shift 2 ;;
    --artifact-catalog) require_value "$@"; ARTIFACT_CATALOG_PATH="$2"; shift 2 ;;
    --destination-root) require_value "$@"; DESTINATION_ROOT="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) usage >&2; die "Unknown argument: $1" ;;
  esac
done

for path in "$STATE_ROOT" "$DESTINATION_ROOT"; do
  [[ "$path" == /* && "$path" != / ]] || die "Unsafe non-absolute path: $path"
done
[[ -f "$STATE_ROOT/contract.json" ]] || die "Missing state contract: $STATE_ROOT"
[[ -f "$STATE_ROOT/coverage.json" ]] || die "Missing state coverage: $STATE_ROOT"
[[ -f "$STATE_ROOT/pair_status.tsv" ]] || die "Missing state pair status: $STATE_ROOT"
[[ -f "$ARTIFACT_CATALOG_PATH" ]] || die "Missing artifact catalogue: $ARTIFACT_CATALOG_PATH"

JOB_TOKEN="${JOB_ID:-manual_$$}"
RUN_TAG="${JOB_TOKEN}_$(date -u +%Y%m%dT%H%M%SZ)"
WORK="/scratch0/cafa3_hydrate_compare_${JOB_TOKEN}"
FRAMEWORK_DIR="$WORK/Protein-Benchmark-Framework-Dissertation"
HYDRATED_CACHE="$WORK/hydrated_cache"
ROUNDTRIP_CACHE="$WORK/roundtrip_cache"
PUBLISHED_ROOT="$WORK/published"
SCRATCH_RESULT="$WORK/result"
FINAL_RESULT="$DESTINATION_ROOT/$RUN_TAG"
FAILED_RESULT="${FINAL_RESULT}.failed"
WORKFLOW_LOG="$WORK/workflow.log"
SUBMISSION_DIR="${SGE_O_WORKDIR:-$PWD}"
FRAMEWORK_REPO_URL="${FRAMEWORK_REPO_URL:-https://github.com/nadroj0-0/Protein-Benchmark-Framework-Dissertation.git}"
FRAMEWORK_COMMIT="${FRAMEWORK_COMMIT:-}"
WORK_OWNED=0
PUBLISHED=0

publish_results() {
  local status="$1" destination="$FINAL_RESULT"
  local staging="${FINAL_RESULT}.staging-${JOB_TOKEN}" copy_status=0
  [[ "$PUBLISHED" == 0 ]] || return 0
  if [[ "$status" != 0 ]]; then
    destination="$FAILED_RESULT"
    staging="${FAILED_RESULT}.staging-${JOB_TOKEN}"
  fi
  [[ ! -e "$destination" && ! -e "$staging" ]] || return 1
  mkdir -p "$staging/logs" || return 1
  [[ ! -d "$SCRATCH_RESULT" ]] || cp -a "$SCRATCH_RESULT/." "$staging/" || copy_status=$?
  [[ ! -f "$WORKFLOW_LOG" ]] || cp -p "$WORKFLOW_LOG" "$staging/logs/workflow.log" || copy_status=$?
  if [[ "$status" == 0 ]]; then
    [[ -f "$staging/HYDRATION_COMPARISON_COMPLETE.json" ]] || copy_status=1
    [[ -f "$staging/artifacts/cafa3_reproduction_hydrated_cache.tar.gz" ]] || copy_status=1
  else
    rm -f "$staging/HYDRATION_COMPARISON_COMPLETE.json"
    printf '{"complete":false,"workflow_exit_status":%s}\n' "$status" \
      > "$staging/WORKFLOW_FAILED.json" || copy_status=$?
  fi
  if [[ "$copy_status" == 0 ]]; then mv "$staging" "$destination" || copy_status=$?; fi
  if [[ "$copy_status" == 0 ]]; then
    PUBLISHED=1
    echo "Published CAFA3 hydration/comparison: $destination"
  elif [[ -d "$staging" && ! -L "$staging" ]]; then
    rm -rf -- "$staging"
  fi
  return "$copy_status"
}

cleanup() {
  local status=$? publish_status=0
  trap - EXIT
  set +e
  publish_results "$status" || publish_status=$?
  if [[ "$WORK_OWNED" == 1 && "$WORK" == /scratch0/cafa3_hydrate_compare_* && ! -L "$WORK" ]]; then
    cd "$HOME"
    rm -rf -- "$WORK"
  else
    echo "Refusing unsafe scratch cleanup: $WORK" >&2
    [[ "$status" != 0 ]] || status=1
  fi
  if [[ "$status" == 0 && "$publish_status" != 0 ]]; then status=$publish_status; fi
  exit "$status"
}
trap cleanup EXIT
trap 'echo "Received termination signal"; exit 130' INT TERM

[[ ! -e "$WORK" ]] || die "Scratch path already exists: $WORK"
mkdir -p "$WORK/tmp" "$SCRATCH_RESULT/reports" "$SCRATCH_RESULT/artifacts" "$DESTINATION_ROOT"
WORK_OWNED=1
export TMPDIR="$WORK/tmp" TMP="$WORK/tmp" TEMP="$WORK/tmp"

if [[ -z "$FRAMEWORK_COMMIT" ]]; then
  [[ -d "$SUBMISSION_DIR/.git" ]] || die "Submit from a clean framework checkout"
  [[ -z "$(git_in_dir "$SUBMISSION_DIR" status --porcelain)" ]] || die "Submission checkout is dirty"
  FRAMEWORK_COMMIT="$(git_in_dir "$SUBMISSION_DIR" rev-parse HEAD)"
fi
[[ "$FRAMEWORK_COMMIT" =~ ^[0-9a-fA-F]{40}$ ]] || die "Invalid FRAMEWORK_COMMIT"

echo "Host             : $(hostname -f 2>/dev/null || hostname)"
echo "Framework commit : $FRAMEWORK_COMMIT"
echo "State            : $STATE_ROOT"
echo "Destination      : $FINAL_RESULT"
echo "Scratch          : $WORK"

git clone --no-checkout "$FRAMEWORK_REPO_URL" "$FRAMEWORK_DIR"
git_in_dir "$FRAMEWORK_DIR" checkout --detach "$FRAMEWORK_COMMIT"

cd "$FRAMEWORK_DIR"
source scripts/reproduction_common.sh
export ARTIFACT_CATALOG="$ARTIFACT_CATALOG_PATH"
load_framework_paths "$FRAMEWORK_DIR"
add_mmfp_singularity_bind "$STATE_ROOT"
add_mmfp_singularity_bind "$DESTINATION_ROOT"
add_mmfp_singularity_bind "$(dirname "$ARTIFACT_CATALOG_PATH")"
activate_or_create_mmfp_env
PYTHON_BIN="$(command -v python)"
CONFIG="$FRAMEWORK_DIR/configs/pfp_benchmark_run.cafa3.json"

state_files=(contract.json coverage.json targets.tsv pair_status.tsv baseline_accepted.tsv)
printf 'file\tsha256_before\n' > "$SCRATCH_RESULT/reports/source_state_before.tsv"
for name in "${state_files[@]}"; do
  [[ -f "$STATE_ROOT/$name" ]] || die "Required state evidence is missing: $name"
  printf '%s\t%s\n' "$name" "$(sha256_file "$STATE_ROOT/$name")" \
    >> "$SCRATCH_RESULT/reports/source_state_before.tsv"
done

echo "==> [1/7] Hydrate every accepted regenerated array"
"$PYTHON_BIN" scripts/embeddings/manage_resumable_embedding_state.py hydrate \
  --state-root "$STATE_ROOT" \
  --output-cache-root "$HYDRATED_CACHE" \
  --preserve-evidence \
  --report "$SCRATCH_RESULT/reports/hydration.json"

echo "==> [2/7] Create one deterministic hydrated-cache archive"
"$PYTHON_BIN" scripts/embeddings/manage_embedding_archive.py create \
  --cache-root "$HYDRATED_CACHE" \
  --archive "$SCRATCH_RESULT/artifacts/cafa3_reproduction_hydrated_cache.tar.gz" \
  --config "$CONFIG" \
  --report "$SCRATCH_RESULT/reports/hydrated_archive_creation.json"

echo "==> [3/7] Extract and validate the archived cache used for comparison"
"$PYTHON_BIN" scripts/embeddings/manage_embedding_archive.py extract \
  --archive "$SCRATCH_RESULT/artifacts/cafa3_reproduction_hydrated_cache.tar.gz" \
  --output-cache-root "$ROUNDTRIP_CACHE" \
  --config "$CONFIG" \
  --report "$SCRATCH_RESULT/reports/hydrated_archive_roundtrip.json"
"$PYTHON_BIN" - \
  "$SCRATCH_RESULT/reports/hydrated_archive_creation.json" \
  "$SCRATCH_RESULT/reports/hydrated_archive_roundtrip.json" <<'PY'
import json
import sys

created = json.load(open(sys.argv[1], encoding="utf-8"))
extracted = json.load(open(sys.argv[2], encoding="utf-8"))
for key in ("archive_sha256", "member_count", "members_by_directory", "member_content_sha256"):
    if created[key] != extracted[key]:
        raise SystemExit(f"Hydrated archive round-trip mismatch: {key}")
PY
rm -rf -- "$HYDRATED_CACHE"

echo "==> [4/7] Authenticate and extract Zijian's published embeddings"
mkdir -p "$PUBLISHED_ROOT"
source scripts/artifact_catalog.sh
declare -A expected_sha=(
  [mmfp_embeddings_prott5.tar.gz]=30dd88fc4acbe3bc267bd8d5ae05e4d967fa7c169a6f063f12d2395fb0ffb00f
  [mmfp_embeddings_struct_ppi.tar.gz]=6d7a243f2c5e2149c162698b4e6a5e297731a4fea835d57b9049bb31f4af32de
  [mmfp_embeddings_text_temporal.tar.gz]=df1bf558fab1c018286a5b389665245917ba951f2c6cd8558546d0b1a3b47e36
)
printf 'archive\tsource\tsha256\n' > "$SCRATCH_RESULT/reports/published_archives.tsv"
for name in \
  mmfp_embeddings_prott5.tar.gz \
  mmfp_embeddings_struct_ppi.tar.gz \
  mmfp_embeddings_text_temporal.tar.gz; do
  source_path="$(resolve_artifact_path "$(zijian_embedding_artifact_id "$name")" "")"
  [[ -f "$source_path" ]] || die "Published embedding archive is missing: $name"
  observed_sha="$(sha256_file "$source_path")"
  [[ "$observed_sha" == "${expected_sha[$name]}" ]] || die "Published archive checksum mismatch: $name"
  tar -tzf "$source_path" >/dev/null
  tar -xzf "$source_path" -C "$PUBLISHED_ROOT"
  printf '%s\t%s\t%s\n' "$name" "$source_path" "$observed_sha" \
    >> "$SCRATCH_RESULT/reports/published_archives.tsv"
done
PUBLISHED_CACHE="$PUBLISHED_ROOT/data/embedding_cache"
declare -A expected_count=(
  [prott5]=69811
  [exp_text_embeddings_temporal]=69517
  [IF1]=67948
  [ppi]=58294
)
printf 'directory\tcount\n' > "$SCRATCH_RESULT/reports/published_cache_counts.tsv"
for directory in prott5 exp_text_embeddings_temporal IF1 ppi; do
  count="$(find "$PUBLISHED_CACHE/$directory" -maxdepth 1 -type f -name '*.npy' -print | wc -l | tr -d ' ')"
  [[ "$count" == "${expected_count[$directory]}" ]] || die "Published cache count mismatch: $directory"
  printf '%s\t%s\n' "$directory" "$count" >> "$SCRATCH_RESULT/reports/published_cache_counts.tsv"
done

echo "==> [5/7] Compare regenerated and published arrays"
"$PYTHON_BIN" scripts/diagnostics/compare_embeddings.py \
  --generated-cache-root "$ROUNDTRIP_CACHE" \
  --published-cache-root "$PUBLISHED_CACHE" \
  --out-csv "$SCRATCH_RESULT/reports/embedding_comparison.csv" \
  --out-json "$SCRATCH_RESULT/reports/embedding_comparison_summary.json" \
  > "$SCRATCH_RESULT/logs_embedding_comparison.txt" 2>&1
gzip -f "$SCRATCH_RESULT/reports/embedding_comparison.csv"

echo "==> [6/7] Prove the resumable source state remained unchanged"
printf 'file\tsha256_after\n' > "$SCRATCH_RESULT/reports/source_state_after.tsv"
for name in "${state_files[@]}"; do
  printf '%s\t%s\n' "$name" "$(sha256_file "$STATE_ROOT/$name")" \
    >> "$SCRATCH_RESULT/reports/source_state_after.tsv"
done
paste "$SCRATCH_RESULT/reports/source_state_before.tsv" \
      "$SCRATCH_RESULT/reports/source_state_after.tsv" \
  | awk -F '\t' 'NR == 1 {next} $1 != $3 || $2 != $4 {exit 1}' \
  || die "The source embedding state changed during hydration/comparison"

echo "==> [7/7] Build the compact completion manifest"
"$PYTHON_BIN" - \
  "$STATE_ROOT/contract.json" \
  "$STATE_ROOT/coverage.json" \
  "$SCRATCH_RESULT/reports/hydration.json" \
  "$SCRATCH_RESULT/reports/hydrated_archive_creation.json" \
  "$SCRATCH_RESULT/reports/embedding_comparison_summary.json" \
  "$SCRATCH_RESULT/HYDRATION_COMPARISON_COMPLETE.json" \
  "$FRAMEWORK_COMMIT" <<'PY'
import hashlib
import json
import sys
from datetime import datetime, timezone

contract_path, coverage_path, hydration_path, archive_path, comparison_path, output_path, commit = sys.argv[1:]
load = lambda path: json.load(open(path, encoding="utf-8"))
comparison_bytes = open(comparison_path, "rb").read()
value = {
    "schema_version": 1,
    "complete": True,
    "completed_at": datetime.now(timezone.utc).isoformat(),
    "framework_commit": commit,
    "state_contract_sha256": load(contract_path)["contract_sha256"],
    "state_coverage": load(coverage_path),
    "hydration": load(hydration_path),
    "hydrated_archive": load(archive_path),
    "comparison_summary_sha256": hashlib.sha256(comparison_bytes).hexdigest(),
    "comparison": load(comparison_path),
    "source_state_unchanged": True,
    "published_embeddings_used_for_training": False,
}
with open(output_path, "w", encoding="utf-8") as handle:
    json.dump(value, handle, indent=2, sort_keys=True)
    handle.write("\n")
PY

"$PYTHON_BIN" - \
  "$SCRATCH_RESULT/reports/embedding_comparison_summary.json" \
  "$SCRATCH_RESULT/COMPARISON_SUMMARY.md" <<'PY'
import json
import sys

summary = json.load(open(sys.argv[1], encoding="utf-8"))
lines = [
    "# CAFA3 Reproduction vs Published Embeddings",
    "",
    "The generated cache contains only arrays accepted by the resumable-state validator.",
    "Missing counts therefore describe unavailable generated arrays, not comparison failures.",
    "",
    "| Modality | Generated | Published | Common | Exact | Numeric | Different | Missing generated | Missing published | Mean cosine |",
    "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
]
for modality, values in sorted(summary.items()):
    statuses = values.get("statuses", {})
    cosine = values.get("cosine_mean")
    cosine_text = "NA" if cosine is None else f"{cosine:.9f}"
    lines.append(
        f"| {modality} | {values.get('generated_count', 0):,} | "
        f"{values.get('published_count', 0):,} | {values.get('common_count', 0):,} | "
        f"{statuses.get('exact_match', 0):,} | {statuses.get('numeric_match', 0):,} | "
        f"{statuses.get('different', 0):,} | {statuses.get('missing_generated', 0):,} | "
        f"{statuses.get('missing_published', 0):,} | {cosine_text} |"
    )
open(sys.argv[2], "w", encoding="utf-8").write("\n".join(lines) + "\n")
PY

printf 'path\tbytes\tsha256\n' > "$SCRATCH_RESULT/output_manifest.tsv"
while IFS= read -r -d '' path; do
  relative="${path#"$SCRATCH_RESULT"/}"
  printf '%s\t%s\t%s\n' "$relative" "$(wc -c < "$path" | tr -d ' ')" \
    "$(sha256_file "$path")" >> "$SCRATCH_RESULT/output_manifest.tsv"
done < <(find "$SCRATCH_RESULT" -type f ! -name output_manifest.tsv -print0 | sort -z)

rm -rf -- "$ROUNDTRIP_CACHE" "$PUBLISHED_ROOT"
publish_results 0
echo "Finished: $(date)"
