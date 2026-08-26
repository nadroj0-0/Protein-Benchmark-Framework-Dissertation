#!/usr/bin/env bash
# Submit the released-artifact control and raw-database CAFA3 reconstruction models.

set -Eeuo pipefail

CAMPAIGN_ROOT="${CAMPAIGN_ROOT:-$HOME/cafa3_reconstruction_end_to_end_20260826}"
SOURCE_CACHE="${SOURCE_CACHE:-/SAN/bioinf/bmpfp/diagnostics/cafa3_embedding_hydration_comparison/7125100_20260730T154052Z/artifacts/cafa3_reproduction_hydrated_cache.tar.gz}"
SOURCE_TARGETS="${SOURCE_TARGETS:-/SAN/bioinf/bmpfp/diagnostics/cafa3_embedding_hydration_comparison/7125100_20260730T154052Z/evidence/canonical_cafa3_published_csvs_b7d1e7d/targets.tsv}"

RAW_ROOT="$CAMPAIGN_ROOT/benchmarks/7061922_raw_database"
ARTIFACT_ROOT="$CAMPAIGN_ROOT/benchmarks/7061973_released_artifacts"
RAW_BENCHMARK="$RAW_ROOT/generated"
ARTIFACT_BENCHMARK="$ARTIFACT_ROOT/generated"
RAW_CHECKSUMS="$RAW_ROOT/builder/output_checksums.sha256"
ARTIFACT_CHECKSUMS="$ARTIFACT_ROOT/builder/output_checksums.sha256"
RAW_MANIFEST="$RAW_ROOT/builder/build_manifest.json"
ARTIFACT_MANIFEST="$ARTIFACT_ROOT/builder/build_manifest.json"
RAW_OBO="$CAMPAIGN_ROOT/ontologies/go-basic-2017-02-01.obo"
ARTIFACT_OBO="$CAMPAIGN_ROOT/ontologies/deepgoplus-cafa3-go.obo"
RAW_CACHE_RESULTS="$CAMPAIGN_ROOT/cache_results"
MODEL_RESULTS="$CAMPAIGN_ROOT/model_results"
RAW_CACHE_TAG=7061922_raw_database_cache
ARTIFACT_MODEL_TAG=7061973_released_artifacts_full_seed42
RAW_MODEL_TAG=7061922_raw_database_full_seed42
CONFIG="$PWD/configs/pfp_benchmark_run.cafa3_reconstructed.json"
EXPECTED="$PWD/configs/pfp_cafa3_published_metrics.json"
CACHE_WRAPPER="$PWD/hpc_jobs/active/hpc_cafa3_reconstruction_cache.sh"
MODEL_WRAPPER="$PWD/hpc_jobs/active/hpc_pfp_benchmark.sh"

die() { echo "ERROR: $*" >&2; exit 2; }
sha256_file() {
  if command -v sha256sum >/dev/null 2>&1; then sha256sum "$1" | awk '{print $1}'; else shasum -a 256 "$1" | awk '{print $1}'; fi
}
verify_benchmark() {
  local benchmark="$1" ledger="$2" name got expected
  for name in bp-training.csv bp-validation.csv bp-test.csv cc-training.csv cc-validation.csv cc-test.csv mf-training.csv mf-validation.csv mf-test.csv; do
    [[ -f "$benchmark/$name" && ! -L "$benchmark/$name" ]] || die "Missing benchmark CSV: $benchmark/$name"
    got="$(sha256_file "$benchmark/$name")"
    expected="$(awk -v n="$name" '{p=$2; sub(/^.*\//,"",p); if (p==n) {print $1; exit}}' "$ledger")"
    [[ -n "$expected" && "$got" == "$expected" ]] || die "Benchmark checksum mismatch: $benchmark/$name"
  done
}

[[ -d .git ]] || die "Run from the framework checkout root"
[[ "$(git symbolic-ref --quiet --short HEAD)" == fixed-evidence-codes ]] || die "Submit from fixed-evidence-codes"
[[ -z "$(git status --porcelain)" ]] || die "Framework checkout must be clean"
FRAMEWORK_COMMIT="$(git rev-parse HEAD)"
[[ "$FRAMEWORK_COMMIT" == "$(git rev-parse origin/fixed-evidence-codes)" ]] || die "Local branch is not at origin/fixed-evidence-codes"
for path in "$RAW_CHECKSUMS" "$ARTIFACT_CHECKSUMS" "$RAW_MANIFEST" "$ARTIFACT_MANIFEST" "$RAW_OBO" "$ARTIFACT_OBO" "$SOURCE_CACHE" "$SOURCE_TARGETS" "$CONFIG" "$EXPECTED" "$CACHE_WRAPPER" "$MODEL_WRAPPER"; do
  [[ -f "$path" && ! -L "$path" ]] || die "Required input is missing or unsafe: $path"
done
[[ "$(sha256_file "$RAW_OBO")" == c560d540174e4f451b49f4f9d45653842bdf0b8195e7cd4ae914a38d3769e182 ]] || die "Raw benchmark ontology checksum mismatch"
[[ "$(sha256_file "$ARTIFACT_OBO")" == e19da74998201c977bc45433a9d37ccb8fdc018ca4b57306ba1c1d4da6e12a9e ]] || die "Artifact benchmark ontology checksum mismatch"
[[ "$(sha256_file "$SOURCE_CACHE")" == c6cdcafd00b0cb871a50beb8cd649ce5c13d5882fcde9c3663a19a86905b9e87 ]] || die "Source cache checksum mismatch"
[[ "$(sha256_file "$SOURCE_TARGETS")" == 34874b6c3f70c6c6ece20019c57a932ad6027ea1bcde113ad07ff86348cfcec6 ]] || die "Source target manifest checksum mismatch"
verify_benchmark "$RAW_BENCHMARK" "$RAW_CHECKSUMS"
verify_benchmark "$ARTIFACT_BENCHMARK" "$ARTIFACT_CHECKSUMS"
for path in "$RAW_CACHE_RESULTS/$RAW_CACHE_TAG" "$RAW_CACHE_RESULTS/${RAW_CACHE_TAG}.failed" "$MODEL_RESULTS/$ARTIFACT_MODEL_TAG" "$MODEL_RESULTS/${ARTIFACT_MODEL_TAG}.failed" "$MODEL_RESULTS/$RAW_MODEL_TAG" "$MODEL_RESULTS/${RAW_MODEL_TAG}.failed"; do
  [[ ! -e "$path" ]] || die "Refusing to overwrite an existing campaign result: $path"
done
mkdir -p "$RAW_CACHE_RESULTS" "$MODEL_RESULTS" "$CAMPAIGN_ROOT/submissions"

RAW_CACHE_JOB="$(qsub -terse -N c3rawemb -v "FRAMEWORK_COMMIT=$FRAMEWORK_COMMIT" "$CACHE_WRAPPER" \
  --benchmark-id cafa3-raw-database-7061922 \
  --benchmark-dir "$RAW_BENCHMARK" \
  --benchmark-checksums "$RAW_CHECKSUMS" \
  --source-cache-archive "$SOURCE_CACHE" \
  --source-targets "$SOURCE_TARGETS" \
  --results-root "$RAW_CACHE_RESULTS" \
  --run-tag "$RAW_CACHE_TAG")"
RAW_CACHE_JOB="${RAW_CACHE_JOB%%.*}"

ARTIFACT_MODEL_JOB="$(qsub -terse -N c3artpfp -v "FRAMEWORK_COMMIT=$FRAMEWORK_COMMIT,RUN_TAG=$ARTIFACT_MODEL_TAG" "$MODEL_WRAPPER" \
  --benchmark-id cafa3-released-artifacts-7061973 \
  --benchmark-dir "$ARTIFACT_BENCHMARK" \
  --embedding-cache-archive "$SOURCE_CACHE" \
  --obo-file "$ARTIFACT_OBO" \
  --results-root "$MODEL_RESULTS" \
  --execution-mode train-eval \
  --config "$CONFIG" \
  --modality-mode full \
  --seed 42 \
  --num-workers 0 \
  --capture-predictions \
  --expected-metrics "$EXPECTED" \
  --benchmark-evidence "$ARTIFACT_MANIFEST" \
  --benchmark-evidence "$ARTIFACT_CHECKSUMS")"
ARTIFACT_MODEL_JOB="${ARTIFACT_MODEL_JOB%%.*}"

RAW_MODEL_JOB="$(qsub -terse -N c3rawpfp -hold_jid "$RAW_CACHE_JOB" -v "FRAMEWORK_COMMIT=$FRAMEWORK_COMMIT,RUN_TAG=$RAW_MODEL_TAG" "$MODEL_WRAPPER" \
  --benchmark-id cafa3-raw-database-7061922 \
  --benchmark-dir "$RAW_BENCHMARK" \
  --embedding-cache-archive "$RAW_CACHE_RESULTS/$RAW_CACHE_TAG/artifacts/cafa3_reconstruction_embedding_cache.tar.gz" \
  --obo-file "$RAW_OBO" \
  --results-root "$MODEL_RESULTS" \
  --execution-mode train-eval \
  --config "$CONFIG" \
  --modality-mode full \
  --seed 42 \
  --num-workers 0 \
  --capture-predictions \
  --expected-metrics "$EXPECTED" \
  --benchmark-evidence "$RAW_MANIFEST" \
  --benchmark-evidence "$RAW_CHECKSUMS")"
RAW_MODEL_JOB="${RAW_MODEL_JOB%%.*}"

SUBMISSION_FILE="$CAMPAIGN_ROOT/submissions/submission_$(date -u +%Y%m%dT%H%M%SZ).tsv"
{
  printf 'role\tjob_id\thold_job_id\trun_tag\tframework_commit\n'
  printf 'raw_cache\t%s\t\t%s\t%s\n' "$RAW_CACHE_JOB" "$RAW_CACHE_TAG" "$FRAMEWORK_COMMIT"
  printf 'artifact_model\t%s\t\t%s\t%s\n' "$ARTIFACT_MODEL_JOB" "$ARTIFACT_MODEL_TAG" "$FRAMEWORK_COMMIT"
  printf 'raw_model\t%s\t%s\t%s\t%s\n' "$RAW_MODEL_JOB" "$RAW_CACHE_JOB" "$RAW_MODEL_TAG" "$FRAMEWORK_COMMIT"
} > "$SUBMISSION_FILE"

printf 'Raw cache job      : %s\n' "$RAW_CACHE_JOB"
printf 'Artifact model job : %s\n' "$ARTIFACT_MODEL_JOB"
printf 'Raw model job      : %s (holds on %s)\n' "$RAW_MODEL_JOB" "$RAW_CACHE_JOB"
printf 'Submission ledger  : %s\n' "$SUBMISSION_FILE"
