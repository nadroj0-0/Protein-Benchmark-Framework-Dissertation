#!/usr/bin/env bash
# Submit the complete 30% coverage / 30% identity clique benchmark campaign.

set -Eeuo pipefail

die() { echo "ERROR: $*" >&2; exit 2; }
git_in_dir() { local directory="$1"; shift; (cd "$directory" && git "$@"); }

FRAMEWORK_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
cd "$FRAMEWORK_ROOT"

CAMPAIGN_TAG="${CAMPAIGN_TAG:-multilinkage30c-prod-$(date -u +%Y%m%dT%H%M%SZ)}"
[[ "$CAMPAIGN_TAG" =~ ^[A-Za-z0-9._-]+$ && "$CAMPAIGN_TAG" =~ [A-Za-z0-9] ]] || \
  die "CAMPAIGN_TAG contains unsafe path characters"
SAN_ROOT="${SAN_ROOT:-/SAN/bioinf/bmpfp}"
CAMPAIGN_ROOT="${CAMPAIGN_ROOT:-$SAN_ROOT/reruns/fixed-evidence-codes/$CAMPAIGN_TAG}"
LOG_ROOT="${LOG_ROOT:-$HOME/pfp_multilinkage_logs/$CAMPAIGN_TAG}"
ARTIFACT_CATALOG="${ARTIFACT_CATALOG:-$SAN_ROOT/manifests/artifact_paths.tsv}"
DRY_RUN="${DRY_RUN:-0}"
MINIMUM_SAN_FREE_GB="${MINIMUM_SAN_FREE_GB:-20}"
[[ "$DRY_RUN" == 0 || "$DRY_RUN" == 1 ]] || die "DRY_RUN must be 0 or 1"
[[ "$MINIMUM_SAN_FREE_GB" =~ ^[1-9][0-9]*$ ]] || \
  die "MINIMUM_SAN_FREE_GB must be a positive integer"
[[ "$CAMPAIGN_ROOT" == /SAN/* && "$CAMPAIGN_ROOT" != /SAN ]] || \
  die "CAMPAIGN_ROOT must be a non-root SAN path"
[[ "$LOG_ROOT" == /* && "$LOG_ROOT" != / ]] || die "LOG_ROOT must be absolute"

MULTILINKAGE_INPUT="${MULTILINKAGE_INPUT:-$SAN_ROOT/frozen_inputs/supervisor_multilinkage/2026-09-02/coverage_30_identity_30/30_30_exp_go_subset.csv.gz}"
MULTILINKAGE_INPUT_SHA256="${MULTILINKAGE_INPUT_SHA256:-0733b8750ecfaf48d532e25b3abe4421411a9338a3b2c3f4e2a5687d1ff36fbc}"
COMMON_CACHE="${COMMON_CACHE:-$SAN_ROOT/reruns/fixed-evidence-codes/evidence12-20260824T072711Z/derived/homology_common_cache}"
UNIREF50_FASTA="${UNIREF50_FASTA:-$SAN_ROOT/frozen_inputs/uniref50/2026_02/uniref50.fasta.gz}"
HOMOLOGY_OBO="${HOMOLOGY_OBO:-$SAN_ROOT/frozen_inputs/ontology/2026-06-19/go-basic.obo}"
HOMOLOGY_ARCHIVE="${HOMOLOGY_ARCHIVE:-$SAN_ROOT/embeddings/homology/2026_02/identity_30/cluster-count-random/pair-resolved-paper-faithful/current-text-ledger_7132992_20260804T051124Z/final/homology_30_embedding_cache.tar.gz}"
HOMOLOGY_ARCHIVE_SHA256="${HOMOLOGY_ARCHIVE_SHA256:-a1cb0cf0fc2e0142a039a146bc86408090632c5eaed001db0f90f235644a188f}"
BENCHMARK_ID="homology-multilinkage-coverage30-identity30-option-c-evidence12"

BRANCH="$(git symbolic-ref --quiet --short HEAD || true)"
[[ "$BRANCH" == fixed-evidence-codes ]] || \
  die "Submit from fixed-evidence-codes, not ${BRANCH:-a detached checkout}"
[[ -z "$(git status --porcelain)" ]] || die "Framework checkout must be clean"
FRAMEWORK_COMMIT="$(git rev-parse HEAD)"
[[ "$FRAMEWORK_COMMIT" =~ ^[0-9a-f]{40}$ ]] || die "Could not resolve framework commit"
SHORT_COMMIT="${FRAMEWORK_COMMIT:0:12}"

command -v sha256sum >/dev/null 2>&1 || die "sha256sum is unavailable"
command -v gzip >/dev/null 2>&1 || die "gzip is unavailable"
for path in "$ARTIFACT_CATALOG" "$MULTILINKAGE_INPUT" "$UNIREF50_FASTA" \
  "$HOMOLOGY_OBO" "$HOMOLOGY_ARCHIVE" "$COMMON_CACHE/CACHE_COMPLETE.json"; do
  [[ -f "$path" && ! -L "$path" ]] || die "Required immutable input is missing or unsafe: $path"
done
gzip -t "$MULTILINKAGE_INPUT" || die "30/30 clique input failed gzip validation"
[[ "$(sha256sum "$MULTILINKAGE_INPUT" | awk '{print $1}')" == \
   "$MULTILINKAGE_INPUT_SHA256" ]] || die "30/30 clique input SHA-256 mismatch"
[[ "$(sha256sum "$HOMOLOGY_ARCHIVE" | awk '{print $1}')" == \
   "$HOMOLOGY_ARCHIVE_SHA256" ]] || die "Homology embedding archive SHA-256 mismatch"

free_kb="$(df -Pk "$SAN_ROOT" | awk 'END {print $4}')"
required_kb=$((MINIMUM_SAN_FREE_GB * 1024 * 1024))
(( free_kb >= required_kb )) || \
  die "SAN has less than ${MINIMUM_SAN_FREE_GB} GB free"
[[ ! -e "$CAMPAIGN_ROOT" ]] || die "Campaign root already exists: $CAMPAIGN_ROOT"

REMOTE_COMMIT="$(git ls-remote origin refs/heads/fixed-evidence-codes | awk 'NR==1 {print $1}')"
[[ "$REMOTE_COMMIT" == "$FRAMEWORK_COMMIT" ]] || \
  die "origin/fixed-evidence-codes does not match $FRAMEWORK_COMMIT"
if [[ "$DRY_RUN" == 0 ]]; then
  command -v qsub >/dev/null 2>&1 || die "qsub is unavailable"
fi

DRY_ROOT=""
if [[ "$DRY_RUN" == 1 ]]; then
  DRY_ROOT="$(mktemp -d "${TMPDIR:-/tmp}/multilinkage30-submit.XXXXXX")"
  trap 'rm -rf -- "$DRY_ROOT"' EXIT
  LEDGER="$DRY_ROOT/jobs.tsv"
else
  mkdir -p "$CAMPAIGN_ROOT/submission" "$LOG_ROOT"
  LEDGER="$CAMPAIGN_ROOT/submission/jobs.tsv"
fi
printf 'stage\tjob_name\tjob_id\thold_jid\texpected_output\tcommand\n' > "$LEDGER"

if [[ "$DRY_RUN" == 0 ]]; then
  cat > "$CAMPAIGN_ROOT/submission/campaign.env" <<EOF
campaign_tag=$CAMPAIGN_TAG
framework_commit=$FRAMEWORK_COMMIT
benchmark_id=$BENCHMARK_ID
multilinkage_coverage=0.3
multilinkage_identity=0.3
multilinkage_input=$MULTILINKAGE_INPUT
multilinkage_input_sha256=$MULTILINKAGE_INPUT_SHA256
embedding_archive=$HOMOLOGY_ARCHIVE
embedding_archive_sha256=$HOMOLOGY_ARCHIVE_SHA256
EOF
fi

NEXT_DRY_JOB=910000
SUBMITTED_JOB_ID=""
submit_job() {
  local stage="$1" name="$2" holds="$3" expected_output="$4"
  shift 4
  local command=(qsub -N "$name" -o "$LOG_ROOT") output job_id command_text
  [[ -z "$holds" ]] || command+=(-hold_jid "$holds")
  command+=("$@")
  printf -v command_text '%q ' "${command[@]}"
  if [[ "$DRY_RUN" == 1 ]]; then
    job_id="$NEXT_DRY_JOB"
    NEXT_DRY_JOB=$((NEXT_DRY_JOB + 1))
    printf 'DRY-RUN: %s\n' "$command_text" >&2
  else
    output="$("${command[@]}")"
    job_id="$(printf '%s\n' "$output" | sed -n \
      -e 's/.*job-array \([0-9][0-9]*\).*/\1/p' \
      -e 's/.*job \([0-9][0-9]*\).*/\1/p' | head -n 1)"
    [[ "$job_id" =~ ^[0-9]+$ ]] || die "Could not parse qsub output for $name: $output"
    echo "Submitted $name as $job_id${holds:+, holding on $holds}" >&2
  fi
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$stage" "$name" "$job_id" "$holds" "$expected_output" "$command_text" >> "$LEDGER"
  SUBMITTED_JOB_ID="$job_id"
}

BENCHMARK_ROOT="$CAMPAIGN_ROOT/benchmarks/homology_option_c"
BENCHMARK_ENV="FRAMEWORK_SOURCE_ROOT=$FRAMEWORK_ROOT,FRAMEWORK_REVISION=$FRAMEWORK_COMMIT,HOMOLOGY_RUNTIME_KIND=array,RUN_ID=$CAMPAIGN_TAG,RESULTS_ROOT=$BENCHMARK_ROOT,ARTIFACT_CATALOG=$ARTIFACT_CATALOG,HOMOLOGY_COMMON_PREPROCESSING_CACHE=$COMMON_CACHE,REQUIRE_HOMOLOGY_COMMON_CACHE=1,REQUIRE_HOMOLOGY_CLUSTER_CACHE=0,UNIREF_LEVEL=50,UNIREF50_FASTA=$UNIREF50_FASTA,MMSEQS_SENSITIVITY=4,MMSEQS_PROFILE=daniel-aligned-defaults,UNIPROT_SOURCE_SCOPE=sprot-and-trembl,SPLIT_POLICY=cluster-count-random,TRAINING_POPULATION=annotated-only,SEED=0,MIN_COUNT=50,MULTILINKAGE_CLUSTER_MEMBERSHIPS=$MULTILINKAGE_INPUT,MULTILINKAGE_CLUSTER_MEMBERSHIPS_SHA256=$MULTILINKAGE_INPUT_SHA256,MULTILINKAGE_COVERAGE=0.3,GO_OBO=$HOMOLOGY_OBO,MINIMUM_SCRATCH_GB=30,SCRATCH_SAFETY_MULTIPLIER=1,MMSEQS_WORK_MULTIPLIER=1,PUBLICATION_SAFETY_MULTIPLIER=1"
submit_job benchmark ml30_bench "" "$BENCHMARK_ROOT" \
  -S /bin/bash -cwd -j y -notify -t 1 -tc 1 -pe smp 4 \
  -l tmem=7G,tscratch=8G,h_rt=48:0:0 \
  -v "$BENCHMARK_ENV" \
  scripts/benchmark_generation/run_homology_cluster_runtime_hpc.sh
BENCHMARK_JOB="$SUBMITTED_JOB_ID"

BENCHMARK_TASK_ROOT="$BENCHMARK_ROOT/runtime_array/source_sprot-and-trembl/uniref50_sensitivity_4/multilinkage_coverage_0p3/framework_${SHORT_COMMIT}/run_${CAMPAIGN_TAG}/job_${BENCHMARK_JOB}/task_1_identity_30"
BENCHMARK_DIR="$BENCHMARK_TASK_ROOT/benchmark/source_sprot-and-trembl/framework_${SHORT_COMMIT}/uniref50_sensitivity_4/mmseqs_daniel-aligned-defaults/multilinkage_coverage_30/identity_30/cluster-count-random/annotated-only/seed_0/min_count_50"
EMBEDDING_EVIDENCE_DIR="$CAMPAIGN_ROOT/embedding_evidence/homology_option_c"
submit_job embedding-bind ml30_bind "$BENCHMARK_JOB" "$EMBEDDING_EVIDENCE_DIR" \
  -v "BENCHMARK_DIR=$BENCHMARK_DIR,BENCHMARK_ID=$BENCHMARK_ID,EMBEDDING_ARCHIVE=$HOMOLOGY_ARCHIVE,EMBEDDING_ARCHIVE_SHA256=$HOMOLOGY_ARCHIVE_SHA256,RUN_CONFIG=configs/pfp_benchmark_run.homology.json,EVIDENCE_OUTPUT=$EMBEDDING_EVIDENCE_DIR,FRAMEWORK_COMMIT=$FRAMEWORK_COMMIT" \
  hpc_jobs/active/hpc_bind_embedding_archive_evidence.sh
EMBEDDING_JOB="$SUBMITTED_JOB_ID"

MODEL_RESULTS_ROOT="$CAMPAIGN_ROOT/model_runs/homology_option_c/full"
MODEL_RUN_DIR="$MODEL_RESULTS_ROOT/run"
submit_job model ml30_full "$EMBEDDING_JOB" "$MODEL_RUN_DIR" \
  -l tmem=32G,tscratch=40G,scratch0free=60G,h_rt=96:0:0 \
  -v "FRAMEWORK_COMMIT=$FRAMEWORK_COMMIT,RUN_TAG=run" \
  hpc_jobs/active/hpc_pfp_benchmark.sh \
  --benchmark-id "$BENCHMARK_ID" \
  --benchmark-dir "$BENCHMARK_DIR" \
  --benchmark-evidence "$BENCHMARK_DIR/validation_report.json" \
  --benchmark-evidence "$BENCHMARK_DIR/output_manifest.json" \
  --benchmark-evidence "$BENCHMARK_DIR/RUN_COMPLETE.json" \
  --embedding-cache-archive "$HOMOLOGY_ARCHIVE" \
  --embedding-evidence "$EMBEDDING_EVIDENCE_DIR/coverage.json" \
  --embedding-evidence "$EMBEDDING_EVIDENCE_DIR/contract.json" \
  --embedding-evidence "$EMBEDDING_EVIDENCE_DIR/targets.tsv" \
  --embedding-evidence "$EMBEDDING_EVIDENCE_DIR/pair_status.tsv" \
  --require-embedding-evidence \
  --obo-file "$HOMOLOGY_OBO" \
  --results-root "$MODEL_RESULTS_ROOT" \
  --config configs/pfp_benchmark_run.homology.json \
  --execution-mode train-eval \
  --modality-mode full \
  --aspect BPO --aspect CCO --aspect MFO \
  --seed 42 --num-workers 0 --capture-predictions
MODEL_JOB="$SUBMITTED_JOB_ID"

ANALYSIS_ROOT="$CAMPAIGN_ROOT/analyses/homology_option_c"
SOURCE_LABEL="multilinkage-coverage30-identity30-full-evidence12"
DIAGNOSTICS_DIR="$ANALYSIS_ROOT/diagnostics"
submit_job diagnostics ml30_diag "$MODEL_JOB" "$DIAGNOSTICS_DIR" \
  -pe smp 1 -l tmem=24G,tscratch=10G,scratch0free=15G,h_rt=24:0:0 \
  -v "FRAMEWORK_COMMIT=$FRAMEWORK_COMMIT" \
  hpc_jobs/active/hpc_contemporary_followup_analysis.sh \
  --analysis diagnostics --source-run "$MODEL_RUN_DIR" --obo "$HOMOLOGY_OBO" \
  --source-label "$SOURCE_LABEL" --output-dir "$DIAGNOSTICS_DIR"
DIAGNOSTICS_JOB="$SUBMITTED_JOB_ID"

PREDICTION_CAPTURE_DIR="$ANALYSIS_ROOT/prediction_capture"
submit_job prediction-capture ml30_cap "$MODEL_JOB" "$PREDICTION_CAPTURE_DIR" \
  -l tmem=32G,tscratch=40G,scratch0free=60G,h_rt=24:0:0 \
  -v "FRAMEWORK_COMMIT=$FRAMEWORK_COMMIT" \
  hpc_jobs/active/hpc_contemporary_followup_prediction_capture.sh \
  --source-run "$MODEL_RUN_DIR" --cache-archive "$HOMOLOGY_ARCHIVE" \
  --obo "$HOMOLOGY_OBO" --benchmark-id "$BENCHMARK_ID" \
  --source-label "$SOURCE_LABEL" --output-dir "$PREDICTION_CAPTURE_DIR"
PREDICTION_CAPTURE_JOB="$SUBMITTED_JOB_ID"

CALIBRATION_DIR="$ANALYSIS_ROOT/calibration"
submit_job calibration ml30_cal "$PREDICTION_CAPTURE_JOB" "$CALIBRATION_DIR" \
  -pe smp 1 -l tmem=24G,tscratch=10G,scratch0free=15G,h_rt=24:0:0 \
  -v "FRAMEWORK_COMMIT=$FRAMEWORK_COMMIT" \
  hpc_jobs/active/hpc_contemporary_followup_analysis.sh \
  --analysis calibration --capture-pair-dir "$PREDICTION_CAPTURE_DIR" \
  --source-run "$MODEL_RUN_DIR" --obo "$HOMOLOGY_OBO" \
  --source-label "$SOURCE_LABEL" --output-dir "$CALIBRATION_DIR"
CALIBRATION_JOB="$SUBMITTED_JOB_ID"

FINAL_HOLDS="$DIAGNOSTICS_JOB,$PREDICTION_CAPTURE_JOB,$CALIBRATION_JOB"
submit_job closure ml30_close "$FINAL_HOLDS" "$CAMPAIGN_ROOT/RUN_COMPLETE.json" \
  -v "CAMPAIGN_ROOT=$CAMPAIGN_ROOT,BENCHMARK_DIR=$BENCHMARK_DIR,EMBEDDING_EVIDENCE_DIR=$EMBEDDING_EVIDENCE_DIR,MODEL_RUN_DIR=$MODEL_RUN_DIR,DIAGNOSTICS_DIR=$DIAGNOSTICS_DIR,PREDICTION_CAPTURE_DIR=$PREDICTION_CAPTURE_DIR,CALIBRATION_DIR=$CALIBRATION_DIR,SUBMISSION_LEDGER=$LEDGER,BENCHMARK_ID=$BENCHMARK_ID,FRAMEWORK_COMMIT=$FRAMEWORK_COMMIT" \
  hpc_jobs/active/hpc_validate_multilinkage_campaign.sh
CLOSURE_JOB="$SUBMITTED_JOB_ID"

if [[ "$DRY_RUN" == 1 ]]; then
  echo "Dry-run submission ledger:"
  cat "$LEDGER"
else
  printf 'submitted_at_utc=%s\nclosure_job=%s\n' \
    "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$CLOSURE_JOB" \
    >> "$CAMPAIGN_ROOT/submission/campaign.env"
  echo "Campaign root: $CAMPAIGN_ROOT"
  echo "Submission ledger: $LEDGER"
  echo "Closure job: $CLOSURE_JOB"
fi
