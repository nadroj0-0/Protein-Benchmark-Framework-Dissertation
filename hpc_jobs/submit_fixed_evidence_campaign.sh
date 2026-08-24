#!/usr/bin/env bash
# Submit the complete corrected-evidence benchmark, model and analysis DAG.

set -Eeuo pipefail

die() { echo "ERROR: $*" >&2; exit 2; }
git_in_dir() { local directory="$1"; shift; (cd "$directory" && git "$@"); }

FRAMEWORK_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
cd "$FRAMEWORK_ROOT"

CAMPAIGN_TAG="${CAMPAIGN_TAG:-evidence12-$(date -u +%Y%m%dT%H%M%SZ)}"
[[ "$CAMPAIGN_TAG" =~ ^[A-Za-z0-9._-]+$ && "$CAMPAIGN_TAG" =~ [A-Za-z0-9] ]] || \
  die "CAMPAIGN_TAG contains unsafe path characters"
SAN_ROOT="${SAN_ROOT:-/SAN/bioinf/bmpfp}"
CAMPAIGN_ROOT="${CAMPAIGN_ROOT:-$SAN_ROOT/reruns/fixed-evidence-codes/$CAMPAIGN_TAG}"
LOG_ROOT="${LOG_ROOT:-$HOME/pfp_fixed_evidence_logs/$CAMPAIGN_TAG}"
ARTIFACT_CATALOG="${ARTIFACT_CATALOG:-$SAN_ROOT/manifests/artifact_paths.tsv}"
DRY_RUN="${DRY_RUN:-0}"
ALLOW_DIRTY_PREVIEW="${ALLOW_DIRTY_PREVIEW:-0}"
MINIMUM_SAN_FREE_GB="${MINIMUM_SAN_FREE_GB:-40}"
for value in DRY_RUN ALLOW_DIRTY_PREVIEW; do
  [[ "${!value}" == 0 || "${!value}" == 1 ]] || die "$value must be 0 or 1"
done
[[ "$CAMPAIGN_ROOT" == /* && "$CAMPAIGN_ROOT" != / ]] || die "CAMPAIGN_ROOT must be absolute"
[[ "$LOG_ROOT" == /* && "$LOG_ROOT" != / ]] || die "LOG_ROOT must be absolute"

BRANCH="$(git branch --show-current)"
[[ "$BRANCH" == fixed-evidence-codes ]] || die "Submit from fixed-evidence-codes, not $BRANCH"
if [[ "$DRY_RUN" == 0 || "$ALLOW_DIRTY_PREVIEW" == 0 ]]; then
  [[ -z "$(git status --porcelain)" ]] || die "Framework checkout must be clean"
fi
FRAMEWORK_COMMIT="$(git rev-parse HEAD)"
SHORT_COMMIT="${FRAMEWORK_COMMIT:0:12}"

EXPECTED_CODES="EXP,IDA,IPI,IMP,IGI,IEP,HTP,HDA,HMP,HGI,HEP,IGC"
PYTHONPATH="$FRAMEWORK_ROOT/benchmark_builders/contemporary_cafa/src:$FRAMEWORK_ROOT/benchmark_builders/homology_cluster/src" \
python3 - "$EXPECTED_CODES" <<'PY'
import sys
from cafa_benchmark_builder.config import SUPERVISOR_EXP_CODES
from homology_cluster_benchmark.config import SUPERVISOR_EVIDENCE_CODES

expected = frozenset(sys.argv[1].split(","))
if SUPERVISOR_EXP_CODES != expected:
    raise SystemExit(f"Contemporary evidence policy differs: {sorted(SUPERVISOR_EXP_CODES)}")
if SUPERVISOR_EVIDENCE_CODES != expected:
    raise SystemExit(f"Homology evidence policy differs: {sorted(SUPERVISOR_EVIDENCE_CODES)}")
PY

if [[ "$DRY_RUN" == 0 ]]; then
  command -v qsub >/dev/null 2>&1 || die "qsub is unavailable"
  REMOTE_COMMIT="$(git ls-remote origin refs/heads/fixed-evidence-codes | awk 'NR==1 {print $1}')"
  [[ "$REMOTE_COMMIT" == "$FRAMEWORK_COMMIT" ]] || \
    die "origin/fixed-evidence-codes does not match $FRAMEWORK_COMMIT"
fi

GLOBAL_ARCHIVE="$SAN_ROOT/embeddings/contemporary/2025_01_to_2026_02_supervisor/variants/text-cutoff-2025-03-08__ppi-paper-faithful/finalized_pfp_cache/contemporary_embedding_cache.tar.gz"
GLOBAL_ARCHIVE_SHA="8c579d492a9e9ee93a3539f722e479da9f917c6aada2f0238893b263066aa70e"
NK_ARCHIVE="$SAN_ROOT/embeddings/contemporary/2025_01_to_2026_02_supervisor_nk_lk/pair-resolved-paper-faithful/temporal-text-ledger_7132994_20260804T051124Z/final/homology_30_embedding_cache.tar.gz"
NK_ARCHIVE_SHA="5584a749e785d7762d9ba9a8d520ddfcf016eace5bb365d199ef737bd4660768"
HOMOLOGY_ARCHIVE="$SAN_ROOT/embeddings/homology/2026_02/identity_30/cluster-count-random/pair-resolved-paper-faithful/current-text-ledger_7132992_20260804T051124Z/final/homology_30_embedding_cache.tar.gz"
HOMOLOGY_ARCHIVE_SHA="a1cb0cf0fc2e0142a039a146bc86408090632c5eaed001db0f90f235644a188f"
TEMPORAL_OBO="$SAN_ROOT/frozen_inputs/ontology/2025-02-06/go-basic.obo"
HOMOLOGY_OBO="$SAN_ROOT/frozen_inputs/ontology/2026-06-19/go-basic.obo"
CLUSTER_CACHE="$SAN_ROOT/derived_inputs/homology/2026_02/mmseqs_cluster_assignments/uniref50_sensitivity_4"
COMMON_CACHE="$CAMPAIGN_ROOT/derived/homology_common_cache"

CLUSTER_BASE="$CLUSTER_CACHE/uniref50_2026_02"
FULL_30="$CLUSTER_BASE/identity_30/contract_e3e44179dd219d65/cluster_assignments.tsv.gz"
FULL_25="$CLUSTER_BASE/identity_25/contract_a6fd69df69a78cb0/cluster_assignments.tsv.gz"
FULL_20="$CLUSTER_BASE/identity_20/contract_fe783953dbe1c275/cluster_assignments.tsv.gz"
FULL_15="$CLUSTER_BASE/identity_15/contract_fb6fcfd60847780e/cluster_assignments.tsv.gz"
FULL_10="$CLUSTER_BASE/identity_10/contract_9206700890e1ccce/cluster_assignments.tsv.gz"
FULL_05="$CLUSTER_BASE/identity_05/contract_abd0c5c4beaae966/cluster_assignments.tsv.gz"

if [[ "$DRY_RUN" == 0 ]]; then
  for path in "$ARTIFACT_CATALOG" "$GLOBAL_ARCHIVE" "$NK_ARCHIVE" "$HOMOLOGY_ARCHIVE" \
    "$TEMPORAL_OBO" "$HOMOLOGY_OBO" "$FULL_30" "$FULL_25" "$FULL_20" \
    "$FULL_15" "$FULL_10" "$FULL_05"; do
    [[ -f "$path" && ! -L "$path" ]] || die "Required immutable input is missing or unsafe: $path"
  done
  [[ -d "$CLUSTER_CACHE" && ! -L "$CLUSTER_CACHE" ]] || die "Cluster cache root is missing"
  free_kb="$(df -Pk "$SAN_ROOT" | awk 'END {print $4}')"
  required_kb=$((MINIMUM_SAN_FREE_GB * 1024 * 1024))
  (( free_kb >= required_kb )) || die "SAN has less than ${MINIMUM_SAN_FREE_GB} GB free"
fi
[[ ! -e "$CAMPAIGN_ROOT" ]] || die "Campaign root already exists: $CAMPAIGN_ROOT"
mkdir -p "$CAMPAIGN_ROOT/submission" "$LOG_ROOT"

LEDGER="$CAMPAIGN_ROOT/submission/jobs.tsv"
printf 'stage\tjob_name\tjob_id\thold_jid\texpected_output\tcommand\n' > "$LEDGER"
printf 'campaign_tag=%s\nframework_commit=%s\nevidence_codes=%s\n' \
  "$CAMPAIGN_TAG" "$FRAMEWORK_COMMIT" "$EXPECTED_CODES" \
  > "$CAMPAIGN_ROOT/submission/campaign.env"

NEXT_DRY_JOB=900000
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
    echo "DRY-RUN: $command_text" >&2
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

GLOBAL_BENCH_ROOT="$CAMPAIGN_ROOT/benchmarks/global_nk"
GLOBAL_BENCH_DIR="$GLOBAL_BENCH_ROOT/benchmark/outputs"
GLOBAL_REPORT_DIR="$GLOBAL_BENCH_ROOT/benchmark/reports"
NK_BENCH_ROOT="$CAMPAIGN_ROOT/benchmarks/nk_lk"
NK_BENCH_DIR="$NK_BENCH_ROOT/benchmark/outputs"
NK_REPORT_DIR="$NK_BENCH_ROOT/benchmark/reports"
HOMOLOGY_BENCH_ROOT="$CAMPAIGN_ROOT/benchmarks/homology"

submit_job benchmark-global e12_gbench "" "$GLOBAL_BENCH_DIR" \
  -v "PROFILE=supervisor,FRAMEWORK_REVISION=$FRAMEWORK_COMMIT,RESULTS_ROOT=$GLOBAL_BENCH_ROOT,RUN_TAG=benchmark" \
  hpc_jobs/active/hpc_contemporary_temporal_benchmark.sh --artifact-catalog "$ARTIFACT_CATALOG"
GLOBAL_BENCH_JOB="$SUBMITTED_JOB_ID"

submit_job benchmark-nk-lk e12_nbench "" "$NK_BENCH_DIR" \
  -v "PROFILE=supervisor-nk-lk,FRAMEWORK_REVISION=$FRAMEWORK_COMMIT,RESULTS_ROOT=$NK_BENCH_ROOT,RUN_TAG=benchmark" \
  hpc_jobs/active/hpc_contemporary_temporal_benchmark.sh --artifact-catalog "$ARTIFACT_CATALOG"
NK_BENCH_JOB="$SUBMITTED_JOB_ID"

submit_job homology-common-cache e12_hcache "" "$COMMON_CACHE" \
  -v "FRAMEWORK_SOURCE_ROOT=$FRAMEWORK_ROOT,SAN_ROOT=$SAN_ROOT,HOMOLOGY_COMMON_PREPROCESSING_CACHE=$COMMON_CACHE" \
  hpc_jobs/active/hpc_homology_uniref50_common_cache.sh
HOMOLOGY_CACHE_JOB="$SUBMITTED_JOB_ID"

submit_job benchmark-homology-array e12_hbench "$HOMOLOGY_CACHE_JOB" "$HOMOLOGY_BENCH_ROOT" \
  -v "FRAMEWORK_SOURCE_ROOT=$FRAMEWORK_ROOT,FRAMEWORK_REVISION=$FRAMEWORK_COMMIT,RUN_ID=$CAMPAIGN_TAG,RESULTS_ROOT=$HOMOLOGY_BENCH_ROOT,HOMOLOGY_COMMON_PREPROCESSING_CACHE=$COMMON_CACHE,HOMOLOGY_CLUSTER_CACHE_ROOT=$CLUSTER_CACHE,REQUIRE_HOMOLOGY_CLUSTER_CACHE=1" \
  hpc_jobs/active/hpc_homology_cluster_runtime_array_12core_uniref50.sh --artifact-catalog "$ARTIFACT_CATALOG"
HOMOLOGY_BENCH_JOB="$SUBMITTED_JOB_ID"

HOMOLOGY_JOB_ROOT="$HOMOLOGY_BENCH_ROOT/runtime_array/source_sprot-and-trembl/uniref50_sensitivity_4/framework_${SHORT_COMMIT}/run_${CAMPAIGN_TAG}/job_${HOMOLOGY_BENCH_JOB}"
GLOBAL_EVIDENCE="$CAMPAIGN_ROOT/embedding_evidence/global_nk"
NK_EVIDENCE="$CAMPAIGN_ROOT/embedding_evidence/nk_lk"

submit_job bind-global e12_gbind "$GLOBAL_BENCH_JOB" "$GLOBAL_EVIDENCE" \
  -v "BENCHMARK_DIR=$GLOBAL_BENCH_DIR,BENCHMARK_ID=contemporary-2025-01-to-2026-02-supervisor,EMBEDDING_ARCHIVE=$GLOBAL_ARCHIVE,EMBEDDING_ARCHIVE_SHA256=$GLOBAL_ARCHIVE_SHA,RUN_CONFIG=configs/pfp_benchmark_run.temporal.json,EVIDENCE_OUTPUT=$GLOBAL_EVIDENCE,FRAMEWORK_COMMIT=$FRAMEWORK_COMMIT" \
  hpc_jobs/active/hpc_bind_embedding_archive_evidence.sh
GLOBAL_BIND_JOB="$SUBMITTED_JOB_ID"

submit_job bind-nk-lk e12_nbind "$NK_BENCH_JOB" "$NK_EVIDENCE" \
  -v "BENCHMARK_DIR=$NK_BENCH_DIR,BENCHMARK_ID=contemporary-2025-01-to-2026-02-supervisor-nk-lk,EMBEDDING_ARCHIVE=$NK_ARCHIVE,EMBEDDING_ARCHIVE_SHA256=$NK_ARCHIVE_SHA,RUN_CONFIG=configs/pfp_benchmark_run.temporal_nk_lk.json,EVIDENCE_OUTPUT=$NK_EVIDENCE,FRAMEWORK_COMMIT=$FRAMEWORK_COMMIT" \
  hpc_jobs/active/hpc_bind_embedding_archive_evidence.sh
NK_BIND_JOB="$SUBMITTED_JOB_ID"

HOMOLOGY_BENCHMARKS=()
HOMOLOGY_BINDINGS=()
HOMOLOGY_BIND_JOBS=()
identities=(30 25 20 15 10 5)
tasks=(1 2 3 4 5 6)
for index in "${!identities[@]}"; do
  identity="${identities[$index]}" task="${tasks[$index]}" identity_dir="$identity"
  [[ "$identity" != 5 ]] || identity_dir=05
  benchmark="$HOMOLOGY_JOB_ROOT/task_${task}_identity_${identity}/benchmark/source_sprot-and-trembl/framework_${SHORT_COMMIT}/uniref50_sensitivity_4/mmseqs_daniel-aligned-defaults/identity_${identity_dir}/cluster-count-random/annotated-only/seed_0/min_count_50"
  binding="$CAMPAIGN_ROOT/embedding_evidence/homology/identity_${identity_dir}"
  HOMOLOGY_BENCHMARKS[$identity]="$benchmark"
  HOMOLOGY_BINDINGS[$identity]="$binding"
  submit_job "bind-homology-$identity" "e12_h${identity}b" "$HOMOLOGY_BENCH_JOB" "$binding" \
    -v "BENCHMARK_DIR=$benchmark,BENCHMARK_ID=homology-uniref50-identity-${identity}-fixed-evidence,EMBEDDING_ARCHIVE=$HOMOLOGY_ARCHIVE,EMBEDDING_ARCHIVE_SHA256=$HOMOLOGY_ARCHIVE_SHA,RUN_CONFIG=configs/pfp_benchmark_run.homology.json,EVIDENCE_OUTPUT=$binding,FRAMEWORK_COMMIT=$FRAMEWORK_COMMIT" \
    hpc_jobs/active/hpc_bind_embedding_archive_evidence.sh
  HOMOLOGY_BIND_JOBS[$identity]="$SUBMITTED_JOB_ID"
done

GLOBAL_MODEL_JOBS=()
GLOBAL_RUNS=()
global_modes=(full sequence-only sequence-text sequence-structure sequence-ppi)
for index in "${!global_modes[@]}"; do
  mode="${global_modes[$index]}"
  case "$mode" in
    full) job_name=e12_gfull ;;
    sequence-only) job_name=e12_gseq ;;
    sequence-text) job_name=e12_gtxt ;;
    sequence-structure) job_name=e12_gstr ;;
    sequence-ppi) job_name=e12_gppi ;;
  esac
  results_root="$CAMPAIGN_ROOT/model_runs/global_nk/$mode"
  run_path="$results_root/run"
  submit_job "model-global-$mode" "$job_name" "$GLOBAL_BIND_JOB" "$run_path" \
    -v "FRAMEWORK_COMMIT=$FRAMEWORK_COMMIT,RUN_TAG=run" \
    hpc_jobs/active/hpc_pfp_benchmark.sh \
    --benchmark-id contemporary-2025-01-to-2026-02-supervisor \
    --benchmark-dir "$GLOBAL_BENCH_DIR" \
    --embedding-cache-archive "$GLOBAL_ARCHIVE" \
    --embedding-evidence "$GLOBAL_EVIDENCE/coverage.json" \
    --embedding-evidence "$GLOBAL_EVIDENCE/contract.json" \
    --embedding-evidence "$GLOBAL_EVIDENCE/targets.tsv" \
    --embedding-evidence "$GLOBAL_EVIDENCE/pair_status.tsv" \
    --require-embedding-evidence --obo-file "$TEMPORAL_OBO" \
    --results-root "$results_root" --config configs/pfp_benchmark_run.temporal.json \
    --execution-mode train-eval --modality-mode "$mode" \
    --aspect BPO --aspect CCO --aspect MFO --seed 42 --num-workers 0 --capture-predictions
  GLOBAL_MODEL_JOBS[$index]="$SUBMITTED_JOB_ID"
  GLOBAL_RUNS[$index]="$run_path"
done

NK_RESULTS_ROOT="$CAMPAIGN_ROOT/model_runs/nk_lk/full"
NK_RUN="$NK_RESULTS_ROOT/run"
submit_job model-nk-lk-full e12_nfull "$NK_BIND_JOB" "$NK_RUN" \
  -v "FRAMEWORK_COMMIT=$FRAMEWORK_COMMIT,RUN_TAG=run" \
  hpc_jobs/active/hpc_pfp_benchmark.sh \
  --benchmark-id contemporary-2025-01-to-2026-02-supervisor-nk-lk \
  --benchmark-dir "$NK_BENCH_DIR" --embedding-cache-archive "$NK_ARCHIVE" \
  --embedding-evidence "$NK_EVIDENCE/coverage.json" \
  --embedding-evidence "$NK_EVIDENCE/contract.json" \
  --embedding-evidence "$NK_EVIDENCE/targets.tsv" \
  --embedding-evidence "$NK_EVIDENCE/pair_status.tsv" \
  --require-embedding-evidence --obo-file "$TEMPORAL_OBO" \
  --results-root "$NK_RESULTS_ROOT" --config configs/pfp_benchmark_run.temporal_nk_lk.json \
  --execution-mode train-eval --modality-mode full \
  --aspect BPO --aspect CCO --aspect MFO --seed 42 --num-workers 0 --capture-predictions
NK_MODEL_JOB="$SUBMITTED_JOB_ID"

HOMOLOGY_MODEL_JOBS=()
HOMOLOGY_RUNS=()
for identity in "${identities[@]}"; do
  identity_dir="$identity"; [[ "$identity" != 5 ]] || identity_dir=05
  benchmark="${HOMOLOGY_BENCHMARKS[$identity]}" binding="${HOMOLOGY_BINDINGS[$identity]}"
  results_root="$CAMPAIGN_ROOT/model_runs/homology/identity_${identity_dir}/full"
  run_path="$results_root/run"
  submit_job "model-homology-$identity" "e12_h${identity}m" "${HOMOLOGY_BIND_JOBS[$identity]}" "$run_path" \
    -v "FRAMEWORK_COMMIT=$FRAMEWORK_COMMIT,RUN_TAG=run" \
    hpc_jobs/active/hpc_pfp_benchmark.sh \
    --benchmark-id "homology-uniref50-identity-${identity}-fixed-evidence" \
    --benchmark-dir "$benchmark" \
    --benchmark-evidence "$benchmark/validation_report.json" \
    --benchmark-evidence "$benchmark/output_manifest.json" \
    --benchmark-evidence "$benchmark/RUN_COMPLETE.json" \
    --embedding-cache-archive "$HOMOLOGY_ARCHIVE" \
    --embedding-evidence "$binding/coverage.json" --embedding-evidence "$binding/contract.json" \
    --embedding-evidence "$binding/targets.tsv" --embedding-evidence "$binding/pair_status.tsv" \
    --require-embedding-evidence --obo-file "$HOMOLOGY_OBO" \
    --results-root "$results_root" --config configs/pfp_benchmark_run.homology.json \
    --execution-mode train-eval --modality-mode full \
    --aspect BPO --aspect CCO --aspect MFO --seed 42 --num-workers 0 --capture-predictions
  HOMOLOGY_MODEL_JOBS[$identity]="$SUBMITTED_JOB_ID"
  HOMOLOGY_RUNS[$identity]="$run_path"
done

GLOBAL_MODEL_HOLDS="$(IFS=,; echo "${GLOBAL_MODEL_JOBS[*]}")"
GLOBAL_PANEL="$CAMPAIGN_ROOT/analyses/global_nk/modality_panel"
panel_args=(hpc_jobs/active/hpc_pfp_modality_panel_analysis.sh --obo-file "$TEMPORAL_OBO" --output-dir "$GLOBAL_PANEL")
for index in "${!global_modes[@]}"; do
  mode="${global_modes[$index]}"
  panel_args+=(--run "$mode=${GLOBAL_RUNS[$index]}" --prediction-run "$mode=${GLOBAL_RUNS[$index]}")
done
submit_job analysis-global-panel e12_gpanel "$GLOBAL_MODEL_HOLDS" "$GLOBAL_PANEL" \
  -v "FRAMEWORK_COMMIT=$FRAMEWORK_COMMIT" "${panel_args[@]}"

GLOBAL_FULL_JOB="${GLOBAL_MODEL_JOBS[0]}"
GLOBAL_FULL_RUN="${GLOBAL_RUNS[0]}"
GLOBAL_DIAGNOSTICS="$CAMPAIGN_ROOT/analyses/global_nk/diagnostics"
submit_job analysis-global-diagnostics e12_gdiag "$GLOBAL_FULL_JOB" "$GLOBAL_DIAGNOSTICS" \
  -v "FRAMEWORK_COMMIT=$FRAMEWORK_COMMIT" hpc_jobs/active/hpc_contemporary_followup_analysis.sh \
  --analysis diagnostics --source-run "$GLOBAL_FULL_RUN" --obo "$TEMPORAL_OBO" \
  --source-label global-nk-full-evidence12 --output-dir "$GLOBAL_DIAGNOSTICS"
GLOBAL_DIAGNOSTICS_JOB="$SUBMITTED_JOB_ID"

GLOBAL_CAPTURE="$CAMPAIGN_ROOT/analyses/global_nk/prediction_capture"
submit_job capture-global e12_gcap "$GLOBAL_FULL_JOB" "$GLOBAL_CAPTURE" \
  -v "FRAMEWORK_COMMIT=$FRAMEWORK_COMMIT" hpc_jobs/active/hpc_contemporary_followup_prediction_capture.sh \
  --source-run "$GLOBAL_FULL_RUN" --cache-archive "$GLOBAL_ARCHIVE" --obo "$TEMPORAL_OBO" \
  --benchmark-id contemporary-2025-01-to-2026-02-supervisor \
  --source-label global-nk-full-evidence12 --output-dir "$GLOBAL_CAPTURE"
GLOBAL_CAPTURE_JOB="$SUBMITTED_JOB_ID"
for analysis in calibration confidence; do
  output="$CAMPAIGN_ROOT/analyses/global_nk/$analysis"
  submit_job "analysis-global-$analysis" "e12_g${analysis:0:3}" "$GLOBAL_CAPTURE_JOB" "$output" \
    -v "FRAMEWORK_COMMIT=$FRAMEWORK_COMMIT" hpc_jobs/active/hpc_contemporary_followup_analysis.sh \
    --analysis "$analysis" --capture-pair-dir "$GLOBAL_CAPTURE" --source-label global-nk-full-evidence12 \
    --output-dir "$output"
done

GLOBAL_CENSUS="$CAMPAIGN_ROOT/analyses/global_nk/knowledge_census"
submit_job analysis-global-census e12_gcensus "$GLOBAL_BENCH_JOB" "$GLOBAL_CENSUS" \
  -v "FRAMEWORK_COMMIT=$FRAMEWORK_COMMIT,BENCHMARK_DIR=$GLOBAL_BENCH_DIR,BENCHMARK_MANIFEST=$GLOBAL_REPORT_DIR/build_manifest.json" \
  hpc_jobs/active/hpc_contemporary_knowledge_cohort_census.sh --output-dir "$GLOBAL_CENSUS"

GLOBAL_FORENSICS="$CAMPAIGN_ROOT/analyses/global_nk/forensics"
submit_job analysis-global-forensics e12_gfor "$GLOBAL_BIND_JOB" "$GLOBAL_FORENSICS" \
  -v "FRAMEWORK_COMMIT=$FRAMEWORK_COMMIT,CONTEMPORARY_BENCHMARK_DIR=$GLOBAL_BENCH_DIR,CONTEMPORARY_MODALITY_STATUS=$GLOBAL_EVIDENCE/pair_status.tsv" \
  hpc_jobs/active/hpc_benchmark_forensics.sh --dataset contemporary --output-dir "$GLOBAL_FORENSICS"

NK_DIAGNOSTICS="$CAMPAIGN_ROOT/analyses/nk_lk/diagnostics"
submit_job analysis-nk-lk-diagnostics e12_ndiag "$NK_MODEL_JOB" "$NK_DIAGNOSTICS" \
  -v "FRAMEWORK_COMMIT=$FRAMEWORK_COMMIT" hpc_jobs/active/hpc_contemporary_followup_analysis.sh \
  --analysis diagnostics --source-run "$NK_RUN" --obo "$TEMPORAL_OBO" \
  --source-label nk-lk-full-evidence12 --output-dir "$NK_DIAGNOSTICS"
NK_DIAGNOSTICS_JOB="$SUBMITTED_JOB_ID"

NK_CAPTURE="$CAMPAIGN_ROOT/analyses/nk_lk/prediction_capture"
submit_job capture-nk-lk e12_ncap "$NK_MODEL_JOB" "$NK_CAPTURE" \
  -v "FRAMEWORK_COMMIT=$FRAMEWORK_COMMIT" hpc_jobs/active/hpc_contemporary_followup_prediction_capture.sh \
  --source-run "$NK_RUN" --cache-archive "$NK_ARCHIVE" --obo "$TEMPORAL_OBO" \
  --benchmark-id contemporary-2025-01-to-2026-02-supervisor-nk-lk \
  --source-label nk-lk-full-evidence12 --output-dir "$NK_CAPTURE"
NK_CAPTURE_JOB="$SUBMITTED_JOB_ID"
for analysis in calibration confidence; do
  output="$CAMPAIGN_ROOT/analyses/nk_lk/$analysis"
  submit_job "analysis-nk-lk-$analysis" "e12_n${analysis:0:3}" "$NK_CAPTURE_JOB" "$output" \
    -v "FRAMEWORK_COMMIT=$FRAMEWORK_COMMIT" hpc_jobs/active/hpc_contemporary_followup_analysis.sh \
    --analysis "$analysis" --capture-pair-dir "$NK_CAPTURE" --source-label nk-lk-full-evidence12 \
    --output-dir "$output"
done

NK_FORENSICS="$CAMPAIGN_ROOT/analyses/nk_lk/forensics"
submit_job analysis-nk-lk-forensics e12_nfor "$NK_BIND_JOB" "$NK_FORENSICS" \
  -v "FRAMEWORK_COMMIT=$FRAMEWORK_COMMIT,CONTEMPORARY_NK_LK_BENCHMARK_DIR=$NK_BENCH_DIR,CONTEMPORARY_NK_LK_MODALITY_STATUS=$NK_EVIDENCE/pair_status.tsv" \
  hpc_jobs/active/hpc_benchmark_forensics.sh --dataset contemporary-nk-lk --output-dir "$NK_FORENSICS"

NK_COMPLETED="$CAMPAIGN_ROOT/analyses/nk_lk/cohort_comparison"
submit_job analysis-nk-lk-cohorts e12_ncohort "$GLOBAL_FULL_JOB,$NK_MODEL_JOB" "$NK_COMPLETED" \
  -v "FRAMEWORK_COMMIT=$FRAMEWORK_COMMIT" hpc_jobs/active/hpc_evaluate_nk_lk_completed_run.sh \
  --prediction-manifest "$NK_RUN/evaluation/prediction_artifacts/prediction_artifact_manifest.json" \
  --membership-tsv "$NK_REPORT_DIR/test_knowledge_cohort_membership.tsv" \
  --nk-evaluation-summary "$GLOBAL_FULL_RUN/evaluation/evaluation_summary.json" \
  --output-dir "$NK_COMPLETED" --bootstrap-replicates 2000 --bootstrap-seed 20260805

HOMOLOGY_DIAGNOSTIC_JOBS=()
HOMOLOGY_DIAGNOSTICS=()
for identity in "${identities[@]}"; do
  identity_dir="$identity"; [[ "$identity" != 5 ]] || identity_dir=05
  run_path="${HOMOLOGY_RUNS[$identity]}" model_job="${HOMOLOGY_MODEL_JOBS[$identity]}"
  base="$CAMPAIGN_ROOT/analyses/homology/identity_${identity_dir}"
  diagnostics="$base/diagnostics"
  submit_job "analysis-homology-${identity}-diagnostics" "e12_h${identity}d" "$model_job" "$diagnostics" \
    -v "FRAMEWORK_COMMIT=$FRAMEWORK_COMMIT" hpc_jobs/active/hpc_contemporary_followup_analysis.sh \
    --analysis diagnostics --source-run "$run_path" --obo "$HOMOLOGY_OBO" \
    --source-label "homology-identity-${identity}-full-evidence12" --output-dir "$diagnostics"
  HOMOLOGY_DIAGNOSTIC_JOBS[$identity]="$SUBMITTED_JOB_ID"
  HOMOLOGY_DIAGNOSTICS[$identity]="$diagnostics"

  capture="$base/prediction_capture"
  submit_job "capture-homology-$identity" "e12_h${identity}c" "$model_job" "$capture" \
    -v "FRAMEWORK_COMMIT=$FRAMEWORK_COMMIT" hpc_jobs/active/hpc_contemporary_followup_prediction_capture.sh \
    --source-run "$run_path" --cache-archive "$HOMOLOGY_ARCHIVE" --obo "$HOMOLOGY_OBO" \
    --benchmark-id "homology-uniref50-identity-${identity}-fixed-evidence" \
    --source-label "homology-identity-${identity}-full-evidence12" --output-dir "$capture"
  capture_job="$SUBMITTED_JOB_ID"
  for analysis in calibration confidence; do
    output="$base/$analysis"
    submit_job "analysis-homology-${identity}-$analysis" "e12_h${identity}${analysis:0:2}" "$capture_job" "$output" \
      -v "FRAMEWORK_COMMIT=$FRAMEWORK_COMMIT" hpc_jobs/active/hpc_contemporary_followup_analysis.sh \
      --analysis "$analysis" --capture-pair-dir "$capture" \
      --source-label "homology-identity-${identity}-full-evidence12" --output-dir "$output"
  done
done

ROOT_ENV="FRAMEWORK_COMMIT=$FRAMEWORK_COMMIT,ROOT_30=${HOMOLOGY_BENCHMARKS[30]},ROOT_25=${HOMOLOGY_BENCHMARKS[25]},ROOT_20=${HOMOLOGY_BENCHMARKS[20]},ROOT_15=${HOMOLOGY_BENCHMARKS[15]},ROOT_10=${HOMOLOGY_BENCHMARKS[10]},ROOT_05=${HOMOLOGY_BENCHMARKS[5]}"
HOMOLOGY_THRESHOLD_AUDIT="$CAMPAIGN_ROOT/analyses/homology/threshold_progression"
submit_job analysis-homology-thresholds e12_hthr "$HOMOLOGY_BENCH_JOB" "$HOMOLOGY_THRESHOLD_AUDIT" \
  -v "$ROOT_ENV,FULL_30=$FULL_30,FULL_25=$FULL_25,FULL_20=$FULL_20,FULL_15=$FULL_15,FULL_10=$FULL_10,FULL_05=$FULL_05,RESULTS_ROOT=$HOMOLOGY_THRESHOLD_AUDIT" \
  hpc_jobs/active/hpc_homology_threshold_progression_audit.sh

HOMOLOGY_EVIDENCE_AUDIT="$CAMPAIGN_ROOT/analyses/homology/evidence_policy"
submit_job analysis-homology-evidence e12_hev "$HOMOLOGY_BENCH_JOB" "$HOMOLOGY_EVIDENCE_AUDIT" \
  -v "$ROOT_ENV,INCLUDE_DANIEL_COMPARATOR=0,RESULTS_ROOT=$HOMOLOGY_EVIDENCE_AUDIT" \
  hpc_jobs/active/hpc_homology_evidence_policy_audit.sh

HOMOLOGY_DIAGNOSTIC_HOLDS="$(IFS=,; echo "${HOMOLOGY_DIAGNOSTIC_JOBS[*]}")"
HOMOLOGY_SPECIFICITY="$CAMPAIGN_ROOT/analyses/homology/specificity_comparison"
homology_spec_args=(hpc_jobs/active/hpc_specificity_comparison.sh)
for identity in "${identities[@]}"; do
  homology_spec_args+=(--source "identity-${identity}=${HOMOLOGY_DIAGNOSTICS[$identity]}/analysis/specificity")
done
homology_spec_args+=(--output-dir "$HOMOLOGY_SPECIFICITY")
submit_job analysis-homology-specificity e12_hspec "$HOMOLOGY_DIAGNOSTIC_HOLDS" "$HOMOLOGY_SPECIFICITY" \
  -v "FRAMEWORK_COMMIT=$FRAMEWORK_COMMIT" "${homology_spec_args[@]}"

CROSS_SPECIFICITY="$CAMPAIGN_ROOT/analyses/cross_benchmark_specificity"
submit_job analysis-cross-specificity e12_xspec "$GLOBAL_DIAGNOSTICS_JOB,$NK_DIAGNOSTICS_JOB,${HOMOLOGY_DIAGNOSTIC_JOBS[30]}" "$CROSS_SPECIFICITY" \
  -v "FRAMEWORK_COMMIT=$FRAMEWORK_COMMIT" hpc_jobs/active/hpc_specificity_comparison.sh \
  --source "global-nk=$GLOBAL_DIAGNOSTICS/analysis/specificity" \
  --source "nk-lk=$NK_DIAGNOSTICS/analysis/specificity" \
  --source "homology-30=${HOMOLOGY_DIAGNOSTICS[30]}/analysis/specificity" \
  --output-dir "$CROSS_SPECIFICITY"

printf 'submitted_at_utc=%s\nledger=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$LEDGER" \
  >> "$CAMPAIGN_ROOT/submission/campaign.env"
echo "Campaign root: $CAMPAIGN_ROOT"
echo "Submission ledger: $LEDGER"
