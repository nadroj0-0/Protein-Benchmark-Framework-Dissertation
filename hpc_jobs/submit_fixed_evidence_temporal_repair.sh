#!/usr/bin/env bash
# Submit the corrected temporal embedding repair, model runs, and analyses.

set -Eeuo pipefail

die() { echo "ERROR: $*" >&2; exit 2; }

FRAMEWORK_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
cd "$FRAMEWORK_ROOT"

: "${CAMPAIGN_ROOT:?Pass the existing corrected-evidence CAMPAIGN_ROOT}"
: "${GLOBAL_LEDGER:?Pass the completed global-NK primary ledger directory}"
: "${NK_LEDGER:?Pass the completed NK+LK primary ledger directory}"

SAN_ROOT="${SAN_ROOT:-/SAN/bioinf/bmpfp}"
ARTIFACT_CATALOG="${ARTIFACT_CATALOG:-$SAN_ROOT/manifests/artifact_paths.tsv}"
TEMPORAL_OBO="${TEMPORAL_OBO:-$SAN_ROOT/frozen_inputs/ontology/2025-02-06/go-basic.obo}"
EMBEDDING_POLICY="configs/contemporary_embedding_generation.json"
TEXT_CUTOFF_DATE="2025-03-08"
CAMPAIGN_TAG="$(basename "$CAMPAIGN_ROOT")"
LOG_ROOT="${LOG_ROOT:-$HOME/pfp_fixed_evidence_logs/${CAMPAIGN_TAG}-temporal-repair}"
REPAIR_ROOT="$CAMPAIGN_ROOT/embedding_repair/generated"
GLOBAL_BATCH="$REPAIR_ROOT/global_nk"
NK_BATCH="$REPAIR_ROOT/nk_lk"
GLOBAL_BENCH_DIR="$CAMPAIGN_ROOT/benchmarks/global_nk/benchmark/outputs"
NK_BENCH_DIR="$CAMPAIGN_ROOT/benchmarks/nk_lk/benchmark/outputs"
GLOBAL_REPORT_DIR="$CAMPAIGN_ROOT/benchmarks/global_nk/benchmark/reports"
NK_REPORT_DIR="$CAMPAIGN_ROOT/benchmarks/nk_lk/benchmark/reports"
GLOBAL_ARCHIVE="$GLOBAL_BATCH/final/homology_30_embedding_cache.tar.gz"
NK_ARCHIVE="$NK_BATCH/final/homology_30_embedding_cache.tar.gz"
EVIDENCE_ROOT="$CAMPAIGN_ROOT/embedding_evidence_repaired"
GLOBAL_EVIDENCE="$EVIDENCE_ROOT/global_nk"
NK_EVIDENCE="$EVIDENCE_ROOT/nk_lk"
MODEL_ROOT="$CAMPAIGN_ROOT/model_runs_repaired"
ANALYSIS_ROOT="$CAMPAIGN_ROOT/analyses_repaired"

[[ "$(git symbolic-ref --quiet --short HEAD || true)" == fixed-evidence-codes ]] || \
  die "Submit from fixed-evidence-codes"
[[ -z "$(git status --porcelain)" ]] || die "Framework checkout must be clean"
FRAMEWORK_COMMIT="$(git rev-parse HEAD)"
REMOTE_COMMIT="$(git ls-remote origin refs/heads/fixed-evidence-codes | awk 'NR==1 {print $1}')"
[[ "$REMOTE_COMMIT" == "$FRAMEWORK_COMMIT" ]] || die "Remote branch does not match HEAD"

for path in "$ARTIFACT_CATALOG" "$TEMPORAL_OBO" \
  "$GLOBAL_LEDGER/output_manifest.json" "$GLOBAL_LEDGER/RUN_COMPLETE.json" \
  "$NK_LEDGER/output_manifest.json" "$NK_LEDGER/RUN_COMPLETE.json"; do
  [[ -f "$path" && ! -L "$path" ]] || die "Required input is missing or unsafe: $path"
done
for path in "$GLOBAL_BENCH_DIR" "$NK_BENCH_DIR"; do
  [[ -d "$path" && ! -L "$path" ]] || die "Benchmark directory is missing: $path"
done
for path in "$GLOBAL_BATCH" "$NK_BATCH" "$EVIDENCE_ROOT" "$MODEL_ROOT" "$ANALYSIS_ROOT"; do
  [[ ! -e "$path" ]] || die "Repair destination already exists: $path"
done
command -v qsub >/dev/null 2>&1 || die "qsub is unavailable"
mkdir -p "$REPAIR_ROOT" "$LOG_ROOT" "$CAMPAIGN_ROOT/submission"

LEDGER="$CAMPAIGN_ROOT/submission/temporal_repair_jobs.tsv"
printf 'stage\tjob_name\tjob_id\thold_jid\texpected_output\tcommand\n' > "$LEDGER"
SUBMITTED_JOB_ID=""
submit_job() {
  local stage="$1" name="$2" holds="$3" expected_output="$4"
  shift 4
  local command=(qsub -N "$name" -o "$LOG_ROOT") output job_id command_text
  [[ -z "$holds" ]] || command+=(-hold_jid "$holds")
  command+=("$@")
  printf -v command_text '%q ' "${command[@]}"
  output="$("${command[@]}")"
  job_id="$(printf '%s\n' "$output" | sed -n \
    -e 's/.*job-array \([0-9][0-9]*\).*/\1/p' \
    -e 's/.*job \([0-9][0-9]*\).*/\1/p' | head -n 1)"
  [[ "$job_id" =~ ^[0-9]+$ ]] || die "Could not parse qsub output for $name: $output"
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$stage" "$name" "$job_id" "$holds" "$expected_output" "$command_text" >> "$LEDGER"
  echo "Submitted $name as $job_id${holds:+, holding on $holds}" >&2
  SUBMITTED_JOB_ID="$job_id"
}

submit_lane_embeddings() {
  local lane="$1" benchmark="$2" source_ledger="$3" batch="$4"
  submit_job "embedding-$lane-gpu" "e12_${lane}eg" "" "$batch/deltas" \
    -v "FRAMEWORK_JOB_ROOT=$FRAMEWORK_ROOT,FRAMEWORK_COMMIT=$FRAMEWORK_COMMIT,EMBEDDING_POLICY=$EMBEDDING_POLICY" \
    hpc_jobs/active/hpc_homology_embedding_gpu_array.sh \
    --benchmark-dir "$benchmark" --ledger-dir "$source_ledger" --batch-root "$batch" \
    --text-cutoff-date "$TEXT_CUTOFF_DATE" --artifact-catalog "$ARTIFACT_CATALOG"
  local gpu_job="$SUBMITTED_JOB_ID"

  submit_job "embedding-$lane-ppi" "e12_${lane}ep" "" "$batch/deltas/ppi" \
    -v "FRAMEWORK_JOB_ROOT=$FRAMEWORK_ROOT,FRAMEWORK_COMMIT=$FRAMEWORK_COMMIT,EMBEDDING_POLICY=$EMBEDDING_POLICY" \
    hpc_jobs/active/hpc_homology_embedding_ppi.sh \
    --benchmark-dir "$benchmark" --ledger-dir "$source_ledger" --batch-root "$batch" \
    --artifact-catalog "$ARTIFACT_CATALOG"
  local ppi_job="$SUBMITTED_JOB_ID"

  submit_job "embedding-$lane-finalize" "e12_${lane}ef" "$gpu_job,$ppi_job" "$batch/final" \
    -v "FRAMEWORK_COMMIT=$FRAMEWORK_COMMIT" \
    hpc_jobs/active/hpc_homology_embedding_finalize.sh \
    --ledger-dir "$source_ledger" --batch-root "$batch" --policy "$EMBEDDING_POLICY"
}

submit_lane_embeddings g "$GLOBAL_BENCH_DIR" "$GLOBAL_LEDGER" "$GLOBAL_BATCH"
GLOBAL_FINAL_JOB="$SUBMITTED_JOB_ID"
submit_lane_embeddings n "$NK_BENCH_DIR" "$NK_LEDGER" "$NK_BATCH"
NK_FINAL_JOB="$SUBMITTED_JOB_ID"

submit_job bind-global e12_rgbind "$GLOBAL_FINAL_JOB" "$GLOBAL_EVIDENCE" \
  -v "BENCHMARK_DIR=$GLOBAL_BENCH_DIR,BENCHMARK_ID=contemporary-2025-01-to-2026-02-supervisor,EMBEDDING_ARCHIVE=$GLOBAL_ARCHIVE,RUN_CONFIG=configs/pfp_benchmark_run.temporal.json,EVIDENCE_OUTPUT=$GLOBAL_EVIDENCE,FRAMEWORK_COMMIT=$FRAMEWORK_COMMIT" \
  hpc_jobs/active/hpc_bind_embedding_archive_evidence.sh
GLOBAL_BIND_JOB="$SUBMITTED_JOB_ID"

submit_job bind-nk-lk e12_rnbind "$NK_FINAL_JOB" "$NK_EVIDENCE" \
  -v "BENCHMARK_DIR=$NK_BENCH_DIR,BENCHMARK_ID=contemporary-2025-01-to-2026-02-supervisor-nk-lk,EMBEDDING_ARCHIVE=$NK_ARCHIVE,RUN_CONFIG=configs/pfp_benchmark_run.temporal_nk_lk.json,EVIDENCE_OUTPUT=$NK_EVIDENCE,FRAMEWORK_COMMIT=$FRAMEWORK_COMMIT" \
  hpc_jobs/active/hpc_bind_embedding_archive_evidence.sh
NK_BIND_JOB="$SUBMITTED_JOB_ID"

GLOBAL_MODEL_JOBS=()
GLOBAL_RUNS=()
GLOBAL_MODES=(full sequence-only sequence-text sequence-structure sequence-ppi)
for index in "${!GLOBAL_MODES[@]}"; do
  mode="${GLOBAL_MODES[$index]}"
  case "$mode" in
    full) name=e12_rgfull ;;
    sequence-only) name=e12_rgseq ;;
    sequence-text) name=e12_rgtxt ;;
    sequence-structure) name=e12_rgstr ;;
    sequence-ppi) name=e12_rgppi ;;
  esac
  results_root="$MODEL_ROOT/global_nk/$mode"
  submit_job "model-global-$mode" "$name" "$GLOBAL_BIND_JOB" "$results_root/run" \
    -v "FRAMEWORK_COMMIT=$FRAMEWORK_COMMIT,RUN_TAG=run" \
    hpc_jobs/active/hpc_pfp_benchmark.sh \
    --benchmark-id contemporary-2025-01-to-2026-02-supervisor \
    --benchmark-dir "$GLOBAL_BENCH_DIR" --embedding-cache-archive "$GLOBAL_ARCHIVE" \
    --embedding-evidence "$GLOBAL_EVIDENCE/coverage.json" \
    --embedding-evidence "$GLOBAL_EVIDENCE/contract.json" \
    --embedding-evidence "$GLOBAL_EVIDENCE/targets.tsv" \
    --embedding-evidence "$GLOBAL_EVIDENCE/pair_status.tsv" \
    --require-embedding-evidence --obo-file "$TEMPORAL_OBO" \
    --results-root "$results_root" --config configs/pfp_benchmark_run.temporal.json \
    --execution-mode train-eval --modality-mode "$mode" \
    --aspect BPO --aspect CCO --aspect MFO --seed 42 --num-workers 0 --capture-predictions
  GLOBAL_MODEL_JOBS[$index]="$SUBMITTED_JOB_ID"
  GLOBAL_RUNS[$index]="$results_root/run"
done

NK_RESULTS_ROOT="$MODEL_ROOT/nk_lk/full"
NK_RUN="$NK_RESULTS_ROOT/run"
submit_job model-nk-lk-full e12_rnfull "$NK_BIND_JOB" "$NK_RUN" \
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

GLOBAL_MODEL_HOLDS="$(IFS=,; echo "${GLOBAL_MODEL_JOBS[*]}")"
GLOBAL_PANEL="$ANALYSIS_ROOT/global_nk/modality_panel"
panel_args=(hpc_jobs/active/hpc_pfp_modality_panel_analysis.sh --obo-file "$TEMPORAL_OBO" --output-dir "$GLOBAL_PANEL")
for index in "${!GLOBAL_MODES[@]}"; do
  mode="${GLOBAL_MODES[$index]}"
  panel_args+=(--run "$mode=${GLOBAL_RUNS[$index]}" --prediction-run "$mode=${GLOBAL_RUNS[$index]}")
done
submit_job analysis-global-panel e12_rgpanel "$GLOBAL_MODEL_HOLDS" "$GLOBAL_PANEL" \
  -v "FRAMEWORK_COMMIT=$FRAMEWORK_COMMIT" "${panel_args[@]}"

GLOBAL_FULL_JOB="${GLOBAL_MODEL_JOBS[0]}"
GLOBAL_FULL_RUN="${GLOBAL_RUNS[0]}"
GLOBAL_DIAGNOSTICS="$ANALYSIS_ROOT/global_nk/diagnostics"
submit_job analysis-global-diagnostics e12_rgdiag "$GLOBAL_FULL_JOB" "$GLOBAL_DIAGNOSTICS" \
  -v "FRAMEWORK_COMMIT=$FRAMEWORK_COMMIT" hpc_jobs/active/hpc_contemporary_followup_analysis.sh \
  --analysis diagnostics --source-run "$GLOBAL_FULL_RUN" --obo "$TEMPORAL_OBO" \
  --source-label global-nk-full-evidence12-repaired --output-dir "$GLOBAL_DIAGNOSTICS"
GLOBAL_DIAGNOSTICS_JOB="$SUBMITTED_JOB_ID"

GLOBAL_CAPTURE="$ANALYSIS_ROOT/global_nk/prediction_capture"
submit_job capture-global e12_rgcap "$GLOBAL_FULL_JOB" "$GLOBAL_CAPTURE" \
  -v "FRAMEWORK_COMMIT=$FRAMEWORK_COMMIT" hpc_jobs/active/hpc_contemporary_followup_prediction_capture.sh \
  --source-run "$GLOBAL_FULL_RUN" --cache-archive "$GLOBAL_ARCHIVE" --obo "$TEMPORAL_OBO" \
  --benchmark-id contemporary-2025-01-to-2026-02-supervisor \
  --source-label global-nk-full-evidence12-repaired --output-dir "$GLOBAL_CAPTURE"
GLOBAL_CAPTURE_JOB="$SUBMITTED_JOB_ID"
for analysis in calibration confidence; do
  output="$ANALYSIS_ROOT/global_nk/$analysis"
  submit_job "analysis-global-$analysis" "e12_rg${analysis:0:3}" "$GLOBAL_CAPTURE_JOB" "$output" \
    -v "FRAMEWORK_COMMIT=$FRAMEWORK_COMMIT" hpc_jobs/active/hpc_contemporary_followup_analysis.sh \
    --analysis "$analysis" --capture-pair-dir "$GLOBAL_CAPTURE" \
    --source-label global-nk-full-evidence12-repaired --output-dir "$output"
done

GLOBAL_FORENSICS="$ANALYSIS_ROOT/global_nk/forensics"
submit_job analysis-global-forensics e12_rgfor "$GLOBAL_BIND_JOB" "$GLOBAL_FORENSICS" \
  -v "FRAMEWORK_COMMIT=$FRAMEWORK_COMMIT,CONTEMPORARY_BENCHMARK_DIR=$GLOBAL_BENCH_DIR,CONTEMPORARY_MODALITY_STATUS=$GLOBAL_EVIDENCE/pair_status.tsv" \
  hpc_jobs/active/hpc_benchmark_forensics.sh --dataset contemporary --output-dir "$GLOBAL_FORENSICS"

NK_DIAGNOSTICS="$ANALYSIS_ROOT/nk_lk/diagnostics"
submit_job analysis-nk-lk-diagnostics e12_rndiag "$NK_MODEL_JOB" "$NK_DIAGNOSTICS" \
  -v "FRAMEWORK_COMMIT=$FRAMEWORK_COMMIT" hpc_jobs/active/hpc_contemporary_followup_analysis.sh \
  --analysis diagnostics --source-run "$NK_RUN" --obo "$TEMPORAL_OBO" \
  --source-label nk-lk-full-evidence12-repaired --output-dir "$NK_DIAGNOSTICS"
NK_DIAGNOSTICS_JOB="$SUBMITTED_JOB_ID"

NK_CAPTURE="$ANALYSIS_ROOT/nk_lk/prediction_capture"
submit_job capture-nk-lk e12_rncap "$NK_MODEL_JOB" "$NK_CAPTURE" \
  -v "FRAMEWORK_COMMIT=$FRAMEWORK_COMMIT" hpc_jobs/active/hpc_contemporary_followup_prediction_capture.sh \
  --source-run "$NK_RUN" --cache-archive "$NK_ARCHIVE" --obo "$TEMPORAL_OBO" \
  --benchmark-id contemporary-2025-01-to-2026-02-supervisor-nk-lk \
  --source-label nk-lk-full-evidence12-repaired --output-dir "$NK_CAPTURE"
NK_CAPTURE_JOB="$SUBMITTED_JOB_ID"
for analysis in calibration confidence; do
  output="$ANALYSIS_ROOT/nk_lk/$analysis"
  submit_job "analysis-nk-lk-$analysis" "e12_rn${analysis:0:3}" "$NK_CAPTURE_JOB" "$output" \
    -v "FRAMEWORK_COMMIT=$FRAMEWORK_COMMIT" hpc_jobs/active/hpc_contemporary_followup_analysis.sh \
    --analysis "$analysis" --capture-pair-dir "$NK_CAPTURE" \
    --source-label nk-lk-full-evidence12-repaired --output-dir "$output"
done

NK_FORENSICS="$ANALYSIS_ROOT/nk_lk/forensics"
submit_job analysis-nk-lk-forensics e12_rnfor "$NK_BIND_JOB" "$NK_FORENSICS" \
  -v "FRAMEWORK_COMMIT=$FRAMEWORK_COMMIT,CONTEMPORARY_NK_LK_BENCHMARK_DIR=$NK_BENCH_DIR,CONTEMPORARY_NK_LK_MODALITY_STATUS=$NK_EVIDENCE/pair_status.tsv" \
  hpc_jobs/active/hpc_benchmark_forensics.sh --dataset contemporary-nk-lk --output-dir "$NK_FORENSICS"

NK_COHORTS="$ANALYSIS_ROOT/nk_lk/cohort_comparison"
submit_job analysis-nk-lk-cohorts e12_rncoh "$GLOBAL_FULL_JOB,$NK_MODEL_JOB" "$NK_COHORTS" \
  -v "FRAMEWORK_COMMIT=$FRAMEWORK_COMMIT" hpc_jobs/active/hpc_evaluate_nk_lk_completed_run.sh \
  --prediction-manifest "$NK_RUN/evaluation/prediction_artifacts/prediction_artifact_manifest.json" \
  --membership-tsv "$NK_REPORT_DIR/test_knowledge_cohort_membership.tsv" \
  --nk-evaluation-summary "$GLOBAL_FULL_RUN/evaluation/evaluation_summary.json" \
  --output-dir "$NK_COHORTS" --bootstrap-replicates 2000 --bootstrap-seed 20260805

HOMOLOGY_30_DIAGNOSTICS_JOB="${HOMOLOGY_30_DIAGNOSTICS_JOB:-7228746}"
HOMOLOGY_30_DIAGNOSTICS="$CAMPAIGN_ROOT/analyses/homology/identity_30/diagnostics/analysis/specificity"
CROSS_SPECIFICITY="$ANALYSIS_ROOT/cross_benchmark_specificity"
submit_job analysis-cross-specificity e12_rxspec \
  "$GLOBAL_DIAGNOSTICS_JOB,$NK_DIAGNOSTICS_JOB,$HOMOLOGY_30_DIAGNOSTICS_JOB" "$CROSS_SPECIFICITY" \
  -v "FRAMEWORK_COMMIT=$FRAMEWORK_COMMIT" hpc_jobs/active/hpc_specificity_comparison.sh \
  --source "global-nk=$GLOBAL_DIAGNOSTICS/analysis/specificity" \
  --source "nk-lk=$NK_DIAGNOSTICS/analysis/specificity" \
  --source "homology-30=$HOMOLOGY_30_DIAGNOSTICS" --output-dir "$CROSS_SPECIFICITY"

printf 'submitted_at_utc=%s\nframework_commit=%s\npolicy=%s\nledger=%s\n' \
  "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$FRAMEWORK_COMMIT" "$EMBEDDING_POLICY" "$LEDGER" \
  > "$CAMPAIGN_ROOT/submission/temporal_repair.env"
echo "Temporal repair DAG submitted. Ledger: $LEDGER"
