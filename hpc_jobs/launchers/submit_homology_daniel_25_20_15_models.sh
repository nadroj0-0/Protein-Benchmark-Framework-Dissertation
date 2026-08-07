#!/usr/bin/env bash
# Bind the accepted hydrated cache to Daniel's 25/20/15 CSVs, then train.

set -Eeuo pipefail

die() { echo "ERROR: $*" >&2; exit 2; }
git_in_dir() { local directory="$1"; shift; (cd "$directory" && git "$@"); }

FRAMEWORK_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$FRAMEWORK_ROOT"
[[ -z "$(git status --porcelain)" ]] || die "Framework checkout must be clean"
FRAMEWORK_COMMIT="$(git rev-parse HEAD)"
SHORT_COMMIT="${FRAMEWORK_COMMIT:0:12}"
RUN_TAG="${SHORT_COMMIT}_$(date -u +%Y%m%dT%H%M%SZ)"

BENCHMARK_JOB_ROOT="/SAN/bioinf/bmpfp/benchmarks/homology/2026_02/supervisor_daniel_buchan/uniref50_sensitivity_4/runtime_array/source_sprot-and-trembl/uniref50_sensitivity_4/framework_686e21246b66/run_runtime-7142238/job_7142238"
BENCHMARK_SUFFIX="benchmark/source_sprot-and-trembl/framework_686e21246b66/uniref50_sensitivity_4/mmseqs_daniel-aligned-defaults"
CACHE_ARCHIVE="/SAN/bioinf/bmpfp/embeddings/homology/2026_02/identity_30/cluster-count-random/pair-resolved-paper-faithful/current-text-ledger_7132992_20260804T051124Z/final/homology_30_embedding_cache.tar.gz"
CACHE_SHA256="a1cb0cf0fc2e0142a039a146bc86408090632c5eaed001db0f90f235644a188f"
OBO_FILE="/SAN/bioinf/bmpfp/frozen_inputs/ontology/2026-06-19/go-basic.obo"
BINDING_ROOT="/SAN/bioinf/bmpfp/embedding_states/homology/2026_02/benchmark_bindings/daniel_stream_cluster_count_random_current_text/$RUN_TAG"
MODEL_ROOT="/SAN/bioinf/bmpfp/model_runs/homology/2026_02"
LOG_ROOT="${LOG_ROOT:-$HOME}"
[[ -d "$LOG_ROOT" ]] || die "Scheduler log directory is missing: $LOG_ROOT"

job_id_from_qsub() {
  local output="$1" job_id
  job_id="$(printf '%s\n' "$output" | sed -n 's/.*job \([0-9][0-9]*\).*/\1/p' | head -n 1)"
  [[ "$job_id" =~ ^[0-9]+$ ]] || die "Could not parse qsub job ID from: $output"
  printf '%s\n' "$job_id"
}

mkdir -p "$BINDING_ROOT"
printf 'identity\tbind_job\tmodel_job\tbenchmark_dir\tbinding_dir\tmodel_results_root\n' > "$BINDING_ROOT/submission_ledger.tsv"

for identity in 25 20 15; do
  case "$identity" in
    25) task=2 ;;
    20) task=3 ;;
    15) task=4 ;;
  esac
  benchmark_dir="$BENCHMARK_JOB_ROOT/task_${task}_identity_${identity}/$BENCHMARK_SUFFIX/identity_${identity}/cluster-count-random/annotated-only/seed_0/min_count_50"
  benchmark_id="homology-daniel-identity-${identity}"
  binding_dir="$BINDING_ROOT/identity_${identity}"
  results_root="$MODEL_ROOT/identity_${identity}/daniel_cluster_count_random/full_corrected_current_text"
  [[ -d "$benchmark_dir" ]] || die "Benchmark is missing: $benchmark_dir"
  for evidence in validation_report.json output_manifest.json RUN_COMPLETE.json; do
    [[ -f "$benchmark_dir/$evidence" ]] || die "Benchmark evidence is missing: $benchmark_dir/$evidence"
  done

  bind_output="$(qsub -N "hd${identity}_bind" -o "$LOG_ROOT" \
    -v "BENCHMARK_DIR=$benchmark_dir,BENCHMARK_ID=$benchmark_id,EMBEDDING_ARCHIVE=$CACHE_ARCHIVE,EMBEDDING_ARCHIVE_SHA256=$CACHE_SHA256,RUN_CONFIG=configs/pfp_benchmark_run.homology.json,EVIDENCE_OUTPUT=$binding_dir,FRAMEWORK_COMMIT=$FRAMEWORK_COMMIT" \
    hpc_jobs/active/hpc_bind_embedding_archive_evidence.sh)"
  bind_job="$(job_id_from_qsub "$bind_output")"
  echo "$bind_output"

  model_output="$(qsub -N "hd${identity}_full" -o "$LOG_ROOT" -hold_jid "$bind_job" \
    -v "FRAMEWORK_COMMIT=$FRAMEWORK_COMMIT" \
    hpc_jobs/active/hpc_pfp_benchmark.sh \
    --benchmark-id "$benchmark_id" \
    --benchmark-dir "$benchmark_dir" \
    --benchmark-evidence "$benchmark_dir/validation_report.json" \
    --benchmark-evidence "$benchmark_dir/output_manifest.json" \
    --benchmark-evidence "$benchmark_dir/RUN_COMPLETE.json" \
    --embedding-cache-archive "$CACHE_ARCHIVE" \
    --embedding-evidence "$binding_dir/coverage.json" \
    --embedding-evidence "$binding_dir/contract.json" \
    --embedding-evidence "$binding_dir/targets.tsv" \
    --embedding-evidence "$binding_dir/pair_status.tsv" \
    --require-embedding-evidence \
    --obo-file "$OBO_FILE" \
    --results-root "$results_root" \
    --config configs/pfp_benchmark_run.homology.json \
    --execution-mode train-eval \
    --modality-mode full \
    --capture-predictions)"
  model_job="$(job_id_from_qsub "$model_output")"
  echo "$model_output"
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$identity" "$bind_job" "$model_job" "$benchmark_dir" "$binding_dir" "$results_root" \
    >> "$BINDING_ROOT/submission_ledger.tsv"
done

echo "Submission ledger: $BINDING_ROOT/submission_ledger.tsv"
