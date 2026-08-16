# Benchmark-Agnostic PFP Execution

`run_pfp_benchmark.sh` is the submitted preparation, training, and evaluation
entry point for CAFA3, contemporary global-NK, contemporary NK+LK, and homology
benchmarks. It does not download inputs or modify the external PFP checkout.

## External Model Contract

PFP is pinned to:

```text
1e04fd6d6d3c40458fd41ec1a881ed6e24de768e
```

The runner requires that checkout's tracked source to be clean and rejects
untracked importable/executable code. Framework Git state is optional
report-only provenance and never blocks execution.

## Operations

```text
prepare-only  validate and stage a nine-CSV benchmark
eval-only     replay supplied checkpoints
train-eval    train and strictly evaluate fresh checkpoints
```

Fresh training invokes upstream `train.py --single` separately for BPO, CCO,
and MFO. Evaluation reloads the upstream dataset/model code and fixes:

```text
norm=cafa
prop=max
no_orphans=false
```

Supported modality modes are `full`, `sequence-only`, `sequence-text`,
`sequence-structure`, and `sequence-ppi`.

## Run

```bash
bash scripts/model_execution/run_pfp_benchmark.sh \
  --benchmark-id IDENTIFIER \
  --benchmark-dir /path/to/nine_csvs \
  --obo-file /path/to/go.obo \
  --pfp-root /path/to/pinned/PFP \
  --work-dir /path/to/new_work_dir \
  --output-dir /path/to/new_output_dir \
  --config configs/pfp_benchmark_run.temporal.json \
  --execution-mode train-eval \
  --embedding-cache-root /path/to/embedding_cache \
  --modality-mode full \
  --seed 42 \
  --num-workers 0 \
  --benchmark-evidence /path/to/build_manifest.json \
  --benchmark-evidence /path/to/output_checksums.sha256 \
  --require-embedding-evidence \
  --embedding-evidence /path/to/coverage.json \
  --embedding-evidence /path/to/contract.json \
  --embedding-evidence /path/to/targets.tsv \
  --embedding-evidence /path/to/pair_status.tsv \
  --capture-predictions
```

Use the matching submitted config:

```text
configs/pfp_benchmark_run.cafa3.json
configs/pfp_benchmark_run.temporal.json
configs/pfp_benchmark_run.temporal_nk_lk.json
configs/pfp_benchmark_run.homology.json
```

The global-NK temporal config requires bound embedding evidence with cache role
`accepted-corrected-2025-03-08`. The intermediate composition base cannot pass
model validation, even if `--require-embedding-evidence` is omitted.

Run `--help` for checkpoint, reference-artifact, tolerance, aspect, and evidence
options.

## Prediction Capture

`--capture-predictions` writes authenticated per-protein probabilities and a
prediction manifest. Select the split explicitly:

```text
--evaluation-split valid
--evaluation-split test
```

Calibration is fitted from a completed validation capture and assessed on a
separate completed test capture from the same benchmark/configuration,
ontology, PFP revision, checkpoint set, and modality mode. Use separate new
output directories because completed publications are immutable.

## Validation

Before training, the runner validates:

- nine CSV schemas and split/sequence overlap policy;
- benchmark identity and build evidence;
- OBO and training term consistency;
- embedding dimensions, coverage, membership, and evidence hashes;
- requested modality availability;
- checkpoints and reference artifacts for replay;
- PFP source and exact revision.

Every run publishes preparation, validation, execution, evaluation, and
provenance reports plus checksums and `RUN_COMPLETE.json`. Output is staged and
atomically published; an existing destination is not overwritten.

Fmax values are observations for the supplied seed and benchmark. The runner
does not convert single-seed differences into significance claims.
