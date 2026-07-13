# Protein Benchmark Framework

A reproducibility and benchmarking wrapper around the PFP/MMFP
protein-function-prediction pipeline. The repository is intended to make
the original CAFA3-style workflow easier to rerun, audit, and eventually
swap onto cleaner or decontaminated benchmarks.

The framework currently supports three levels of reproduction:

1.  **Evaluation only** using downloaded author artefacts.
2.  **Retrain and evaluate** using downloaded data/embeddings.
3.  **Regenerate embeddings, retrain, and evaluate** from the prepared
    benchmark pipeline.

The long-term design goal is to separate benchmark-specific values from
reusable logic. Verification scripts enforce data and embedding
contracts before expensive training is launched.

## Repository layout

``` text
.
├── benchmark_builders/
│   └── contemporary_cafa/                 # 2025→2026 CAFA-style benchmark builder
├── embedding_inventory/                    # CSV-native embedding inventory/reuse planner
├── scripts/
│   ├── reproduction/                      # Main PFP reproduction entrypoints
│   ├── embeddings/                        # Embedding wrappers and FASTA builder
│   ├── verification/                      # Shell and Python verification gates
│   ├── data_acquisition/                  # HPC/raw database download and inspection helpers
│   ├── diagnostics/                       # Environment probes and comparison diagnostics
│   ├── benchmark_generation/              # Reusable contemporary benchmark runner
│   └── validation/                        # Historical benchmark validation workflows
├── hpc_jobs/
│   ├── active/                            # qsub wrappers for current cluster workflows
│   ├── examples/                          # Scheduler examples/templates
│   └── archive/                           # Historical reproduction attempts
├── configs/
│   ├── cafa3.json                         # Default CAFA3 verification config
│   ├── embedding_inventory.cafa3_published.json
│   ├── embedding_inventory.future.example.json
│   └── paths.example.sh                   # Example local/HPC path configuration
```

Generated data, model checkpoints, cloned upstream repositories, and
embedding caches are intentionally not committed. See `.gitignore` for
the excluded paths.

## Embedding inventory and reuse planning

`embedding_inventory/` provides a benchmark-agnostic planner that reads the
nine PFP-compatible CSVs directly, validates existing arrays, applies
modality-specific scientific reuse rules, and emits reuse/generation/masking
manifests without copying or generating embeddings.

It supports both `paper-faithful` and `maximize-coverage` action policies,
requires an explicit source benchmark and explicit text directory, permits
cross-ID ProtT5 reuse only for identical complete sequence SHA-256 values, and
uses only explicit alias mappings. Structure, text, and PPI provenance
uncertainty is routed to manual review.

Run it from the repository root:

```bash
python scripts/verification/inventory_embeddings.py \
  --benchmark-dir /path/to/nine-csv-benchmark \
  --source-benchmark-dir /path/to/cache-source-benchmark \
  --embedding-cache /path/to/embedding_cache \
  --config configs/embedding_inventory.future.example.json \
  --policy maximize-coverage \
  --output-dir results/embedding_inventory/my_benchmark
```

See [`embedding_inventory/README.md`](embedding_inventory/README.md) for the
decision model, canonical CAFA3 configuration, output contract, alias format,
integration tests, future temporal/homology commands, and scientific
limitations.

## Local path configuration

Machine-specific paths should live outside committed scripts. Use
`configs/paths.example.sh` as a template:

``` bash
cp configs/paths.example.sh configs/paths.local.sh
```

Then edit `configs/paths.local.sh` for the current machine. Typical
variables include `PFP_DIR`, `CAFA_ASSESSMENT_DIR`, `CAFA3_RAW_DIR`,
`PROTEIN_DATABASES_DIR`, `STRING_H5_FILE`, `STRING_ALIAS_FILE`,
`CONDA_EXE`, `MMFP_ENV`, `MMFP_ENV_DIR`, `PFP_GIT_URL`,
`PFP_CLONE_DIR`, `PFP_EXTERNAL_DIR`, `PFP_DATA_DIR`,
`DEPENDENCY_ENV`, and `VERIFY_CSV_WORKDIR`.

## Contemporary benchmark builder

The 2025→2026 CAFA-style temporal benchmark builder lives in:

``` text
benchmark_builders/contemporary_cafa/
```

It is kept as a self-contained subproject so that the immutable PFP
reproduction wrappers and the new benchmark-generation code remain
separate.

The production runner and UCL/SGE entrypoint are:

```bash
bash scripts/benchmark_generation/run_contemporary_temporal_benchmark.sh
qsub hpc_jobs/active/hpc_contemporary_temporal_benchmark.sh
```

See `benchmark_builders/contemporary_cafa/README.md` for the temporal contract,
named profiles, required frozen inputs, reports and QC gates.

## Quick start

Run the lightest reproduction path first:

``` bash
bash scripts/reproduction/reproduce_eval_only.sh
```

Retrain using downloaded artefacts:

``` bash
bash scripts/reproduction/reproduce_retrain_eval.sh
```

Run the full embedding-generation workflow:

``` bash
bash scripts/reproduction/reproduce_embeddings_retrain_eval.sh
```

The final route clones/builds the upstream PFP environment before
invoking `scripts/embeddings/generate_embeddings_run_all.sh`.

## HPC jobs

Cluster submission wrappers live in `hpc_jobs/active/` and are intended
to be submitted with `qsub`. They contain Sun Grid Engine resource
directives, scratch-space setup, result-copy logic, and then invoke the
normal framework entrypoints under `scripts/reproduction/`.

Reusable implementation scripts stay under `scripts/`; the HPC job
wrappers should remain thin scheduler-facing launchers.

## Embedding-generation workflow

`scripts/embeddings/generate_embeddings_run_all.sh` is a sub-orchestrator. It assumes:

-   the upstream `PFP` repository has already been cloned;
-   the required Python environment is active;
-   the current working directory is the upstream `PFP` repository root.

Pipeline stages:

``` text
[0/8] External dependencies
[1/8] Data preparation
[2/8] Verify split contract
[3/8] Build proteins.fasta
[4/8] Sequence embeddings
[5/8] Text embeddings
[6/8] Structure embeddings
[7/8] PPI embeddings
[8/8] Verify generated embeddings
```

Both verification scripts are executed in **strict mode** as part of the
embedding-generation pipeline. Any contract violation stops the pipeline
before expensive downstream computation is launched.

The embedding verifier is called with an explicit configuration path:

``` bash
python "${REPO_ROOT}/scripts/verification/verify_embeddings.py" \
  --data-dir data \
  --config "${REPO_ROOT}/configs/cafa3.json" \
  --strict
```

Passing the configuration explicitly avoids relying on the current
working directory.

## Verification gates

### Split verification

`verify_splits.py` validates the internal data contract produced by the
data-preparation stage.

For each configured aspect/split it checks:

-   `{aspect}_{split}_names.npy` exists and is non-empty;
-   protein IDs are unique within each split;
-   `{aspect}_{split}_labels.npz` has one row per protein;
-   the label matrix column count matches `{aspect}_go_terms.json`;
-   every protein has a non-empty sequence;
-   `sequences.json` has no missing entries relative to `names.npy`;
-   train/valid/test protein IDs are disjoint within each aspect.

Usage:

``` bash
python scripts/verification/verify_splits.py --data-dir data --strict
```

For alternative aspect or split names:

``` bash
python scripts/verification/verify_splits.py \
  --data-dir data \
  --aspects BPO CCO MFO \
  --splits train valid test \
  --strict
```

### Embedding verification

`verify_embeddings.py` checks that generated embeddings are complete and
numerically valid before training.

For every protein ID in the configured dataset it verifies:

-   embedding presence;
-   coverage against the configured threshold;
-   sampled embedding dimension;
-   sampled numerical finiteness (`NaN`/`Inf` detection).

Usage:

``` bash
python scripts/verification/verify_embeddings.py \
  --data-dir data \
  --config configs/cafa3.json \
  --strict
```

The verifier is benchmark-agnostic. Benchmark-specific values live
entirely in JSON configuration files.

### Configuration

`configs/cafa3.json` fully specifies a benchmark configuration for
embedding verification. The verification logic is benchmark-agnostic;
adapting to a new benchmark requires creating a new configuration file
rather than modifying the Python source.

The configuration defines:

-   aspects;
-   splits;
-   cache directory;
-   catastrophic coverage threshold;
-   sample size;
-   modality cache directories;
-   expected embedding dimensions;
-   minimum coverage thresholds.

## Device selection

GPU embedding scripts default to CUDA:

``` bash
DEVICE=cuda
```

For local CPU testing:

``` bash
DEVICE=cpu bash scripts/embeddings/generate_embeddings_sequence.sh
DEVICE=cpu bash scripts/embeddings/generate_embeddings_structure.sh
```

CPU execution may be very slow for model-heavy stages such as ESM-IF1.

## Suggested development order

Before launching expensive cluster runs:

1.  Run the data-preparation stage.
2.  Run `scripts/verification/verify_splits.py` in strict mode.
3.  Build `data/proteins.fasta` with `scripts/embeddings/generate_embeddings_fasta.py`.
4.  Generate embeddings.
5.  Run `scripts/verification/verify_embeddings.py` in strict mode with an explicit
    configuration path.
6.  Start retraining and evaluation.

This verification layer is designed to catch malformed datasets,
incomplete embedding caches, incorrect embedding dimensions and
non-finite values before GPU-intensive training begins, reducing wasted
compute and improving reproducibility.
