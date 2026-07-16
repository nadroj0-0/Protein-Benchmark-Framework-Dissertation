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
├── benchmark_reuse_planner/                # Exact CSV-to-CSV reuse/regenerate partition
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
│   ├── embedding_inventory.contemporary.json
│   ├── embedding_inventory.homology.example.json
│   └── paths.example.sh                   # Example local/HPC path configuration
```

Generated data, model checkpoints, cloned upstream repositories, and
embedding caches are intentionally not committed. See `.gitignore` for
the excluded paths.

## Persistent frozen inputs

The persistent UCL project store is populated by the idempotent SAN acquisition
workflow rather than the older `$HOME/protein_databases` downloader:

```bash
bash scripts/data_acquisition/populate_san_frozen_inputs.sh --dry-run
bash scripts/data_acquisition/populate_san_frozen_inputs.sh --profile all
```

It freezes the temporal and homology database inputs, STRING v12.0, MMseqs2,
canonical CAFA3/DeepGOPlus references, and Zijian's published MMFP artefacts
under `/SAN/bioinf/bmpfp`. Downloads are resumable, release-guarded,
checksum-checked where a trusted checksum is known, structurally validated,
and accompanied by SHA-256/provenance sidecars. See
[`scripts/data_acquisition/README.md`](scripts/data_acquisition/README.md) for
profiles, storage estimates, verification modes, and the distinction from the
legacy home-directory script.

## Embedding inventory and reuse planning

`embedding_inventory/` provides a benchmark-agnostic planner that reads the
nine PFP-compatible CSVs directly, validates existing arrays, applies
modality-specific scientific reuse rules, and emits an operational two-bucket
plan (`reuse` or `regenerate`) without copying or generating
embeddings.

It supports `paper-faithful` and `maximize-coverage`, separate strict target
and permissive source contracts, exact-sequence-only cross-ID ProtT5 reuse, and
explicit diagnostic aliases. Exact canonical text/structure/PPI reuse additionally
requires deterministic benchmark, cache, archive, and PFP-reference proof;
cross-benchmark uncertainty is retained as a detailed reason and receives
action `regenerate`.

Run it from the repository root:

```bash
python scripts/verification/inventory_embeddings.py \
  --benchmark-dir /path/to/nine-csv-benchmark \
  --source-benchmark-dir /path/to/cache-source-benchmark \
  --embedding-cache /path/to/embedding_cache \
  --config configs/embedding_inventory.contemporary.json \
  --policy maximize-coverage \
  --report-level compact \
  --output-dir /path/to/non-repository-results/my_benchmark
```

See [`embedding_inventory/README.md`](embedding_inventory/README.md) for the
decision model, canonical CAFA3 configuration, output contract, alias format,
integration tests, future temporal/homology commands, and scientific
limitations.

To exercise the planner against a completed contemporary benchmark and
Zijian's real published cache on UCL Grid Engine:

```bash
qsub -v BENCHMARK_DIR=/path/to/contemporary/run/outputs \
  hpc_jobs/active/hpc_contemporary_embedding_inventory.sh
```

This downloads the published embedding archives only into scratch and returns
the two binary manifests, per-modality ID lists, and the ProtT5 generation
FASTA. It does not generate embeddings. The scheduler-neutral implementation is
`scripts/verification/run_contemporary_embedding_inventory.sh`; both a single
benchmark directory and per-CSV path overrides are supported.

The overnight generation-and-assembly workflow consumes the stricter CSV-only
reuse plan, regenerates only its `regenerate` partition, extracts only its
`reuse` partition from authenticated published caches, and packages the merged
cache for PFP:

```bash
qsub hpc_jobs/active/hpc_contemporary_embedding_generation.sh \
  --target-benchmark-dir /absolute/path/to/contemporary/outputs \
  --reuse-plan-dir /absolute/path/to/reuse/plan
```

The four modalities run concurrently using the existing three-GPU plus CPU
layout. Large dependencies, models, PDBs, and unpacked arrays remain in
job-owned scratch. Only compressed generated-cache archives, the final merged
cache archive, logs, provenance, coverage reports, and assembly reports are
published home. See [`scripts/embeddings/README.md`](scripts/embeddings/README.md)
and [`hpc_jobs/README.md`](hpc_jobs/README.md) for the contract and overrides.

## Local path configuration

Machine-specific paths should live outside committed scripts. Use
`configs/paths.example.sh` as a template:

``` bash
cp configs/paths.example.sh configs/paths.local.sh
```

Then edit `configs/paths.local.sh` for the current machine. Typical
variables include `PFP_DIR`, `CAFA_ASSESSMENT_DIR`, `CAFA3_RAW_DIR`,
`PROTEIN_DATABASES_DIR`, `STRING_H5_FILE`, `STRING_ALIAS_FILE`,
`CONDA_EXE`, `MMFP_ENV`, `MMFP_ENV_DIR`, `MMFP_PYTHON`,
`MMFP_TORCH_INDEX_URL`, `MMFP_PYG_WHEEL_BASE`,
`MMFP_SINGULARITY_DIR`, `MMFP_SINGULARITY_IMAGE`,
`MMFP_SINGULARITY_VENV`, `MMFP_SINGULARITY_IMAGE_URI`,
`PFP_GIT_URL`, `PFP_CLONE_DIR`, `PFP_EXTERNAL_DIR`, `PFP_DATA_DIR`,
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

Run the hardened full CAFA3 audit on the UCL cluster:

```bash
qsub hpc_jobs/active/hpc_cafa3_full_from_scratch_reproduction.sh
```

This newer route regenerates all modalities in parallel, compares them against
authenticated published embeddings without training on those published arrays,
trains fresh checkpoints, evaluates them against the paper values, and copies
back a compact provenance-rich report. See `hpc_jobs/README.md` for the exact
contract and output layout.

The final route clones/builds the upstream PFP environment before
invoking `scripts/embeddings/generate_embeddings_run_all.sh`.

New `mmfp` environments use Python `3.9.23` and the package versions supplied
by Zijian for the MMFP/PFP paper environment. Those versions are pinned
directly in `scripts/reproduction_common.sh`. Dependencies for which no version
was supplied, including `fair-esm` and the PyTorch-Geometric stack, remain
unpinned. PyTorch defaults to the official CUDA 12.6 wheel index, and the
compiled PyTorch-Geometric extensions are selected from the wheel index that
matches the installed PyTorch/CUDA build. An existing `mmfp` environment is
reused without modification, but every active workflow validates the exact
Python and supplied package versions, required imports, PyG binary compatibility,
and `pip check` before proceeding.

The official PyTorch 2.8 Linux wheels require glibc 2.28 or newer. Environment
creation checks this before creating a partial environment and fails with a
container-runtime instruction on older hosts such as CentOS 7; it never silently
substitutes a different Python or PyTorch version.

On the UCL CentOS 7 cluster, stop jobs that use `mmfp` and rebuild the compatible
Singularity-backed entrypoint with:

```bash
REBUILD_MMFP=YES bash scripts/environment/rebuild_mmfp_singularity.sh
```

See `scripts/environment/README.md` for the runtime layout and validation notes.

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
