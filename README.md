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
├── reproduce_eval_only.sh                 # Download artefacts and reproduce evaluation only
├── reproduce_retrain_eval.sh              # Download artefacts, retrain, then evaluate
├── reproduce_embeddings_retrain_eval.sh   # Full pipeline entry point: embeddings -> retrain -> eval
├── generate_embeddings_run_all.sh         # Sub-orchestrator for all embedding-generation stages
├── generate_embeddings_dependencies.sh    # External dependency setup
├── generate_embeddings_prepare_data.sh    # Raw CSVs -> splits, labels, sequences
├── generate_embeddings_sequence.sh        # ProtT5 sequence embeddings
├── generate_embeddings_text.sh            # PubMedBERT / UniProt text embeddings
├── generate_embeddings_structure.sh       # ESM-IF1 structure embeddings
├── generate_embeddings_ppi.sh             # STRING/SPACE PPI embeddings
├── verify_embeddings.py                   # Embedding completeness/correctness gate
├── verify_splits.py                       # Split-contract verification gate
├── verify_csv.sh                          # CSV sanity checks
├── benchmark_builders/
│   └── contemporary_cafa/                 # 2025→2026 CAFA-style benchmark builder
├── scripts/
│   ├── data_acquisition/                  # HPC/raw database download and inspection helpers
│   └── hpc/                               # Cluster environment probes and utilities
├── configs/
│   ├── cafa3.json                         # Default CAFA3 verification config
│   └── paths.example.sh                   # Example local/HPC path configuration
└── HPC Cluster/                           # Example cluster submission scripts/notes
```

Generated data, model checkpoints, cloned upstream repositories, and
embedding caches are intentionally not committed. See `.gitignore` for
the excluded paths.

## Local path configuration

Machine-specific paths should live outside committed scripts. Use
`configs/paths.example.sh` as a template:

``` bash
cp configs/paths.example.sh configs/paths.local.sh
```

Then edit `configs/paths.local.sh` for the current machine. Typical
variables include `PFP_DIR`, `CAFA_ASSESSMENT_DIR`, `CAFA3_RAW_DIR`,
`PROTEIN_DATABASES_DIR`, `MMFP_ENV`, and `VERIFY_CSV_WORKDIR`.

## Contemporary benchmark builder

The 2025→2026 CAFA-style temporal benchmark builder lives in:

``` text
benchmark_builders/contemporary_cafa/
```

It is kept as a self-contained subproject so that the immutable PFP
reproduction wrappers and the new benchmark-generation code remain
separate.

## Quick start

Run the lightest reproduction path first:

``` bash
bash reproduce_eval_only.sh
```

Retrain using downloaded artefacts:

``` bash
bash reproduce_retrain_eval.sh
```

Run the full embedding-generation workflow:

``` bash
bash reproduce_embeddings_retrain_eval.sh
```

The final route clones/builds the upstream PFP environment before
invoking `generate_embeddings_run_all.sh`.

## Embedding-generation workflow

`generate_embeddings_run_all.sh` is a sub-orchestrator. It assumes:

-   the upstream `PFP` repository has already been cloned;
-   the required Python environment is active;
-   the current working directory is the upstream `PFP` repository root.

Pipeline stages:

``` text
[0/7] External dependencies
[1/7] Data preparation
[2/7] Verify split contract
[3/7] Sequence embeddings
[4/7] Text embeddings
[5/7] Structure embeddings
[6/7] PPI embeddings
[7/7] Verify generated embeddings
```

Both verification scripts are executed in **strict mode** as part of the
embedding-generation pipeline. Any contract violation stops the pipeline
before expensive downstream computation is launched.

The embedding verifier is called with an explicit configuration path:

``` bash
python "${HERE}/verify_embeddings.py" \
  --data-dir data \
  --config "${HERE}/configs/cafa3.json" \
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
python verify_splits.py --data-dir data --strict
```

For alternative aspect or split names:

``` bash
python verify_splits.py \
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
python verify_embeddings.py \
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
DEVICE=cpu bash generate_embeddings_sequence.sh
DEVICE=cpu bash generate_embeddings_structure.sh
```

CPU execution may be very slow for model-heavy stages such as ESM-IF1.

## Current known gap

The sequence embedding stage expects:

``` text
data/proteins.fasta
```

The wrapper does not yet guarantee generation of this FASTA before the
ProtT5 stage. This has intentionally been left unresolved until the full
pipeline has been exercised and the exact failure mode confirmed.

## Suggested development order

Before launching expensive cluster runs:

1.  Run the data-preparation stage.
2.  Run `verify_splits.py` in strict mode.
3.  Generate embeddings.
4.  Run `verify_embeddings.py` in strict mode with an explicit
    configuration path.
5.  Start retraining and evaluation.

This verification layer is designed to catch malformed datasets,
incomplete embedding caches, incorrect embedding dimensions and
non-finite values before GPU-intensive training begins, reducing wasted
compute and improving reproducibility.
