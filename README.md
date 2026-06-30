# Protein Benchmark Framework

A reproducibility and benchmarking wrapper around the PFP/MMFP protein-function-prediction pipeline. The repository is intended to make the original CAFA3-style workflow easier to rerun, audit, and eventually swap onto cleaner or decontaminated benchmarks.

The framework currently supports three levels of reproduction:

1. **Evaluation only** using downloaded author artefacts.
2. **Retrain and evaluate** using downloaded data/embeddings.
3. **Regenerate embeddings, retrain, and evaluate** from the prepared benchmark pipeline.

The long-term design goal is to separate benchmark-specific values from reusable logic. Verification scripts now enforce data and embedding contracts before expensive training is launched.

## Repository layout

```text
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
├── verify_csv.sh                          # CSV sanity checks
├── configs/
│   ├── cafa3.json                         # Default CAFA3 verification config
│   └── verify_splits.py                   # Split-contract verification gate
└── HPC Cluster/                           # Example cluster submission scripts/notes
```

Generated data, model checkpoints, cloned upstream repositories, and embedding caches are intentionally not committed. See `.gitignore` for the excluded paths.

## Quick start

Run the lightest reproduction path first:

```bash
bash reproduce_eval_only.sh
```

Retrain using downloaded artefacts:

```bash
bash reproduce_retrain_eval.sh
```

Run the fuller embedding-generation workflow:

```bash
bash reproduce_embeddings_retrain_eval.sh
```

That final route calls `generate_embeddings_run_all.sh` after cloning/building the upstream PFP environment and changing into the upstream repository root.

## Embedding-generation workflow

`generate_embeddings_run_all.sh` is a sub-orchestrator, not the usual top-level entry point. It assumes:

- the upstream `PFP` repository has already been cloned;
- the Python environment required by MMFP/PFP is active;
- the current working directory is the upstream `PFP` repository root;
- this wrapper repository is available via the script path used by the caller.

The current stage order is:

```text
[0/7] External dependencies
[1/7] Data preparation
[2/7] Verify split contract
[3/7] Sequence embeddings
[4/7] Text embeddings
[5/7] Structure embeddings
[6/7] PPI embeddings
[7/7] Verify generated embeddings
```

The split verifier runs immediately after data preparation. This catches malformed split outputs before launching expensive embedding jobs.

The embedding verifier runs after all modality-generation scripts. It is called with an explicit config path:

```bash
python "${HERE}/verify_embeddings.py" \
  --data-dir data \
  --config "${HERE}/configs/cafa3.json" \
  --strict
```

Passing the config explicitly avoids relying on the current working directory.

## Verification gates

### Split verification

`configs/verify_splits.py` validates the internal data contract produced by the data-preparation stage.

For each configured aspect/split it checks that:

- `{aspect}_{split}_names.npy` exists and is non-empty;
- protein IDs within a split are unique;
- `{aspect}_{split}_labels.npz` has one row per protein;
- the label matrix column count matches `{aspect}_go_terms.json`;
- every protein in `names.npy` has a non-empty sequence;
- `sequences.json` has no missing entries relative to `names.npy`;
- train/valid/test protein IDs are disjoint within each aspect.

Current usage:

```bash
python configs/verify_splits.py --data-dir data --strict
```

For non-CAFA aspect or split names, use:

```bash
python configs/verify_splits.py \
  --data-dir data \
  --aspects BPO CCO MFO \
  --splits train valid test \
  --strict
```

A future cleanup should make this script read the shared JSON config used by `verify_embeddings.py`.

### Embedding verification

`verify_embeddings.py` checks that generated embeddings are complete and numerically usable before training.

For every protein ID appearing in the configured splits, it checks each configured modality for:

- `{protein_id}.npy` presence;
- coverage against the expected minimum threshold;
- sampled embedding dimension;
- sampled numerical finiteness (`NaN`/`Inf` detection).

Run directly with:

```bash
python verify_embeddings.py \
  --data-dir data \
  --config configs/cafa3.json \
  --strict
```

The verifier is benchmark-agnostic: benchmark-specific values live in JSON config files. The default CAFA3 config is:

```text
configs/cafa3.json
```

This config defines:

- `aspects`
- `splits`
- `cache_dir`
- `catastrophic_factor`
- `sample_size`
- modality cache directories, expected dimensions, and minimum coverage thresholds

To adapt the embedding verifier to a new benchmark, create a new JSON config and pass it via `--config`.

## Config design

`configs/cafa3.json` fully specifies the embedding-verification run. There is no partial merge with in-script defaults. This is deliberate: a config file should be explicit enough that reviewers can see exactly which aspects, splits, modalities, dimensions, and thresholds are being validated.

Example skeleton:

```json
{
  "aspects": ["BPO", "CCO", "MFO"],
  "splits": ["train", "valid", "test"],
  "cache_dir": "embedding_cache",
  "catastrophic_factor": 0.6,
  "sample_size": 25,
  "modalities": {
    "sequence": {
      "dirs": ["prott5"],
      "dim": 1024,
      "min_coverage": 0.99
    }
  }
}
```

If a benchmark has a different modality set, the config should list exactly the modalities that are expected for that benchmark.

## Device selection

GPU embedding scripts default to CUDA:

```bash
DEVICE=cuda
```

For local CPU testing, override the environment variable:

```bash
DEVICE=cpu bash generate_embeddings_sequence.sh
DEVICE=cpu bash generate_embeddings_structure.sh
```

CPU execution may be very slow for model-heavy stages such as ESM-IF1.

## Current known gaps

The sequence embedding stage expects:

```text
data/proteins.fasta
```

At present, the pipeline comments this as a known gap: no current wrapper step guarantees that this FASTA is generated before ProtT5 runs. This has deliberately not been patched yet; confirm the failure mode during a real run before adding the FASTA-generation step.

`verify_splits.py` and `verify_embeddings.py` also currently use different configuration mechanisms. `verify_embeddings.py` reads the JSON config; `verify_splits.py` still uses CLI `--aspects` and `--splits`. A later refactor should make both verifiers consume the same benchmark config.

## Suggested development order

Before launching expensive cluster runs:

1. Run the data-preparation stage.
2. Run `configs/verify_splits.py` in strict mode.
3. Generate embeddings.
4. Run `verify_embeddings.py` in strict mode with an explicit config path.
5. Only then start retraining/evaluation.

This avoids wasting GPU/cluster time on malformed splits, partial embedding caches, wrong dimensions, or non-finite embeddings.
