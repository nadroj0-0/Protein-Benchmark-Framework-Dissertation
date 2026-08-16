# Embedding Workflows

The retained workflows generate and authenticate the four PFP modalities:
ProtT5 sequence, temporal/current UniProt text, ESM-IF1 structure, and STRING
PPI. Upstream PFP is pinned to
`1e04fd6d6d3c40458fd41ec1a881ed6e24de768e`.

## CAFA3

The fresh CAFA3 driver prepares the benchmark, generates all modalities,
validates coverage, publishes resumable state, and compares regenerated arrays
with authenticated published arrays:

```bash
bash scripts/reproduction/run_cafa3_full_from_scratch_reproduction.sh --help
```

If the state is incomplete, continue one modality with:

```bash
bash scripts/reproduction/run_cafa3_embedding_retry.sh --help
```

Retry output is merged only after subset-equivalence and state-contract checks.
The full driver stops after embedding comparison; model training belongs to
`scripts/model_execution/run_pfp_benchmark.sh`.

## Contemporary Global-NK

The submitted route is:

1. build a paper-faithful inventory and exact reuse plan;
2. generate the coarse-plan cache, published at
   `archive/contemporary_embedding_cache.tar.gz`;
3. initialize its archive-backed state and run any required missing-pair retries;
4. finalize that coarse base after its embedding gate passes;
5. run `run_contemporary_embedding_retry.sh --modality text
   --generate-all-text-only --text-cutoff-date 2025-03-08` against the retained
   state contract;
6. compose the finalized base and complete corrected text layer with
   `compose_contemporary_text_replacement.py`;
7. initialize a new state from the composed root; and
8. finalize that new state as the accepted corrected cache.

The all-target replacement is part of the accepted from-scratch contract. It
prevents a reused published text layer from silently retaining an older cutoff;
it is not included to reproduce the earlier mistaken result.

Entry points:

```bash
bash scripts/verification/run_contemporary_embedding_inventory.sh --help
bash scripts/verification/run_benchmark_reuse_plan.sh --help
bash scripts/embeddings/run_contemporary_embedding_generation.sh --help
bash scripts/embeddings/initialize_contemporary_embedding_state.sh --help
bash scripts/embeddings/run_contemporary_embedding_retry.sh --help
python scripts/embeddings/compose_contemporary_text_replacement.py --help
python scripts/embeddings/finalize_embedding_state.py --help
```

The finalizer is therefore invoked twice with different roots: first to
authenticate the coarse non-text base used by composition, then to authenticate
the composed corrected cache. The earlier text layer is never a reported model
input.

`configs/contemporary_embedding_resume.json` defines acceptance dimensions,
coverage floors, cutoff, and paper-faithful PPI policy.

Evidence-only hash upgrades remain directly available without retiring source
bytes:

```bash
python scripts/embeddings/manage_resumable_embedding_state.py \
  upgrade-evidence-hashes --help
```

The finalizer is a separate destructive consolidation step. Its full required
arguments must be supplied, and `--retire-source-embeddings` explicitly permits
retiring baseline/retry embedding bytes after the final archive is verified.

## Pair-Resolved NK+LK and Homology

The coarse reuse planner proves target ID/sequence eligibility. Resolve actual
source archives separately:

```bash
python scripts/embeddings/resolve_embedding_reuse_sources.py \
  --coarse-plan-dir /path/to/coarse_plan \
  --output-dir /path/to/new_ledger \
  --cache-source contemporary=contemporary-global-nk=/path/to/cache.tar.gz=SHA256 \
  --source-text-policy contemporary=same-role \
  --cache-source cafa3=cafa3-regenerated=/path/to/cache.tar.gz=SHA256 \
  --source-text-policy cafa3=never
```

For homology, both authenticated sources use `source-current` text policy.
Source priority is the order of `--cache-source` arguments: contemporary before
CAFA3.

Generate each required modality in a new output directory:

```bash
bash scripts/embeddings/run_homology_embedding_modality.sh \
  --pfp-root /path/to/pinned/PFP \
  --work-dir /path/to/new_work \
  --output-dir /path/to/new_modality_run \
  --benchmark-dir /path/to/nine_csvs \
  --ledger-dir /path/to/ledger \
  --modality sequence \
  --policy configs/contemporary_nk_lk_embedding_generation.json
```

For NK+LK text, that policy requires cutoff `2025-03-08`; both the modality
runner and final pair assembler reject a missing or different cutoff.

Repeat for `text`, `structure`, and `ppi`; temporal NK+LK text additionally
uses `--text-cutoff-date 2025-03-08`. Homology uses
`configs/homology_embedding_generation.json` and current text.

Assemble only authenticated completed modality runs:

```bash
python scripts/embeddings/assemble_pair_resolved_embedding_cache.py \
  --ledger-dir /path/to/ledger \
  --generated-run sequence=/path/to/sequence_run \
  --generated-run text=/path/to/text_run \
  --generated-run structure=/path/to/structure_run \
  --generated-run ppi=/path/to/ppi_run \
  --policy configs/contemporary_nk_lk_embedding_generation.json \
  --output-archive /path/to/new_embedding_cache.tar.gz \
  --report-dir /path/to/new_assembly_report
```

The assembler verifies every completion marker's modality, ledger hash, policy
hash, PFP revision, archive location, and archive SHA-256 before reading it.
Use `bind_embedding_archive_evidence.py` to bind the final cache to its benchmark
and policy before model execution.

## Runtime Notes

All output paths must be outside the source checkout. `PYTHON_BIN` may select
the interpreter. ESM-IF1 uses a standard `PYTHONPATH` overlay where required;
no container-specific environment variable is necessary. AlphaFold downloads
are bounded and recorded, and model caches may be placed through `HF_HOME` and
`TORCH_HOME`.
