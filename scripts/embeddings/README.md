# Embedding Workflows

The scripts in this directory are scheduler-neutral wrappers around the
upstream PFP embedding scripts. Run them from a PFP repository root with the
framework paths supplied by the caller. PFP is treated as immutable upstream
code.

## Original CAFA3 path

`generate_embeddings_run_all.sh` follows the published PFP README route. Its
parallel mode runs ProtT5, text, and ESM-IF1 on three separate GPUs while STRING
PPI extraction runs on CPU. The modality wrappers remain the reference path for
the contemporary workflow:

```text
generate_embeddings_sequence.sh
generate_embeddings_text.sh
generate_embeddings_structure.sh
generate_embeddings_ppi.sh
```

## Contemporary reuse and regeneration

`run_contemporary_embedding_generation.sh` consumes two immutable inputs:

1. Nine PFP-compatible contemporary CSVs.
2. A completed `benchmark_reuse_planner` plan containing exact binary
   `reuse` and `regenerate` partitions.

It binds the plan to the nine CSV SHA-256 values, creates minimal PFP split
views for only the regenerate proteins, runs a small parallel preflight, and
then resumes the same caches for the full parallel run.

The workflow deliberately keeps action and source coupled:

| Planner action | Permitted cache source |
|---|---|
| `reuse` | Authenticated arrays from Zijian's published archive |
| `regenerate` | Arrays created during this job |

A regenerate protein never falls back to a published array. Generated files
outside the regenerate partition are rejected. Every selected array must have
the exact PFP dimension, a numeric dtype, and finite values.

ProtT5 must exist for every regenerate protein. Text, structure, and PPI can be
absent when their external source is unavailable; each absence is listed by ID
and the final PFP cache omits that array so the existing PFP loader applies its
normal mask and zero vector.

### Temporal text

PFP's text extraction functions are called directly by
`run_pfp_temporal_text.py`, with the benchmark's t0 date supplied at runtime.
The PFP source file is not edited. `generate_embeddings_text_temporal_cls.sh`
runs PFP's PubMedBERT embedder and simultaneously invokes
`reduce_text_embeddings_to_cls.py`. Each completed `(1, L, 768)` hidden-state
array is atomically replaced by `array[0, 0, :].astype(float32)`, the exact
`(768,)` CLS vector selected by PFP's flatten-and-truncate loader behavior and
confirmed by Zijian. This streaming reduction is required to keep scratch use
bounded.

### PPI compatibility

PFP's PPI extractor divides a diagnostic mapped-CAFA count by the number of
CAFA IDs. Contemporary IDs are UniProt accessions, so that diagnostic set is
empty. `build_pfp_ppi_compat_copy.py` validates the exact upstream expression
and writes a separate scratch copy with a denominator guard. The generated PPI
values and mapping logic are unchanged; `generate_embeddings_ppi.sh` keeps its
original source path by default and accepts the copied path through
`PPI_EXTRACT_SCRIPT` only for this workflow.

## Final artifact

`assemble_contemporary_embedding_cache.py` writes the standard PFP layout:

```text
data/embedding_cache/
├── prott5/
├── exp_text_embeddings_temporal/
├── IF1/
└── ppi/
```

The HPC wrapper packages that directory as
`contemporary_embedding_cache.tar.gz`; unpacked arrays remain in scratch. The
assembly report records source, action, availability, dtype, dimension, and
transfer method for every protein and modality. The compact reuse-plan
completion, output, run, and summary manifests are copied beside the assembly
reports, while the acquisition table records their hashes and original paths.

## Lightweight validation

No network, model, or full-data access is needed for the focused contract
tests:

```bash
python -m unittest discover -s scripts/embeddings/tests -v
bash -n scripts/embeddings/run_contemporary_embedding_generation.sh
bash -n hpc_jobs/active/hpc_contemporary_embedding_generation.sh
```

The HPC wrapper performs the real bounded preflight before the full run. It
must be submitted from a clean committed framework checkout so the scratch
clone executes the reviewed revision.
