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

For the hardened end-to-end historical experiment, use
`scripts/reproduction/run_cafa3_full_from_scratch_reproduction.sh` through its
HPC wrapper. That workflow still calls these modality wrappers, but adds the
current IF1/PPI compatibility copies, temporal CLS reduction, a reversible
bounded preflight, exhaustive output validation, and comparison against the
authenticated published cache before fresh training.

The preflight helper `prepare_cafa3_embedding_preflight.py` backs up all 27
prepared split view files plus the exact full `proteins.fasta`, subsets every
ontology/split, and restores each original byte only after authenticating its
SHA-256. This preserves the full FASTA ordering used for ProtT5 batching. It
never changes the canonical CSVs or upstream PFP source.

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

### ESM-IF1 compatibility

Zijian's stated NumPy 2.0.2 and Biotite 0.38.0 pins cannot import the compiled
`biotite.structure` extension from the published Biotite wheel. The main MMFP
environment keeps those author-supplied pins. During this workflow only, an
isolated NumPy 1.26.4 overlay is installed into job-owned scratch and exposed
only to the IF1 process. The overlay is reported in
`reports/if1_environment.json` and is removed with the rest of scratch.

fair-esm 2.0.0's `get_encoder_output()` creates its coordinate, confidence, and
mask tensors on CPU even when PFP moves the IF1 model to CUDA.
`build_pfp_if1_compat_copy.py` validates the exact upstream source blocks and
writes a separate scratch copy that creates the same encoder inputs explicitly
on the model device. It also exits non-zero when dependencies cannot import or
all PDB files fail. The source PFP checkout is not modified, and the source and
compatibility-copy hashes are recorded in `reports/pfp_if1_compatibility.json`.

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
bash -n scripts/embeddings/generate_embeddings_structure.sh
bash -n hpc_jobs/active/hpc_contemporary_embedding_generation.sh
```

The HPC wrapper performs the real bounded preflight before the full run. It
must be submitted from a clean committed framework checkout so the scratch
clone executes the reviewed revision.

## Resumable CAFA3 generation

The full CAFA3 reproduction no longer discards valid arrays when one external
source has incomplete coverage. `manage_resumable_embedding_state.py` owns one
persistent, benchmark-bound cache outside PFP. The default HPC location is:

```text
/SAN/bioinf/bmpfp/embedding_states/cafa3_full_reproduction/
```

Its contract includes the nine CSV hashes, every protein sequence SHA-256, the
PFP and framework commits, the author environment report, policy hash, text
cutoff, GO/STRING inputs, and exact upstream/compatibility script hashes. An
existing state is rejected if any contract field changes. Persistent arrays are
accepted only after safe-ID, numeric dtype, finite-value, exact-dimension, and
SHA-256 checks, then copied atomically under a filesystem lock.

There are only two pair states:

| State | Meaning |
|---|---|
| `accepted` | A validated array exists in the one cumulative cache. |
| `needs_retry` | No accepted array exists; retry regardless of the diagnostic reason. |

`failure_ledger.tsv` retains one cumulative row per currently failed pair with
attempt count and latest reason. It does not create per-attempt embedding cache
copies. `needs_retry.tsv`, `coverage.json`, and one of
`GENERATION_INCOMPLETE.json` or `EMBEDDING_GATE_PASSED.json` are regenerated
atomically after every merge.

The historical gate is tied to the published CAFA3 cache counts, not the older
generic lower bounds:

```text
ProtT5      69,811 / 69,811
text        69,517 / 69,811
ESM-IF1     67,948 / 69,811
PPI         58,294 / 69,811
```

Pairs can remain in `needs_retry` after the gate passes; the marker means the
historical published coverage floor has been reached, not that every external
source contains every protein.

### AlphaFold acquisition

PFP's checked-in `check_alphafold_coverage.py` and its hard-coded 1,000-worker
main path remain unchanged. The resumable workflow selects the opt-in
`framework-bounded` mode in `generate_embeddings_structure.sh`. The framework
imports PFP's CAFA-to-UniProt mapping and AlphaFold interpretation functions,
calls them with eight workers for only uncached IDs, downloads PDBs atomically,
records source URLs, versions and SHA-256 values, and copies the requested PDB
view into job scratch for IF1 inference. Authenticated PDBs live once under the
state `source_cache`; a killed retry does not force their acquisition again.

### Retry and resume sequence

1. Submit the full workflow in `initial` mode. If coverage is insufficient, it
   publishes a normal `.incomplete` report and exits zero after saving every
   valid array. It does not train or evaluate.
2. Submit `hpc_cafa3_embedding_retry.sh` once per modality that still has
   missing pairs. Each job builds a PFP view containing only that modality's
   missing IDs plus 20 accepted controls.
3. At least five controls must regenerate and match the accepted arrays within
   `rtol=1e-5`, `atol=1e-6` before retry outputs can be merged. Source-unavailable
   controls are reported separately; numerical differences fail loudly.
4. Once `EMBEDDING_GATE_PASSED.json` exists, rerun the full wrapper with
   `--embedding-mode resume`. It hydrates only validated arrays, then performs
   the published-cache comparison, training, and evaluation.

Failure reason never decides retry eligibility. A 404, mapping absence, API
error, invalid array, process failure, or source absence all remain
`needs_retry`; the reason exists for analysis only.

The persistent staging directory should be removed only after a final embedding
release has been authenticated and published elsewhere. No cleanup command is
automated because deleting persistent SAN state is intentionally a manual,
reviewed operation.
