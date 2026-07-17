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

The completion marker no longer invokes `git` from the Python runtime. The HPC
wrapper already resolves and verifies both commits before entering the workflow;
it passes those values as environment metadata to the final JSON writer. A
missing `git` executable inside the MMFP container therefore cannot invalidate
an otherwise complete embedding run.

## Archive-backed contemporary retries

The successful 2025_01 to 2026_02 supervisor-profile cache contains hundreds
of thousands of `.npy` files. Extracting all of them on SAN would exceed the
project store's file quota. The contemporary retry path therefore treats the
authenticated `contemporary_embedding_cache.tar.gz` as an immutable baseline
and stores only newly recovered arrays under `retry_state/cache/`. Together the
archive and retry delta form one logical cache.

`initialize_contemporary_embedding_state.sh` binds this state to:

- all nine contemporary CSV hashes;
- every target protein and sequence SHA-256 from the exact reuse plan;
- the benchmark build and reuse-plan manifests;
- PFP/framework commits, environment report, compatibility scripts and runtime
  policy;
- the baseline archive and its complete assembly report.

It verifies that the assembly report covers every target/modality pair and that
the archive contains exactly the arrays reported as available. It indexes the
archive without extracting it persistently. Contract drift is rejected.
The HPC initializer and retry wrappers explicitly bind their caller-selected
benchmark, reuse-plan, baseline and state directories into the MMFP Singularity
runtime; SAN is not visible inside that container unless it is bound.

The gate in `configs/contemporary_embedding_resume.json` scales Zijian's
published CAFA3 coverage proportions to the contemporary target count. For the
current 156,421 proteins this means:

```text
ProtT5      ceil(156,421 * 69,811 / 69,811) = 156,421
text        ceil(156,421 * 69,517 / 69,811) = 155,763
ESM-IF1     ceil(156,421 * 67,948 / 69,811) = 152,247
PPI         ceil(156,421 * 58,294 / 69,811) = 130,616
```

These are scaled coverage floors, not the raw CAFA3 counts. Retaining the raw
counts would make the larger contemporary benchmark pass with much poorer
relative coverage.

`run_contemporary_embedding_retry.sh` selects only currently missing pairs for
one modality. Twenty accepted controls are materialized temporarily from the
baseline archive or retry delta and regenerated with the requested subset. At
least five must match within `rtol=1e-5`, `atol=1e-6` before valid outputs are
atomically merged. Failures remain in the retry ledger regardless of reason.
The wrappers never edit PFP and always remove job-owned scratch.

The state records the framework commit used at initialization. During active
development, a later retry may run from a newer framework commit without
invalidating the scientific inputs. The default retry behavior therefore warns
on a framework-commit difference while still enforcing the exact PFP commit,
text cutoff, environment fingerprint, and hashes of all PFP extraction and
framework compatibility scripts recorded by the state. Use
`--strict-framework-commit` when a frozen release requires whole-repository
revision equality.

### Same-node text and structure diagnostic

`run_contemporary_embedding_reproducibility.sh` investigates a failed
subset-equivalence gate without changing the cumulative state. For either
`text` or `structure` it:

1. deterministically selects 20 accepted controls from the existing state;
2. materializes their accepted arrays into disposable scratch;
3. builds one shared control-only PFP input view;
4. freezes and hashes one exact text TSV or set of AlphaFold PDB files;
5. generates the controls twice as separate model invocations on the same GPU;
6. compares repeat 1 with repeat 2 and both repeats with the baseline;
7. records maximum/mean absolute difference, RMSE, L2 difference, relative
   difference, cosine similarity, existing `allclose` status and exact equality.

The corresponding HPC wrapper is pinned to `animal-206-2.local`, requests one
GPU, verifies the assigned host at runtime, and publishes the reference and two
small generated control caches with JSON, TSV and Markdown reports. It records
GPU/driver/PyTorch settings and hashes the exact sources and inputs. Numerical
differences are observations, not automatic failures; missing, malformed or
non-finite control arrays remain hard integrity failures.

There is deliberately no `merge` operation in this workflow. The experiment
does not silently relax the production tolerance or alter accepted SAN arrays.
Its repeat-to-repeat distribution provides the evidence needed for a later,
explicit tolerance decision. Model weights or authenticated AlphaFold PDBs may
be added to the state's non-scientific `source_cache`; the accepted-array
ledger, coverage and benchmark contract are not changed.

### Reusing frozen dependencies

Embedding workflows resolve static inputs in the same order as the rest of the
framework: explicit path, `ARTIFACT_CATALOG`, then network fallback. The
catalogued embedding inputs are the nine canonical CAFA3 CSVs, STRING v12.0
aliases/network embeddings, and Zijian's published embedding bundles. On UCL
HPC, pass:

```bash
--artifact-catalog /SAN/bioinf/bmpfp/manifests/artifact_paths.tsv
```

The catalogue does not replace dynamic per-protein acquisition. UniProt text,
AlphaFold structures not already present in resumable source state, and model
caches remain workflow-managed. This keeps the code portable while avoiding
repeat downloads of large frozen files.

## Lightweight validation

No network, model, or full-data access is needed for the focused contract
tests:

```bash
python -m unittest discover -s scripts/embeddings/tests -v
bash -n scripts/embeddings/run_contemporary_embedding_generation.sh
bash -n scripts/embeddings/initialize_contemporary_embedding_state.sh
bash -n scripts/embeddings/run_contemporary_embedding_retry.sh
bash -n scripts/embeddings/run_contemporary_embedding_reproducibility.sh
bash -n scripts/embeddings/generate_embeddings_structure.sh
bash -n hpc_jobs/active/hpc_contemporary_embedding_generation.sh
bash -n hpc_jobs/active/hpc_contemporary_embedding_state_initialize.sh
bash -n hpc_jobs/active/hpc_contemporary_embedding_retry.sh
bash -n hpc_jobs/active/hpc_contemporary_embedding_reproducibility.sh
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
