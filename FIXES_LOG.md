# Compatibility Fixes Log

This log records compatibility changes made around the immutable upstream PFP
code. Changes tagged `[compat]` must preserve the original model and embedding
semantics; behavior-changing corrections belong on a later `[fix]` branch.

## 2026-07-17 - `[compat]` Freeze homology GOA input to archived release 234

### Observed failure

Homology pilot `7071299.1` stopped before input staging because its runtime
guard discovered that EBI's mutable GOA current endpoint had advanced from
UniProt-GOA release 234 (2026-06-17) to release 235 (2026-07-08). The SAN
population attempt had left only an unpublished `.partial` whose size did not
match the pinned release-234 size.

### Compatibility change

- Replace the mutable GOA current URL with EBI's immutable historical
  `goa_uniprot_all.gaf.234.gz` URL in both the SAN catalogue and scratch-first
  homology runtime.
- Remove dependencies on live `current_release_numbers.txt` and the unavailable
  archived MD5 sidecar.
- Authenticate release 234 using its pinned byte size, pinned SHA-256, embedded
  `!date-generated: 2026-06-17`, and embedded GO version `2026-06-15`.
- Keep UniProt release `2026_02`, GOA release `234`, and the ontology snapshot
  unchanged; no experiment is silently upgraded to release 235.

### Scientific behavior

This corrects input addressing only. It preserves the previously selected
homology benchmark snapshot and changes no clustering, splitting, labelling,
term-filtering, model, or evaluation policy.

### Validation

- EBI returned HTTP 200 and content length `11,664,243,116` for the archived
  release-234 GAF, exactly matching the committed catalogue.
- Focused tests assert that production acquisition and runtime paths use only
  the immutable URL and pinned SHA-256.

## 2026-07-17 - `[compat]` Preserve and resume the contemporary embedding cache

### Observed failure

Contemporary job `7070469` completed all four preflight and full-generation
processes, assembled the combined cache, and wrote authenticated archives. Its
final completion-marker writer then tried to invoke `git` from inside the MMFP
Python container, where `git` is unavailable, so the wrapper labelled the
otherwise complete result `.failed`. The combined cache retained complete
ProtT5 coverage but incomplete text, structure, and PPI source coverage.

The project SAN also has a 250,000-file quota. Extracting the combined cache's
hundreds of thousands of per-protein arrays into a conventional persistent
state would exceed that quota.

### Compatibility change

- Remove the final Python runtime `git` subprocess. The host HPC wrapper passes
  the commits it already resolved and verified into the completion metadata.
- Preserve the authenticated combined archive as an immutable baseline and
  store only newly recovered arrays in one retry delta.
- Extend the resumable state manager through opt-in baseline arguments; CAFA3
  flat-cache behavior remains unchanged.
- Bind the contemporary state to the exact nine CSV hashes, protein sequence
  hashes, benchmark/reuse-plan manifests, environment, commits, scripts,
  runtime values, baseline archive, and assembly report.
- Add archive indexing, archive-member verification, scratch materialization,
  hydration, one-modality retry selection, control equivalence, and atomic
  delta merging.
- Scale Zijian's published CAFA3 coverage proportions to the contemporary
  target count instead of reusing CAFA3's smaller absolute counts.

### Scientific behavior

Models, weights, text cutoff, pooling, PPI source, IF1 extraction, benchmark
membership, and PFP missing-modality behavior are unchanged. The changes affect
completion reporting, persistent representation, validation, and retry scope.

### Validation

- All 24 focused embedding tests pass, including baseline archive verification,
  archive-plus-delta coverage, control materialization, retry merge, hydration,
  contemporary workspace selection, and no-runtime-git behavior.
- Python compilation, shell syntax checks, and `git diff --check` pass.
- The preserved SAN archive and all nine benchmark CSVs/five pickles pass
  SHA-256 verification.
- The state initializer must complete on the cluster before retry jobs are
  submitted.

## 2026-07-16 - `[compat]` Preserve partial embedding work and retry by modality

### Observed failure

Full CAFA3 reproduction job `7070501` generated complete ProtT5 coverage and
large valid text/PPI/IF1 caches, but IF1 reached only `41,748 / 69,811` because
the AlphaFold acquisition phase recorded `27,161` API/connection errors. The
independent 85% structure gate correctly blocked training, but unconditional
scratch cleanup then discarded every valid generated array.

### Compatibility change

- Add one provenance-bound persistent cache keyed by `(protein_id, modality)`.
- Validate and atomically merge every successful array before deciding whether
  the overall coverage gate passed.
- Retain every absent or invalid pair as `needs_retry`; reasons are diagnostic
  and never make a pair terminal.
- Add a one-modality retry workflow and Grid Engine wrapper. Accepted controls
  must demonstrate subset equivalence before a retry merge.
- Add full-workflow resume mode, which hydrates only authenticated arrays and
  refuses to train until the historical published-cache count floor passes.
- Persist one authenticated AlphaFold PDB source cache and acquire only missing
  PDBs through a bounded framework path that reuses PFP's mapping/API logic.
- Leave the pinned PFP checkout and its 1,000-worker main function unchanged.

### Scientific behavior

Model IDs, weights, pooling, PFP data preparation, text cutoff, STRING source,
IF1 compatibility logic, training, and evaluation are unchanged. The change is
operational: acquisition concurrency, durable checkpointing, independent array
validation, retry granularity, and failure signaling.

### Validation

- Fixture tests cover contract drift, partial merge, invalid arrays, pair-level
  retry, persistent gate state, retry workspace selection, subset equivalence,
  and authenticated PDB reuse.
- Existing embedding tests continue to pass.
- Shell syntax and `git diff --check` pass before publication.
- A one-modality CUDA retry remains required before the SAN state is treated as
  production validated.

## 2026-07-16 - `[compat]` Restore ESM-IF1 extraction in the frozen MMFP runtime

### Observed failure

Contemporary embedding job `7070127` produced no structure embeddings. The
job's zero-output guard stopped publication, but the upstream extractor had
already converted its dependency import failure into a successful process exit.

### Verified causes

1. The frozen MMFP runtime uses NumPy `2.0.2` and Biotite `0.38.0`. The installed
   Biotite wheel contains extensions compiled against the NumPy 1.x ABI, so
   importing `biotite.structure` fails under NumPy 2.0.2. A disposable
   NumPy `1.26.4` overlay was verified to import both `biotite.structure` and
   `esm.inverse_folding` successfully without changing the main environment.
2. PFP moves ESM-IF1 to CUDA and then calls fair-esm `get_encoder_output()`.
   fair-esm `2.0.0` calls its coordinate batch converter without a device, so
   the generated coordinates, confidence values, padding masks, and tokens stay
   on CPU while the model is on CUDA. This reproduces the reported
   `cuda:0`/`cpu` `index_select` mismatch.

### Compatibility change

- Keep the main MMFP environment on Zijian's supplied NumPy `2.0.2` pin.
- Install NumPy `1.26.4` into a job-local scratch overlay and expose it only to
  the ESM-IF1 extraction process.
- Validate the overlay by importing NumPy, Biotite structure, and fair-esm
  inverse-folding modules before extraction.
- Generate a checksummed runtime copy of PFP's IF1 extractor. The copy performs
  the same fair-esm encoder call after passing the model device explicitly to
  `CoordBatchConverter`; it also raises on dependency failure and when every PDB
  fails. The checked-out upstream PFP source is never edited.
- Write `if1_environment.json` and `pfp_if1_compatibility.json` into the run
  reports for provenance.

### Scientific behavior

The extraction model, weights, coordinates, tokenization, confidence values,
padding mask, encoder call, and output mean pooling are unchanged. The runtime
copy changes device placement and failure signaling only.

### Validation

- Disposable cluster overlay: NumPy `1.26.4`, Biotite `0.38.0`, and fair-esm
  inverse folding imported successfully.
- Compatibility-copy generation against pinned PFP source: passed.
- Generated extractor Python compilation: passed.
- Embedding workflow unit tests: passed (`7` tests).
- Shell syntax checks for all changed wrappers: passed.
- `git diff --check`: passed.
- A real CUDA extraction smoke test remains required on the cluster before the
  full embedding job is treated as recovered.

## 2026-07-16 - `[compat]` Separate the MMseqs2 release and binary identities

### Observed failure

Homology pilot `7070128` downloaded every frozen input and the official
MMseqs2 `18-8cc5c` release archive, then rejected the executable because
`mmseqs version` returned the full Git commit
`8cc5ce367b5638c4306c2d7cfc652dd099a4643f` rather than the release tag.

### Verified storage evidence

The pilot used `158.0 GiB` at its measured peak and the scratch filesystem
still had approximately `1.06 TiB` free. Scratch exhaustion did not cause the
temporary suspension or final failure.

### Compatibility change

- Retain `18-8cc5c` as the pinned release asset and minimum-feature version.
- Pin the complete binary-reported Git commit
  `8cc5ce367b5638c4306c2d7cfc652dd099a4643f`; the builder additionally checks
  that it begins with the commit prefix embedded in the release tag (`8cc5c`).
- Record the release tag, expected release version and observed binary identity
  separately.
- Preserve exact-match behavior for MMseqs2 builds that report the release tag
  directly.

This changes version interpretation and provenance only. It does not change
MMseqs2 commands, clustering parameters, input data, or benchmark policy.

## 2026-07-16 - `[compat]` Support the CAFA3 full audit on Morecambe's Git

### Observed failure

Full CAFA3 reproduction job `7070493` validated the frozen MMFP environment,
then failed before canonical input acquisition with `Unknown option: -C`.
Morecambe's installed Git does not support `git -C`. The substantive workflow
and its dependency downloader still contained four reachable uses even though
the HPC wrapper already used a portable directory helper.

### Compatibility change

- Replace all four reachable `git -C` calls with a `git_in_dir` helper that
  changes directory inside a subshell and runs ordinary Git.
- Preserve the same clean-checkout validation, detached commit checkout and
  exact commit recording.
- Leave the caller's working directory unchanged after each Git operation.
- Add a regression contract covering the HPC wrapper, substantive workflow and
  dependency downloader so `git -C` cannot be reintroduced into this execution
  chain.

### Scientific behavior

This changes command portability only. It does not change any downloaded input,
embedding algorithm, model configuration, training operation or evaluation
policy.

### Validation

- Shell syntax checks for all three execution-chain scripts: passed.
- Focused full-reproduction tests: passed (`6` tests).
- Reachable execution-chain `git -C` scan: zero matches.
- `git diff --check`: passed.
## 2026-07-17 - `[compat]` Scope the contemporary retry revision gate

### Observed failure

Contemporary text, structure and PPI retry jobs `7073862`-`7073864` stopped
before generation because the state was initialized at framework commit
`36e01c5`, while later storage/catalogue fixes moved the submitted checkout to
`d633a73`. The blanket repository revision equality treated unrelated active
development as scientific contract drift.

### Compatibility change

- Keep the state initialization commit as provenance.
- Warn, rather than fail, when a development retry uses a different framework
  commit.
- Continue to fail on a PFP commit mismatch, text-cutoff mismatch, environment
  fingerprint mismatch, or a changed hash for any recorded PFP extraction or
  framework IF1/PPI compatibility script.
- Add `--strict-framework-commit` for frozen-release audits that require exact
  whole-repository equality.

### Scientific behavior

No embedding algorithm, source model, pooling operation, benchmark population,
or acceptance threshold changes. The fix narrows the provenance gate from the
entire evolving repository to the runtime sources that can affect generated
arrays.
