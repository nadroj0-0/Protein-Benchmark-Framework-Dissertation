# Compatibility Fixes Log

This log records compatibility changes made around the immutable upstream PFP
code. Changes tagged `[compat]` must preserve the original model and embedding
semantics; behavior-changing corrections belong on a later `[fix]` branch.

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
