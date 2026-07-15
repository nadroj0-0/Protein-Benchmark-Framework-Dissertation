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
