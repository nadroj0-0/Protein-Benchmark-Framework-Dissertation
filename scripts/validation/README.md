# CAFA3 Validation

The retained validation routes answer different questions and should not be
collapsed into one claim.

## Released Pickles to Nine CSVs

```bash
bash scripts/validation/run_cafa3_deepgoplus_validation.sh --help
```

This validates the transformation from released DeepGOPlus-shaped pickles to
the nine CSV files consumed by PFP.

## Official Sources to Pickles

```bash
bash scripts/validation/run_cafa3_deepgoplus_pickle_generation_validation.sh --help
```

This reconstructs and compares the official-source pickle layer before CSV
conversion. `compare_deepgoplus_pickles.py` imports shared comparison logic
from `compare_cafa3_outputs.py`.

## Historical Reconstruction

```bash
bash scripts/validation/run_cafa3_historical_validation.sh
```

The primary submitted defaults are:

```text
training snapshot:       September 2016
target universe:         official CAFA3 targets
test source:             official ground truth
endpoint policy:         assigned-date proxy
backfill policy:         exclude pre-t0
benchmark ontology:      packaged DeepGOPlus ontology
```

The official-ground-truth route preserves released CAFA IDs and target FASTA
sequences and gates the expected 3,328 proteins plus ontology-specific export
populations. It validates the artifact path consumed by PFP; it does not claim
to reconstruct an unavailable private organizer snapshot.

The optional raw-GOA route compares public GAF snapshots. Its endpoint,
backfill, and ontology controls are explicit environment variables documented
by `--help`. Because the later public GAF was generated after the CAFA3
organizer snapshot, raw-GOA output is a sensitivity analysis rather than exact
organizer reconstruction.

## Performance and Inputs

Plain and gzip inputs are streamed. If `pigz` is available it may be used for
decompression; Python gzip remains the fallback. Temporary extraction defaults
to `${TMPDIR:-/tmp}` and can be redirected through the runner's explicit work
paths.

Optional local source overrides and `ARTIFACT_CATALOG` avoid repeated downloads.
Every accepted source is hashed and recorded. Validation output includes
machine-readable comparisons, overlap and population diagnostics, logs,
checksums, and completion markers.

These checks establish artifact/schema/population equivalence within their
declared inputs. They do not establish that public historical snapshots are
biologically identical to unavailable private freezes.
