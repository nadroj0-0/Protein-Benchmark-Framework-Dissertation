# Contemporary CAFA Benchmark Builder

This package builds deterministic temporal protein-function benchmarks and
exports the nine wide CSV files consumed by PFP. It also retains the validated
historical DeepGOPlus and official CAFA-file conversion modes.

## Submitted Profiles

| Profile | Test eligibility at t0 | Training records | Submitted use |
|---|---|---|---|
| `supervisor` | no qualifying annotation in any ontology | Swiss-Prot + TrEMBL | global-NK |
| `supervisor-nk-lk` | no qualifying annotation in the evaluated ontology | Swiss-Prot + TrEMBL | NK+LK |

Both profiles use CAFA3 target taxa, the supervisor evidence-code set,
snapshot-membership endpoints, backfill allowance, split 0.9, seed 0, and
`min_count=50`. The NK+LK profile changes only the knowledge-eligibility and
per-ontology overlap contract.

The qualifying evidence codes are:

```text
EXP IDA IPI IMP IGI IEP HTP HDA HMP HGI HEP TAS NAS IGC RCA ND IC
```

The older `cafa3-reconstructed` and `contemporary-cafa3-style` profiles remain
available as explicit historical capabilities but are not the submitted
contemporary results.

## Temporal Contract

A test protein must:

1. exist in the target population at t0;
2. map unambiguously to a t1 UniProt record;
3. satisfy the selected profile's t0 knowledge rule;
4. acquire a qualifying post-cutoff t1 annotation;
5. retain its t0 sequence in exported test data;
6. pass ontology, sequence, binding, and term-universe validation.

Proteins first appearing at t1 are excluded. GO terms are canonicalised in the
source ontology and projected into the frozen t0 prediction graph. `alt_id` and
unambiguous `replaced_by` mappings are supported; `consider` is never chosen
automatically.

The submitted snapshot inputs are:

```text
UniProt 2025_01 and 2026_02
GOA UniProt releases 225 and 234
GO prediction graph 2025-02-06
t0 source-resolution graph 2025-03-16
t1 source-resolution graph 2026-06-19
```

Rows assigned on or before `2025-03-08` are not annotation gains. Rows after
the frozen GOA 234 endpoint are rejected.

The propagated t0 development population defines the ontology term universe.
Test data never contributes terms. CSV inputs are sorted before the seeded
90/10 development split, and exact-sequence isolation is validated within each
ontology.

## Run

Install the package:

```bash
python -m pip install -e benchmark_builders/contemporary_cafa
```

The repository runner resolves explicit inputs and then `ARTIFACT_CATALOG`:

```bash
DB_ROOT=/absolute/path/to/pfp_inputs \
ARTIFACT_CATALOG=/absolute/path/to/artifact_paths.tsv \
PROFILE=supervisor \
bash scripts/benchmark_generation/run_contemporary_temporal_benchmark.sh
```

Use `PROFILE=supervisor-nk-lk` for the NK+LK benchmark. The direct package CLI
is available with:

```bash
python -m cafa_benchmark_builder --help
```

`--uniprot-t0` and `--uniprot-t1` are repeatable so Swiss-Prot and TrEMBL can
be supplied together. Paths are never fixed in Python.

## Outputs

The PFP bundle is:

```text
bp-training.csv      bp-validation.csv      bp-test.csv
cc-training.csv      cc-validation.csv      cc-test.csv
mf-training.csv      mf-validation.csv      mf-test.csv
```

The build also publishes DeepGOPlus-shaped pickles, checksums, a
`build_manifest.json`, protein-flow and exclusion tables, annotation-gain and
knowledge-cohort reports, source-resolution diagnostics, and validation
reports. Diagnostic files are written before strict QC raises so a failed run
retains its exact exclusion evidence.

Production fails on missing or malformed inputs, unresolved source ontology
IDs, invalid dates, forbidden train/test overlap, duplicate IDs, non-binary
labels, within-ontology protein or sequence leakage, or an empty required
split. Expected temporal exclusions are counted rather than silently dropped.

## Historical Routes

These remain separate from the contemporary snapshot front end:

```text
--source-mode deepgoplus
--source-mode cafa3-files
```

The repository-level CAFA3 historical runner defaults to September-2016
training, official CAFA3 targets and ground truth, assigned-date proxy,
pre-t0-backfill exclusion, and the packaged ontology:

```bash
bash scripts/validation/run_cafa3_historical_validation.sh
```

The raw-GOA route is a sensitivity analysis of public snapshots, not a claim
to reconstruct the unavailable organizer snapshot exactly.
