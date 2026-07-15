# PFP benchmark reuse planner

This package is an independent, benchmark-level planner. It reads only the nine
PFP-compatible CSV files in each named benchmark and partitions every unique
target protein into exactly one of two actions:

- `reuse`: at least one embedded benchmark contains the same case-sensitive
  protein ID and the identical complete sequence.
- `regenerate`: the ID is absent from the embedded union, or the union contains
  that ID with a different sequence. All four modalities (`prott5`, `text`,
  `structure`, and `ppi`) are scheduled.

GO labels and ontology/split placement are validated and reported as
memberships, but never affect the action. Sequence-only matching is forbidden.
An ID conflict anywhere in the embedded union or within the target is fatal.
`reuse` preserves the earlier run's files and missing/masked-modality behaviour;
it does not claim that all four modality files exist. `regenerate` only schedules
generation. Later generation failures and unavailable external sources belong to
generation logs and are outside this planner's scope.
The package does not import the existing `embedding_inventory` planner and does
not locate, open, validate, hash, or otherwise inspect embedding files.

## Install and run

From the repository root:

```bash
python -m pip install -e benchmark_reuse_planner
```

One embedded benchmark:

```bash
python -m pfp_benchmark_reuse plan \
  --embedded-benchmark 'cafa3=/path/to/cafa3 csvs' \
  --target-benchmark 'homology=/path/to/homology csvs' \
  --output-dir '/path/to/reuse plan'
```

Multiple embedded benchmarks:

```bash
python -m pfp_benchmark_reuse plan \
  --embedded-benchmark 'cafa3=/path/to/cafa3 csvs' \
  --embedded-benchmark 'contemporary=/path/to/contemporary csvs' \
  --target-benchmark 'homology=/path/to/homology csvs' \
  --output-dir '/path/to/reuse plan'
```

For a checkout-only invocation without installation, prefix the command with
`PYTHONPATH=benchmark_reuse_planner/src`.

`--embedded-benchmark NAME=PATH` is repeatable and required. Exactly one
`--target-benchmark NAME=PATH` is accepted. Names must match
`[A-Za-z0-9][A-Za-z0-9._-]*` and must be unique across all inputs. The first
`=` separates the name, so paths may contain spaces or additional `=` signs.

## Input contract

Each input directory must provide:

```text
bp-training.csv     bp-validation.csv     bp-test.csv
cc-training.csv     cc-validation.csv     cc-test.csv
mf-training.csv     mf-validation.csv     mf-test.csv
```

Every CSV begins with `proteins,sequences` followed by one or more unique `GO:`
columns. Protein IDs and sequences are non-empty, rows match the header width,
and every GO label is the literal `0` or `1`. IDs are opaque and case-sensitive;
whitespace, controls, path separators, `.` and `..` are rejected. Sequences are
not stripped, uppercased, or otherwise normalized.

Identical occurrences across CSVs are deduplicated and all membership filenames
are retained. A repeated ID with a different sequence is rejected. Exact duplicate
rows inside one CSV are accepted; contradictory duplicate rows are rejected.
Unrelated non-required files in an input directory are ignored.

## Outputs

The destination must not already exist. Reports are written to a unique sibling
staging directory, checked, and renamed only after success. Cleanup handles
ordinary failures and Python interrupts. `RUN_COMPLETE.json` is the final file
written in staging.

The package publishes:

```text
reuse_proteins.tsv
regenerate_proteins.tsv
reuse_proteins.txt
regenerate_proteins.txt
regenerate_proteins.fasta
known_embedded_proteins.tsv
summary.json
summary.md
run_manifest.json
output_manifest.json
RUN_COMPLETE.json
```

Both action TSVs use the same columns:

| Column | Meaning |
|---|---|
| `protein_id` | Exact target protein ID |
| `sequence` | Complete target sequence |
| `sequence_sha256` | SHA-256 of the exact UTF-8 sequence |
| `action` | Exactly `reuse` or `regenerate` |
| `reason` | Stable decision token |
| `matching_embedded_benchmarks` | Sorted JSON list of exact-pair matches |
| `embedded_benchmark_memberships` | Sorted JSON list of `name:csv-file` memberships for the ID |
| `target_memberships` | Sorted JSON list of target CSV filenames |
| `regenerate_modalities` | Empty for reuse; sorted JSON list of all four modalities otherwise |

Reason tokens are `exact-id-sequence-match`, `protein-id-absent`, and
`sequence-mismatch`.

`known_embedded_proteins.tsv` contains the complete deduplicated embedded union,
including proteins absent from the target. TXT files contain one sorted bare ID
per line. FASTA contains every regenerate protein once with a bare `>protein_id`
header and the exact target sequence on one line.

The known-protein TSV columns are `protein_id`, `sequence`, `sequence_sha256`,
`embedded_benchmarks`, and `embedded_benchmark_memberships`. `summary.json`
records the comparison policy, action semantics, benchmark paths, regenerate
modalities, and counts; `summary.md` presents the same decision model for readers.

`run_manifest.json` records canonicalized command arguments, benchmark names,
resolved paths, every input CSV size/SHA-256 identity, and counts. Embedded
arguments are sorted by name, so their CLI order cannot change the result.
`output_manifest.json` records the relative path, byte size, and SHA-256 of the
nine payload reports written before it. It excludes itself and the subsequently
written completion marker; `RUN_COMPLETE.json` records the output-manifest
identity.

## Development checks

```bash
cd benchmark_reuse_planner
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src python -m unittest discover -s tests -v
PYTHONPYCACHEPREFIX=/tmp/pfp-reuse-pycache python -m compileall -q src tests
ruff check --no-cache src tests
```
