# Implementation Report

## Summary

Implemented a separate benchmark-building codebase for a contemporary
2025->2026 CAFA-style temporal benchmark. It does not modify PFP. Its final
target is the nine wide CSVs consumed by PFP's `prepare_cafa3_data.py`.

The code is located under:

```text
src/cafa_benchmark_builder/
```

## Files Created

```text
pyproject.toml
requirements.txt
README.md
REPORT.md

src/cafa_benchmark_builder/__init__.py
src/cafa_benchmark_builder/__main__.py
src/cafa_benchmark_builder/config.py
src/cafa_benchmark_builder/io_utils.py
src/cafa_benchmark_builder/models.py
src/cafa_benchmark_builder/parsers.py
src/cafa_benchmark_builder/goa.py
src/cafa_benchmark_builder/ontology.py
src/cafa_benchmark_builder/builder.py
src/cafa_benchmark_builder/cli.py

tests/test_builder_smoke.py
tests/fixtures/go-mini.obo
tests/fixtures/uniprot-t0.fasta
tests/fixtures/uniprot-t1.fasta
tests/fixtures/goa-t0.gaf
tests/fixtures/goa-t1.gaf
```

## Module Roles

### `config.py`

Defines benchmark configuration and the final CAFA3 evidence policy:

```text
EXP, IDA, IPI, IMP, IGI, IEP, TAS, IC
```

Also contains ontology prefix/namespace mappings:

```text
bp -> biological_process
cc -> cellular_component
mf -> molecular_function
```

### `io_utils.py`

Small helper for opening plain text or `.gz` files. Used by all parsers so GOA
and UniProt can be streamed directly from compressed snapshots.

### `models.py`

Dataclasses for parsed records:

```text
ProteinRecord
GafRecord
```

### `parsers.py`

Parses UniProt sequence files.

Supported:

```text
FASTA / FASTA.gz
DAT / DAT.gz
```

For FASTA it extracts:

```text
accession
sequence
OX taxon
reviewed/unreviewed if header begins sp| or tr|
entry name
```

For DAT it extracts:

```text
primary accession
all accessions
sequence
NCBI taxon
reviewed status
entry name
```

### `goa.py`

Streams GAF rows. It does not load a full 20+ GB GOA file into memory row by
row.

Filters:

```text
DB == UniProtKB
evidence in final CAFA3 policy
aspect in P/C/F
optional target taxa
NOT removed
```

### `ontology.py`

DeepGOPlus-compatible OBO parser, closely based on `utils.Ontology`.

Preserved behavior:

```text
Ontology(..., with_rels=True)
get_anchestors()
alt_id handling
obsolete term removal
relationship parents included when with_rels=True
```

The misspelling `get_anchestors` is preserved because DeepGOPlus calls that
method name.

### `builder.py`

Main execution logic.

Stages:

1. Load GO.
2. Load t0/t1 UniProt sequence maps.
3. Stream t0/t1 GOA annotations.
4. Build DeepGOPlus-style `train_data.pkl` table.
5. Build DeepGOPlus-style `terms.pkl` using propagated training annotations and
   `min_count`.
6. Split train/valid using DeepGOPlus `np.random.seed(0)` and 90/10 split.
7. Build temporal test set from gained annotations.
8. Export TEMPROT-style ontology-specific wide CSVs.

### `cli.py`

Command-line interface. All major paths are configurable:

```text
--uniprot-t0
--uniprot-t1
--goa-t0
--goa-t1
--go-obo
--output-dir
--target-taxon
--target-taxa-file
--min-count
--split
--seed
--reviewed-only
--max-gaf-records
```

## Faithfulness To Historical Sources

### Directly Mirrored

DeepGOPlus `cafa3_data.py`:

- `Ontology(..., with_rels=True)`
- propagated annotations via `go.get_anchestors`
- `train_data.pkl`, `test_data.pkl`, `terms.pkl`
- term frequency threshold `min_count >= 50`

DeepGOPlus `deepgoplus.py`:

- `np.random.seed(seed=0)`
- `np.random.shuffle(index)`
- `split=0.9`

TEMPROT `dataset.py`:

- ontology-specific term selection
- `proteins`, `sequences`, binary GO columns
- ancestor propagation during CSV export
- duplicate sequence removal across splits

PFP `prepare_cafa3_data.py`:

- exact final filenames use `training`, `validation`, `test`
- GO columns are columns beginning with `GO:`

### Intentional Changes

1. The historical scripts assume fixed paths. This implementation uses a CLI.
2. Historical DeepGOPlus reads prepared CAFA files, not raw modern GOA. This
   implementation streams raw GOA because that is required for 2025->2026.
3. TEMPROT writes `*-testing.csv`; PFP expects `*-test.csv`, so this builder
   writes PFP-compatible filenames directly.
4. TEMPROT uses repeated `list_terms.index(term)` lookups. This builder uses a
   dictionary index for equivalent output with better scaling.

## Important Methodological Caveat

The implemented default test policy follows the dissertation/supervisor
temporal-holdout requirement:

```text
test proteins = proteins that gained experimental annotations between t0 and t1
and are not present in the t0 training set
```

Historical CAFA also distinguishes NK and LK targets. The public
`create_benchmark.py` implements NK/LK discovery, but PFP's wide CSV interface
does not require preserving the NK/LK category labels. The builder leaves room
to add benchmark variants later.

## Missing Historical Details

The official CAFA3 benchmark README says:

- evidence codes include `TAS` and `IC`;
- CGD backfilled annotations assigned before t0 were removed;
- MFO benchmark proteins whose only MFO annotation was `GO:0005515` were
  removed.

The public Python `create_benchmark.py` does not implement all of those final
rules. Older CAFA2 MATLAB code implements related protein-binding removal, but
the exact CAFA3 postprocessing script was not found locally.

Therefore this implementation is faithful to the recoverable execution graph,
but the CGD/protein-binding special cases remain documented methodology gaps.

## Tests Run

Command:

```bash
python3 -m unittest discover -s tests
```

Result:

```text
Ran 2 tests
OK
```

The tests verify:

- UniProt FASTA parsing.
- GOA GAF parsing.
- final CAFA3 `TAS` and `IC` evidence kept.
- `NOT` annotations removed.
- gained t1 annotation creates a test protein.
- nine PFP-compatible CSV filenames are written.

## Next Recommended Steps

1. Add a real CAFA3 target-taxon file.
2. Run a tiny real-data smoke test with `--max-gaf-records`.
3. Verify output CSV headers and counts.
4. Run a larger t0-only parser count check.
5. Only then run the full 21 GB GOA build on HPC.
