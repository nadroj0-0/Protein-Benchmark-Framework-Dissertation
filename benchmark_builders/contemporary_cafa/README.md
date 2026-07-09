# Contemporary CAFA Benchmark Builder

This is a small, separate codebase for generating 2025->2026 CAFA-style
benchmark CSVs that can be consumed by immutable PFP via
`scripts/prepare_cafa3_data.py`.

It recreates the historical execution chain we reverse engineered:

```text
Raw UniProt / GOA / GO ontology
    -> CAFA-style temporal annotation comparison
    -> DeepGOPlus-style proteins/sequences/annotations tables
    -> DeepGOPlus min_count term filtering
    -> DeepGOPlus 90/10 train/valid split
    -> TEMPROT-style ontology-specific wide CSV export
    -> PFP-compatible bp/cc/mf training/validation/test CSVs
```

PFP itself is not modified.

## Provenance

This builder intentionally follows the historical sources as closely as
possible:

- DeepGOPlus `cafa3_data.py`
  - loads GO with `Ontology(..., with_rels=True)`
  - propagates annotations with `go.get_anchestors(go_id)`
  - writes `train_data.pkl`, `test_data.pkl`, `terms.pkl`
  - filters terms by `min_count >= 50`
- DeepGOPlus `deepgoplus.py`
  - splits train/valid using `np.random.seed(seed=0)`,
    `np.random.shuffle(index)`, and `split=0.9`
- TEMPROT `src/<ontology>/dataset.py`
  - reads `train_data_train.pkl`, `train_data_valid.pkl`, `test_data.pkl`,
    `terms.pkl`
  - keeps ontology-specific terms
  - writes `proteins`, `sequences`, then one binary column per GO term
  - removes duplicated sequences across train/test/valid
- PFP `scripts/prepare_cafa3_data.py`
  - consumes `bp-training.csv`, `bp-validation.csv`, `bp-test.csv`
  - same for `cc` and `mf`

The final CAFA3 evidence policy is taken from the official CAFA3 benchmark
README:

```text
EXP, IDA, IPI, IMP, IGI, IEP, TAS, IC
```

This deliberately differs from the public `CAFA_benchmark/create_benchmark.py`,
which only keeps:

```text
EXP, IDA, IPI, IMP, IGI, IEP
```

The README is treated as the final official policy.

## What It Produces

The final output directory contains:

```text
bp-training.csv
bp-validation.csv
bp-test.csv

cc-training.csv
cc-validation.csv
cc-test.csv

mf-training.csv
mf-validation.csv
mf-test.csv
```

Each CSV has the PFP/TEMPROT shape:

```text
proteins,sequences,GO:...,GO:...,GO:...
P12345,MSEQ...,1,0,1
```

If intermediates are enabled, it also writes DeepGOPlus-style pickles:

```text
train_data.pkl
train_data_train.pkl
train_data_valid.pkl
test_data.pkl
terms.pkl
```

## Installation

From this directory:

```bash
python3 -m venv .venv
. .venv/bin/activate
pip install -e .
```

The code only depends on NumPy and pandas.

## Smoke Tests

Run the lightweight tests:

```bash
python3 -m unittest discover -s tests
```

The tests use tiny local fixtures only. They do not touch real GOA or UniProt
files.

## Example Command

The builder has three source modes.

### Historical CAFA3 File-to-Pickle Validation Mode

Use this mode to validate the recovered DeepGOPlus `cafa3_data.py` layer:

```bash
python3 -m cafa_benchmark_builder \
  --source-mode cafa3-files \
  --go-obo "/Users/jordansydney-darlin/CAFA3 Supplementary material/external_repos/TEMPROT/data-cafa/go.obo" \
  --train-sequences-file "/Users/jordansydney-darlin/CAFA3 Supplementary material/external_repos/TEMPROT/data-cafa/CAFA3_training_data/uniprot_sprot_exp.fasta" \
  --train-annotations-file "/Users/jordansydney-darlin/CAFA3 Supplementary material/external_repos/TEMPROT/data-cafa/CAFA3_training_data/uniprot_sprot_exp.txt" \
  --test-sequences-file "/Users/jordansydney-darlin/CAFA3 Supplementary material/external_repos/TEMPROT/data-cafa/CAFA3_targets/targets_all.fasta" \
  --test-annotations-file "/Users/jordansydney-darlin/CAFA3 Supplementary material/external_repos/TEMPROT/data-cafa/benchmark20171115/groundtruth/leafonly_all.txt" \
  --output-dir /tmp/cafa3_deepgoplus_pickles \
  --min-count 50
```

This writes:

```text
train_data.pkl
test_data.pkl
terms.pkl
```

It does not create the 90/10 train/validation split or the nine PFP CSVs. Those
belong to the next historical layer.

### Historical DeepGOPlus Validation Mode

Use this mode to validate the recovered CAFA3 -> DeepGOPlus -> TEMPROT -> PFP
path without reparsing raw GOA:

```bash
python3 -m cafa_benchmark_builder \
  --source-mode deepgoplus \
  --deepgoplus-dir "/Users/jordansydney-darlin/CAFA3 Supplementary material/external_repos/TEMPROT/data-cafa" \
  --go-obo "/Users/jordansydney-darlin/CAFA3 Supplementary material/external_repos/TEMPROT/data-cafa/go.obo" \
  --output-dir /tmp/cafa3_deepgoplus_export
```

This reads:

```text
train_data.pkl
train_data_train.pkl
train_data_valid.pkl
test_data.pkl
terms.pkl
```

and exports the nine PFP-compatible CSVs. This is the best mode for proving the
final CSV export layer is faithful to the historical benchmark.

### Raw Snapshot Mode

Example structure for the real contemporary benchmark:

```bash
python3 -m cafa_benchmark_builder \
  --source-mode snapshots \
  --uniprot-t0 /home/jsydneyd/protein_databases/uniprot/release_2025_01/uniprot_sprot.dat.gz \
  --uniprot-t1 /home/jsydneyd/protein_databases/uniprot/release_2026_02/uniprot_sprot.dat.gz \
  --goa-t0 /home/jsydneyd/protein_databases/goa/release_2025_01/goa_uniprot_all.gaf.225.gz \
  --goa-t1 /home/jsydneyd/protein_databases/goa/release_2026_02/goa_uniprot_all.gaf.234.gz \
  --go-obo /home/jsydneyd/protein_databases/go/2026_06_15/go.obo \
  --reviewed-only \
  --output-dir /home/jsydneyd/contemporary_cafa_2025_2026/csvs
```

By default this uses the official CAFA3-style broad training scope: all loaded
Swiss-Prot proteins with qualifying annotations. To switch to your supervisor's
"same model organisms as CAFA3" variant, add:

```bash
--taxon-policy cafa3-targets
```

To use a custom taxon list instead:

```bash
--taxon-policy custom --target-taxa-file /path/to/taxa.txt
```

Evidence policies are also named presets:

```bash
--evidence-policy cafa3-final
--evidence-policy cafa3-public-python
--evidence-policy supervisor
```

For strict CAFA3 validation use `cafa3-final`. For the dissertation benchmark,
switch to `supervisor` if that is the agreed policy.

For a tiny smoke run on real files, use:

```bash
python3 -m cafa_benchmark_builder \
  ...same paths... \
  --min-count 1 \
  --max-gaf-records 10000 \
  --output-dir /tmp/cafa_smoke
```

Do not use `--max-gaf-records` for the real benchmark.

## Method Implemented

1. Load GO from OBO.
2. Stream UniProt t0 and t1 sequences from FASTA or DAT.
3. Stream GOA t0 and t1 GAF rows, filtering early to the loaded UniProt
   accession universe so irrelevant GOA rows are skipped before expensive
   object construction.
4. Keep only:
   - `DB == UniProtKB`
   - evidence in final CAFA3 policy
   - aspect in `P`, `C`, `F`
   - optionally configured target taxa
   - no `NOT` qualifier
5. Training annotations are experimental t0 annotations for proteins with t0
   sequences.
6. Test annotations are gained t1 annotations for proteins not present in the
   t0 training set.
7. Propagate annotations through GO using DeepGOPlus-style ancestor expansion.
8. Count propagated training annotations.
9. Keep GO terms with `min_count >= 50`.
10. Split train/valid with DeepGOPlus `seed=0`, `split=0.9`.
11. Export ontology-specific TEMPROT/PFP CSVs.

## Known Caveats

The official CAFA3 README documents two postprocessing rules that are not fully
recoverable from the public Python code alone:

1. CGD backfill removal for Candida annotations assigned before t0.
2. Removal of MFO benchmark proteins whose only MFO annotation was
   `GO:0005515`.

Older CAFA2 MATLAB code contains protein-binding removal logic, but the exact
CAFA3 Python execution path for those final postprocessing rules is not present
locally. This builder therefore records the gap rather than pretending it is
fully solved.

This temporal benchmark also follows TEMPROT's exact-sequence duplicate removal
only. It is temporally split, but it is not a homology-decontaminated benchmark:
near-identical homologs are intentionally left to a separate homology-aware
benchmark design.

## PFP Handoff

After generating the nine CSVs, run PFP's existing preparer:

```bash
python scripts/prepare_cafa3_data.py \
  --cafa3-dir /path/to/generated/csvs \
  --output-dir data
```

Then continue with PFP's existing embedding/training workflow.
