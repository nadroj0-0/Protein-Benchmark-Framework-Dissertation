# Contemporary CAFA Benchmark Builder

This package constructs a deterministic temporal protein-function benchmark and
exports the nine wide CSV files consumed by immutable PFP. It also retains the
validated historical DeepGOPlus and CAFA-file conversion modes.

## Production contract

For the snapshot mode, a test protein must:

1. exist in the selected target population at `t0`;
2. map unambiguously to a `t1` UniProt record through primary or secondary
   accessions;
3. satisfy the profile's knowledge policy at `t0`;
4. have a qualifying annotation at `t1` dated after the configured `t0` cutoff;
5. use its `t0` sequence in every exported test file;
6. survive the sequence-change, ontology, protein-binding and training-term
   policies.

Proteins first appearing at `t1` are never test candidates. The CAFA-style
profiles apply eligibility per ontology: a protein may be an LK/type2 benchmark
in an ontology where it had no qualifying `t0` annotation even if another
ontology was already annotated. The supervisor profile uses the stricter global
rule from the dissertation specification, so any qualifying `t0` annotation
places the protein in training and excludes it from every test ontology.

## Named profiles

| Profile | Training taxa | Target taxa | Evidence | Test knowledge rule | Training records |
|---|---|---|---|---|---|
| `cafa3-reconstructed` | all taxa | CAFA3 targets | final CAFA3 eight-code set | ontology-specific NK/LK | reviewed |
| `contemporary-cafa3-style` | all taxa | CAFA3 targets | final CAFA3 eight-code set | ontology-specific NK/LK | reviewed |
| `supervisor` | CAFA3 targets | CAFA3 targets | supervisor-specified set | globally unannotated at t0 | reviewed + unreviewed |

All profiles include reviewed and unreviewed proteins in the target population
when those UniProt records are supplied. This preserves the difference between
the broad Swiss-Prot CAFA-style training population and the target-organism
candidate population. CLI switches can override each policy without editing
source code.

The final CAFA3 evidence set is:

```text
EXP IDA IPI IMP IGI IEP TAS IC
```

The supervisor set is:

```text
EXP IDA IPI IMP IGI IEP HTP HDA HMP HGI HEP TAS NAS IGC RCA ND IC
```

## Frozen inputs

The production 2025 to 2026 run expects:

```text
UniProt 2025_01
UniProt 2026_02
GOA UniProt release 225
GOA UniProt release 234
GO product release 2025-02-06 (frozen prediction graph)
GO product release 2025-03-16 (t0 source-ID resolution)
GO product release 2026-06-19, data-version 2026-06-15 (t1 source-ID resolution)
```

The GAF date window is bounded explicitly: rows assigned on or before
`2025-03-08` are not annotation gains, and rows after the frozen GOA 234
endpoint (`2026-06-17`) are rejected. The historical CAFA3 profile uses the
documented `2017-02-13` to `2017-11-15` growth period.

GOA 225 declares GO version 2025-03-07, but the public GO product archive does
not retain a standalone 2025-03-07 release directory. The builder therefore
freezes the prediction graph to the last preceding public release (2025-02-06)
and uses the first following public release (2025-03-16) only to resolve source
IDs. GOA IDs are canonicalised in the source graph, including primary IDs,
`alt_id` values and unambiguous `replaced_by` values, then mapped into the frozen
graph. Terms absent from the frozen graph are reported and excluded.

The nearest retained t0 source product (2025-03-16) postdates the frozen graph.
If a raw GOA 225 ID is obsolete or absent there but the exact raw ID is live in
the frozen 2025-02-06 graph, the default policy classifies this as a source
snapshot mismatch and uses that frozen term. This is not a `consider` mapping.
`consider` suggestions are never selected automatically. Rows absent from both
graphs remain unresolved and fail strict QC. Use
`--no-frozen-source-fallback` to disable this narrow compatibility rule.

To represent all proteins from the target organisms, supply both Swiss-Prot and
TrEMBL DAT files. A Swiss-Prot-only build is a diagnostic variant, not the full
supervisor benchmark.

## Pipeline

```text
UniProt DAT/FASTA t0 + t1
        |
        +-- canonical records and secondary-accession crosswalk
        |
GOA GAF t0 + t1 -- evidence / NOT / date filtering
        |
GO OBO t0 + t1 -- GO-ID canonicalisation into frozen t0 graph
        |
        +-- t0 qualifying annotations -> training population
        +-- profile-eligible targets with post-t0 annotations -> test
        |
DeepGOPlus-shaped train_data.pkl / test_data.pkl / terms.pkl
        |
np.random.seed(0), np.random.shuffle, 90/10 train/validation
        |
TEMPROT-shaped ontology separation and exact-sequence de-duplication
        |
bp/cc/mf x training/validation/test CSVs for PFP
```

The training-defined term universe uses propagated `t0` annotations and
`min_count >= 50`. Snapshot rows and GO columns are sorted before the seeded
split so Python hash randomisation cannot change the outputs.

## Direct CLI use

Install the builder and its compatible dependencies for direct CLI use:

```bash
python -m pip install -e benchmark_builders/contemporary_cafa
```

Then run with paths supplied explicitly:

```bash
python -m cafa_benchmark_builder \
  --profile contemporary-cafa3-style \
  --uniprot-t0 /path/2025_01/uniprot_sprot.dat.gz \
  --uniprot-t0 /path/2025_01/uniprot_trembl.dat.gz \
  --uniprot-t1 /path/2026_02/uniprot_sprot.dat.gz \
  --uniprot-t1 /path/2026_02/uniprot_trembl.dat.gz \
  --goa-t0 /path/goa_uniprot_all.gaf.225.gz \
  --goa-t1 /path/goa_uniprot_all.gaf.234.gz \
  --go-obo /path/go-2025-02-06/go-basic.obo \
  --go-obo-t0 /path/go-2025-03-16/go-basic.obo \
  --go-obo-t1 /path/go-2026-06-19/go-basic.obo \
  --output-dir /path/run/outputs \
  --report-dir /path/run/reports
```

Paths are never hard-coded in Python. `--uniprot-t0` and `--uniprot-t1` are
repeatable, allowing Swiss-Prot and TrEMBL inputs to be combined.

The repository-level runner resolves the dissertation database layout and calls
the same CLI. Existing inputs are reused; missing inputs are downloaded into
`DB_ROOT`:

```bash
DB_ROOT="$HOME/protein_databases" \
PROFILE=contemporary-cafa3-style \
bash scripts/benchmark_generation/run_contemporary_temporal_benchmark.sh
```

For a deliberately incomplete Swiss-Prot-only diagnostic:

```bash
ALLOW_SPROT_ONLY=1 bash scripts/benchmark_generation/run_contemporary_temporal_benchmark.sh
```

The historical UniProt 2025_01 TrEMBL data are available only inside the 181 GB
`knowledgebase2025_01.tar.gz` archive, while the current 2026_02 TrEMBL DAT is
about 110 GB compressed. The runner streams these sources through
`filter_uniprot_dat.py` and stores only records from the CAFA3 taxa. It does not
retain either full TrEMBL source unless a persistent local copy was explicitly
provided.

## Cluster run

From the cloned framework on Morecambe:

```bash
qsub hpc_jobs/active/hpc_contemporary_temporal_benchmark.sh
```

Choose the supervisor profile at submission time:

```bash
qsub -v PROFILE=supervisor hpc_jobs/active/hpc_contemporary_temporal_benchmark.sh
```

The wrapper requests 64 GB RAM and 200 GB scratch. It copies any available
frozen inputs from `${PROTEIN_DATABASE_ROOT:-$HOME/protein_databases}`, downloads
anything missing directly into scratch, activates and uses `mmfp`, invokes the
normal runner, and copies outputs/reports/logs to
`$HOME/contemporary_cafa_benchmark_results`, and removes its scratch directory
on success, error, interruption or `qdel` termination.

To use a different persistent database root:

```bash
qsub -v PROTEIN_DATABASE_ROOT=/path/to/databases,PROFILE=supervisor \
  hpc_jobs/active/hpc_contemporary_temporal_benchmark.sh
```

Advanced callers may supply colon-separated `UNIPROT_T0_INPUTS` and
`UNIPROT_T1_INPUTS`; the wrapper copies those files to scratch before running.

The `.sh` suffix is intentional. `qsub` reads the `#$` directives inside a shell
script and does not require a `.qsub` filename extension.

## Outputs

PFP-compatible CSVs:

```text
bp-training.csv      bp-validation.csv      bp-test.csv
cc-training.csv      cc-validation.csv      cc-test.csv
mf-training.csv      mf-validation.csv      mf-test.csv
```

DeepGOPlus-shaped intermediates:

```text
train_data.pkl
train_data_train.pkl
train_data_valid.pkl
test_data.pkl
terms.pkl
```

Audit material:

```text
benchmark_build_report.md
benchmark_statistics.json
build_manifest.json
protein_flow.tsv
exclusion_reasons.tsv
annotation_gain_summary.tsv
taxon_summary.tsv
evidence_summary.tsv
input_checksums.sha256
output_checksums.sha256
input_acquisition.tsv
unresolved_source_go_annotations.tsv
outside_frozen_go_annotations.tsv
```

The GO diagnostic TSVs are written before strict QC raises, so the exact
protein, GO ID, evidence, date, source row, ontology metadata, classification
and action survive a failed production job. Valid t1 terms outside the frozen
t0 graph are expected ontology drift: they are reported and excluded from the
label space without failing the build.

## Failure gates

A production build exits non-zero when it finds:

- an input file missing;
- an invalid GAF date where backfill filtering is required;
- an unresolvable GO ID in its matching source ontology;
- a test protein absent at `t0`;
- a qualifying `t0` annotation on a selected test protein;
- train/test protein overlap;
- malformed CSV schemas, duplicate IDs or non-binary labels;
- exact sequence overlap between exported splits;
- an empty ontology/split under strict QC.

Expected exclusions such as a valid post-`t0` GO term outside the frozen graph,
an accession merge, a sequence change, or a `t1`-only protein are counted in the
reports rather than silently discarded.

## Tests

```bash
cd benchmark_builders/contemporary_cafa
PYTHONDONTWRITEBYTECODE=1 python -m unittest discover -s tests -v
```

The tests cover evidence/NOT filtering, backfill removal, `t0`-presence rules,
`t0` sequence use, secondary-accession mapping, ambiguous merges, sequence
changes, GO alternate/replacement IDs, all nine PFP CSVs, historical modes and
byte-identical repeated fixture outputs.

## Historical modes

These validated paths remain available and are intentionally separate from the
new snapshot front end:

```text
--source-mode deepgoplus
--source-mode cafa3-files
```

They reproduce the released pickle-to-CSV and official-file-to-pickle stages.
They do not apply the contemporary temporal candidate policy.

## Raw CAFA3 historical validation

The separate repository-level historical runner supports two independent
controls:

```text
HISTORICAL_TRAINING_SNAPSHOT=september-2016 | february-2017-legacy
TARGET_UNIVERSE_POLICY=official-cafa3-targets | reconstructed-all-qualifying
```

The primary validation defaults to `september-2016` plus
`official-cafa3-targets`. It uses UniProtKB release `2016_08` (07-Sep-2016), the
last public monthly Swiss-Prot release before the official CAFA3 training
package date of 24-Sep-2016. The runner verifies the archive against UniProt's
published byte size and MD5 and records a SHA-256 inventory. This is the closest
defensible public snapshot, not a claim that the organiser's private freeze is
bit-identical.

The released `data-cafa.tar.gz` supplies the official training annotation table,
aggregate CAFA3 target FASTA, target mapping files and DeepGOPlus reference
pickles. Official-target mode preserves released CAFA IDs and exact FASTA
sequences. UniProt mappings are conservative and fully reported; unmapped,
ambiguous and custom-source targets remain in the candidate catalogue. The
February-2017/reconstructed combination remains a named legacy comparison.

```bash
HISTORICAL_TRAINING_SNAPSHOT=september-2016 \
TARGET_UNIVERSE_POLICY=official-cafa3-targets \
bash scripts/validation/run_cafa3_historical_validation.sh
```

Large inputs are acquired into scratch when no override is supplied. Optional
local overrides are `HISTORICAL_TRAINING_UNIPROT_ARCHIVE` and
`OFFICIAL_CAFA3_ARCHIVE_INPUT`.
