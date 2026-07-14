# Homology-cluster benchmark builder

This isolated package builds Daniel's Part B dissertation benchmark from a frozen UniRef90
sequence population. It runs MMseqs2 independently at 30%, 25%, 20%, 15%, 10%, or 5%
identity with 80% longest-sequence coverage, retains clusters connected to qualifying GOA
annotations, assigns whole clusters globally to training/validation/test, and publishes the nine
CSV files consumed by immutable PFP plus DeepGOPlus-shaped pickle intermediates.

It does **not** alter the contemporary temporal benchmark, immutable PFP, or the embedding
inventory. It does not generate embeddings or train/evaluate PFP.

## Scientific data flow

```text
frozen GOA 234 annotation
    │ exact/explicit accession handling; own GO labels only
    ▼
UniProtKB accession ── frozen idmapping_selected column 9 ──► UniRef90 identifier
    │                                                        │
    │ frozen UniProt sequence                                │ frozen UniRef90 FASTA
    │                                                        ▼
    └──────────────────────────────────────────── MMseqs2 cluster at identity X
                                                             │ whole-cluster assignment
                                                             ▼
                                                   development / test
                                                             │
                              own labels + GO propagation + development-defined term universe
                                                             │
                                      development-only 90:10 cluster split
                                                             ▼
                                               training / validation / test
                                                             │
                                         PFP CSVs and DeepGOPlus-shaped pickles
```

These three identities never collapse into one field:

1. UniProtKB accession: GOA key, benchmark protein ID, sequence/embedding key.
2. UniRef90 identifier: one frozen clustering-scaffold entry.
3. MMseqs2 cluster representative: lower-identity cluster key for splitting.

Cluster membership determines retention and split only. It never transfers GO labels between
homologues.

## Frozen endpoint

Production builds lock and record:

- UniProt / UniRef release `2026_02`;
- GOA release `234`, whose established local header is generated `2026-06-17` against GO
  `2026-06-15`;
- ontology `data-version: releases/2026-06-15`, distributed from the established
  `2026-06-19` GO release directory.

Every input is declared as a local path, URL, or both. Local paths always win and are never
downloaded again, while a supplied source URL remains in provenance. Resolved path, URL, release
label, byte size, expected hash, observed SHA-256, and acquisition action are recorded. A
production build requires both an expected SHA-256 for all five inputs and a reviewed
`--frozen-input-manifest` with exactly five release/URL/filename/size/hash/acquisition/embedded-
metadata/notes entries. Fixture mode may generate a clearly synthetic manifest. The manifest
separates byte reproducibility, recorded provenance, and an explicit authoritative-origin review;
a caller computing hashes over arbitrary local bytes cannot bypass the missing review contract.
This is an auditable attestation, not a cryptographic signature or independent proof of origin.

## Methodological locks

- Identities: exactly `30, 25, 20, 15, 10, 5` percent. Other values fail.
- Coverage: exactly `0.8` with MMseqs2 `--cov-mode 0`.
- MMseqs2, not CD-HIT/DIAMOND or post-hoc homologue removal.
- Frozen UniRef90 FASTA, not a CAFA or temporal benchmark population.
- Daniel's exact 17-code set:
  `EXP, IDA, IPI, IMP, IGI, IEP, HTP, HDA, HMP, HGI, HEP, TAS, NAS, IGC, RCA, ND, IC`.
- A locked global 80:20 development/test cluster split, followed by a locked 90:10
  training/validation cluster split inside development.
- One global split before BP/CC/MF separation.
- Production `min_count >= 50`, derived from the complete development population before its
  90:10 cluster split for each threshold;
  lower fixture values require explicit `--fixture-mode`.
- Ontology roots `GO:0008150`, `GO:0005575`, and `GO:0003674` are retained in propagated labels,
  `terms.pkl`, and PFP CSV columns when they meet development support, matching PFP/TEMPROT shape.
- `annotated-only` is the complete production population.
- `all-cluster-members` is represented but hard-fails before input checks or clustering because absence from GOA
  is not a negative label and no label-transfer policy was authorized.

## Exact MMseqs2 workflow

For identity `X`, the builder constructs argument arrays equivalent to:

```bash
mmseqs createdb FROZEN_UNIREF90_FASTA WORK/uniref90_db \
  --dbtype 1 --shuffle 0

mmseqs cluster WORK/uniref90_db WORK/clusters_db WORK/mmseqs_tmp \
  --min-seq-id X \
  -c 0.8 \
  --cov-mode 0 \
  --cluster-mode 0 \
  --alignment-mode 3 \
  --seq-id-mode 0 \
  --cluster-reassign 1 \
  -s 7.5 \
  -e 1e-4 \
  --threads THREADS

mmseqs createtsv WORK/uniref90_db WORK/uniref90_db \
  WORK/clusters_db WORK/uniref90_clusters.tsv
```

`--alignment-mode 3 --seq-id-mode 0` makes the identity fraction literal identical aligned
residues divided by alignment columns, including internal gaps. `--shuffle 0` preserves the
frozen FASTA order in `createdb`. `-c 0.8 --cov-mode 0` requires aligned residues divided
by the longer of query/target lengths to be at least 0.8, matching the documented UniRef
"80% overlap with the longest seed" interpretation. Cluster mode 0 is greedy set cover; it is not
connected components or an all-pairs guarantee. `-s 7.5` is prefilter sensitivity, not identity.
The reviewed implementation locks `-e 1e-4`; whether that value remains scientifically preferable
at 5–15% identity is still a documented future sensitivity decision, not a runtime override.
Production requires a non-placeholder exact
`--expected-mmseqs-version`, compares it to a parseable successful `mmseq version` probe, resolves
the executable once, records its path/SHA-256 when readable, and rechecks the executable before
publication. Fixture precomputed assignments are explicitly relaxed because MMseqs is not run. Successful MMseqs logs are
published; failed-run logs and bounded semantic failure reports are retained outside any completed
run. See [METHODOLOGY.md](METHODOLOGY.md).

## GOA, accession, and ontology policy

The GAF parser streams plain or gzip text, requires `!gaf-version: 2.2`, requires exactly 17 data
columns and every mandatory GAF 2.2 field, and records the embedded GAF/date/GO metadata. It
retains only:

- database `UniProtKB`;
- object type `protein`;
- one of Daniel's exact evidence codes;
- aspect `P`, `C`, or `F` consistent with the resolved term namespace;
- annotations without an exact pipe-delimited `NOT` token.

Reference, With/From, assigned date (field 14), Assigned By, annotation extension, and gene-product
form are preserved in annotation manifests; assigned date is never used to reconstruct snapshot
membership. Field-17 form-specific rows and accessions with an explicit `-N` isoform suffix are
excluded unless a future independently verified isoform mapping policy is added. A suffix is never
stripped. UniProt DAT secondary accessions
can be explicitly collapsed to their primary accession; FASTA supplies exact accessions only.
Ambiguous aliases, mapping conflicts, blank/absent mappings, mappings missing from the frozen
FASTA, and missing canonical sequences remain separate statuses.

GO primary IDs and alt IDs resolve to the live canonical ID. Obsolete IDs resolve only through a
single `replaced_by`; `consider` is never selected. Unknown, unresolved obsolete, and
namespace-mismatched terms are excluded with reason. The established DeepGOPlus-style ancestor
policy includes OBO relationship targets when enabled; the exact relationship types observed are
recorded. Daniel did not specify GO-root handling. Retaining the three ontology roots is inherited
from the established DeepGOPlus/TEMPROT/PFP interface, and they remain after propagation.

## Split algorithms

Both algorithms sort IDs before any seeded operation and assign indivisible MMseqs2 clusters.

### `cluster-count-random`

1. Sort cluster IDs.
2. Shuffle using Python's isolated `random.Random(seed)`.
3. Allocate `round(cluster_count * 0.8)` clusters to development, bounded so both sides are nonempty.
4. Repeat on development with `seed + 1` and fraction `0.9` for training vs validation.

### `sequence-balanced`

1. Weight each retained cluster by its frozen UniRef90 member count—not PFP protein count.
2. Construct the smaller side (20% test, then 10% validation) to avoid systematic large-side
   overshoot, using four fixed deterministic candidate orders.
3. Refine candidates with strict-improvement single moves, bisect-assisted one-for-one exchanges,
   and a fixed-size two-move pool; passes and pools are bounded for production scale.
4. Select by absolute target deviation plus a SHA-256 membership tie-break, enforce nonempty sides,
   and repeat on development with `seed + 1`.

Reports include requested/achieved ratios by cluster count, UniRef90 member count, and qualifying
UniProt protein count, plus the largest indivisible cluster. Sequence-balanced deviation above two
percentage points warns; production fails if either stage exceeds five points. Fixture runs report
otherwise-impossible giant-cluster cases without pretending they are production-valid. The bounded
heuristic does not claim globally optimal subset-sum allocation.

## Label and PFP contracts

Each protein receives only its own canonicalized GO annotations. Propagation occurs before term
support counting and retains all three roots. The complete development population defines
`terms.pkl` at `min_count` before its 90:10 cluster partition; test never contributes. Every split's
PFP labels are intersected with that universe. Each threshold has
an independent term universe. A six-threshold invocation—or the later `summarize` command over six
independent HPC outputs—publishes pairwise overlap, gained/lost terms, and identical-clustering
warnings.

The nine CSVs are:

```text
bp-training.csv    bp-validation.csv    bp-test.csv
cc-training.csv    cc-validation.csv    cc-test.csv
mf-training.csv    mf-validation.csv    mf-test.csv
```

Every file begins `proteins,sequences`, followed by lexicographically ordered ontology-specific
GO IDs from that threshold's development-defined term universe. The same header/order is used across all three
splits. Rows are sorted by UniProtKB accession and labels are binary. Proteins without an evaluable
term in an ontology are absent from that ontology's CSV and counted; they are never emitted as
all-zero negatives.

DeepGOPlus-shaped intermediates preserve each protein's complete propagated annotation
tuple, matching the established historical chain. `terms.pkl` defines the evaluable universe and
the nine CSVs are the explicitly restricted PFP views. This retains out-of-universe labels for
audit without exposing them as PFP target columns.

The intermediates are:

- `train_data.pkl`: sorted training + validation supervised rows;
- `train_data_train.pkl`, `train_data_valid.pkl`, `test_data.pkl`: dataframes with exactly
  `proteins`, `sequences`, `annotations`;
- `terms.pkl`: one-column dataframe named `terms`.

## CLI

Run from the package root or set `PYTHONPATH=benchmark_builders/homology_cluster/src`.

Command preview (no path resolution, download, MMseqs2 execution, or publication):

```bash
python3 -m homology_cluster_benchmark build \
  --identity 30 \
  --split-policy sequence-balanced \
  --training-population annotated-only \
  --uniref90-fasta /frozen/uniref90.fasta.gz \
  --idmapping /frozen/idmapping_selected.tab.gz \
  --uniprot-sequences /frozen/uniprot.fasta.gz \
  --goa /frozen/goa_uniprot_all.gaf.234.gz \
  --go-obo /frozen/go-basic.obo \
  --frozen-input-manifest /frozen/homology-inputs.reviewed.json \
  --mmseqs-bin /path/to/mmseqs \
  --expected-mmseqs-version 15-REPLACE_WITH_EXACT_BUILD \
  --output-dir /results/homology \
  --temp-dir /scratch/homology \
  --threads 16 \
  --seed 0 \
  --dry-run
```

Remove `--dry-run` for a real future build. Use `--identity all` only when all six independent runs
are intentionally requested. The production command refuses existing final run directories; an
unsafe resume/overwrite mode is not implemented.

Production execution additionally requires all five `--*-sha256` values, the reviewed frozen-input
manifest, an exact MMseqs version token, and a clean,
commit-addressable framework checkout. `--identity all` requires the five inputs to be staged as
local paths so multi-gigabyte sources are not downloaded once per threshold. Precomputed cluster
assignments, `min_count < 50`, disabled strict QC, and empty ontology/split outputs require explicit
`--fixture-mode`. Fixture artifacts declare `production_eligible: false` and
`benchmark_scope: fixture-only` in their summaries and completion markers.

Each run publishes beneath:

```text
OUTPUT_ROOT/identity_30/sequence-balanced/annotated-only/seed_0/min_count_50/
```

Publication validation checks marker/manifest hashes and the saved strict validation result:

```bash
python3 -m homology_cluster_benchmark validate \
  --run-dir OUTPUT_ROOT/identity_30/sequence-balanced/annotated-only/seed_0/min_count_50
```

After six independently scheduled runs, validate their shared frozen-input, release, evidence,
split/population/seed/min-count, root/relationship, builder/git, and MMseqs2 scientific fingerprint
and publish scoped cross-threshold reports. Python/pandas/NumPy versions remain provenance fields;
the aggregator does not overclaim that they are part of this scientific fingerprint:

```bash
python3 -m homology_cluster_benchmark summarize \
  --output-dir AGGREGATE_OUTPUT_ROOT \
  --run-dir /results/job-30 \
  --run-dir /results/job-25 \
  --run-dir /results/job-20 \
  --run-dir /results/job-15 \
  --run-dir /results/job-10 \
  --run-dir /results/job-05
```

Each `--run-dir` may be a publication leaf or a timestamped parent containing exactly one leaf.
Identity and shared selectors come only from hashed child publication metadata. The aggregate stores
small child hash/path references; it does not copy CSVs, pickles, or membership manifests.

### URL and hash behavior

For each input, use a local option, `--OPTION-url URL`, or both, and supply
`--OPTION-sha256 HEX` for production. A missing local path never silently falls back to an
unrelated current file.
Use `--no-downloads` to require all paths locally. The known stable ontology URL may be supplied;
exact archived 2026_02 UniRef/idmapping and GOA 234 URLs remain configurable rather than invented.

## Output schemas

Large complete audits use deterministic gzip (`mtime=0`, no embedded filename). Rejection decisions
are complete counts; detailed rejected rows are a deterministic first-N sample per reason controlled
by `--excluded-sample-per-reason` (default 1000). Accepted annotations and retained members remain
complete compressed streams.

| File | Row key and fields |
|---|---|
| `frozen_input_manifest.json` | exact reviewed five-source declaration copied byte-for-byte (synthetic and non-authoritative in fixture mode) |
| `input_manifest.json` | resolved inputs, byte/provenance/origin eligibility dimensions, frozen-manifest fingerprint/hash, embedded metadata, capacity observations |
| `publication_metadata.json` | hashed fixture/eligibility/scope/identity/policy/population/seed/min-count/input/MMseq/repository contract mirrored exactly by the marker |
| `parameters.json` | every split, label, MMseqs2, seed, evidence, release, and root-policy parameter |
| `disk_preflight.json` | scratch/publication free bytes, declared input bytes, safety multiplier, and required-space estimates |
| `run_provenance.json` | command, git state, Python/NumPy/pandas/MMseqs2/OS versions, inputs, determinism boundary |
| `mmseqs_commands.tsv` | stage, argument index, literal argument, shell-quoted display command |
| `uniprot_to_uniref90.tsv` | raw/canonical UniProt accession, accession and lifecycle action/status, UniRef90 ID, mapping detail, FASTA/sequence availability |
| `mmseqs_cluster_membership.tsv.gz` | the one complete canonical MMseqs representative/member publication; raw `createtsv` remains scratch-only |
| `protein_cluster_assignments.tsv` | complete UniProt→UniRef90→MMseqs2→split chain and mapping diagnostics |
| `cluster_split_assignments.tsv` | MMseqs2 cluster, split, UniRef90 weight, qualifying-protein count, assignment stage |
| `retained_clusters.tsv` | retained cluster, split, counts, singleton/giant flags |
| `retained_cluster_members.tsv.gz` | complete retained cluster/split/member, sequence SHA-256/length, qualifying-annotation connection flag |
| `qualifying_annotations.tsv.gz` | complete accepted raw/canonical protein and GO provenance plus mapping chain/split |
| `annotation_decision_counts.tsv` | complete accepted/rejected counts by disposition/reason/evidence/database/object type/aspect |
| `excluded_annotations_sample.tsv.gz` | bounded detailed rejection sample, never a complete rejection table |
| `attrition_summary.json/.tsv` | mutually exclusive row and protein-candidate terminal buckets with explicit denominators |
| `mapping_summary.json` | qualifying/canonical sequence counts, ambiguous aliases, mapping-status counts |
| `evidence_summary.tsv` | evidence, accepted/rejected disposition, reason, row count |
| `go_term_summary.tsv` | term/namespace, direct rows, development-universe membership, unrestricted/evaluable support per split |
| `taxonomy_summary.tsv` | qualifying annotation or supervised-protein stage, split, taxon, count |
| `split_summary.tsv` | split and achieved cluster/member/qualifying, label-intermediate, and evaluable-PFP counts and ratios |
| `cluster_size_summary.tsv` | all/per-split count, total, singleton/giant count, min/median/mean/p95/max |
| `benchmark_summary.json/.md` | population/result totals, label losses, limitations, human-readable summary |
| `validation_report.json/.md` | named pass/fail checks, warnings, metrics, scientific boundary |
| `logs/mmseqs/*` | successful command stdout/stderr; fixture runs contain an explicit `NOT_EXECUTED.json` instead |
| `output_manifest.json` | payload relative path, final byte size, SHA-256, deterministic-payload flag |
| `RUN_COMPLETE.json` | written last; binds output manifest and publication metadata hashes and mirrors all publication-policy fields exactly |

Failure diagnostics live outside successful run directories under
`OUTPUT_ROOT/_failed_runs/<identity>-<uuid>/{FAILURE.json,logs/...}`. They never contain
`RUN_COMPLETE.json` and are not benchmark publications. A six-threshold aggregate separately writes
`term_universe_overlap.tsv`, `term_universe_changes.tsv`, `identical_cluster_assignments.tsv`,
`all_thresholds_summary.json`, its validation reports/manifest, and its own completion marker.

## Atomicity and validation

The builder writes a sibling staging directory, validates the complete staged benchmark, hashes
the payload, renames the payload atomically, re-verifies every hash at the final path, then writes
`RUN_COMPLETE.json` atomically. Existing outputs are never overwritten. A failed run cannot publish
a marker: the owned stage is removed while bounded diagnostics and MMseqs2 logs are preserved under
the output root's `_failed_runs/` directory.

Strict validation covers input/hash/release metadata, MMseqs completion and exact member-assignment
content, the identifier chain, MMseq/UniRef/protein disjointness, and one combined disk-backed exact-
sequence check spanning retained UniRef entries and qualifying UniProt sequences. It also replays the
split, independently recounts the development-defined term universe from serialized training and
validation pickles, verifies exact CSV and pickle contents (not schemas alone), checks binary and
namespace labels, reconciles complete versus sampled audits and the output file set, and verifies
final hashes plus marker/publication-metadata agreement.
Warnings cover giant/imbalanced clusters, poor indivisible balance, mapping and annotation-chain
loss, ontology/taxonomy/evidence imbalance, low evaluable support, and identical threshold outputs.

Internal validation proves software and benchmark-contract properties. It cannot prove that
MMseqs2 discovered every biologically meaningful low-identity edge, that the split is biologically
optimal, that frozen source content is externally correct beyond recorded metadata/hashes, or that
PFP training will succeed on the full population.

## Shell and HPC use

The real shell entrypoint is:

```bash
scripts/benchmark_generation/run_homology_cluster_benchmark.sh
```

It uses `set -euo pipefail`, resolves the repository root from `BASH_SOURCE`, prefers supplied local
inputs, prints the exact command, and delegates all benchmark logic to Python.

The Grid Engine wrapper is:

```bash
qsub -v IDENTITY=30,SPLIT_POLICY=sequence-balanced,TRAINING_POPULATION=annotated-only,\
FRAMEWORK_REVISION=<reviewed-commit>,\
RESULTS_ROOT=/persistent/homology_cluster_benchmark_results,NO_DOWNLOADS=1,\
UNIREF90_FASTA=/persistent/uniref90.fasta.gz,\
UNIREF90_FASTA_SHA256=<64-hex>,\
IDMAPPING=/persistent/idmapping_selected.tab.gz,\
IDMAPPING_SHA256=<64-hex>,\
UNIPROT_SEQUENCES=/persistent/uniprot.fasta.gz,\
UNIPROT_SEQUENCES_SHA256=<64-hex>,\
GOA=/persistent/goa_uniprot_all.gaf.234.gz,\
GOA_SHA256=<64-hex>,\
GO_OBO=/persistent/go-basic.obo,\
GO_OBO_SHA256=<64-hex>,\
FROZEN_INPUT_MANIFEST=/persistent/homology-inputs.reviewed.json,\
MMSEQS_BIN=/path/on/compute/node/mmseqs,\
EXPECTED_MMSEQS_VERSION=<exact-version-token> \
  hpc_jobs/active/hpc_homology_cluster_benchmark.sh
```

The wrapper activates the exact existing `MMFP_ENV_DIR`, verifies the canonical `CONDA_PREFIX`, and
refuses to create/install anything. It clones a clean commit-addressable framework checkout, checks
MMseqs2 on the compute node, stages inputs without changing basenames, and refuses `THREADS` above
allocated `NSLOTS`. Before staging it records local/manifest bytes and checks scratch/MMseqs plus
persistent-publication estimates; before copy it records rounded allocated scratch usage and applies
a configurable copy-safety multiplier plus reserve. Neither value is described as exact.
Success is scratch validate → partial copy → copied validate → atomic final rename → cleanup.
Ordinary failure or forwarded INT/TERM publishes only a marker-free `.failed` diagnostic directory;
copy failure preserves and prints the exact scratch path. No SIGKILL handling claim is made. The
inherited `64G / 200G / 72h` directives are known repository syntax, **not** a claim
that they are sufficient for full UniRef90; a site-specific parallel-environment directive is also
not guessed. Measure the 30% smoke and preflight input, MMseqs database/index, temporary alignment,
SQLite, retained-manifest, and copy-back space before production. No job is submitted by this
package.

## Future execution examples — not executed in this implementation task

The following tiny fixture smoke exercises parsing, mapping, splitting, labels, export, validation,
and publication without executing MMseqs2. It is **NOT EXECUTED as written here**:

```bash
FIXTURE_MODE=1 NO_DOWNLOADS=1 IDENTITY=30 MIN_COUNT=1 THREADS=1 \
UNIREF90_FASTA=benchmark_builders/homology_cluster/tests/fixtures/uniref90.fasta \
IDMAPPING=benchmark_builders/homology_cluster/tests/fixtures/idmapping_selected.tab \
UNIPROT_SEQUENCES=benchmark_builders/homology_cluster/tests/fixtures/uniprot.fasta \
GOA=benchmark_builders/homology_cluster/tests/fixtures/goa.gaf \
GO_OBO=benchmark_builders/homology_cluster/tests/fixtures/go-mini.obo \
CLUSTER_ASSIGNMENTS=benchmark_builders/homology_cluster/tests/fixtures/clusters.tsv \
OUTPUT_ROOT=/tmp/homology-cluster-fixture-smoke \
TEMP_DIR=/tmp/homology-cluster-fixture-scratch \
bash scripts/benchmark_generation/run_homology_cluster_benchmark.sh
```

The following representative-subset smoke executes MMseqs2 but marks every result fixture-only and
non-production. The five subset files must preserve the frozen GOA/OBO metadata and contain at least
three retained clusters. Replace every placeholder before use. It is also **NOT EXECUTED** here:

```bash
FIXTURE_MODE=1 NO_DOWNLOADS=1 IDENTITY=30 MIN_COUNT=1 THREADS=4 \
UNIREF90_FASTA=/path/to/smoke/uniref90.fasta \
IDMAPPING=/path/to/smoke/idmapping_selected.tab \
UNIPROT_SEQUENCES=/path/to/smoke/uniprot.fasta \
GOA=/path/to/smoke/goa.gaf \
GO_OBO=/path/to/smoke/go-basic.obo \
MMSEQS_BIN=/path/to/mmseqs \
OUTPUT_ROOT=/path/to/smoke-output \
TEMP_DIR=/path/to/smoke-scratch \
bash scripts/benchmark_generation/run_homology_cluster_benchmark.sh \
  --allow-empty-fixture-outputs
```

The `qsub` command above is the future one-threshold production form and was **NOT EXECUTED**.
Replace the five `<64-hex>` values and other placeholders, then use validated persistent paths.

Once the 30% smoke above has been reviewed, the exact six-job submission pattern is the following
loop. It is documentation only and was **NOT EXECUTED**. Export every common variable shown in the
one-threshold command first; each iteration changes only `IDENTITY`:

```bash
for IDENTITY in 30 25 20 15 10 5; do
  export IDENTITY
  qsub -v IDENTITY,SPLIT_POLICY,TRAINING_POPULATION,FRAMEWORK_REVISION,RESULTS_ROOT,\
NO_DOWNLOADS,UNIREF90_FASTA,UNIREF90_FASTA_SHA256,IDMAPPING,IDMAPPING_SHA256,\
UNIPROT_SEQUENCES,UNIPROT_SEQUENCES_SHA256,GOA,GOA_SHA256,GO_OBO,GO_OBO_SHA256,\
FROZEN_INPUT_MANIFEST,MMSEQS_BIN,EXPECTED_MMSEQS_VERSION \
    hpc_jobs/active/hpc_homology_cluster_benchmark.sh
done
```

The first cluster-storage smoke sequence is deliberately ordered:

1. Inspect the compute-node executable with `mmseqs version`, resolve its path, hash it if readable,
   and replace the version placeholder only after review.
2. Fill and review all five frozen-manifest entries and verify their filenames, sizes, hashes,
   embedded GOA/OBO metadata, and authoritative URLs against the staged bytes.
3. Run the representative 30% fixture-only MMseqs command above against cluster scratch and
   persistent result storage; do not submit the six production jobs yet.
4. Validate that direct fixture publication and inspect its `disk_preflight.json`, MMseq logs,
   cluster balance, retained population, audit sizes, peak storage, memory, and walltime.
5. Export exactly the representative-subset variables from the command above and run the same 30%
   fixture through the HPC adapter with
   `qsub -v FIXTURE_MODE,NO_DOWNLOADS,IDENTITY,MIN_COUNT,THREADS,UNIREF90_FASTA,IDMAPPING,UNIPROT_SEQUENCES,GOA,GO_OBO,MMSEQS_BIN,RESULTS_ROOT hpc_jobs/active/hpc_homology_cluster_benchmark.sh`.
   Validate the copied leaf and inspect `hpc_capacity_preflight.json` and
   `hpc_capacity_precopy.json`.
6. Adjust the provisional scheduler resources from those observations, repeat the 30% smoke if
   needed, and obtain review before using the six-job loop.

## Embedding inventory and immutable PFP

After a threshold publishes:

1. point `scripts/verification/inventory_embeddings.py` at the directory containing the nine CSVs;
2. use `configs/embedding_inventory.homology.example.json`;
3. review per-protein/per-modality `reusable`, `missing`, `invalid`, sequence mismatch, and
   provenance-unknown/manual-review results;
4. generate only embeddings not safely reusable;
5. prepare/train immutable PFP.

Split reassignment alone does not invalidate sequence/protein embeddings. Existing modality-specific
provenance rules still apply: exact-sequence ProtT5 reuse is distinct from text/structure/PPI
provenance.

## Local tests

```bash
cd benchmark_builders/homology_cluster
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s tests -v
PYTHONDONTWRITEBYTECODE=1 python3 -m pytest -q -p no:cacheprovider
PYTHONPYCACHEPREFIX=/tmp/homology-pyc python3 -m compileall -q src tests
python3 -m ruff check --no-cache .

cd ../..
bash -n \
  scripts/benchmark_generation/run_homology_cluster_benchmark.sh \
  hpc_jobs/active/hpc_homology_cluster_benchmark.sh
```

The real-MMseqs integration test skips explicitly when MMseqs2 is absent. Tiny tests use a
synthetic, complete `createtsv` fixture and production `min_count` remains 50.

## Limitations and supervisor decisions

- Confirm whether `-e 1e-4` should remain fixed after measured sensitivity analysis at 5–15%.
- Root retention is locked for PFP/TEMPROT compatibility; change it only through a reviewed
  benchmark-contract revision.
- Confirm which OBO relationship types should remain in the DeepGOPlus-compatible propagation
  policy; version 0.1 records the observed types.
- Provide/pin exact archived 2026_02 UniRef90/idmapping and GOA 234 URLs plus expected hashes when
  available.
- `idmapping_selected` can prove present, blank, ambiguous, and absent mappings, but absence alone
  cannot prove that a UniProt accession is obsolete. Such rows are explicitly
  `obsolete-status-unknown-absent-from-idmapping`; an independently frozen deleted-accession source
  is required for a positive obsolete classification.
- `all-cluster-members` needs a supervisor-approved supervised-label policy before it may export
  PFP rows. No zero-negative or homology-transfer policy is assumed.
- Full resource directives need a measured smoke run; none was performed locally.
- Greedy set cover and finite prefilter sensitivity do not prove that every biologically meaningful
  low-identity edge was discovered, especially at 5–15% identity.
- UniRef/MMseq membership is disk-backed and GOA scratch records use low-compression streaming
  spools, but the qualifying accession set, qualifying UniProt sequences, mapping decisions, and
  final supervised pandas frames remain memory-resident because the required pickle contract is
  dataframe-shaped. Peak memory must be measured on the representative smoke before production.
