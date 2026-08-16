# Homology-Cluster Benchmark Builder

This package builds the framework-generated homology benchmarks reported in
the dissertation. It clusters a frozen UniRef50 population independently at
30, 25, 20, 15, 10, and 5 percent identity, assigns whole clusters to splits,
and publishes PFP-compatible CSVs plus DeepGOPlus-shaped pickles.

The submitted workflow does not consume external cluster assignments. Generic
partition comparison remains available for diagnostics, but imported
comparison partitions and their execution routes are not part of this branch.

## Submitted Contract

The canonical profile is `framework-uniref50-s4-defaults`:

```text
UniRef scaffold:             UniRef50 2026_02
MMseqs release:              18-8cc5c
MMseqs full commit:          8cc5ce367b5638c4306c2d7cfc652dd099a4643f
identity thresholds:         30, 25, 20, 15, 10, 5 percent
sensitivity:                 4
coverage:                    0.8
cov-mode:                    0
cluster-mode:                0
alignment-mode:              3
cluster reassignment:        disabled
input shuffle option:        omitted
E-value option:              omitted
UniProt source population:   Swiss-Prot + TrEMBL
split policy:                cluster-count-random
development/test:            0.8 / 0.2
training/validation:         0.9 / 0.1 within development
training population:         annotated-only
seed:                        0
minimum term count:          50
```

`sequence-balanced` remains dormant implementation capability. It is not used
or reported as a submitted result.

The qualifying evidence-code set is:

```text
EXP IDA IPI IMP IGI IEP HTP HDA HMP HGI HEP TAS NAS IGC RCA ND IC
```

Exact pipe-delimited `NOT` annotations are excluded. Cluster membership
controls retention and splitting only: every supervised protein keeps only its
own GO annotations. No representative or neighbour transfers labels.

## Inputs and Caches

The run requires:

- frozen UniRef50 FASTA;
- 22-column UniProt `idmapping_selected`;
- Swiss-Prot and TrEMBL DAT files;
- GOA GAF 2.2;
- GO OBO;
- authenticated common-preprocessing cache;
- authenticated per-threshold MMseqs cluster cache, built when absent unless
  reuse-only operation is explicitly requested.

The common cache binds source bytes and preprocessing implementation. The
cluster cache binds UniRef level/release/hash, exact MMseqs identity and
executable hash, threshold, coverage, modes, sensitivity, and database flags.
Runtime threads and optional framework Git metadata are not scientific cache
identity.

`--require-cluster-cache` prevents an expensive or unintended reclustering.
Matching cached membership is fully revalidated before downstream work.

## Run

Install the package:

```bash
python -m pip install -e benchmark_builders/homology_cluster
```

The scheduler-neutral launcher defaults to the complete submitted profile:

```bash
ARTIFACT_CATALOG=/absolute/path/to/artifact_paths.tsv \
OUTPUT_ROOT=/absolute/path/to/benchmark_outputs \
TEMP_DIR=/absolute/path/to/work \
bash scripts/benchmark_generation/run_homology_cluster_benchmark.sh
```

The launcher resolves authenticated `mmseqs2_executable` from
`ARTIFACT_CATALOG` first, then falls back to explicit `MMSEQS_BIN` and finally
`mmseqs` on `PATH`. Every selected binary must report the exact full commit.

Preview resolved commands without running the build:

```bash
ARTIFACT_CATALOG=/absolute/path/to/artifact_paths.tsv \
DRY_RUN=1 bash scripts/benchmark_generation/run_homology_cluster_benchmark.sh
```

The direct CLI supports one threshold, a subset, or the complete sweep:

```bash
python -m homology_cluster_benchmark build --identity 20 --help
python -m homology_cluster_benchmark build --identity all --help
```

The launcher reads the common cache's frozen manifest and generates the
runtime attrition policy with `homology_cluster_benchmark.runtime_contract`.
Framework Git state is optional report-only provenance and never blocks a run.

## Mapping and Splitting

UniProt accessions remain the supervised protein and embedding keys. UniRef50
IDs define the clustering scaffold, while MMseqs representatives identify the
lower-identity clusters used for retention and splitting.

DAT secondary accessions are resolved conservatively. Ambiguous aliases and
conflicting primary records are reported and excluded or rejected rather than
assigned arbitrarily. GO `alt_id` values canonicalise; obsolete IDs resolve
only through one `replaced_by`; `consider` is never selected automatically.

Random splitting shuffles deterministically sorted cluster IDs with isolated
seeded RNG state. The complete development population defines the term
universe before training/validation division, so test data cannot introduce a
term.

## Publication and Validation

Each threshold publishes the nine CSVs, five pickles, scientific parameters,
source-aware manifests, mapping/annotation/cluster/split/term summaries,
attrition evidence, hashes, validation reports, and `RUN_COMPLETE.json`.
Publication is staged and atomic; an existing destination is not overwritten.

Outputs are namespaced by source scope, UniRef/sensitivity, MMseqs profile,
benchmark scope, identity, split, population, seed, and `min_count`. Framework
commit and scheduler resources do not enter paths or scientific fingerprints.

The `summarize` command validates and aggregates exactly one completed child
for each threshold. Legacy framework fields are report-only; frozen inputs,
normalised attrition policy, clustering profile, and scientific fingerprint
must agree.

Internal validation proves declared byte identities, mapping decisions,
complete cluster membership, whole-cluster and exact-sequence separation,
own-protein labels, development-defined terms, output schemas, and publication
hashes. It does not prove external biological correctness or downstream model
performance.

See `METHODOLOGY.md` for the scientific interpretation boundary.
