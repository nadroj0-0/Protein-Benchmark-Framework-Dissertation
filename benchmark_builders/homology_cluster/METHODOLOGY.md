# Homology Benchmark Methodology

## Question

The homology experiment asks whether PFP performance changes when complete
sequence clusters, rather than individual proteins, are kept on one side of
the development/test boundary. It tests robustness across the six implemented
clustering thresholds; it does not prove that every cross-split protein pair
is below the named percentage identity.

## Submitted Population

The framework starts from frozen UniRef50 2026_02 and maps eligible Swiss-Prot
and TrEMBL accessions through `idmapping_selected`. Qualifying GOA annotations
use the declared evidence-code set and exclude `NOT` assertions. GO labels are
canonicalised and propagated through the frozen ontology.

MMseqs2 release `18-8cc5c`, full commit
`8cc5ce367b5638c4306c2d7cfc652dd099a4643f`, clusters the selected UniRef50
scaffold independently at 30, 25, 20, 15, 10, and 5 percent identity with:

```text
sensitivity 4
coverage 0.8
cov-mode 0
cluster-mode 0
alignment-mode 3
cluster reassignment disabled
no explicit input-shuffle option
no explicit E-value option
```

This is the neutral `framework-uniref50-s4-defaults` profile. Framework code
generates every submitted partition. External comparison partitions are not
required to reproduce the reported framework series.

## Labels

Cluster membership determines retention and split assignment, never label
transfer. Each UniProt accession receives only its own qualifying annotations.
An unannotated member is not treated as a negative example and cannot donate or
receive labels.

The complete development population defines terms with `min_count >= 50`.
Training and validation are then separated without allowing test labels to
affect the term universe.

## Splits

The submitted policy is `cluster-count-random` with seed 0:

1. sort cluster identifiers deterministically;
2. shuffle them with isolated seeded RNG state;
3. assign complete clusters 80:20 to development/test;
4. assign complete development clusters 90:10 to training/validation.

Only annotated proteins become supervised rows. The dormant
`sequence-balanced` implementation is retained as generic capability but is
not a submitted result.

## Integrity

Frozen input manifests bind release metadata and observed bytes. Common and
cluster caches are authenticated and revalidated before reuse. Exact MMseqs
binary identity, clustering flags, source scope, thresholds, split policy,
seed, and term policy enter the scientific contract.

Framework Git revision, dirty state, machine paths, scheduler resources, and
runtime timestamps are optional provenance only. They neither block execution
nor alter scientific identity.

Every threshold is staged, validated, hashed, and atomically published with a
completion marker. Validation checks complete scaffold assignment,
whole-cluster separation, exact-sequence separation, own-protein labels,
schema consistency, nonempty required views, and manifest integrity.

## Interpretation Boundary

The threshold is an MMseqs clustering parameter applied to a UniRef50
scaffold. It is not an exhaustive all-pairs maximum-identity guarantee. Similar
performance across thresholds supports stability to the tested partitioning
regimes. It does not establish unrestricted generalisation below 10 or 5
percent identity.

The 10 and 5 percent framework assignment files may be compared directly. If
they are identical in a fresh run, they must be reported as one effective
partitioning outcome rather than two independent experiments.

Internal validation establishes implementation and contract correctness. It
cannot establish external biological truth, annotation completeness, or
statistical significance from a single model seed.
