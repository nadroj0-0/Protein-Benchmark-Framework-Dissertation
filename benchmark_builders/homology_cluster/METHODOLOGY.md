# Methodology and evidence

## Decision classes

This implementation separates three sources of decisions.

### Directly specified by Daniel

- Cluster the frozen UniRef90 sequence population with MMseqs2.
- Run six independent thresholds: 30%, 25%, 20%, 15%, 10%, and 5% identity.
- Preserve 80% sequence coverage.
- Retain clusters connected to at least one latest-GOA annotation using his literal 17-code set.
- Create a new 80:20 whole-cluster development/test split.
- Support both random cluster-count allocation and approximately sequence-balanced allocation.
- Represent the choice between annotated-only training and all cluster members.

The complete operative Part B wording is preserved here:

> B) Homology cluster split
>
> Here we want to cluster uniprot at various degrees of sequence identity at a fixed sequence
> coverage. You would take the uniref90 sequences and then use the mmseqs tool to cluster the
> sequences. With 30, 25, 20, 15, 10 and 5% sequence identity but preserve the coverage that uniref
> use (80%).
>
> Once you have the clusters we only care about the ones that contain a sequence that has an
> experimental evidence code from the latest GOA. Keep those cluster and the sequences assigned to
> them. Then randomly split the clusters 80:20 to train-test set, and this then tells you which
> sequences assigned to which cluster are in the train and test set. This ensures any sequence in a
> given cluster can only
>
> Alternatively, you could try and ensure approx 80% sequencs are in the training set such that if a
> sequence is assigned to a train or test set all its cluster members must join it i.e sequence from
> the same cluster can not appear in both the training and test set.
>
> You then have a further choice should your training data be comprised only of the sequences with
> an experimental GO annotation or will you use all the cluster members in training. This might be
> constrained by how long it takes to make the extra embedding vecotrs you don't already have
> calculated

The exact evidence list Daniel supplied for these benchmark instructions is preserved literally:
`EXP, IDA, IPI, IMP, IGI, IEP, HTP, HDA, HMP, HGI, HEP, TAS, NAS, IGC, RCA, ND, IC`.
Although it is represented as an explicit configuration field for provenance and tests, production
and fixture validation reject any changed set because the experiment forbids additions or removals.

Daniel supplied [James Lingford's MMseqs2 clustering guide](https://www.jameslingford.com/blog/mmseqs-cluster/).
The blog motivates the workflow and `-e 1e-4`; this repair locks that reviewed value while official
MMseqs2 sources determine parameter semantics. Its future biological review is not a runtime option.

### Required by immutable PFP / established benchmark shape

- Nine ontology/split CSVs beginning with plural `proteins,sequences` and development-defined GO columns.
- A separate validation split; version 0.1 takes it only from development clusters.
- DeepGOPlus-shaped dataframe pickles and one-column `terms.pkl`.
- Each protein carries its own labels; unannotated proteins cannot be supervised rows without an
  additional policy.
- `min_count >= 50` in production and deterministic term ordering.
- Daniel did not specify GO-root handling. Retaining the BP/CC/MF roots is inherited from the
  established DeepGOPlus/TEMPROT/PFP interface.

### Implementation choices

- Disk-backed SQLite indexes for the full UniRef90 and MMseqs2 member maps.
- Streaming GAF, FASTA/DAT, and idmapping parsers; only the qualifying annotated subset is held in
  Python annotation/mapping structures.
- Exact local-or-URL input declarations and SHA-256 verification.
- A required reviewed five-source production manifest that distinguishes reproducible bytes,
  recorded acquisition provenance, and an explicit authoritative-origin review. This is an
  auditable declaration rather than a signature or independent origin proof.
- Required production hashes, frozen release checks, and a clean commit-addressable checkout.
- Disk preflight on scratch and publication filesystems plus streamed, bounded downloads.
- Staged publication, final-path hash recheck, and completion marker last.
- Deterministic split tie-breaking through an isolated RNG or SHA-256 seed keys, with bounded
  multi-candidate, single-move, one-for-one, and two-move refinement for sequence balance.
- Disk-backed combined UniRef/UniProt exact-sequence leakage checks and memoized immutable GO
  ancestor closures.
- Strict failure on malformed GAF rows, 22-column idmapping violations, incomplete/duplicate
  MMseqs assignments, mapping conflicts, schemas, or leakage.

## Authoritative sources and how they support implementation

### MMseqs2

- [Official repository](https://github.com/soedinglab/MMseqs2) establishes the maintained tool and
  workflows.
- [Official wiki](https://github.com/soedinglab/MMseqs2/wiki) documents clustering workflows,
  coverage modes, sensitivity, and `createtsv` representative/member output.
- [Latest user guide PDF](https://mmseqs.com/latest/userguide.pdf) documents literal versus
  score-derived identity, alignment modes, cluster modes, and effective parameters.
- [Official parameter definitions](https://github.com/soedinglab/MMseqs2/blob/master/src/commons/Parameters.cpp)
  identify cluster mode 0 as greedy set cover and expose defaults/ranges.
- [Official releases](https://github.com/soedinglab/MMseqs2/releases) show that release 12 fixed
  `--cluster-reassign` coverage-mode handling; the exact modern runtime version is therefore recorded.

Consequences:

- `--shuffle 0` makes `createdb` preserve the frozen FASTA order where MMseqs2 permits it.
- `--min-seq-id X --alignment-mode 3 --seq-id-mode 0` expresses literal aligned identity as identical residues
  divided by alignment columns including internal gaps. Without mode 3, identity is normally
  score-derived. Literal identity is an auditability choice required to give Daniel's percentages
  their direct operational meaning; it is not a claim that literal identity is always the best
  biological similarity measure.
- `-c 0.8 --cov-mode 0` applies `alignment residues / max(query length, target length) >= 0.8`.
- `--cluster-mode 0` produces greedy set-cover clusters. Members meet the representative policy
  after reassign, but arbitrary member-member pairs need not meet it; separate clusters do not prove
  the absence of every qualifying edge.
- `--cluster-reassign 1` corrects cascaded assignments that no longer meet the representative
  criterion. Production pins an exact expected version token, requires one parseable successful
  runtime token, resolves and hashes the executable when readable, and rechecks it before
  publication. A mere arbitrary release >=12 is insufficient; precomputed fixtures are relaxed.
- `-s 7.5` is the documented maximum clustering sensitivity and is deliberately higher than the
  automatic 6.0 used at identities <=30%. It does not guarantee prefilter recall.
- `-e 1e-4` is locked and recorded. The blog uses it, but official sources do not establish that it
  is biologically optimal at 5–15%; E-value also depends on frozen database size.
- The blog's `--max-seqs 20` is intentionally not copied because official documentation warns that
  lowering the prefilter-result cap reduces sensitivity.

The headerless `mmseqs createtsv DB DB CLU out.tsv` contract is
`cluster_representative<TAB>cluster_member`, normally including a representative self-row. The
builder verifies exactly two fields, known identifiers, unique members, complete input coverage,
representative existence, and self-rows.

### UniRef and UniProt

- [UniRef overview](https://www.uniprot.org/help/uniref) and
  [UniRef90 help](https://www.uniprot.org/help/uniref90) establish at least 90% identity and 80%
  overlap with the longest seed for UniRef90 construction.
- [UniRef seed versus representative](https://www.uniprot.org/help/uniref_seed) explains that the
  seed is the longest sequence while the displayed representative is selected with additional
  reviewed/annotation/organism/length criteria.
- [FASTA header specification](https://www.uniprot.org/help/fasta-headers) supports using the first
  whitespace-delimited token such as `UniRef90_A0...` as the MMseqs2 key; `RepID=` is metadata.
- [Official UniProt FTP](https://ftp.uniprot.org/pub/databases/uniprot/) identifies current and
  previous-release structures.
- [ID-mapping README](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/README)
  defines `idmapping_selected.tab.gz` as headerless with 22 fields and UniRef90 at field 9
  (zero-based index 8).

Consequences:

- The frozen UniRef90 FASTA—not UniProt proteins or old benchmark rows—is the MMseqs population.
- Every FASTA ID must start `UniRef90_`; sequence IDs/hashes/lengths are indexed without retaining
  all sequence strings in Python memory.
- The idmapping parser preserves empty tab fields and trailing columns. It never whitespace-splits.
- Canonical UniProt accession, UniRef90 ID, and MMseqs cluster remain separate.
- UniRef IDs can change between releases. A mutable current URL is not frozen evidence; local path,
  release label, size, and SHA-256 bind the supplied snapshot.

### GOA, GAF, and GO ontology

- [GOA downloads](https://www.ebi.ac.uk/GOA/downloads.html) and
  [official GOA FTP](https://ftp.ebi.ac.uk/pub/databases/GO/goa/) identify the annotation source and
  archive/current structures.
- [GAF 2.2 specification](https://geneontology.org/docs/go-annotation-file-gaf-format-2.2/)
  defines 17 tab-separated fields, comment lines, object ID, pipe-delimited relation/qualifier,
  GO ID, evidence, aspect, object type, taxon, assigned date, and gene-product form.
- [GO annotation guidance](https://geneontology.org/docs/go-annotations/) explains `NOT` as an
  explicit negative assertion; it is excluded from a positive-label benchmark.
- [GO evidence guide](https://geneontology.org/docs/guide-go-evidence-codes/) shows why Daniel's
  literal set must not be described as experimental-only: it contains author, curator, and
  computational codes as well as experiments.
- [GO ontology documentation](https://geneontology.org/docs/ontology-documentation/) supports
  namespace, alt ID, obsolete/replacement, and graph traversal handling.
- [GO ontology downloads](https://geneontology.org/docs/download-ontology/) supports pinning a
  release-specific OBO product.

Consequences:

- Exact `DB == UniProtKB` and `DB_Object_Type == protein` filters are applied.
- Qualifiers split on `|`; only exact `NOT` rejects. Positive relations are preserved.
- Field 14 is provenance, never a t0/t1 reconstruction mechanism in this single-snapshot benchmark.
- Reference, With/From, field 14, Assigned By, annotation extension, and field 17 are retained in
  annotation audit records. Field-17 form-specific rows are excluded rather than stripping an
  isoform suffix.
- Alt IDs canonicalize; an obsolete term needs one unambiguous `replaced_by`; `consider` never
  becomes an automatic label.
- Aspect/namespace mismatch is excluded and reported.
- A raw GOA row contributes labels only if that same raw accession has a unique retained mapping;
  an ambiguous secondary alias cannot leak its annotation onto a mapped primary accession.
- Lack of GOA annotation is open-world missing knowledge, not a negative function label.

## Frozen input policy and unresolved URLs

The local contemporary workflow establishes UniProt/UniRef-era release `2026_02`, GOA `234`
generated `2026-06-17`, and GO `releases/2026-06-15` from the `2026-06-19` distribution directory.
Version 0.1 verifies those embedded GOA/OBO conventions.

While 2026_02 remains current, a stable historical UniRef90/idmapping URL may not yet exist. GOA's
current publication state can also lag release announcements. Therefore the implementation accepts
explicit URLs but does not invent frozen archive URLs. Before a production run, provide local frozen
files or checksum-pinned URLs and record official MD5/metalink evidence alongside the computed
SHA-256 where available.

Production requires expected SHA-256 values for every input and a clean, commit-addressable builder
checkout. Hashes are rechecked immediately before publication to detect input replacement during a
run. A configured release label alone is not evidence that bytes came from that release.

The 22-column idmapping table does not contain a deleted-accession lifecycle field. Presence, blank
mapping, ambiguity, secondary-to-primary mapping, and absence are distinct statuses, but an absent
accession is reported as `obsolete-status-unknown-absent-from-idmapping`, not asserted obsolete. A
positive obsolete classification requires a separately frozen authoritative deleted-accession
source that version 0.1 does not accept.

## Cluster retention and label ownership

An MMseqs cluster is retained when at least one qualifying GOA UniProtKB accession maps uniquely to
a UniRef90 member assigned to that cluster. Every UniRef90 member of a retained cluster remains in
the retained-member manifest. Only the qualifying annotated UniProt proteins are supervised PFP
rows in production.

Retention is a population decision, not label transfer. Protein A's annotation never becomes
protein B's annotation merely because their UniRef/MMseqs records share a cluster. Unannotated
members are quantified as additional sequence/embedding burden but never emitted as zero-label
negatives.

## Global cluster splitting and leakage

The 80:20 split is performed once at MMseqs-cluster level. Complete development annotations are
then propagated and counted to freeze the term universe; only afterward is development split at the
same level using the locked 90:10 fraction. BP, CC, and MF are views of this global assignment.
Sequence-balanced weights are retained UniRef90 members; qualifying labelled UniProt counts are
reporting-only secondary weights.

The validator rejects:

- one MMseq cluster assigned to more than one split;
- a UniRef90 member with missing/duplicate MMseq membership;
- one UniProt ID in multiple splits;
- conflicting sequences for one UniProt accession;
- identical UniProt sequences across any pair of train/validation/test;
- identical sequences crossing splits between any retained UniRef90 member and any qualifying
  UniProt protein, using one combined disk-backed query;
- ontology CSV rows inconsistent with global assignment.

Whole-cluster assignment proves separation under the actual MMseqs output. Because greedy set cover
is not exhaustive connected components, it does not prove no undetected qualifying edge exists
between output clusters.

## GO propagation and term universe

Each protein's canonical direct terms propagate through the frozen graph using the established
DeepGOPlus-compatible parent behavior. Relationship targets are included when configured and the
observed types are recorded. The BP, CC, and MF roots remain in propagated annotations and the term
universe when supported, matching the PFP/TEMPROT contract.

The complete 80% development population determines support before its 90:10 whole-cluster split.
A term enters the threshold-specific universe when its propagated support is at least 50; test never
contributes. The builder independently recounts this universe from serialized training plus
validation artifacts, catching missing eligible, below-threshold, and test-only terms. No
temporal/CAFA `terms.pkl` is reused, and no cross-threshold shared universe is forced.

The DeepGOPlus-shaped protein pickles retain each protein's full propagated annotation tuple.
`terms.pkl` and the nine PFP CSVs are the evaluable views: all CSV bits are restricted to the
development-defined universe, and out-of-universe label losses are counted explicitly.

## Determinism and publication

All non-MMseq transformations use sorted IDs, explicit seed 0 by default, stable term/row/manifest
ordering, isolated RNG state, deterministic SHA-256 tie keys, and no Python hash-order dependence.
The tiny fixture is rerun and deterministic payload hashes are compared. Timestamp-bearing
`run_provenance.json`, `output_manifest.json`, and `RUN_COMPLETE.json` are explicitly variable.

MMseqs2 byte identity is not claimed across versions or environments. Exact version, parameters,
input hashes, OS, and runtime versions are recorded instead.

The staged payload is fully validated, hashed, atomically renamed, rehashed at its final path, and
only then marked complete. A marker-free directory is never a successful run. Successful MMseqs2
logs are payload files; failed-run logs and bounded validation diagnostics are retained separately.
Precomputed assignment fixtures have an explicit MMseqs `NOT_EXECUTED` record. All fixture and
representative-subset smoke artifacts declare `production_eligible: false` and cannot be confused
with dissertation production output.

`publication_metadata.json` is written before and hashed by `output_manifest.json`. It binds
fixture/eligibility/scope, identity, split policy, population, seed, min-count, frozen/run input
manifest hashes, exact MMseqs identity, and repository commit. `RUN_COMPLETE.json` repeats those
fields and both metadata/manifest hashes; revalidation requires exact agreement. This prevents a
marker-only relabelling but is not a digital signature against an actor who reauthors every file.

Large audit streams are deterministic gzip. Accepted annotations, retained members, and the one
canonical complete MMseqs membership remain complete; detailed exclusions are bounded first-N
samples per reason while a complete cross-dimensional decision-count table preserves every reject.
Raw MMseqs `createtsv` output remains scratch-only. Six-run aggregation takes exactly six repeated
run roots, validates their hashed common fingerprint, and publishes only small child references and
cross-threshold reports.

The HPC state machine permits final publication only after scratch validation, partial copy, copied
validation, and an atomic rename. Failure and forwarded INT/TERM publish marker-free `.failed`
diagnostics. Copy failure preserves the exact scratch path. Capacity is checked before local input
staging/MMseqs work and again from rounded allocated scratch usage before persistent copy;
configured multipliers, reserves, observations, and `estimates_exact: false` are recorded.

## Remaining decisions for Daniel

1. Confirm `-e 1e-4` after sensitivity/resource measurement at 5–15% identity.
2. Treat root retention as a compatibility lock unless a reviewed contract revision changes it.
3. Confirm whether all OBO relationship targets should remain in propagation or whether a named
   relationship subset is scientifically preferred.
4. Resolve a supervised label policy if `all-cluster-members` is ever to produce PFP rows.
5. Confirm production scratch/memory/walltime after a measured 30% smoke run.
6. Pin exact frozen URLs/hashes once archived 2026_02 UniRef/idmapping and GOA 234 endpoints are
   published/proven.
7. Decide whether an authoritative frozen deleted-accession source should be added to classify
   obsolete UniProt accessions rather than leaving absent cases explicitly unknown.

These are reported uncertainties, not silent substitutions for Daniel's fixed experiment.
