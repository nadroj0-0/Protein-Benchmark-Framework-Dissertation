# Methodology and evidence

## Decision classes

### Directly specified by Daniel

Daniel requires MMseqs2 clustering of UniRef90 at 30%, 25%, 20%, 15%, 10%, and 5% identity with
80% coverage. Clusters connected to a qualifying latest-GOA annotation are retained. Whole
clusters—not individual proteins—are split 80:20. He described both random cluster-count and
approximately sequence-balanced allocation, and raised annotated-only versus all-cluster-member
training as a choice.

The exact supplied qualifying evidence-code set is:
`EXP, IDA, IPI, IMP, IGI, IEP, HTP, HDA, HMP, HGI, HEP, TAS, NAS, IGC, RCA, ND, IC`.
It includes more than strictly experimental evidence and is not relabelled as experimental-only.

These instructions rule out CD-HIT, DIAMOND, BLAST, Foldseek, post-hoc homologue removal, changing
coverage, changing the six identity thresholds, splitting members of one cluster, or transferring
annotations from a representative/neighbor.

### Required by established PFP and benchmark contracts

- The nine BP/CC/MF training/validation/test CSVs begin with `proteins,sequences`.
- Five DeepGOPlus-shaped pickles retain their established names and dataframe shapes.
- Complete development defines the term universe before its development-only 90:10 cluster split;
  test never defines trainable terms.
- Each protein carries only its own propagated annotations.
- Production uses `min_count >= 50` and deterministic term ordering.
- BP/CC/MF roots remain a documented DeepGOPlus/TEMPROT/PFP compatibility lock.

`annotated-only` is implemented. `all-cluster-members` fails before data work because unannotated
proteins are not negatives and no supervised-label policy exists.

### Unresolved supervisor policy

Daniel said “uniprot” but did not select Swiss-Prot, TrEMBL, or complete UniProtKB. The code
therefore requires one of `sprot-only`, `trembl-only`, or `sprot-and-trembl`; none is claimed to be
scientifically preferred. Arrays for different scopes must be separate experiments with separate
manifests, result roots/run IDs, attrition evidence, and approval.

Final attrition thresholds are also not a software choice. The implementation defines and measures
the metrics, but a reviewer must set limits from evidence such as the 30% pilot.

## Source populations and frozen provenance

UniRef90 FASTA is always the MMseqs2 population. The selected UniProt source controls which GOA
accessions may supply supervised rows and authorize cluster retention.

Production input roles are:

- always: UniRef90 FASTA, `idmapping_selected`, GOA GAF, GO OBO;
- `sprot-only`: Swiss-Prot DAT, for five total entries;
- `trembl-only`: TrEMBL DAT, for five total entries;
- `sprot-and-trembl`: both DAT files, for six total entries.

Schema-v2 manifests bind role, source population, release, local filename, authoritative URL,
expected/observed SHA-256, size, acquisition action, embedded metadata, notes, and an explicit
authoritative-origin review. The manifest scope must exactly match configured roles. UniRef90,
idmapping, and selected UniProt inputs must share the same frozen release label. A release label or
self-computed hash alone does not prove origin.

DAT is authoritative in production because it includes secondary accessions. FASTA remains useful
for tiny fixtures but cannot reproduce that alias information and is reported as diagnostic-only.
Swiss-Prot DAT records must declare `Reviewed` and TrEMBL DAT records must declare `Unreviewed` on
their `ID` lines; this content check prevents a product from satisfying the other product's role by
filename and manifest declaration alone.
Combined DAT scanning is fixed-order and every record retains its source. A SQLite collision audit
reports per-source primary/secondary counts. Conflicting sequences, duplicate primaries even with
identical sequences, and ambiguous aliases fail production rather than depending on incidental
input order. A GOA accession absent from the selected source cannot retain a cluster merely because
idmapping links it to a UniRef90 entry.

Per-source mapping counts preserve distinct stages: selected-UniProt resolution by primary or
secondary accession, mapping to a present UniRef90 member, and assignment to an MMseqs2 cluster.
Selected populations are reported even when their observed count is zero.

The locked endpoint labels are UniProt/UniRef `2026_02`, GOA `234`, and ontology
`releases/2026-06-15`. Exact archived URLs and hashes must be reviewed when production bytes are
staged; mutable current URLs are not frozen evidence.

## Authoritative references

### MMseqs2

- [Official repository](https://github.com/soedinglab/MMseqs2)
- [Official wiki](https://github.com/soedinglab/MMseqs2/wiki)
- [User guide](https://mmseqs.com/latest/userguide.pdf)
- [Daniel's supplied guide](https://www.jameslingford.com/blog/mmseqs-cluster/)

The implementation uses literal aligned identity with `--alignment-mode 3 --seq-id-mode 0`,
`-c 0.8 --cov-mode 0` for aligned residues divided by the longer sequence, greedy set-cover
cluster mode 0, `--cluster-reassign 1`, `-s 7.5`, and `-e 1e-4`. Greedy set cover does not imply
every member pair meets the threshold, nor does finite prefilter sensitivity prove every remote
homology edge was discovered. Production pins one exact successful MMseqs2 version token,
executable path/hash when readable, and the literal command arguments.

### UniRef and UniProt

- [UniRef overview](https://www.uniprot.org/help/uniref)
- [UniRef90 help](https://www.uniprot.org/help/uniref90)
- [UniProt FASTA headers](https://www.uniprot.org/help/fasta-headers)
- [Official UniProt FTP](https://ftp.uniprot.org/pub/databases/uniprot/)

The headerless idmapping table is parsed as exactly 22 tab-separated columns with UniRef90 at field
9. Empty fields are preserved. UniProt accession, UniRef90 ID, and MMseqs representative never
collapse into one identifier. Absent idmapping rows remain lifecycle-unknown unless a separately
frozen deleted-accession source is approved.

### GOA, GAF, and GO

- [GOA downloads](https://www.ebi.ac.uk/GOA/downloads.html)
- [GOA FTP](https://ftp.ebi.ac.uk/pub/databases/GO/goa/)
- [GAF 2.2 specification](https://geneontology.org/docs/go-annotation-file-gaf-format-2.2/)
- [GO evidence guide](https://geneontology.org/docs/guide-go-evidence-codes/)
- [GO ontology downloads](https://geneontology.org/docs/download-ontology/)

GAF input must have 17 columns and embedded version 2.2. Exact `DB == UniProtKB`, protein object,
P/C/F namespace consistency, Daniel's evidence set, and absence of exact `NOT` are enforced.
Field-17 form-specific and explicit isoform-suffix rows are excluded rather than stripped. GO alt
IDs canonicalize; obsolete terms require one `replaced_by`; `consider` is never promoted. Lack of
GOA annotation is open-world missing knowledge, not a negative function label.

## Retention, splitting, and label ownership

A cluster is retained only when a qualifying raw GOA accession resolves to a selected-source
canonical UniProt sequence, maps uniquely to a present UniRef90 member, and that member has an
MMseqs2 cluster assignment. All UniRef90 members of retained clusters remain in the retained-member
audit; only qualifying annotated UniProt proteins become supervised rows.

The global development/test split is made once at cluster level. Complete development annotations
are propagated and counted to define the term universe; development clusters are then split 90:10
for training/validation. BP, CC, and MF are views of the same assignments. Sequence-balanced
weights are retained UniRef90 member counts. Deterministic candidate ordering and tie hashes do not
claim global subset-sum optimality.

The validator rejects cluster, protein, or exact-sequence leakage; incomplete/duplicate MMseqs
membership; conflicting accession sequences; output rows inconsistent with the global split; and
terms inconsistent with an independent development recount. Whole-cluster separation proves the
contract for the actual MMseqs output, not absence of every possible biological edge.

## Reviewed attrition policy

Every production run loads a structured policy bound to source scope, frozen releases, full
framework commit, and frozen-manifest hash before hashing the large input files or beginning
ontology/MMseqs2 work. The actual input hashes are then checked against the manifest, and the policy
hash is rechecked before evaluation and publication.
It evaluates:

- GOA-to-selected-UniProt and selected-UniProt-to-UniRef90 mapping;
- qualifying annotation-row retention;
- retained UniRef90-cluster member coverage;
- raw qualifying protein evaluability;
- propagated-term and BP/CC/MF evaluability;
- development/test and training/validation member-ratio deviations.

Every observation records an exact numerator, denominator, ratio, definition, bound type, allowed
limit, and outcome. Production eligibility requires structural provenance plus policy authorization.
Diagnostic pilot mode can complete outside reviewed limits for evidence collection, but is always
diagnostic, reports `production_authorized: false`, remains non-production, and is barred from
downstream PFP experiments. Template placeholder text cannot satisfy reviewed policy or override
fields.

A reviewed exception is a separate JSON document bound to the exact failed metric names and
observed values, run input-manifest hash, scope, commit, reviewer/date, justification, and pilot/run
identifier. There is no boolean bypass. One array-wide override is rejected because different
identity tasks can have different failures.

## Pilot approval and resource evidence

Task 1 is locked to the 30% diagnostic pilot. The full array requires a human-completed approval
whose hashes bind the validated pilot publication, completion marker, attrition report, worker task
context, reviewed production attrition policy, and a separate measurement-evidence document. The
task context proves job ID, task 1, 30%, run ID, scope, full commit, requested `smp 2`, `NSLOTS=2`,
and MMseqs threads 2. The completion marker must bind the same run ID.

Authorization reconstructs the pilot's registered observations from their saved numerators and
denominators, verifies each ratio and definition, and re-evaluates them against the reviewed policy.
The array remains blocked if any reviewed limit rejects the 30% pilot. Approval and measurement
placeholder text is rejected rather than treated as human review.

Runtime, peak memory, peak scratch, and output size must be finite positive measurements. Their
measurement-evidence JSON names a reviewed source for each value. Grid Engine accounting may
supply runtime/peak memory; explicit scratch monitoring and filesystem measurement must supply
their corresponding evidence. End-of-run scratch usage is not silently described as a peak. The
pilot never approves itself.

## Grid Engine and publication safety

The array mapping is `1→30`, `2→25`, `3→20`, `4→15`, `5→10`, `6→5`. Each task requests
`-pe smp 2`; `NSLOTS` is authoritative and must equal MMseqs threads. Array concurrency and
within-task threads are different: all six runnable tasks can request up to 12 slots.

Production uses one shared local, checksum-verified input collection with `NO_DOWNLOADS=1`. It
requires an exact 40-character lowercase framework SHA, detached checkout, exact HEAD equality,
and a clean tree before input staging or MMseqs2. The full launcher validates approval and policy;
the queued worker rechecks their launcher-time hashes and reruns authorization from the detached
revision before expensive work.

Each task's scratch path includes job ID, array task, identity, scope, and collision-resistant run
ID. A task atomically creates its scratch and persistent claim; it never overwrites an existing
path or ownership marker. Success, command failure, either validation failure, INT, TERM, and copy
failure attempt marker-free diagnostics and always attempt safe owned-scratch cleanup. Unsafe,
empty, root, relative, symlinked, pre-existing, outside-base, or ownership-mismatched paths are
refused. Persistent ownership prevents a colliding failure from moving or deleting an earlier
successful publication.

## Determinism and validation boundary

Non-MMseq transformations use sorted IDs, isolated seeded RNG state, deterministic hashes, stable
serialization, deterministic gzip, and disk-backed joins. MMseq byte identity across versions or
environments is not claimed; exact runtime provenance is recorded instead. Publications are staged,
strictly validated, hashed, atomically renamed, rehashed at the final path, and marked complete
last. Aggregation validates exactly six common-scope/common-commit/common-manifest publications.

Internal QC validates declared inputs and methodology, mappings, schemas, leakage, term-universe
construction, attrition authorization, and publication integrity. It does not prove biological
optimality, external source correctness beyond reviewed metadata/hash evidence, exhaustive
low-identity recall, scheduler resource sufficiency, or successful PFP training.

## Remaining decisions for Daniel

1. Select Swiss-Prot, TrEMBL, or both for the production eligible UniProtKB population.
2. Review pilot observations and approve explicit limits for every attrition metric.
3. Confirm provisional memory, scratch, and walltime after the real 30% pilot.
4. Revisit `-e 1e-4` after low-identity sensitivity/resource evidence if required.
5. Confirm the established GO relationship and root-retention policies or authorize a revision.
6. Approve a scientifically valid policy before any all-cluster-member supervision.
7. Pin authoritative archived inputs and decide whether a frozen deleted-accession source is needed.

No full dataset, full MMseqs2 run, production benchmark, or real Grid Engine job was used to claim
these software-level results.
