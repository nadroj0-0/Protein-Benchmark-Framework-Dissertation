# Homology-cluster benchmark builder

This isolated package implements Daniel's Part B dissertation benchmark. It clusters the frozen
UniRef90 FASTA independently at 30%, 25%, 20%, 15%, 10%, and 5% sequence identity, with 80%
longest-sequence coverage, retains clusters connected to a qualifying GOA annotation, and assigns
whole clusters to development/test and then training/validation. It publishes the nine CSV files
consumed by immutable PFP and five DeepGOPlus-shaped pickle files.

It does not modify the contemporary or historical benchmark implementations, embedding inventory,
embedding generation, PFP training, or immutable PFP.

## Scientific contract

The three identifiers remain separate:

1. A UniProtKB accession is the GOA key, supervised protein ID, and sequence/embedding key.
2. A UniRef90 identifier is one member of the frozen clustering scaffold.
3. An MMseqs2 representative identifies a lower-identity cluster used for retention and splitting.

Cluster membership controls retention and split assignment only. Each protein retains only its own
GO annotations; no representative, neighbour, or unannotated member supplies labels to another
protein.

The locked experiment is:

- identities `30, 25, 20, 15, 10, 5` percent;
- coverage `0.8` with `--cov-mode 0`;
- MMseqs2 cluster mode 0, alignment mode 3, literal identity mode 0, cluster reassign enabled,
  sensitivity `7.5`, and E-value `1e-4`;
- Daniel's exact supplied qualifying evidence-code set:
  `EXP, IDA, IPI, IMP, IGI, IEP, HTP, HDA, HMP, HGI, HEP, TAS, NAS, IGC, RCA, ND, IC`;
- removal of exact pipe-delimited `NOT` annotations;
- one global 80:20 development/test whole-cluster split, followed by a 90:10
  training/validation whole-cluster split within development;
- either `cluster-count-random` or `sequence-balanced` splitting;
- production `min_count >= 50`, defined from complete development before its 90:10 split;
- `annotated-only` supervised rows.

The supplied evidence set is broader than strictly experimental evidence and is therefore always
described as Daniel's qualifying set. `all-cluster-members` remains an explicit hard failure: lack
of a GOA annotation is not a negative label, and no homology-transfer policy is authorized.

The generated MMseqs2 commands are equivalent to:

```bash
mmseqs createdb FROZEN_UNIREF90_FASTA WORK/uniref90_db --dbtype 1 --shuffle 0
mmseqs cluster WORK/uniref90_db WORK/clusters_db WORK/mmseqs_tmp \
  --min-seq-id X -c 0.8 --cov-mode 0 --cluster-mode 0 \
  --alignment-mode 3 --seq-id-mode 0 --cluster-reassign 1 \
  -s 7.5 -e 1e-4 --threads THREADS
mmseqs createtsv WORK/uniref90_db WORK/uniref90_db \
  WORK/clusters_db WORK/uniref90_clusters.tsv
```

`-s 7.5` is prefilter sensitivity, not one of the six array values. See
[METHODOLOGY.md](METHODOLOGY.md) for parameter semantics and source references.

## Explicit UniProt source scope

Daniel did not specify whether eligible UniProtKB proteins should come from Swiss-Prot, TrEMBL, or
both. Production therefore has no default and requires exactly one explicit scope:

| Scope | Required authoritative sequence source | Forbidden source | Manifest entries |
|---|---|---|---:|
| `sprot-only` | frozen Swiss-Prot DAT | TrEMBL | 5 |
| `trembl-only` | frozen TrEMBL DAT | Swiss-Prot | 5 |
| `sprot-and-trembl` | both frozen DAT files | neither | 6 |

The four always-required inputs are frozen UniRef90 FASTA, 22-column
`idmapping_selected`, GOA GAF, and GO OBO. The selected source scope controls accession eligibility
and supervised rows; it never replaces UniRef90 as the MMseqs2 population.

Production uses separate `--uniprot-sprot-sequences` and
`--uniprot-trembl-sequences` declarations. DAT is required because it records primary and secondary
accessions. FASTA is allowed only for synthetic fixtures and clearly diagnostic use; reports state
that it cannot reproduce DAT secondary-accession handling. DAT content is also checked against its
declared role: Swiss-Prot records must say `Reviewed` on the `ID` line and TrEMBL records must say
`Unreviewed`; renaming or misdeclaring one official product as the other fails.

Combined mode scans sources in a fixed order and retains each record's source population. A
disk-backed collision audit reports primary/secondary and collision counts per involved source. Conflicting sequences,
unexpected duplicate primary accessions (including identical sequences), and ambiguous secondary
aliases fail production. Results do not depend on caller dictionary or traversal order. GOA
accessions outside the selected source may be visible in `idmapping_selected`, but cannot authorize
cluster retention.

The chosen scope is recorded in parameters, the frozen and run input manifests, mapping and
benchmark summaries, attrition policy/report, publication metadata, completion marker, validation
metrics, aggregate metadata, and source-prefixed output paths.

## Scope-aware frozen manifest

`frozen_input_manifest.template.json` is a deliberately non-authorizing Swiss-Prot example. A
production manifest uses schema version 2 and records for every required role:

- logical role and source population;
- release, local filename, authoritative URL, expected/observed SHA-256, and byte size;
- acquisition action, embedded metadata, notes, and authoritative-origin review evidence.

The manifest's `uniprot_source_scope` determines the exact 5/5/6 role set. Missing, extra,
mislabelled, cross-release, placeholder, or hash-mismatched entries fail. UniRef90, idmapping, and
selected UniProt DAT files must share the frozen UniProt/UniRef release. Local paths win over URLs,
but the reviewed URL remains provenance. A self-hash over arbitrary bytes does not establish
authoritative origin.

The currently locked labels are UniProt/UniRef `2026_02`, GOA `234`, and ontology
`releases/2026-06-15`. Exact archived paths and hashes remain inputs to review, not values invented
by the code.

## Mapping, labels, and term universe

The parsers stream plain/gzip FASTA, DAT, GAF, OBO, and idmapping sources. UniRef90 and cluster
membership are indexed in SQLite. GAF must be version 2.2 with 17 columns; only `UniProtKB`
proteins, aspects P/C/F, Daniel's qualifying codes, and non-`NOT` assertions enter the mapping
chain. Primary, secondary-to-primary, ambiguous, blank, absent, missing-from-FASTA, and
out-of-selected-scope outcomes remain distinct.

Source-specific reports separate selected-UniProt resolution by primary/secondary accession from
successful UniRef90 mapping and from final MMseqs2 assignment. Every selected source receives an
explicit row, including a zero-count row, so a later mapping-stage failure cannot erase earlier
accession-resolution evidence.

GO alt IDs canonicalize. Obsolete IDs resolve only through one `replaced_by`; `consider` is never
chosen automatically. Propagated labels retain the BP/CC/MF roots inherited from the established
DeepGOPlus/TEMPROT/PFP contract. The complete development population defines each threshold's term
universe; test never contributes. A term is not removed merely because one downstream split has no
positive examples.

The two split policies assign indivisible, deterministically ordered MMseqs2 clusters:

- `cluster-count-random` shuffles sorted IDs with isolated seeded RNG state and targets cluster
  counts;
- `sequence-balanced` targets UniRef90 member counts with deterministic bounded candidate, move,
  and swap refinement.

Achieved ratios and the largest indivisible cluster are reported. Split deviations are governed by
the reviewed attrition policy; there is no hidden hard-coded biological tolerance.

## Production attrition gate

Production requires a reviewed JSON policy based on evidence such as the 30% diagnostic pilot.
`attrition_policy.template.json` contains placeholders and cannot authorize a run. It binds the
scope, releases, full framework commit, frozen-manifest hash, author/reviewer/date, rationale, and
an explicit minimum or maximum for every registered metric:

| Metric | Numerator / denominator | Bound |
|---|---|---|
| `goa_to_selected_uniprot_mapping_ratio` | qualifying GOA accessions resolving in selected UniProt / all qualifying GOA accessions | minimum |
| `selected_uniprot_to_uniref90_mapping_ratio` | selected accessions uniquely mapped to a present UniRef90 member / selected accessions | minimum |
| `qualifying_annotation_retention_ratio` | qualifying annotation rows completing the retained-cluster chain / qualifying annotation rows entering it | minimum |
| `retained_cluster_member_ratio` | members of retained clusters / all frozen UniRef90 members | minimum |
| `evaluable_protein_ratio` | raw qualifying GOA proteins with an evaluable propagated term / raw qualifying GOA proteins | minimum |
| `propagated_term_evaluable_ratio` | labelled rows retaining a development-universe term / labelled rows with propagated annotations | minimum |
| `bp_evaluable_ratio` | labelled rows with an evaluable BP term / labelled rows with propagated annotations | minimum |
| `cc_evaluable_ratio` | labelled rows with an evaluable CC term / labelled rows with propagated annotations | minimum |
| `mf_evaluable_ratio` | labelled rows with an evaluable MF term / labelled rows with propagated annotations | minimum |
| `development_split_deviation` | absolute achieved-minus-requested development member fraction / one ratio unit | maximum |
| `training_split_deviation` | absolute achieved-minus-requested training-within-development member fraction / one ratio unit | maximum |

The code does not choose final biological thresholds. A failed limit prevents production
publication. A deliberate exception requires `attrition_override.template.json`'s structured
document bound to the exact failed metric names and observed ratios, run input-manifest hash,
scope, commit, reviewer/date, justification, and run identifier. There is no
`--ignore-attrition` flag. Because failures can differ by identity, the six-task launcher rejects
one array-wide override; any exception must be reviewed for the specific failed task before a
deliberate rerun.

Diagnostic pilot mode records all measurements but is always `diagnostic-pilot`,
`production_eligible: false`, `production_authorized: false`, and unsuitable for downstream PFP
experiments even if its measurement-only limits pass. Reviewed policy/override placeholder text is
rejected as data. The reviewed manifest and production policy contracts are parsed and mutually
bound before any large input is hashed, before ontology parsing, and before MMseqs2. The actual
input bytes are then hashed and bound to that manifest, and the policy hash is rechecked before
evaluation and publication.

## Output contract and isolation

Every threshold continues to publish:

```text
train_data.pkl              bp-training.csv    bp-validation.csv    bp-test.csv
train_data_train.pkl        cc-training.csv    cc-validation.csv    cc-test.csv
train_data_valid.pkl        mf-training.csv    mf-validation.csv    mf-test.csv
test_data.pkl
terms.pkl
```

CSV rows begin with `proteins,sequences`; remaining columns are sorted ontology-specific GO IDs.
Rows are sorted by accession and labels are binary. Proteins with no evaluable term in an ontology
are counted and omitted from that ontology view, never emitted as all-zero negatives.

Each Python publication is isolated beneath:

```text
OUTPUT_ROOT/source_SCOPE/framework_SHA12/identity_XX/SPLIT/annotated-only/seed_N/min_count_N/
```

It also contains source-aware frozen/input manifests, `parameters.json`,
`run_provenance.json`, mapping/annotation/cluster/split/term summaries,
`attrition_policy.json`, `attrition_report.json`, an optional validated override, strict validation
reports, hashes, and `RUN_COMPLETE.json`. The marker is written only after staged validation,
atomic rename, and final-path hash verification. Existing outputs are never overwritten.

`summarize` requires exactly one valid publication for each identity and rejects missing,
duplicate, mixed-scope, mixed-commit, mixed-manifest, or mixed-method runs. It writes threshold,
split, and attrition metrics plus term-universe overlaps/changes and identical-cluster warnings; it
does not copy the large child payloads.

## Grid Engine contract

The normal shell entrypoint delegates scientific work to Python:

```bash
scripts/benchmark_generation/run_homology_cluster_benchmark.sh
```

The array worker is `hpc_jobs/active/hpc_homology_cluster_benchmark.sh`. It requests
`-pe smp 8`, uses scheduler-provided `NSLOTS`, and requires `THREADS == NSLOTS` (exactly eight in
production). Six simultaneously runnable tasks can therefore consume up to 48 CPU slots; Grid
Engine decides actual concurrency.

The locked mapping is:

```text
SGE_TASK_ID 1 -> 30%     SGE_TASK_ID 2 -> 25%     SGE_TASK_ID 3 -> 20%
SGE_TASK_ID 4 -> 15%     SGE_TASK_ID 5 -> 10%     SGE_TASK_ID 6 -> 5%
```

Missing/invalid task IDs and conflicting `IDENTITY` fail. Production requires `JOB_ID`, an exact
40-character lowercase `FRAMEWORK_REVISION`, detached checkout, exact `git rev-parse HEAD`, and a
clean tree. It requires shared local checksum-pinned inputs with `NO_DOWNLOADS=1`; six tasks never
download the source collection independently. The launcher exports a whitelist rather than the
submission shell's full environment.

Scratch and final paths include job ID, task ID, identity, source scope, run ID, split policy, and
framework revision. Scratch is atomically claimed with an ownership marker. Success, command or
validation failure, INT, TERM, and copy failure all attempt marker-free diagnostics and then remove
only task-owned scratch. Empty, root, relative, symlinked, pre-existing, out-of-base, or
marker-mismatched deletion targets are refused. Persistent paths are atomically claimed so a
colliding task cannot alter a prior final or another task's partial output.

The full launcher re-hashes pilot evidence when called; the queued worker checks the exported
hashes and re-runs authorization from its detached framework checkout before input staging. This
closes the queue-time mutation window.

## Final scratch-first HPC workflow

The guarded workflow below remains available when a persistent, manually reviewed frozen-input
collection exists. For the current limited-home-quota setup, use the newer final runtime wrappers:

```bash
# Recommended 30% pilot plus automatic review; does not submit the full array.
qsub hpc_jobs/active/hpc_homology_cluster_runtime_pilot.sh

# All six identities. A pilot is recommended but is not mechanically required.
qsub hpc_jobs/active/hpc_homology_cluster_runtime_array.sh
```

Both are thin Grid Engine wrappers around
`scripts/benchmark_generation/run_homology_cluster_runtime_hpc.sh`. Inputs supplied through
`UNIREF90_FASTA`, `IDMAPPING`, `UNIPROT_SPROT_SEQUENCES`, `UNIPROT_TREMBL_SEQUENCES`, `GOA`, and
`GO_OBO` are copied to node-local scratch. Any omitted input is downloaded from its reviewed
official source at runtime. The driver checks that the moving UniProt/GOA endpoints still identify
UniProt `2026_02` and GOA `234` before accepting their bytes, validates GOA and GO metadata,
computes complete SHA-256 manifests, and installs pinned MMseqs2 `18-8cc5c` into scratch when no
executable path is supplied.

The default source scope is `sprot-and-trembl`; set `UNIPROT_SOURCE_SCOPE=sprot-only` or
`trembl-only` to run another declared population. Every array task has independent node-local
scratch, so missing inputs are downloaded once per active task. The 30% pilot requests approximately
300 GB total: UCL Grid Engine treats consumable resources as per-slot requests, so eight SMP slots
use `tscratch=38G` (304 GB total) and `tmem=8G` (64 GB total), while `scratch0free=300G` is only a
host free-space threshold. It runs as a measurement job: the old speculative multiplier defaults
are reduced globally to neutral `1x` values, and the pilot records rather than enforces the
resulting estimate. The full wrapper's resource request must be recalibrated from the pilot before
submission; it is not evidence that the old 1200 GB estimate was necessary.

The runtime driver records job-owned allocated bytes every 120 seconds and at explicit checkpoints
covering input staging, MMseqs2 installation, builder execution, and validation. Copied results
include `logs/disk_usage.tsv`, a per-directory `logs/disk_usage_by_path.tsv`, and
`logs/disk_usage_summary.tsv`. The summary's `peak_work_bytes` is the measured pilot high-water
mark; apply deliberate headroom before converting it into a future Grid Engine request.

Each task writes the five DeepGOPlus-shaped pickles, nine PFP CSVs, provenance, validation,
attrition observations, automatic review, and logs below
`$HOME/homology_cluster_benchmark_results` (or `RESULTS_ROOT`). Inputs and disposable clustering
state are never copied home. The final copy is staged and atomically renamed. Copy failure is
reported as a failed job, any incomplete home copy is removed, and job-owned scratch is still
deleted as required by cluster etiquette.

The runtime policy intentionally uses non-blocking bounds while preserving all attrition
measurements for review. This allows the six-task array to run without treating a pilot as an
authorization dependency. It validates software and output contracts, not biological optimality.

## Pilot-first production workflow (future commands only)

These commands are documentation for Jordan to run after review. They were not executed while
implementing this change.

### 1. Declare one shared frozen input collection

This example deliberately selects Swiss-Prot only; use the corresponding TrEMBL variables or both
sources for a separately named array after Daniel's scope decision.

```bash
cd '/Users/jordansydney-darlin/Documents/University/Data Science with Machine Learning/Project/Protein-Benchmark-Framework'
export PYTHONPATH="$PWD/benchmark_builders/homology_cluster/src"
export UNIPROT_SOURCE_SCOPE=sprot-only
export UNIPROT_RELEASE='2026_02'
export GOA_RELEASE='234'
export ONTOLOGY_RELEASE='releases/2026-06-15'
export FRAMEWORK_REVISION='<40-lowercase-hex-reviewed-commit>'
export RUN_ID='<collision-resistant-pilot-run-id>'
export RESULTS_ROOT='/persistent/homology-results'
export FROZEN_INPUT_MANIFEST='/persistent/inputs/homology-sprot.reviewed.json'
export UNIREF90_FASTA='/persistent/inputs/uniref90.fasta.gz'
export UNIREF90_FASTA_SHA256='<reviewed-64-hex>'
export IDMAPPING='/persistent/inputs/idmapping_selected.tab.gz'
export IDMAPPING_SHA256='<reviewed-64-hex>'
export UNIPROT_SPROT_SEQUENCES='/persistent/inputs/uniprot_sprot.dat.gz'
export UNIPROT_SPROT_SEQUENCES_SHA256='<reviewed-64-hex>'
unset UNIPROT_TREMBL_SEQUENCES UNIPROT_TREMBL_SEQUENCES_SHA256
export GOA='/persistent/inputs/goa_uniprot_all.gaf.234.gz'
export GOA_SHA256='<reviewed-64-hex>'
export GO_OBO='/persistent/inputs/go-basic.obo'
export GO_OBO_SHA256='<reviewed-64-hex>'
export MMSEQS_BIN='/path/on/compute/node/mmseqs'
export EXPECTED_MMSEQS_VERSION='<exact-mmseqs-version-token>'
export SPLIT_POLICY=sequence-balanced
export TRAINING_POPULATION=annotated-only
export SEED=0
export MIN_COUNT=50
python3 -m json.tool "$FROZEN_INPUT_MANIFEST" >/dev/null
```

### 2. Preview the one-item pilot without contacting Grid Engine

```bash
DRY_RUN=1 bash hpc_jobs/launchers/submit_homology_cluster_pilot.sh
```

The preview prints mode, exact `qsub -t 1 -pe smp 8` command, task range, scope, full commit,
manifest, result root, methodology fields, input paths/hashes, and every exported value.

### 3. Manually submit the 30% diagnostic pilot later

```bash
DRY_RUN=0 bash hpc_jobs/launchers/submit_homology_cluster_pilot.sh
```

### 4. Validate and review the copied pilot

```bash
export PILOT_RUN_DIR='<published diagnostic-pilot benchmark leaf>'
export PILOT_COMPLETION_MARKER="$PILOT_RUN_DIR/RUN_COMPLETE.json"
export PILOT_ATTRITION_REPORT="$PILOT_RUN_DIR/attrition_report.json"
export PILOT_TASK_CONTEXT='<outer task result>/logs/hpc_task_context.json'
python3 -m homology_cluster_benchmark validate --run-dir "$PILOT_RUN_DIR"
python3 -m json.tool "$PILOT_COMPLETION_MARKER"
python3 -m json.tool "$PILOT_ATTRITION_REPORT"
python3 -m json.tool "$PILOT_TASK_CONTEXT"
```

Review all mapping/attrition ratios, five pickles, nine CSVs, validation report, MMseqs2 version,
copy-back, logs, and source-specific counts. Obtain runtime and peak memory from reviewed Grid
Engine accounting, peak scratch from explicit task scratch monitoring, and output bytes from a
reviewed filesystem measurement. Do not relabel end-of-run scratch size as peak usage.

### 5. Create reviewed policy, measurement evidence, and approval

Copy the three non-authorizing templates. Edit the reviewed policy first, using the pilot
observations to set every limit. Then calculate evidence hashes, complete the measurement file,
calculate its hash, and only then complete the approval. Keep `approved: false` until a human has
reviewed all evidence; change it to `true` only as the final manual approval action.

```bash
cp benchmark_builders/homology_cluster/attrition_policy.template.json \
  '/persistent/review/homology-attrition-policy.json'
cp benchmark_builders/homology_cluster/pilot_measurement_evidence.template.json \
  '/persistent/review/homology-pilot-measurements.json'
cp benchmark_builders/homology_cluster/pilot_approval.template.json \
  '/persistent/review/homology-pilot-approval.json'
export ATTRITION_POLICY='/persistent/review/homology-attrition-policy.json'
export PILOT_MEASUREMENT_EVIDENCE='/persistent/review/homology-pilot-measurements.json'
export PILOT_APPROVAL='/persistent/review/homology-pilot-approval.json'

vi "$ATTRITION_POLICY"
python3 -m json.tool "$ATTRITION_POLICY" >/dev/null
export REVIEWED_ATTRITION_POLICY_SHA256="$(sha256sum "$ATTRITION_POLICY" | awk '{print $1}')"
export PILOT_COMPLETION_MARKER_SHA256="$(sha256sum "$PILOT_COMPLETION_MARKER" | awk '{print $1}')"
export PILOT_ATTRITION_REPORT_SHA256="$(sha256sum "$PILOT_ATTRITION_REPORT" | awk '{print $1}')"
export PILOT_TASK_CONTEXT_SHA256="$(sha256sum "$PILOT_TASK_CONTEXT" | awk '{print $1}')"
export FROZEN_INPUT_MANIFEST_SHA256="$(sha256sum "$FROZEN_INPUT_MANIFEST" | awk '{print $1}')"

vi "$PILOT_MEASUREMENT_EVIDENCE"
python3 -m json.tool "$PILOT_MEASUREMENT_EVIDENCE" >/dev/null
export PILOT_MEASUREMENT_EVIDENCE_SHA256="$(
  sha256sum "$PILOT_MEASUREMENT_EVIDENCE" | awk '{print $1}'
)"

vi "$PILOT_APPROVAL"
python3 -m json.tool "$PILOT_APPROVAL" >/dev/null

python3 -m homology_cluster_benchmark authorize-array \
  --attrition-policy "$ATTRITION_POLICY" \
  --pilot-approval "$PILOT_APPROVAL" \
  --pilot-completion-marker "$PILOT_COMPLETION_MARKER" \
  --pilot-attrition-report "$PILOT_ATTRITION_REPORT" \
  --pilot-run-dir "$PILOT_RUN_DIR" \
  --pilot-task-context "$PILOT_TASK_CONTEXT" \
  --pilot-measurement-evidence "$PILOT_MEASUREMENT_EVIDENCE" \
  --frozen-input-manifest "$FROZEN_INPUT_MANIFEST" \
  --framework-revision "$FRAMEWORK_REVISION" \
  --uniprot-source-scope "$UNIPROT_SOURCE_SCOPE" \
  --split-policy "$SPLIT_POLICY" \
  --training-population "$TRAINING_POPULATION" \
  --expected-mmseqs-version "$EXPECTED_MMSEQS_VERSION" \
  --uniprot-release "$UNIPROT_RELEASE" \
  --goa-release "$GOA_RELEASE" \
  --ontology-release "$ONTOLOGY_RELEASE"
```

Populate the JSON fields with the calculated values above. In particular,
`reviewed_attrition_policy_sha256` is `$REVIEWED_ATTRITION_POLICY_SHA256`; it is not the pilot
report's `policy_sha256`. The latter identifies the software-generated, non-production measurement
policy used to let the diagnostic pilot finish and cannot authorize the full array. The approval's
other evidence hashes must use the corresponding calculated variables, and `approved` is changed
to `true` only after manual review.

The approval binds task 1/30%, successful marker hash, complete validated pilot publication,
attrition report hash and scope/commit bindings, reviewed attrition-policy hash, task-context hash,
measurement-evidence hash, pilot job/run IDs, commit, frozen manifest, scope, methodology, MMseqs2
version, positive finite measurements, reviewer/date, and notes. Authorization reconstructs every
pilot metric from its numerator and denominator, checks the recorded definitions/ratio, and
requires the reviewed production limits to accept those observations. The marker, task context,
and measurement evidence must all bind the same run ID. Template placeholders are rejected. The
pilot never approves itself.

### 6. Preview the authorized six-task array

Set a new collision-resistant `RUN_ID` for the full array, then preview:

```bash
export RUN_ID='<collision-resistant-full-array-run-id>'
DRY_RUN=1 bash hpc_jobs/launchers/submit_homology_cluster_array.sh
```

Authorization still runs; only `qsub` is suppressed.

### 7. Manually submit the full array later

```bash
DRY_RUN=0 bash hpc_jobs/launchers/submit_homology_cluster_array.sh
```

### 8. Monitor the future array on the cluster

```bash
qstat -j '<array-job-id>'
```

This monitoring command is documentation only and was not run during implementation.

### 9. Validate and aggregate all six completed publications

```bash
python3 -m homology_cluster_benchmark summarize \
  --output-dir '/persistent/homology-aggregate' \
  --run-dir '<30-percent-publication-leaf>' \
  --run-dir '<25-percent-publication-leaf>' \
  --run-dir '<20-percent-publication-leaf>' \
  --run-dir '<15-percent-publication-leaf>' \
  --run-dir '<10-percent-publication-leaf>' \
  --run-dir '<5-percent-publication-leaf>'
```

## Local validation

```bash
cd benchmark_builders/homology_cluster
PYTHONPYCACHEPREFIX=/tmp/homology-cluster-pycache PYTHONPATH=src \
  python3 -m unittest discover -s tests
PYTHONPYCACHEPREFIX=/tmp/homology-cluster-pycache python3 -m compileall -q src tests
ruff check --no-cache src tests

cd ../..
bash -n scripts/benchmark_generation/run_homology_cluster_benchmark.sh \
  hpc_jobs/active/hpc_homology_cluster_benchmark.sh \
  hpc_jobs/launchers/_homology_cluster_common.sh \
  hpc_jobs/launchers/submit_homology_cluster_pilot.sh \
  hpc_jobs/launchers/submit_homology_cluster_array.sh
```

Launcher tests use only a temporary fake `qsub` that records arguments and controlled return
statuses. The real-MMseqs2 integration test skips explicitly when MMseqs2 is absent.

## Validation boundary and unresolved decisions

Internal validation proves declared hashes/metadata, parser/filter behavior, mapping and source
scope, complete MMseq membership, whole-cluster and exact-sequence separation, own-protein labels,
development-defined terms, output schemas/content, attrition authorization, and publication hashes.
It cannot prove external biological correctness, exhaustive low-identity edge discovery, resource
sufficiency, or downstream model performance.

Daniel still needs to decide the production UniProt source scope. Reviewed pilot evidence must set
the attrition limits and confirm memory, scratch, walltime, and provisional `64G / 200G / 72h`
resources. The `-e 1e-4`, GO relationship policy, root compatibility lock, archived input URLs, and
any future all-cluster-member supervision remain explicit review points.

No real scheduler command, full source download, full MMseqs2 clustering, production benchmark,
embedding generation, or PFP training/evaluation was performed while implementing this workflow.
