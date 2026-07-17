# Implementation report

Status: software and miniature-fixture hardening completed 2026-07-14. No full-data or real-HPC
result is claimed.

## Existing behavior preserved

The implementation remains an isolated `homology_cluster_benchmark` package. UniRef90 is the
MMseqs2 clustering population; the identities remain 30%, 25%, 20%, 15%, 10%, and 5%; coverage is
0.8; the existing MMseqs2 clustering, alignment, sensitivity, E-value, and reassignment parameters
are unchanged. Retained clusters are still selected by Daniel's supplied qualifying evidence-code
set, split as whole clusters, and exported as the same five DeepGOPlus-shaped pickles and nine
PFP-compatible CSVs. Labels remain protein-owned; no homology label transfer or unannotated-negative
policy was introduced.

The contemporary temporal builder, historical CAFA3 validation implementation, embedding
inventory/reuse planner, immutable PFP, embedding generation, and PFP training/evaluation code were
not modified.

## Source-scope and frozen-input contract

Production now requires an explicit `sprot-only`, `trembl-only`, or `sprot-and-trembl` choice.
Scope controls eligible UniProtKB annotation/sequence records, not the UniRef90 clustering
scaffold. Separate Swiss-Prot and TrEMBL declarations are used. The schema-v2 frozen manifest
requires exactly five entries for either one-source scope and six for the combined scope.

Every entry binds logical role, source population, release, local filename and reviewed origin URL,
expected/observed SHA-256, byte size, acquisition action, embedded metadata, origin review, and
notes. Required, irrelevant, mislabelled, cross-release, placeholder, and hash-mismatched inputs
fail. Production uses DAT; Swiss-Prot records must declare `Reviewed` and TrEMBL records must
declare `Unreviewed` on their `ID` lines. Missing-ID records and source-role swaps fail.

Combined-source traversal is fixed-order and source-tagged. A disk-backed accession audit detects
duplicate primary accessions, conflicting sequences, and ambiguous secondary aliases without
depending on caller order. Reports measure primary and secondary identifiers and collision counts
for every involved source. Mapping reports separately count selected-UniProt primary/secondary
resolution, present-UniRef90 mappings, and MMseqs2 assignments, including explicit zero rows for a
selected source. An out-of-scope idmapping hit cannot authorize cluster retention.

## Reviewed attrition gate

Production requires a structured reviewed policy bound to source scope, releases, the full
framework commit, and the frozen-manifest hash. Its eleven metrics cover GOA-to-selected-UniProt and
selected-UniProt-to-UniRef90 mapping, qualifying annotation retention, retained-cluster member
coverage, raw/evaluable proteins, propagated and BP/CC/MF evaluability, and both split deviations.
Each observation records numerator, denominator, definition, ratio, bound, and outcome.

The reviewed manifest/policy contract is parsed before large inputs are hashed. Actual bytes are
then bound to the manifest, and the policy hash is rechecked before evaluation and publication.
Exceeding a reviewed limit prevents production publication unless a separate structured override
exactly binds the failed metrics and run evidence. There is no boolean bypass, and the array
launcher rejects one array-wide override.

A diagnostic pilot may finish to collect observations, but its attrition report always records
`diagnostic: true` and `production_authorized: false`; its publication remains non-production.
Policy, override, approval, and measurement template placeholders are rejected as reviewed data.

## Pilot approval and array authorization

The pilot is locked to array task 1 / 30%. Full-array authorization requires a validated pilot
publication, completion marker, attrition report, worker task context, separately reviewed resource
measurements, reviewed production attrition policy, and human approval.

Approval hashes bind the marker, report, task context, measurement evidence, frozen manifest, and
reviewed policy. The marker, task context, and measurements must agree on pilot job/run identity,
scope, full commit, task 1 / 30%, `smp 2`, `NSLOTS=2`, MMseqs2 threads, methodology, and exact
MMseqs2 version. The pilot attrition report is also tied to the publication's actual run-input
manifest hash. Authorization reconstructs every pilot metric from its numerator/denominator,
verifies the saved ratio and definitions, and requires the reviewed production limits to accept
the pilot observations. The pilot cannot approve itself and no skip flag exists.

## Grid Engine and scratch safety

The deterministic task map is `1→30`, `2→25`, `3→20`, `4→15`, `5→10`, `6→5`. The worker and both
launchers request Grid Engine `smp 2`; production requires `NSLOTS=2`, and MMseqs2 threads equal
`NSLOTS`. Array concurrency remains separate from within-task threading.

Production requires an exact 40-character lowercase commit, detached checkout, exact HEAD match,
and a clean tree before input work. Shared frozen inputs are local and checksum-verified with
downloads disabled; six tasks do not independently download databases. Pilot and full launchers
construct `qsub -t 1` and `qsub -t 1-6` respectively, export an explicit whitelist rather than the
ambient environment, and support scheduler-free dry run. Tests use only a recording fake `qsub`.

Scratch and persistent paths include scope, framework revision, run ID, job ID, task ID, and
identity. Atomic claims and ownership markers prevent collisions. Success, builder failure,
validation failure, INT, TERM, and publication/copy failure attempt marker-free diagnostics and
then remove only owned scratch. Cleanup failure propagates nonzero status; unsafe, broad,
pre-existing, symlinked, or ownership-mismatched targets are refused.

## Publication and aggregation

Each threshold preserves the five-pickle/nine-CSV output contract plus scope-aware manifests,
mapping/annotation/cluster/split/term/attrition summaries, runtime provenance, scientific
fingerprint, strict validation report, output manifest, and last-written completion marker.
Publication remains staged, atomic, non-overwriting, and rehashed at its final path.

Aggregation requires exactly six distinct completed identities. It verifies common source scope,
framework revision, frozen-manifest hash, methodology, scientific fingerprint, and other locked
fields; mixed or incomplete sets fail. The aggregate compares counts, split and attrition metrics,
term universes, and cross-threshold cluster information without copying the large child payloads.

## Local validation results

Observed final commands and results:

```text
PYTHONPATH=src python3 -m unittest discover -s tests -q
  Ran 102 tests in 30.700s — OK (skipped=1)

PYTHONPATH=src python3 -m pytest -q -p no:cacheprovider
  101 passed, 1 skipped in 30.52s

benchmark_builders/contemporary_cafa:
  PYTHONPATH=src python3 -m unittest discover -s tests -q
  Ran 24 tests in 0.607s — OK

benchmark_builders/contemporary_cafa historical-focused module:
  PYTHONPATH=src python3 -m unittest -q tests.test_historical_and_diagnostics
  Ran 6 tests in 0.065s — OK

embedding_inventory:
  PYTHONPATH=. python3 -m unittest discover -s tests -q
  Ran 55 tests in 2.964s — OK (skipped=1)
```

The homology skip is the explicit real-MMseqs2 smoke because `mmseqs` is not installed locally.
The embedding skip is the opt-in real CAFA3 cache integration because
`PFP_RUN_REAL_INTEGRATION=1` was not requested. The homology suite's immutable-PFP ingestion test
did run and consumed all nine generated fixture CSVs from the unchanged local PFP checkout.

Python compilation with a `/tmp` bytecode cache, Ruff, `git diff --check`, and `bash -n` over all
changed shell files passed. ShellCheck was not installed and therefore was not run.

## Explicitly not run

- No `qsub`, `qstat`, `qdel`, SSH, or real scheduler command.
- No full UniRef90, UniProt, GOA, idmapping, or ontology download.
- No full benchmark or full MMseqs2 clustering run on the Mac.
- No embedding generation, PFP training, or PFP evaluation.
- No real production publication or capacity validation.

## Final runtime wrapper added after the original report

Two thin direct-`qsub` wrappers now support the limited-home-quota cluster workflow. The pilot
wrapper runs only task 1 / 30% and automatic review. The array wrapper runs all six fixed identity
tasks and does not require prior pilot evidence. Both delegate to one scratch-first orchestration
script; no scientific transformation was duplicated in the HPC layer.

When local input variables are absent, every task downloads UniRef90, idmapping, the selected
Swiss-Prot/TrEMBL DAT source(s), GOA 234, and GO OBO into node-local scratch after checking the
moving endpoint release markers. Pinned MMseqs2 `18-8cc5c` is also installed into scratch when no
binary path is supplied. Complete hashes, frozen manifest, automatic attrition observations,
strict output validation, five pickles, nine CSVs, and review reports are copied home; inputs and
temporary indexes are discarded.

The full array limits concurrency to two and requests 1200 GB scratch per active task. Copy-back is
atomic. A full home filesystem or quota error makes the task fail but never changes the requirement
to delete owned scratch. These behaviors are fixture-tested; no full download, real MMseqs2 run,
or Grid Engine submission was performed while adding the wrappers.

## Remaining decisions and risks

Daniel still needs to select the production UniProt source scope. A reviewer must set every
attrition limit from the real 30% pilot, approve the frozen input origins/hashes and exact MMseqs2
version, and confirm memory, scratch, walltime, output size, and copy-back behavior before the full
array. The provisional resource requests are syntax, not capacity evidence. Low-identity recall and
biological suitability remain outside internal software validation.

Shell traps cannot clean after SIGKILL or node loss. Claim cleanup intentionally fails closed and
may leave an orphan claim if unexpected files appear. Scratch capacity estimation is conservative
and counts shared input bytes even though they are not copied per task. The worker is larger than a
minimal submission shim because revision, evidence, collision, ownership, publication, and cleanup
checks are operationally enforced there; scientific transformations remain in the Python package.

The exact future nine-step operator workflow, including dry-run, manual pilot/full submission,
review-document hashing, monitoring, and aggregation commands, is in `README.md`. Git commit and
push evidence belongs in the task handoff and repository history; it is reported only after the
operations visibly succeed.
