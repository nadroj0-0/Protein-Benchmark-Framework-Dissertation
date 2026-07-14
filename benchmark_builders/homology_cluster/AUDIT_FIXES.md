# Homology Cluster Builder Audit Fixes

Status: identified during the full implementation audit on 2026-07-14. These
items are recorded for follow-up and have not yet been implemented.

## Required before production

### 1. Represent the complete UniProtKB sequence input faithfully

The builder currently accepts one `uniprot_sequences` input, while the official
UniProt release publishes Swiss-Prot and TrEMBL separately. Production must
either:

- accept both frozen inputs separately; or
- create a deterministic combined derivative whose two source files, hashes,
  releases, construction command, output hash, and review are recorded.

The reviewed solution must also state whether FASTA or DAT is authoritative.
DAT preserves secondary accession aliases; FASTA normally exposes only the
primary accession.

Acceptance criteria:

- both reviewed and unreviewed UniProtKB populations are represented;
- every upstream source and derived byte stream is hash-bound;
- manifests and reports describe the derivation without relying on an
  undocumented local file;
- fixtures cover primary and secondary accession mapping.

### 2. Make severe attrition a production gate

Low mapping coverage, proteins without evaluable terms, and annotation mapping
exclusions are currently warnings. A structurally valid but severely truncated
benchmark can therefore still be marked production-eligible.

Add an explicit, reviewed attrition policy. The policy should fail production
when its limits are exceeded, unless an override with a recorded justification
is deliberately supplied. Do not choose biological thresholds silently; use
pilot-run evidence and record the approved values.

Acceptance criteria:

- production eligibility incorporates the attrition decision;
- reports retain all numerator, denominator, reason, and ratio fields;
- an override is explicit and auditable rather than automatic;
- boundary and failure tests cover every gated metric.

### 3. Always clear HPC scratch

The current wrapper retains scratch after copy/publication failure. This
conflicts with the established UCL cluster policy for this project. The wrapper
must attempt to publish diagnostics, report any publication failure, and then
remove its scratch run directory on every normal EXIT path, including failed
copy-back.

Acceptance criteria:

- success, command failure, validation failure, INT, TERM, and copy failure all
  remove the job-owned scratch directory;
- diagnostics are copied when possible, but failed persistent storage never
  causes scratch retention;
- wrapper tests, README, and methodology text describe the same behavior.

### 4. Require an exact framework revision for production HPC jobs

`FRAMEWORK_REVISION` is currently optional. Require a full commit SHA for
production, check out that revision in detached mode, and verify that `HEAD`
matches it before starting expensive work.

Acceptance criteria:

- production fails before input staging when the revision is absent or invalid;
- all six threshold jobs can be proven to use the same commit;
- fixture/test mode remains explicitly separated.

### 5. Request and use an appropriate HPC parallel environment

Without a Grid Engine parallel-environment allocation, `NSLOTS` defaults to one
and MMseqs2 runs single-threaded. Confirm the correct UCL parallel-environment
name and slot policy, then add or document the exact submission contract.

Acceptance criteria:

- requested slots and `THREADS` agree;
- the wrapper still refuses oversubscription;
- the selected slot count and provisional disk/runtime requests are justified
  by a measured 30% identity pilot before the six production jobs.

### 6. Publish the implementation before production use

At audit time the homology package and wrappers were untracked. Production
already rejects dirty/non-commit-addressable source, so the reviewed
implementation must be committed and pushed before an HPC clone can run it.

## Recommended order

1. Fix unconditional scratch cleanup and its tests/documentation.
2. Require and verify `FRAMEWORK_REVISION`.
3. Resolve the complete UniProtKB input and provenance contract.
4. Add the reviewed attrition gate.
5. Confirm UCL parallel-environment settings and run a 30% pilot.
6. Review pilot mapping, attrition, runtime, memory, and disk reports.
7. Only then submit the six production thresholds.

## Current validation boundary

The fixture implementation has strong automated coverage, but no full UniRef90
MMseqs2 clustering run or production HPC build has yet validated real-data
mapping coverage, runtime, memory, scratch estimates, or copy-back capacity.
