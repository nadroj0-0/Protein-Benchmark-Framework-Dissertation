# Homology-cluster builder audit fixes

Status updated 2026-07-14 after local implementation and fixture validation. “Resolved” below means
the software contract and automated miniature tests pass; it does not mean a full UniRef90 or real
Grid Engine production run has validated the biological population or resource estimates.

## 1. UniProtKB population and provenance — resolved at interface/fixture level

The undifferentiated `uniprot_sequences` input was replaced with explicit
`sprot-only`, `trembl-only`, and `sprot-and-trembl` scopes plus separate Swiss-Prot and TrEMBL
declarations. Production has no default. UniRef90 remains the clustering scaffold; source scope
controls accession eligibility and supervised rows.

Frozen manifest schema v2 requires exactly five Swiss-Prot roles, five TrEMBL roles, or six
combined roles. It binds role, source population, release, URL, local filename, expected/observed
hash, size, acquisition action, embedded metadata, and origin review. Production sequence sources
must be DAT. FASTA is fixture/diagnostic-only and its missing secondary-accession information is
reported. DAT `ID` content is checked too: Swiss-Prot requires `Reviewed` and TrEMBL requires
`Unreviewed`, so a renamed or misdeclared product fails its source role.

Combined-mode ingestion is fixed-order and source-tagged. A disk-backed collision audit fails
duplicate primaries, conflicting sequences, and ambiguous secondary aliases in strict production
mode. Source-specific reports separate primary/secondary selected-UniProt resolution, present
UniRef90 mapping, and MMseqs2 assignment, with explicit zero rows for selected sources. Tests cover all three
scopes, required/irrelevant sources, 5/5/6 manifests, role/scope/release/hash mismatches,
determinism, primary/secondary mappings, collisions, and out-of-scope retention denial.

Remaining boundary: no full official Swiss-Prot/TrEMBL bytes were processed. Daniel still needs to
choose the production scope.

## 2. Severe attrition gate — resolved at policy/fixture level

Production now requires a reviewed structured policy bound to source scope, releases, full commit,
and frozen-manifest hash. Eleven registered metrics cover selected-source and UniRef90 mapping,
qualifying annotation retention, retained-cluster members, raw/evaluable proteins, propagated and
ontology-specific evaluability, and both split deviations. Reports preserve definitions,
numerators, denominators, ratios, limits, and outcomes.

Exceeding a limit prevents production publication unless an exact structured override binds the
observed failures and run evidence. There is no boolean bypass, and the full-array launcher rejects
one array-wide override. Diagnostic pilot output always reports both non-production eligibility and
`production_authorized: false`. Policy validation occurs before expensive scientific work, and the
manifest/policy contract fails before large-input hashing; the policy hash is rechecked before
publication. Boundary, just-failing, wrong-binding, placeholder,
malformed-policy, malformed-override, and valid-override tests pass.

Remaining boundary: no biological limit values have been approved. They must come from reviewed
pilot evidence rather than software defaults.

## 3. Always clear HPC scratch — resolved for tested normal exit paths

The worker atomically creates a job/task-owned scratch directory and persistent claim. Success,
builder failure, both validation failures, INT, TERM, pre-copy failure, and persistent copy failure
attempt marker-free diagnostic publication and then safe cleanup. Copy failure no longer preserves
scratch. Cleanup failures propagate non-zero status and are not masked.

Empty, root, relative, symlinked, pre-existing, outside-base, or ownership-mismatched paths are
refused. A colliding invocation cannot overwrite an ownership marker, delete another task's
scratch, or move a prior successful publication to `.failed`. Tests cover successful/failing
cleanup, both signals, copy/validation failures, unsafe markers/paths, symlinks, pre-existing work,
persistent collisions, and idempotent absence handling.

Remaining boundary: SIGKILL and node loss cannot execute shell traps; no such guarantee is made.

## 4. Exact framework revision — resolved

Production and diagnostic pilots require exactly 40 lowercase hexadecimal characters. The worker
checks out that commit detached, compares exact `git rev-parse HEAD`, rejects a symbolic branch or
dirty tree, and records the requested/observed commit in task context, scientific fingerprint,
publication metadata, and completion marker. Tests cover missing, branch, tag, abbreviated,
uppercase/malformed, checkout failure, HEAD mismatch, dirty checkout, and an accepted full detached
SHA before later work.

## 5. Grid Engine parallel environment and array workflow — implemented; capacity pending pilot

The worker requests `-pe smp 8`. `NSLOTS` is authoritative, production requires eight, and MMseqs2
threads must match. Task IDs map exactly `1→30`, `2→25`, `3→20`, `4→15`, `5→10`, `6→5`.
The pilot launcher uses `-t 1`; the full launcher uses `-t 1-6`. Both export a whitelist and shared
local checksum-pinned inputs with `NO_DOWNLOADS=1`. Six concurrently runnable tasks may request up
to 48 slots, but Grid Engine controls scheduling.

The full launcher requires a validated 30% diagnostic publication, reviewed attrition policy,
human approval, task context, and separately sourced resource measurements. Launcher-time evidence
hashes are rechecked by the queued worker, which reruns authorization from its detached revision.
Approval binds the reviewed-policy hash and one pilot run ID across marker, task context, and
measurement evidence. Authorization reproduces every saved pilot metric and requires the reviewed
production policy to accept the pilot observations; templates cannot pass as reviewed text.
Fake-`qsub` tests prove dry run makes no scheduler call, exact task/PE arguments are used, approval
failure blocks submission, and a fake scheduler failure propagates.

Remaining boundary: the provisional `64G / 200G / 72h` resources are valid syntax, not evidence of
sufficiency. Runtime, memory, scratch peak, output size, and copy-back require a real reviewed 30%
pilot.

## 6. Publish before production — resolved by the repository commit containing this file

Production rejects dirty or non-commit-addressable framework state. The implementation must be
pushed before a cluster clone can use it. Push success is operational evidence outside this file
and must only be reported after visible Git output; this document does not assume it.

## Current validation boundary

Local tests validate parser behavior, source/manifests, mapping and collision contracts, attrition
and approval bindings, whole-cluster/sequence leakage, output files, aggregation, scheduler command
construction, revision checks, publication isolation, and scratch cleanup using miniature fixtures
and fake executables. Immutable PFP ingestion is exercised when the local checkout is available.

No real `qsub`, `qstat`, `qdel`, SSH, full source download, full MMseqs2 clustering, production HPC
build, embedding generation, or PFP training/evaluation was performed. Real-data mapping coverage,
biological correctness, runtime, peak memory, peak scratch, and copy-back capacity remain unproven.
