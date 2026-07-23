# PFP benchmark forensics

This package produces a strict, reusable characterization of one or more
PFP-compatible nine-CSV benchmarks. It is intended for published CAFA3,
contemporary temporal, and future homology-controlled benchmark outputs.

It reports:

- label, root-only, annotation-count, sequence-length, and term-support profiles;
- exact root-only provenance when pre-projection PFP pickles are supplied;
- complete organism/taxonomy distributions for every aspect and split;
- auditable resolution of legitimate taxonomy changes across ordered snapshots;
- optional real protein-family or other category distributions from explicit maps;
- per-modality coverage for every aspect and split;
- modality co-availability patterns;
- direct descriptive deltas between configured benchmarks;
- SHA-256 input/output manifests and an atomic completion marker.

The tool never calls taxa "protein families". True protein-family analysis is
available through `category_sources` only when an InterPro, Pfam, or equivalent
protein-to-family map is explicitly supplied.

## Install and run

```bash
python -m pip install -e benchmark_forensics
pfp-benchmark-forensics \
  --config configs/benchmark_forensics.example.json \
  --output-dir /path/to/benchmark_forensics
```

Existing output is refused by default. Use `--replace` to publish a complete
replacement only after the new staged report has been generated successfully.

This is a CPU-only analysis. A GPU does not accelerate CSV parsing, joins,
counting, OBO traversal, or manifest hashing. Runtime is normally governed by
the total CSV/DAT bytes and input hashing rather than model computation.

## Taxonomy source resolution

Each taxonomy source may declare a unique `name` and integer `priority`.
Higher priority wins, independently of configuration order. The supplied HPC
configuration ranks release recency first and reviewed Swiss-Prot over
unreviewed TrEMBL within a release. This retains older releases as coverage
fallbacks without allowing them to overwrite newer mappings.

If lower-priority sources assign a different taxon to the same accession, the
selected and alternative assignments are written to `taxonomy_conflicts.tsv`,
and the selected source is recorded in `protein_membership.tsv`. This is an
audited historical resolution, not a silent overwrite.

UniProt DAT inputs distinguish each record's primary accession from secondary
historical aliases. When one alias occurs on records with different taxa, the
benchmark sequence is compared with those records. A unique matching taxon is
selected and audited. Zero-match, multi-taxon exact-match, and equal-priority
cross-source conflicts are written as `unresolved`, left taxonomically unmapped,
and do not abort the remaining descriptive analysis. This avoids guessing while
preventing one historical alias from discarding an otherwise valid benchmark
report.

When `name` and `priority` are omitted, the source path is used as the name and
priority defaults to zero. Equal-priority disagreements are retained as
unresolved rather than resolved by file or configuration order.

## Root-only provenance

For a root-only CSV row, the tool distinguishes:

- `source_root_only`: the pre-projection source contained only the aspect root;
- `projection_created`: the source contained informative terms in the aspect,
  but none survived the configured label-universe projection;
- `source_no_aspect_annotation`: a source row exists but has no resolvable term
  in the aspect;
- `source_unresolved`: no matching source annotation row was supplied.

The tool only attributes `projection_created` to `min_count=50` when that policy
is declared in the configuration and the supplied pickles are known to be the
pre-projection inputs to that CSV export.

## Modality coverage meanings

For standard `embedding_inventory.tsv.gz` inputs, four separate states are
reported rather than collapsed into one ambiguous percentage:

- `artifact_exists`;
- `artifact_valid`;
- `scientifically_eligible`;
- `planned_reuse`.

These are not interchangeable. For example, an array can physically exist but
lack sufficient provenance for scientific reuse.
