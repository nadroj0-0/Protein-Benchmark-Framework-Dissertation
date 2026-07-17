# Data Acquisition

This directory contains two intentionally separate acquisition workflows.

## Persistent SAN acquisition

`populate_san_frozen_inputs.sh` is the authoritative persistent-input loader for
the dissertation project store:

```text
/SAN/bioinf/bmpfp
```

Its catalogue is committed as `san_frozen_inputs.tsv`. The catalogue records
the logical role, release, SAN-relative path, authoritative URL, expected byte
count where known, published or previously authenticated checksum where known,
and structural validator for every object.

The six large 2026 UniProt/UniRef/GOA/GO SHA-256 values were recorded by the
authenticated homology runtime acquisition on 16 July 2026. The MMseqs2 hash
was verified against the pinned `18-8cc5c` release archive. Zenodo file sizes
and MD5 values come from the immutable record metadata for records 7409660 and
19498341; the three published embedding archives also retain the SHA-256 values
already enforced by the framework's reproduction workflows. For immutable
dated sources that publish no trusted digest, the first successful HTTPS
acquisition records a local SHA-256 which all later full verifications enforce.

The profiles are:

| Profile | Contents |
|---|---|
| `temporal` | UniProt 2025_01 and 2026_02 data required by the contemporary temporal builder; GOA 225 and 234; the frozen t0/t1 GO products. |
| `homology` | UniProtKB 2026_02 Swiss-Prot and TrEMBL DAT files, UniRef90, `idmapping_selected`, GOA 234, and GO 2026-06-19. |
| `embedding-inputs` | STRING v12.0 network embeddings used by PFP PPI extraction. |
| `references` | Canonical CAFA3 CSVs and GO OBO from Zenodo 7409660, DeepGOPlus CAFA intermediates, and Zijian's published MMFP embeddings/checkpoints/splits from Zenodo 19498341. |
| `tools` | Pinned MMseqs2 `18-8cc5c` Linux AVX2 archive. |
| `all` | Every row above, without duplicate downloads for files shared by profiles. |

Inspect the plan before starting the large transfer:

```bash
bash scripts/data_acquisition/populate_san_frozen_inputs.sh --dry-run
```

At the time the catalogue was frozen, `all` represented approximately 382 GiB
of files with known sizes plus a conservative 11 GiB allowance for small files
whose servers did not expose a pinned size. The default preflight also requires
40 GiB to remain free. The script sums only files that are currently missing,
so later profile runs do not reserve space for data already present.

Populate everything:

```bash
bash scripts/data_acquisition/populate_san_frozen_inputs.sh --profile all
```

Or populate only what one workflow needs:

```bash
bash scripts/data_acquisition/populate_san_frozen_inputs.sh \
  --profile homology --profile tools
```

Because `--profile all` transfers roughly 400 GB, run the real acquisition as
a scheduled cluster job rather than as a long process on the login node.

### Idempotency and integrity

- Existing authenticated files are skipped; they are never downloaded again.
- Interrupted transfers remain as `<filename>.partial` and are resumed.
- Downloads are validated before an atomic rename publishes the final path.
- Known file sizes and checksums are pinned in the committed catalogue.
- Every downloaded file receives a SHA-256 sidecar and provenance TSV.
- `/SAN/bioinf/bmpfp/manifests/frozen_input_catalog.tsv` is rebuilt
  deterministically from the per-file provenance records.
- Mutable UniProt endpoints are checked before and after transfer. GOA 234 is
  downloaded from EBI's immutable historical release URL and verified using
  its pinned size, SHA-256 and embedded GAF release metadata. A later GOA
  release therefore cannot be stored under the frozen release-234 path.
- Normal reruns perform quick metadata checks. `--verify-only` or
  `--full-verify` re-read and structurally validate every selected file, which
  is intentionally slow for the largest archives.
- A process lock prevents two acquisitions writing the same SAN paths at once.

The script refuses to trust an already-present, unpinned file unless it has the
SHA-256 sidecar produced by an earlier successful acquisition. This prevents a
partially copied or manually substituted file from being silently adopted.

### Deliberate exclusions

This input loader does not generate benchmarks, embeddings, PDBs, checkpoints,
or model results. Those are derived outputs and belong under the corresponding
`benchmarks/`, `embeddings/`, `models/`, and `runs/` SAN directories. Per-run
AlphaFold downloads and model caches remain workflow-managed rather than being
treated as frozen global inputs.

## Historical home-directory acquisition

`protein_database_download.sh` is retained unchanged as evidence of the older
`$HOME/protein_databases` workflow used during initial investigation. It is not
the SAN population entrypoint and should not be repurposed: changing it would
erase the behavior that earlier logs and diary entries refer to.

`inspect_protein_databases.sh` likewise describes that older home-directory
layout.
