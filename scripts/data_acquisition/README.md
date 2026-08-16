# Data Acquisition

`populate_san_frozen_inputs.sh` is the canonical frozen-input acquisition
program. Its historical filename is retained for compatibility; it accepts any
absolute persistent-store path and has no dependency on a particular storage
system. Its default store is `$HOME/pfp_inputs`.

## Input Specification

`san_frozen_inputs.tsv` declares logical role, release, URL, relative path,
expected byte size, trusted checksum when available, validation rule, and
workflow profile. Mutable upstream release URLs are guarded by the declared
release and observed bytes. Every acquired file receives SHA-256 and provenance
sidecars even where an authoritative checksum was unavailable.

The profiles are:

```text
temporal
homology
embedding-inputs
references
tools
all
```

## Run

Inspect the plan first:

```bash
bash scripts/data_acquisition/populate_san_frozen_inputs.sh \
  --root /absolute/path/to/pfp_inputs \
  --profile all \
  --dry-run
```

Populate or resume:

```bash
bash scripts/data_acquisition/populate_san_frozen_inputs.sh \
  --root /absolute/path/to/pfp_inputs \
  --profile all
```

Audit an existing store without downloading:

```bash
bash scripts/data_acquisition/populate_san_frozen_inputs.sh \
  --root /absolute/path/to/pfp_inputs \
  --profile all \
  --verify-only
```

Use `--list` to inspect selected entries and `--full-verify` to rehash existing
files. `--skip-derived` omits temporal derivatives and homology cache products.

## Integrity and Publication

Downloads use partial files and atomic replacement. Release metadata is checked
before and after mutable-source transfers. Existing files are accepted only
after the configured size, checksum, metadata, and structural checks pass.

The workflow also creates the filtered temporal TrEMBL products, authenticated
UniRef50 common-preprocessing cache, MMseqs cluster-cache root, and an extracted
MMseqs2 executable from the pinned Linux AVX2 archive. The executable must report
the exact accepted commit before it is catalogued. These products are exposed in
`manifests/artifact_paths.tsv`; point later workflows at that catalogue with:

```bash
export ARTIFACT_CATALOG=/absolute/path/to/pfp_inputs/manifests/artifact_paths.tsv
```

Large files are deliberately not tracked by Git. The repository contains only
their acquisition contracts and processing source.
