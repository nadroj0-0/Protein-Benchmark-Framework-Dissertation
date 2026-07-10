#!/usr/bin/env python3
"""Stream UniProt DAT records for selected NCBI taxonomy IDs."""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path


TAXON_RE = re.compile(r"NCBI_TaxID=(\d+)")


def read_taxa(path: Path) -> frozenset[str]:
    taxa = set()
    with path.open() as handle:
        for line in handle:
            value = line.strip()
            if value and not value.startswith("#"):
                taxa.add(value.removeprefix("taxon:"))
    if not taxa:
        raise ValueError(f"No taxonomy IDs found in {path}")
    return frozenset(taxa)


def filter_records(taxa: frozenset[str]) -> tuple[int, int]:
    records = 0
    kept = 0
    record: list[str] = []
    record_taxa: set[str] = set()

    for line in sys.stdin:
        record.append(line)
        if line.startswith("OX"):
            record_taxa.update(TAXON_RE.findall(line))
        if line.rstrip("\n") != "//":
            continue

        records += 1
        if record_taxa & taxa:
            sys.stdout.writelines(record)
            kept += 1
        record = []
        record_taxa = set()

    if record:
        raise ValueError("Truncated UniProt DAT stream: final record has no // terminator")
    return records, kept


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--taxa-file", type=Path, required=True)
    args = parser.parse_args()
    records, kept = filter_records(read_taxa(args.taxa_file))
    print(f"UniProt DAT filter: processed={records} kept={kept}", file=sys.stderr)
    if kept == 0:
        raise SystemExit("UniProt DAT filter retained zero records")


if __name__ == "__main__":
    main()
