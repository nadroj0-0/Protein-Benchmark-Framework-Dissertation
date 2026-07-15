#!/usr/bin/env python3
"""Create a provenance-recorded PPI extractor copy for a UniProt-only benchmark."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path


OLD = "mapped_cafa3/len(cafa3_ids)*100:.1f"
NEW = "mapped_cafa3/max(1, len(cafa3_ids))*100:.1f"


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if not args.source.is_file():
        raise SystemExit(f"Missing PFP PPI extractor: {args.source}")
    if args.output.exists():
        raise SystemExit(f"Refusing to overwrite compatibility copy: {args.output}")

    source_text = args.source.read_text(encoding="utf-8")
    occurrences = source_text.count(OLD)
    if occurrences != 1:
        raise SystemExit(
            f"Expected exactly one validated PPI denominator expression, found {occurrences}"
        )
    output_text = source_text.replace(OLD, NEW)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(output_text, encoding="utf-8")
    shutil.copymode(args.source, args.output)

    report = {
        "schema_version": 1,
        "source": str(args.source.resolve()),
        "source_sha256": digest(args.source),
        "compatibility_copy": str(args.output.resolve()),
        "compatibility_copy_sha256": digest(args.output),
        "change": "Guard CAFA3-only coverage logging against an empty CAFA ID set.",
        "scientific_output_change": False,
        "upstream_source_modified": False,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
