#!/usr/bin/env python3
"""Stream one suffix-matched member from a .tar.gz archive to stdout."""

from __future__ import annotations

import argparse
import shutil
import sys
import tarfile
from pathlib import Path


def extract(archive: Path | None, suffix: str) -> str:
    source = archive.open("rb") if archive is not None else sys.stdin.buffer
    found = None
    try:
        with tarfile.open(fileobj=source, mode="r|gz") as tar:
            for member in tar:
                if not member.isfile() or not member.name.endswith(suffix):
                    continue
                if found is not None:
                    raise ValueError(f"Archive contains multiple members ending in {suffix}")
                extracted = tar.extractfile(member)
                if extracted is None:
                    raise ValueError(f"Could not read archive member {member.name}")
                shutil.copyfileobj(extracted, sys.stdout.buffer, length=1024 * 1024)
                found = member.name
        if found is None:
            raise ValueError(f"Archive has no member ending in {suffix}")
        return found
    finally:
        if archive is not None:
            source.close()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, help="Archive path; default reads stdin")
    parser.add_argument("--suffix", required=True)
    args = parser.parse_args()
    member = extract(args.archive, args.suffix)
    print(f"Extracted archive member: {member}", file=sys.stderr)


if __name__ == "__main__":
    main()
