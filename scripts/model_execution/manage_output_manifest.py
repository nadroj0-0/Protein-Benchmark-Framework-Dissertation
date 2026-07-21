#!/usr/bin/env python3
"""Write or verify the hash manifest for a completed model-execution tree."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from common import atomic_write_json, sha256_file


EXCLUDED = {"output_manifest.json", "WORKFLOW_COMPLETE.json"}


def payload(root: Path, include_nested_control_files: bool = False) -> dict:
    files = []
    for path in sorted(root.rglob("*")):
        if not path.is_file():
            continue
        relative = path.relative_to(root).as_posix()
        if (
            relative in EXCLUDED
            or (not include_nested_control_files and path.name in EXCLUDED)
        ):
            continue
        files.append(
            {"path": relative, "bytes": path.stat().st_size, "sha256": sha256_file(path)}
        )
    return {
        "schema_version": 1,
        "payload_file_count": len(files),
        "payload_bytes": sum(item["bytes"] for item in files),
        "files": files,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=("write", "verify"))
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument(
        "--include-nested-control-files",
        action="store_true",
        help="Bind nested manifests/completion markers while excluding only root controls",
    )
    args = parser.parse_args()
    root = args.root.resolve()
    manifest = root / "output_manifest.json"
    observed = payload(root, args.include_nested_control_files)
    if args.action == "write":
        if manifest.exists():
            raise ValueError(f"Output manifest already exists: {manifest}")
        atomic_write_json(manifest, observed)
        print(manifest)
        return 0
    if not manifest.is_file():
        raise FileNotFoundError(f"Output manifest is missing: {manifest}")
    expected = json.loads(manifest.read_text(encoding="utf-8"))
    if expected != observed:
        raise ValueError("Output manifest does not match the published payload")
    print(f"Verified {observed['payload_file_count']} files under {root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
