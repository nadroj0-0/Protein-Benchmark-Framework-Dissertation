#!/usr/bin/env python3
"""Reduce PFP PubMedBERT hidden-state arrays to the consumed CLS vector."""

from __future__ import annotations

import argparse
import json
import os
import time
from pathlib import Path

import numpy as np


def process_file(path: Path) -> str:
    array = np.load(path, allow_pickle=False)
    if array.shape == (768,):
        return "already_cls"
    if array.ndim != 3 or array.shape[0] != 1 or array.shape[2] != 768:
        raise ValueError(f"Unexpected text embedding shape for {path.name}: {array.shape}")
    cls = np.asarray(array[0, 0, :], dtype=np.float32)
    temporary = path.with_suffix(".cls.tmp.npy")
    np.save(temporary, cls)
    os.replace(temporary, path)
    return "reduced"


def embedding_paths(directory: Path):
    for path in sorted(directory.glob("*.npy")):
        if path.name.endswith((".tmp.npy", ".cls.tmp.npy")):
            continue
        yield path


def sweep(directory: Path, processed: set[str]) -> dict[str, int]:
    counts = {"reduced": 0, "already_cls": 0}
    for path in embedding_paths(directory):
        if path.name in processed:
            continue
        try:
            counts[process_file(path)] += 1
        except FileNotFoundError:
            # The producer may have atomically renamed a file between the
            # directory scan and load. Leave it for the next sweep.
            continue
        processed.add(path.name)
    return counts


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--watch-pid", type=int)
    parser.add_argument("--poll-seconds", type=float, default=2.0)
    parser.add_argument("--report", type=Path)
    return parser.parse_args()


def process_exists(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def main() -> int:
    args = parse_args()
    args.directory.mkdir(parents=True, exist_ok=True)
    totals = {"reduced": 0, "already_cls": 0, "sweeps": 0}
    processed: set[str] = set()

    while True:
        counts = sweep(args.directory, processed)
        totals["reduced"] += counts["reduced"]
        totals["already_cls"] += counts["already_cls"]
        totals["sweeps"] += 1
        if args.watch_pid is None or not process_exists(args.watch_pid):
            final = sweep(args.directory, processed)
            totals["reduced"] += final["reduced"]
            totals["already_cls"] += final["already_cls"]
            totals["sweeps"] += 1
            break
        time.sleep(args.poll_seconds)

    final_paths = list(embedding_paths(args.directory))
    final_count = len(final_paths)
    report = {
        "schema_version": 1,
        "directory": str(args.directory.resolve()),
        "final_file_count": final_count,
        "newly_reduced": totals["reduced"],
        "sweeps": totals["sweeps"],
        "all_cls": all(
            np.load(path, mmap_mode="r", allow_pickle=False).shape == (768,)
            for path in final_paths
        ),
    }
    if not report["all_cls"]:
        raise SystemExit("One or more text arrays are not 768-dimensional CLS vectors")
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(
            json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
