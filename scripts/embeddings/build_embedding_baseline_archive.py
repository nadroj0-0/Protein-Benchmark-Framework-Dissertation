#!/usr/bin/env python3
"""Validate a generated PFP cache and publish an archive baseline pair."""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import os
import tempfile
from collections import Counter
from pathlib import Path

from manage_embedding_archive import create_archive_from_directories, sha256_file
from manage_resumable_embedding_state import load_policy, load_targets, validate_array


REPORT_MODALITIES = {
    "sequence": "prott5",
    "text": "text",
    "structure": "structure",
    "ppi": "ppi",
}


def atomic_gzip_report(
    destination: Path,
    targets: dict[str, str],
    policies: dict[str, dict],
    available: dict[str, set[str]],
) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{destination.name}.", suffix=".tmp", dir=str(destination.parent)
    )
    try:
        with os.fdopen(descriptor, "wb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
                with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as text:
                    writer = csv.DictWriter(
                        text,
                        fieldnames=["protein_id", "modality", "status", "dimension"],
                        delimiter="\t",
                    )
                    writer.writeheader()
                    for protein_id in sorted(targets):
                        for modality in sorted(policies):
                            writer.writerow({
                                "protein_id": protein_id,
                                "modality": REPORT_MODALITIES[modality],
                                "status": (
                                    "available"
                                    if protein_id in available[modality]
                                    else "missing"
                                ),
                                "dimension": policies[modality]["dimension"],
                            })
        os.replace(temporary_name, destination)
    finally:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass


def build_baseline(
    generated_cache_root: Path,
    data_dir: Path,
    policy_path: Path,
    archive_path: Path,
    assembly_report: Path,
) -> dict:
    if archive_path.exists() or assembly_report.exists():
        raise ValueError("Baseline archive outputs must not already exist")
    targets = load_targets(data_dir)
    policies = load_policy(policy_path)["modalities"]
    available: dict[str, set[str]] = {}
    counts: Counter[str] = Counter()

    for modality, specification in policies.items():
        directory = generated_cache_root / specification["cache_directory"]
        directory.mkdir(parents=True, exist_ok=True)
        accepted: set[str] = set()
        for path in sorted(directory.glob("*.npy"), key=lambda item: item.name):
            protein_id = path.stem
            if protein_id not in targets:
                raise ValueError(f"Generated cache contains unknown target: {path}")
            validate_array(path, int(specification["dimension"]))
            accepted.add(protein_id)
        available[modality] = accepted
        counts[modality] = len(accepted)

    assembly_report.parent.mkdir(parents=True, exist_ok=True)
    archive_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_report = assembly_report.parent / f".{assembly_report.name}.building"
    temporary_archive = archive_path.parent / f".{archive_path.name}.building"
    for path in (temporary_report, temporary_archive):
        path.unlink(missing_ok=True)
    try:
        atomic_gzip_report(temporary_report, targets, policies, available)
        archive_summary = create_archive_from_directories(
            generated_cache_root,
            temporary_archive,
            [specification["cache_directory"] for specification in policies.values()],
            allow_unexpected_top_level=True,
        )
        os.replace(temporary_report, assembly_report)
        os.replace(temporary_archive, archive_path)
    finally:
        temporary_report.unlink(missing_ok=True)
        temporary_archive.unlink(missing_ok=True)

    return {
        "schema_version": 1,
        "target_count": len(targets),
        "available_by_modality": dict(sorted(counts.items())),
        "available_pairs": sum(counts.values()),
        "missing_pairs": len(targets) * len(policies) - sum(counts.values()),
        "archive": str(archive_path.resolve()),
        "archive_sha256": sha256_file(archive_path),
        "archive_size_bytes": archive_path.stat().st_size,
        "archive_member_count": archive_summary["member_count"],
        "assembly_report": str(assembly_report.resolve()),
        "assembly_report_sha256": sha256_file(assembly_report),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--generated-cache-root", type=Path, required=True)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--policy", type=Path, required=True)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--assembly-report", type=Path, required=True)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    try:
        result = build_baseline(
            args.generated_cache_root,
            args.data_dir,
            args.policy,
            args.archive,
            args.assembly_report,
        )
        if args.report:
            args.report.parent.mkdir(parents=True, exist_ok=True)
            args.report.write_text(
                json.dumps(result, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
    except (OSError, ValueError, json.JSONDecodeError) as error:
        raise SystemExit(f"ERROR: {error}") from error
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
