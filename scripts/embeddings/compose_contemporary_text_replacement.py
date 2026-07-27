#!/usr/bin/env python3
"""Replace one finalized contemporary text layer without retaining old text."""

from __future__ import annotations

import argparse
import json
import os
import shutil
import tempfile
import uuid
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from build_embedding_baseline_archive import atomic_gzip_report
from manage_embedding_archive import create_archive, extract_archive, sha256_file
from manage_resumable_embedding_state import (
    load_policy,
    load_target_tables,
    validate_array,
)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def load_json(path: Path) -> dict:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return value


def atomic_write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent)
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(value, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, path)
    finally:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass


def resolve_declared_file(root: Path, relative_name: object, label: str) -> Path:
    if not isinstance(relative_name, str) or not relative_name:
        raise ValueError(f"{label} does not declare a file")
    relative = Path(relative_name)
    if relative.is_absolute() or ".." in relative.parts:
        raise ValueError(f"{label} declares an unsafe path: {relative_name}")
    resolved_root = root.resolve()
    resolved = (root / relative).resolve()
    if resolved != resolved_root and resolved_root not in resolved.parents:
        raise ValueError(f"{label} escapes its run root: {relative_name}")
    if not resolved.is_file() or resolved.is_symlink():
        raise ValueError(f"{label} is missing or unsafe: {resolved}")
    return resolved


def modality_directories(config_path: Path) -> dict[str, str]:
    config = load_json(config_path)
    modalities = config.get("modalities")
    if not isinstance(modalities, dict) or set(modalities) != {
        "sequence",
        "text",
        "structure",
        "ppi",
    }:
        raise ValueError("Run config must define exactly four PFP modalities")
    result = {}
    for modality, specification in modalities.items():
        if not isinstance(specification, dict):
            raise ValueError(f"Invalid run-config modality: {modality}")
        directory = specification.get("directory")
        if not isinstance(directory, str) or Path(directory).name != directory:
            raise ValueError(f"Unsafe run-config directory for {modality}: {directory}")
        result[modality] = directory
    if len(set(result.values())) != len(result):
        raise ValueError("Run config repeats an embedding directory")
    return result


def validate_cache(
    cache_root: Path,
    targets: dict[str, str],
    policy: dict,
) -> tuple[dict[str, set[str]], dict[str, int]]:
    available: dict[str, set[str]] = {}
    counts: Counter[str] = Counter()
    for modality, specification in policy["modalities"].items():
        directory = cache_root / specification["cache_directory"]
        if not directory.is_dir() or directory.is_symlink():
            raise ValueError(f"Missing or unsafe {modality} cache: {directory}")
        accepted: set[str] = set()
        for path in sorted(directory.iterdir(), key=lambda item: item.name):
            if not path.is_file() or path.is_symlink() or path.suffix != ".npy":
                raise ValueError(f"Unexpected cache entry: {path}")
            protein_id = path.stem
            if protein_id not in targets:
                raise ValueError(f"Cache contains an unknown target: {path}")
            validate_array(path, int(specification["dimension"]))
            accepted.add(protein_id)
        available[modality] = accepted
        counts[modality] = len(accepted)
    return available, dict(sorted(counts.items()))


def require_archive_hash(path: Path, expected: object, label: str) -> str:
    if not isinstance(expected, str) or len(expected) != 64:
        raise ValueError(f"{label} has no valid SHA-256 declaration")
    observed = sha256_file(path)
    if observed != expected:
        raise ValueError(f"{label} SHA-256 mismatch: {expected} != {observed}")
    return observed


def compose(args: argparse.Namespace) -> dict:
    if args.output_root.exists():
        raise ValueError(f"Output root already exists: {args.output_root}")
    if args.work_dir.exists() and any(args.work_dir.iterdir()):
        raise ValueError(f"Work directory is not empty: {args.work_dir}")
    for path in (
        args.base_final_root,
        args.replacement_run_root,
        args.plan_dir,
    ):
        if not path.is_dir():
            raise ValueError(f"Required directory is missing: {path}")
    for path in (args.policy, args.config, args.input_acquisition):
        if not path.is_file():
            raise ValueError(f"Required file is missing: {path}")

    policy = load_policy(args.policy)
    directories = modality_directories(args.config)
    for modality, specification in policy["modalities"].items():
        if directories[modality] != specification["cache_directory"]:
            raise ValueError(
                f"Policy/config directory mismatch for {modality}: "
                f"{specification['cache_directory']} != {directories[modality]}"
            )
    text_directory = directories["text"]

    base_marker_path = args.base_final_root / "FINAL_CACHE_COMPLETE.json"
    base_marker = load_json(base_marker_path)
    if base_marker.get("complete") is not True or base_marker.get("validated") is not True:
        raise ValueError("Base finalized cache is not complete and validated")
    archive_name = base_marker.get("archive_name")
    if not isinstance(archive_name, str) or Path(archive_name).name != archive_name:
        raise ValueError("Base marker has an unsafe archive name")
    base_archive = args.base_final_root / archive_name
    if not base_archive.is_file() or base_archive.is_symlink():
        raise ValueError(f"Base archive is missing or unsafe: {base_archive}")
    base_archive_sha256 = require_archive_hash(
        base_archive, base_marker.get("archive_sha256"), "Base archive"
    )

    replacement_marker_path = args.replacement_run_root / "TEXT_GENERATION_COMPLETE.json"
    replacement_marker = load_json(replacement_marker_path)
    required_replacement_values = {
        "complete": True,
        "mode": "full-text-generation-only",
        "requested_cutoff": args.expected_cutoff,
        "effective_cutoff": args.expected_cutoff,
        "old_text_carried_forward": False,
        "hydration_performed": False,
        "state_modified": False,
    }
    for key, expected in required_replacement_values.items():
        if replacement_marker.get(key) != expected:
            raise ValueError(
                f"Replacement marker mismatch for {key}: "
                f"{replacement_marker.get(key)!r} != {expected!r}"
            )
    replacement_archive = resolve_declared_file(
        args.replacement_run_root,
        replacement_marker.get("archive"),
        "Replacement text archive",
    )
    replacement_archive_sha256 = require_archive_hash(
        replacement_archive,
        replacement_marker.get("archive_sha256"),
        "Replacement text archive",
    )

    target_tables = (
        args.plan_dir / "reuse_proteins.tsv",
        args.plan_dir / "regenerate_proteins.tsv",
    )
    targets = load_target_tables(target_tables)
    if len(targets) != replacement_marker.get("target_count"):
        raise ValueError(
            "Replacement target count differs from the bound planner population: "
            f"{replacement_marker.get('target_count')} != {len(targets)}"
        )

    args.work_dir.mkdir(parents=True, exist_ok=True)
    base_cache = args.work_dir / "base_cache"
    replacement_cache = args.work_dir / "replacement_cache"
    roundtrip_cache = args.work_dir / "roundtrip_cache"
    base_extract = extract_archive(base_archive, base_cache, args.config)
    replacement_extract = extract_archive(
        replacement_archive, replacement_cache, args.config
    )
    replacement_counts = replacement_extract["members_by_directory"]
    for modality, directory in directories.items():
        observed = int(replacement_counts.get(directory, 0))
        if modality == "text":
            if observed != replacement_marker.get("text_available"):
                raise ValueError(
                    f"Replacement text count mismatch: {observed} != "
                    f"{replacement_marker.get('text_available')}"
                )
        elif observed != 0:
            raise ValueError(
                f"Replacement archive unexpectedly contains {observed} {modality} arrays"
            )

    old_text = base_cache / text_directory
    new_text = replacement_cache / text_directory
    old_text_count = sum(1 for _ in old_text.glob("*.npy"))
    if old_text_count <= 0:
        raise ValueError("Base cache contains no old text layer to replace")
    shutil.rmtree(old_text)
    os.replace(new_text, old_text)
    for modality, directory in directories.items():
        if modality != "text" and any((replacement_cache / directory).iterdir()):
            raise ValueError("Replacement extraction contains unconsumed cache entries")

    available, combined_counts = validate_cache(base_cache, targets, policy)
    if combined_counts["text"] != replacement_marker.get("text_available"):
        raise ValueError("Combined cache did not retain the complete replacement text layer")
    for modality, directory in directories.items():
        if modality == "text":
            continue
        expected = int(base_extract["members_by_directory"].get(directory, 0))
        if combined_counts[modality] != expected:
            raise ValueError(
                f"Base {modality} membership changed: "
                f"{expected} != {combined_counts[modality]}"
            )

    args.output_root.parent.mkdir(parents=True, exist_ok=True)
    staging = args.output_root.parent / (
        f".{args.output_root.name}.staging-{uuid.uuid4().hex}"
    )
    published = False
    try:
        archive_output = staging / "archive" / "contemporary_embedding_cache.tar.gz"
        assembly_output = staging / "reports" / "assembly" / "embedding_assembly.tsv.gz"
        archive_output.parent.mkdir(parents=True)
        assembly_output.parent.mkdir(parents=True)
        archive_report = create_archive(base_cache, archive_output, args.config)
        atomic_gzip_report(assembly_output, targets, policy["modalities"], available)
        shutil.copy2(
            args.input_acquisition,
            staging / "reports" / "input_acquisition.tsv",
        )

        provenance = staging / "provenance"
        provenance.mkdir()
        shutil.copy2(base_marker_path, provenance / "BASE_FINAL_CACHE_COMPLETE.json")
        shutil.copy2(
            replacement_marker_path,
            provenance / "TEXT_GENERATION_COMPLETE.json",
        )
        for name in (
            "CACHE_ARCHIVE_VALIDATED.json",
            "evidence/contract.json",
            "evidence/coverage.json",
        ):
            source = args.base_final_root / name
            if source.is_file():
                destination = provenance / Path(name).name
                shutil.copy2(source, destination)

        roundtrip_report = extract_archive(
            archive_output, roundtrip_cache, args.config
        )
        for key in (
            "archive_sha256",
            "member_count",
            "members_by_directory",
            "member_content_sha256",
        ):
            if roundtrip_report.get(key) != archive_report.get(key):
                raise ValueError(f"Replacement archive round trip differs for {key}")
        _, roundtrip_counts = validate_cache(roundtrip_cache, targets, policy)
        if roundtrip_counts != combined_counts:
            raise ValueError("Replacement archive round trip changed modality counts")

        manifest = {
            "schema_version": 1,
            "complete": True,
            "created_at": utc_now(),
            "operation": "replace-finalized-text-layer",
            "variant": args.variant_name,
            "expected_cutoff": args.expected_cutoff,
            "old_text_carried_forward": False,
            "replacement_semantics": "old text directory removed before corrected text install",
            "base_final_root": str(args.base_final_root.resolve()),
            "base_archive_sha256": base_archive_sha256,
            "base_member_content_sha256": base_extract["member_content_sha256"],
            "base_counts": {
                modality: int(base_extract["members_by_directory"].get(directory, 0))
                for modality, directory in sorted(directories.items())
            },
            "removed_old_text_count": old_text_count,
            "replacement_run_root": str(args.replacement_run_root.resolve()),
            "replacement_archive_sha256": replacement_archive_sha256,
            "combined_archive": "archive/contemporary_embedding_cache.tar.gz",
            "combined_archive_sha256": archive_report["archive_sha256"],
            "combined_member_content_sha256": archive_report["member_content_sha256"],
            "combined_counts": combined_counts,
            "target_count": len(targets),
            "assembly_report": "reports/assembly/embedding_assembly.tsv.gz",
            "assembly_report_sha256": sha256_file(assembly_output),
            "roundtrip_validated": True,
        }
        atomic_write_json(staging / "COMPOSITION_COMPLETE.json", manifest)
        os.replace(staging, args.output_root)
        published = True
        return manifest
    finally:
        if staging.exists() and not published:
            shutil.rmtree(staging, ignore_errors=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-final-root", type=Path, required=True)
    parser.add_argument("--replacement-run-root", type=Path, required=True)
    parser.add_argument("--plan-dir", type=Path, required=True)
    parser.add_argument("--policy", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--input-acquisition", type=Path, required=True)
    parser.add_argument("--expected-cutoff", required=True)
    parser.add_argument("--variant-name", required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    try:
        result = compose(args)
        if args.report:
            atomic_write_json(args.report, result)
    except (OSError, ValueError, KeyError, json.JSONDecodeError) as error:
        raise SystemExit(f"ERROR: {error}") from error
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
