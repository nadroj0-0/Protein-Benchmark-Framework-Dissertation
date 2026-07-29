#!/usr/bin/env python3
"""Add a validated PPI-only delta to a finalized contemporary cache."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import shutil
import tarfile
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


SAFE_ID = re.compile(r"^[A-Za-z0-9_.-]+$")


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


def require_archive_hash(path: Path, expected: object, label: str) -> str:
    if not isinstance(expected, str) or not re.fullmatch(r"[0-9a-f]{64}", expected):
        raise ValueError(f"{label} has no valid SHA-256 declaration")
    observed = sha256_file(path)
    if observed != expected:
        raise ValueError(f"{label} SHA-256 mismatch: {expected} != {observed}")
    return observed


def resolve_declared_file(root: Path, relative_name: object, label: str) -> Path:
    if not isinstance(relative_name, str) or not relative_name:
        raise ValueError(f"{label} does not declare a file")
    relative = Path(relative_name)
    if relative.is_absolute() or ".." in relative.parts:
        raise ValueError(f"{label} declares an unsafe path: {relative_name}")
    resolved_root = root.resolve()
    resolved = (root / relative).resolve()
    if resolved_root not in resolved.parents:
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
    result: dict[str, str] = {}
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


def load_delta_mapping(path: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"protein_id", "npy_sha256"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError("Delta mapping report has an unexpected schema")
        for row in reader:
            protein_id = row["protein_id"]
            digest = row["npy_sha256"]
            if not SAFE_ID.fullmatch(protein_id):
                raise ValueError(f"Unsafe protein ID in delta mapping: {protein_id!r}")
            if protein_id in result:
                raise ValueError(f"Duplicate protein in delta mapping: {protein_id}")
            if not re.fullmatch(r"[0-9a-f]{64}", digest):
                raise ValueError(f"Invalid array SHA-256 for {protein_id}")
            result[protein_id] = digest
    if not result:
        raise ValueError("Delta mapping report is empty")
    return result


def extract_delta(
    archive: Path,
    destination: Path,
    expected_hashes: dict[str, str],
    dimension: int,
) -> set[str]:
    if destination.exists():
        raise ValueError(f"Delta extraction destination already exists: {destination}")
    destination.mkdir(parents=True)
    seen: set[str] = set()
    try:
        with tarfile.open(archive, "r:gz") as handle:
            for member in handle:
                parts = Path(member.name).parts
                if (
                    not member.isfile()
                    or len(parts) != 2
                    or parts[0] != "ppi"
                    or not parts[1].endswith(".npy")
                ):
                    raise ValueError(f"Unexpected delta archive member: {member.name}")
                protein_id = Path(parts[1]).stem
                if not SAFE_ID.fullmatch(protein_id):
                    raise ValueError(f"Unsafe delta archive member: {member.name}")
                if protein_id in seen:
                    raise ValueError(f"Repeated delta archive member: {member.name}")
                if protein_id not in expected_hashes:
                    raise ValueError(f"Delta archive member is absent from mapping: {protein_id}")
                source = handle.extractfile(member)
                if source is None:
                    raise ValueError(f"Could not read delta member: {member.name}")
                output = destination / parts[1]
                descriptor, temporary_name = tempfile.mkstemp(
                    prefix=f".{parts[1]}.", suffix=".tmp", dir=str(destination)
                )
                try:
                    digest = hashlib.sha256()
                    size = 0
                    with source, os.fdopen(descriptor, "wb") as target:
                        for block in iter(lambda: source.read(1024 * 1024), b""):
                            target.write(block)
                            digest.update(block)
                            size += len(block)
                        target.flush()
                        os.fsync(target.fileno())
                    if size != member.size:
                        raise ValueError(f"Delta member size changed: {member.name}")
                    if digest.hexdigest() != expected_hashes[protein_id]:
                        raise ValueError(f"Delta member SHA-256 mismatch: {protein_id}")
                    os.replace(temporary_name, output)
                finally:
                    try:
                        os.unlink(temporary_name)
                    except FileNotFoundError:
                        pass
                validate_array(output, dimension)
                seen.add(protein_id)
    except BaseException:
        shutil.rmtree(destination, ignore_errors=True)
        raise
    if seen != set(expected_hashes):
        missing = sorted(set(expected_hashes) - seen)
        raise ValueError(f"Delta archive omitted mapped arrays: {missing[:5]}")
    return seen


def compose(args: argparse.Namespace) -> dict:
    if args.output_root.exists():
        raise ValueError(f"Output root already exists: {args.output_root}")
    if args.work_dir.exists() and any(args.work_dir.iterdir()):
        raise ValueError(f"Work directory is not empty: {args.work_dir}")
    for path in (args.base_final_root, args.delta_root, args.plan_dir):
        if not path.is_dir():
            raise ValueError(f"Required directory is missing: {path}")
    for path in (args.policy, args.config, args.input_acquisition):
        if not path.is_file():
            raise ValueError(f"Required file is missing: {path}")

    policy = load_policy(args.policy)
    directories = modality_directories(args.config)
    for modality, specification in policy["modalities"].items():
        if directories[modality] != specification["cache_directory"]:
            raise ValueError(f"Policy/config directory mismatch for {modality}")

    target_tables = (
        args.plan_dir / "reuse_proteins.tsv",
        args.plan_dir / "regenerate_proteins.tsv",
    )
    targets = load_target_tables(target_tables)
    if len(targets) != args.expected_target_count:
        raise ValueError(
            f"Target count differs: {len(targets)} != {args.expected_target_count}"
        )

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

    delta_marker_path = args.delta_root / "DELTA_COMPLETE.json"
    delta_marker = load_json(delta_marker_path)
    required_delta_values = {
        "schema_name": "validated-unambiguous-ppi-delta",
        "complete": True,
        "policy": "widened_unambiguous",
        "roundtrip_validated": True,
        "target_count": args.expected_target_count,
        "base_accepted_count": args.expected_base_ppi_count,
        "delta_count": args.expected_delta_count,
        "final_union_count": args.expected_final_ppi_count,
        "base_overlap_count": 0,
    }
    for key, expected in required_delta_values.items():
        if delta_marker.get(key) != expected:
            raise ValueError(
                f"Delta marker mismatch for {key}: {delta_marker.get(key)!r} != {expected!r}"
            )
    delta_archive = resolve_declared_file(
        args.delta_root, delta_marker.get("archive"), "PPI delta archive"
    )
    delta_archive_sha256 = require_archive_hash(
        delta_archive, delta_marker.get("archive_sha256"), "PPI delta archive"
    )
    mapping_report = resolve_declared_file(
        args.delta_root, delta_marker.get("mapping_report"), "PPI delta mapping"
    )
    require_archive_hash(
        mapping_report, delta_marker.get("mapping_report_sha256"), "PPI delta mapping"
    )
    delta_hashes = load_delta_mapping(mapping_report)
    if len(delta_hashes) != args.expected_delta_count:
        raise ValueError(
            f"Delta mapping count differs: {len(delta_hashes)} != {args.expected_delta_count}"
        )
    if not set(delta_hashes).issubset(targets):
        raise ValueError("PPI delta contains proteins outside the benchmark target set")

    args.work_dir.mkdir(parents=True, exist_ok=True)
    base_cache = args.work_dir / "base_cache"
    delta_cache = args.work_dir / "delta_cache"
    roundtrip_cache = args.work_dir / "roundtrip_cache"
    base_extract = extract_archive(base_archive, base_cache, args.config)
    available_before, base_counts = validate_cache(base_cache, targets, policy)
    if base_counts["ppi"] != args.expected_base_ppi_count:
        raise ValueError(
            f"Base PPI count differs: {base_counts['ppi']} != {args.expected_base_ppi_count}"
        )
    if set(delta_hashes) & available_before["ppi"]:
        raise ValueError("PPI delta overlaps the paper-faithful cache")

    delta_ids = extract_delta(
        delta_archive,
        delta_cache,
        delta_hashes,
        int(policy["modalities"]["ppi"]["dimension"]),
    )
    ppi_directory = base_cache / directories["ppi"]
    for protein_id in sorted(delta_ids):
        destination = ppi_directory / f"{protein_id}.npy"
        if destination.exists():
            raise ValueError(f"Refusing to replace existing PPI array: {protein_id}")
        os.replace(delta_cache / f"{protein_id}.npy", destination)

    available_after, combined_counts = validate_cache(base_cache, targets, policy)
    if combined_counts["ppi"] != args.expected_final_ppi_count:
        raise ValueError(
            f"Combined PPI count differs: {combined_counts['ppi']} != "
            f"{args.expected_final_ppi_count}"
        )
    for modality in ("sequence", "text", "structure"):
        if available_after[modality] != available_before[modality]:
            raise ValueError(f"Base {modality} membership changed while adding PPI")

    args.output_root.parent.mkdir(parents=True, exist_ok=True)
    staging = args.output_root.parent / (
        f".{args.output_root.name}.staging-{uuid.uuid4().hex}"
    )
    published = False
    try:
        archive_output = staging / "archive/contemporary_embedding_cache.tar.gz"
        assembly_output = staging / "reports/assembly/embedding_assembly.tsv.gz"
        archive_output.parent.mkdir(parents=True)
        assembly_output.parent.mkdir(parents=True)
        archive_report = create_archive(base_cache, archive_output, args.config)
        atomic_gzip_report(assembly_output, targets, policy["modalities"], available_after)
        shutil.copy2(args.input_acquisition, staging / "reports/input_acquisition.tsv")

        provenance = staging / "provenance"
        provenance.mkdir()
        shutil.copy2(base_marker_path, provenance / "BASE_FINAL_CACHE_COMPLETE.json")
        shutil.copy2(delta_marker_path, provenance / "DELTA_COMPLETE.json")
        shutil.copy2(mapping_report, provenance / "delta_mapping.tsv.gz")
        base_contract = args.base_final_root / "evidence/contract.json"
        if base_contract.is_file():
            shutil.copy2(base_contract, provenance / "BASE_CONTRACT.json")
        base_variant_root = args.base_final_root.parent
        base_variant_marker = base_variant_root / "VARIANT_COMPLETE.json"
        if base_variant_marker.is_file():
            shutil.copy2(
                base_variant_marker, provenance / "BASE_VARIANT_COMPLETE.json"
            )
        corrected_text_marker = (
            base_variant_root
            / "source_baseline/provenance/TEXT_GENERATION_COMPLETE.json"
        )
        if corrected_text_marker.is_file():
            shutil.copy2(
                corrected_text_marker, provenance / "TEXT_GENERATION_COMPLETE.json"
            )

        roundtrip_report = extract_archive(archive_output, roundtrip_cache, args.config)
        for key in (
            "archive_sha256",
            "member_count",
            "members_by_directory",
            "member_content_sha256",
        ):
            if roundtrip_report.get(key) != archive_report.get(key):
                raise ValueError(f"Combined archive round trip differs for {key}")
        roundtrip_available, roundtrip_counts = validate_cache(
            roundtrip_cache, targets, policy
        )
        if roundtrip_counts != combined_counts:
            raise ValueError("Combined archive round trip changed modality counts")
        if roundtrip_available != available_after:
            raise ValueError("Combined archive round trip changed cache membership")

        manifest = {
            "schema_version": 1,
            "complete": True,
            "created_at": utc_now(),
            "operation": "add-validated-ppi-delta",
            "variant": args.variant_name,
            "ppi_policy": "widened-unambiguous",
            "base_final_root": str(args.base_final_root.resolve()),
            "base_archive_sha256": base_archive_sha256,
            "base_member_content_sha256": base_extract["member_content_sha256"],
            "base_counts": base_counts,
            "delta_root": str(args.delta_root.resolve()),
            "delta_archive_sha256": delta_archive_sha256,
            "delta_count": len(delta_ids),
            "replacement_count": 0,
            "combined_archive": "archive/contemporary_embedding_cache.tar.gz",
            "combined_archive_sha256": archive_report["archive_sha256"],
            "combined_member_content_sha256": archive_report[
                "member_content_sha256"
            ],
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
    parser.add_argument("--delta-root", type=Path, required=True)
    parser.add_argument("--plan-dir", type=Path, required=True)
    parser.add_argument("--policy", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--input-acquisition", type=Path, required=True)
    parser.add_argument("--variant-name", required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--expected-target-count", type=int, required=True)
    parser.add_argument("--expected-base-ppi-count", type=int, required=True)
    parser.add_argument("--expected-delta-count", type=int, required=True)
    parser.add_argument("--expected-final-ppi-count", type=int, required=True)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    try:
        result = compose(args)
        if args.report:
            atomic_write_json(args.report, result)
    except (OSError, ValueError, KeyError, json.JSONDecodeError, tarfile.TarError) as error:
        raise SystemExit(f"ERROR: {error}") from error
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
