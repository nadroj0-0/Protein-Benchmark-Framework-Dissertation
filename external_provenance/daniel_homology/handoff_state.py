#!/usr/bin/env python3
"""Streaming validation and provenance state for the MMseqs handoff bundle."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_sha256(value: Any) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def write_json_atomic(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".partial.{os.getpid()}")
    with temporary.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def load_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def parse_thresholds(value: str) -> list[int]:
    thresholds = [int(item) for item in value.split(",") if item]
    if not thresholds or len(thresholds) != len(set(thresholds)):
        raise ValueError("threshold list must be nonempty and contain no duplicates")
    return thresholds


def contract_from_args(args: argparse.Namespace) -> dict[str, Any]:
    scientific = {
        "input_sha256": args.input_sha256,
        "input_release": args.input_release,
        "thresholds_percent": parse_thresholds(args.thresholds),
        "coverage": 0.8,
        "coverage_mode": 0,
        "cluster_mode": 0,
        "alignment_mode": 3,
        "sequence_identity_mode": 0,
        "sensitivity": 7.5,
        "createdb_shuffle": "MMseqs default",
        "evalue": "MMseqs default",
        "cluster_reassignment": bool(args.cluster_reassign),
        "clustering_backtrace": False,
        "alignment_statistics_backtrace_when_exported": bool(
            args.export_alignment_statistics
        ),
        "export_alignment_statistics": bool(args.export_alignment_statistics),
        "export_cluster_fasta": bool(args.export_cluster_fasta),
    }
    execution = {
        "mmseqs_executable_sha256": args.mmseqs_sha256,
        "mmseqs_version": args.mmseqs_version,
        "threads": args.threads,
    }
    return {
        "schema_version": 1,
        "scientific_contract": scientific,
        "execution_contract": execution,
        "contract_sha256": canonical_sha256({"scientific": scientific, "execution": execution}),
    }


def initialize_run(args: argparse.Namespace) -> int:
    path = Path(args.contract)
    candidate = contract_from_args(args)
    candidate.update(
        {
            "created_at_utc": utc_now(),
            "host": platform.node(),
            "platform": platform.platform(),
            "input_path": str(Path(args.input_path).resolve()),
            "mmseqs_path": str(Path(args.mmseqs_path).resolve()),
        }
    )
    if path.exists():
        existing = load_json(path)
        if existing.get("contract_sha256") != candidate["contract_sha256"]:
            print(
                "Existing output directory is bound to a different input, MMseqs binary, "
                "thread count, threshold set, or scientific configuration.",
                file=sys.stderr,
            )
            print(f"Existing: {existing.get('contract_sha256')}", file=sys.stderr)
            print(f"Requested: {candidate['contract_sha256']}", file=sys.stderr)
            return 2
        print(f"Reusing matching run contract: {candidate['contract_sha256']}")
        return 0
    write_json_atomic(path, candidate)
    print(f"Created run contract: {candidate['contract_sha256']}")
    return 0


def db_files(db_base: Path) -> list[Path]:
    return sorted(path for path in db_base.parent.glob(db_base.name + "*") if path.is_file())


def write_createdb_marker(args: argparse.Namespace) -> int:
    db_base = Path(args.db_base)
    files = db_files(db_base)
    if not files or not any(path.stat().st_size > 0 for path in files):
        print(f"No nonempty MMseqs database files found for {db_base}", file=sys.stderr)
        return 2
    contract = load_json(Path(args.contract))
    marker = {
        "schema_version": 1,
        "stage": "createdb",
        "completed_at_utc": utc_now(),
        "run_contract_sha256": contract["contract_sha256"],
        "database_base": str(db_base.resolve()),
        "database_files": [
            {"name": path.name, "size_bytes": path.stat().st_size} for path in files
        ],
    }
    write_json_atomic(Path(args.marker), marker)
    return 0


def check_createdb_marker(args: argparse.Namespace) -> int:
    marker_path = Path(args.marker)
    if not marker_path.exists():
        return 1
    try:
        marker = load_json(marker_path)
        contract = load_json(Path(args.contract))
        if marker.get("run_contract_sha256") != contract.get("contract_sha256"):
            return 1
        files = db_files(Path(args.db_base))
        if not files or not any(path.stat().st_size > 0 for path in files):
            return 1
    except (OSError, ValueError, KeyError, json.JSONDecodeError):
        return 1
    return 0


def write_cluster_marker(args: argparse.Namespace) -> int:
    cluster_base = Path(args.cluster_base)
    files = db_files(cluster_base)
    if not files or not any(path.stat().st_size > 0 for path in files):
        print(f"No nonempty MMseqs cluster database files found for {cluster_base}", file=sys.stderr)
        return 2
    contract = load_json(Path(args.contract))
    marker = {
        "schema_version": 1,
        "stage": "cluster",
        "completed_at_utc": utc_now(),
        "run_contract_sha256": contract["contract_sha256"],
        "identity_percentage": int(args.threshold),
        "cluster_database_base": str(cluster_base.resolve()),
        "cluster_database_files": [
            {"name": path.name, "size_bytes": path.stat().st_size} for path in files
        ],
    }
    write_json_atomic(Path(args.marker), marker)
    return 0


def check_cluster_marker(args: argparse.Namespace) -> int:
    marker_path = Path(args.marker)
    if not marker_path.exists():
        return 1
    try:
        marker = load_json(marker_path)
        contract = load_json(Path(args.contract))
        if marker.get("run_contract_sha256") != contract.get("contract_sha256"):
            return 1
        if int(marker.get("identity_percentage")) != int(args.threshold):
            return 1
        files = db_files(Path(args.cluster_base))
        if not files or not any(path.stat().st_size > 0 for path in files):
            return 1
    except (OSError, ValueError, KeyError, TypeError, json.JSONDecodeError):
        return 1
    return 0


def validate_tsv(args: argparse.Namespace) -> int:
    input_path = Path(args.input)
    expected_rows = int(args.expected_rows)
    rows = 0
    self_assignments = 0
    malformed = 0
    first_error: dict[str, Any] | None = None
    with input_path.open("rb") as handle:
        for line_number, raw_line in enumerate(handle, 1):
            line = raw_line.rstrip(b"\r\n")
            fields = line.split(b"\t")
            if len(fields) != 2 or not fields[0] or not fields[1]:
                malformed += 1
                if first_error is None:
                    first_error = {
                        "line": line_number,
                        "preview": raw_line[:160].decode("utf-8", errors="replace"),
                    }
                continue
            rows += 1
            if fields[0] == fields[1]:
                self_assignments += 1

    passed = malformed == 0 and rows == expected_rows and self_assignments > 0
    report = {
        "schema_version": 1,
        "validated_at_utc": utc_now(),
        "input": str(input_path.resolve()),
        "size_bytes": input_path.stat().st_size,
        "expected_rows": expected_rows,
        "valid_rows": rows,
        "malformed_rows": malformed,
        "representative_self_assignments": self_assignments,
        "first_error": first_error,
        "passed": passed,
        "scope_note": (
            "Streaming handoff validation checks syntax, exact assignment-row count, and "
            "the presence of representative self-assignments. The benchmark framework performs "
            "the authoritative full membership and cache-contract validation after import."
        ),
    }
    write_json_atomic(Path(args.report), report)
    if not passed:
        print(json.dumps(report, indent=2, sort_keys=True), file=sys.stderr)
        return 2
    print(
        f"Validated {rows:,} cluster assignments; "
        f"{self_assignments:,} representative self-assignments."
    )
    return 0


def write_threshold_marker(args: argparse.Namespace) -> int:
    assignments = Path(args.assignments)
    report = load_json(Path(args.validation_report))
    if not report.get("passed"):
        print("Validation report is not passing", file=sys.stderr)
        return 2
    with gzip.open(assignments, "rb") as handle:
        while handle.read(8 * 1024 * 1024):
            pass
    contract = load_json(Path(args.contract))
    optional_exports: dict[str, Any] = {}
    for name, supplied_path in (
        ("alignment_statistics", args.alignment_statistics),
        ("cluster_fasta", args.cluster_fasta),
    ):
        if not supplied_path:
            continue
        export_path = Path(supplied_path)
        if not export_path.is_file() or export_path.stat().st_size == 0:
            print(f"Missing or empty optional export: {export_path}", file=sys.stderr)
            return 2
        with gzip.open(export_path, "rb") as handle:
            while handle.read(8 * 1024 * 1024):
                pass
        optional_exports[name] = {
            "file": export_path.name,
            "size_bytes": export_path.stat().st_size,
            "sha256": sha256_file(export_path),
        }

    marker = {
        "schema_version": 1,
        "stage": "cluster_assignments",
        "completed_at_utc": utc_now(),
        "run_contract_sha256": contract["contract_sha256"],
        "identity_percentage": int(args.threshold),
        "assignments_file": assignments.name,
        "assignments_size_bytes": assignments.stat().st_size,
        "assignments_sha256": sha256_file(assignments),
        "assignment_rows": report["valid_rows"],
        "representative_self_assignments": report["representative_self_assignments"],
        "validation_report": Path(args.validation_report).name,
        "optional_exports": optional_exports,
    }
    write_json_atomic(Path(args.marker), marker)
    return 0


def check_threshold_marker(args: argparse.Namespace) -> int:
    marker_path = Path(args.marker)
    if not marker_path.exists():
        return 1
    try:
        marker = load_json(marker_path)
        contract = load_json(Path(args.contract))
        assignments = Path(args.assignments)
        if marker.get("run_contract_sha256") != contract.get("contract_sha256"):
            return 1
        if int(marker.get("identity_percentage")) != int(args.threshold):
            return 1
        if not assignments.is_file():
            return 1
        if assignments.stat().st_size != int(marker.get("assignments_size_bytes")):
            return 1
        if sha256_file(assignments) != marker.get("assignments_sha256"):
            return 1
        for export in marker.get("optional_exports", {}).values():
            export_path = marker_path.parent / export["file"]
            if not export_path.is_file():
                return 1
            if export_path.stat().st_size != int(export["size_bytes"]):
                return 1
            if sha256_file(export_path) != export["sha256"]:
                return 1
    except (OSError, ValueError, KeyError, TypeError, json.JSONDecodeError):
        return 1
    return 0


def write_run_manifest(args: argparse.Namespace) -> int:
    output_dir = Path(args.output_dir)
    thresholds = parse_thresholds(args.thresholds)
    contract = load_json(Path(args.contract))
    results = []
    for threshold in thresholds:
        result_dir = output_dir / f"threshold_{threshold:02d}"
        marker_path = result_dir / "COMPLETE.json"
        if not marker_path.is_file():
            print(f"Missing completion marker: {marker_path}", file=sys.stderr)
            return 2
        marker = load_json(marker_path)
        assignments = result_dir / marker["assignments_file"]
        if sha256_file(assignments) != marker["assignments_sha256"]:
            print(f"Assignment checksum mismatch: {assignments}", file=sys.stderr)
            return 2
        results.append(marker)
        for export in marker.get("optional_exports", {}).values():
            export_path = result_dir / export["file"]
            if sha256_file(export_path) != export["sha256"]:
                print(f"Optional export checksum mismatch: {export_path}", file=sys.stderr)
                return 2
    manifest = {
        "schema_version": 1,
        "status": "COMPLETE",
        "completed_at_utc": utc_now(),
        "run_contract": contract,
        "results": results,
        "handoff_note": (
            "These are representative/member assignments from the isolated MMseqs stage. "
            "Benchmark retention, splitting, labels, and final publication occur after import."
        ),
    }
    write_json_atomic(output_dir / "RUN_MANIFEST.json", manifest)
    return 0


def verify_run(args: argparse.Namespace) -> int:
    output_dir = Path(args.output_dir)
    manifest_path = output_dir / "RUN_MANIFEST.json"
    manifest = load_json(manifest_path)
    if manifest.get("status") != "COMPLETE":
        print("Run manifest is not complete", file=sys.stderr)
        return 2
    for result in manifest.get("results", []):
        threshold = int(result["identity_percentage"])
        result_dir = output_dir / f"threshold_{threshold:02d}"
        assignments = result_dir / result["assignments_file"]
        observed = sha256_file(assignments)
        if observed != result["assignments_sha256"]:
            print(f"Checksum mismatch: {assignments}", file=sys.stderr)
            return 2
        with gzip.open(assignments, "rb") as handle:
            while handle.read(8 * 1024 * 1024):
                pass
        for export in result.get("optional_exports", {}).values():
            export_path = result_dir / export["file"]
            if sha256_file(export_path) != export["sha256"]:
                print(f"Checksum mismatch: {export_path}", file=sys.stderr)
                return 2
            with gzip.open(export_path, "rb") as handle:
                while handle.read(8 * 1024 * 1024):
                    pass
    print(f"Verified {len(manifest['results'])} completed threshold result(s).")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)

    init = subparsers.add_parser("initialize-run")
    init.add_argument("--contract", required=True)
    init.add_argument("--input-path", required=True)
    init.add_argument("--input-sha256", required=True)
    init.add_argument("--input-release", required=True)
    init.add_argument("--mmseqs-path", required=True)
    init.add_argument("--mmseqs-sha256", required=True)
    init.add_argument("--mmseqs-version", required=True)
    init.add_argument("--threads", required=True, type=int)
    init.add_argument("--thresholds", required=True)
    init.add_argument("--cluster-reassign", action="store_true")
    init.add_argument("--export-alignment-statistics", action="store_true")
    init.add_argument("--export-cluster-fasta", action="store_true")
    init.set_defaults(func=initialize_run)

    write_db = subparsers.add_parser("write-createdb-marker")
    write_db.add_argument("--contract", required=True)
    write_db.add_argument("--db-base", required=True)
    write_db.add_argument("--marker", required=True)
    write_db.set_defaults(func=write_createdb_marker)

    check_db = subparsers.add_parser("check-createdb-marker")
    check_db.add_argument("--contract", required=True)
    check_db.add_argument("--db-base", required=True)
    check_db.add_argument("--marker", required=True)
    check_db.set_defaults(func=check_createdb_marker)

    write_cluster = subparsers.add_parser("write-cluster-marker")
    write_cluster.add_argument("--contract", required=True)
    write_cluster.add_argument("--threshold", required=True, type=int)
    write_cluster.add_argument("--cluster-base", required=True)
    write_cluster.add_argument("--marker", required=True)
    write_cluster.set_defaults(func=write_cluster_marker)

    check_cluster = subparsers.add_parser("check-cluster-marker")
    check_cluster.add_argument("--contract", required=True)
    check_cluster.add_argument("--threshold", required=True, type=int)
    check_cluster.add_argument("--cluster-base", required=True)
    check_cluster.add_argument("--marker", required=True)
    check_cluster.set_defaults(func=check_cluster_marker)

    validate = subparsers.add_parser("validate-tsv")
    validate.add_argument("--input", required=True)
    validate.add_argument("--expected-rows", required=True, type=int)
    validate.add_argument("--report", required=True)
    validate.set_defaults(func=validate_tsv)

    write_threshold = subparsers.add_parser("write-threshold-marker")
    write_threshold.add_argument("--contract", required=True)
    write_threshold.add_argument("--threshold", required=True, type=int)
    write_threshold.add_argument("--assignments", required=True)
    write_threshold.add_argument("--validation-report", required=True)
    write_threshold.add_argument("--marker", required=True)
    write_threshold.add_argument("--alignment-statistics")
    write_threshold.add_argument("--cluster-fasta")
    write_threshold.set_defaults(func=write_threshold_marker)

    check_threshold = subparsers.add_parser("check-threshold-marker")
    check_threshold.add_argument("--contract", required=True)
    check_threshold.add_argument("--threshold", required=True, type=int)
    check_threshold.add_argument("--assignments", required=True)
    check_threshold.add_argument("--marker", required=True)
    check_threshold.set_defaults(func=check_threshold_marker)

    manifest = subparsers.add_parser("write-run-manifest")
    manifest.add_argument("--output-dir", required=True)
    manifest.add_argument("--contract", required=True)
    manifest.add_argument("--thresholds", required=True)
    manifest.set_defaults(func=write_run_manifest)

    verify = subparsers.add_parser("verify-run")
    verify.add_argument("--output-dir", required=True)
    verify.set_defaults(func=verify_run)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
