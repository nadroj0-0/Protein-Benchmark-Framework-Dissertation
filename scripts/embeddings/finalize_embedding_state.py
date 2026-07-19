#!/usr/bin/env python3
"""Transactionally consolidate an archive-backed embedding retry state."""

from __future__ import annotations

import argparse
import csv
import fcntl
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import uuid
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterator


FRAMEWORK_ROOT = Path(__file__).resolve().parents[2]
STATE_MANAGER = FRAMEWORK_ROOT / "scripts/embeddings/manage_resumable_embedding_state.py"
ARCHIVE_MANAGER = FRAMEWORK_ROOT / "scripts/embeddings/manage_embedding_archive.py"
PREPARER = FRAMEWORK_ROOT / "scripts/model_execution/prepare_pfp_benchmark.py"
VALIDATOR = FRAMEWORK_ROOT / "scripts/model_execution/validate_pfp_embedding_cache.py"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_sha256(value: object) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


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


def tail(path: Path, lines: int = 30) -> str:
    try:
        values = path.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError:
        return ""
    return "\n".join(values[-lines:])


def run_logged(label: str, command: list[str], log_path: Path) -> None:
    print(f"==> {label}", flush=True)
    print("Command: " + " ".join(command), flush=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w", encoding="utf-8") as log:
        result = subprocess.run(
            command,
            stdout=log,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )
    if result.returncode != 0:
        excerpt = tail(log_path)
        raise RuntimeError(
            f"{label} failed with status {result.returncode}; see {log_path}\n{excerpt}"
        )


def load_json(path: Path) -> dict:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return value


def source_snapshot(state_root: Path) -> dict[str, str]:
    names = (
        "contract.json",
        "coverage.json",
        "targets.tsv",
        "pair_status.tsv",
        "baseline_accepted.tsv",
        "EVIDENCE_HASHES_COMPLETE.json",
    )
    result: dict[str, str] = {}
    for name in names:
        path = state_root / name
        if not path.is_file():
            raise ValueError(f"Required embedding-state evidence is missing: {path}")
        result[name] = sha256_file(path)
    return result


def accepted_counts(pair_status_path: Path) -> dict[str, int]:
    counts: dict[str, int] = {}
    with pair_status_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"protein_id", "modality", "state", "embedding_sha256"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError("pair_status.tsv has not been upgraded with embedding hashes")
        for row in reader:
            if row["state"] != "accepted":
                continue
            digest = row["embedding_sha256"]
            if len(digest) != 64 or any(char not in "0123456789abcdef" for char in digest):
                raise ValueError(
                    f"Accepted pair lacks a valid embedding hash: "
                    f"{row['protein_id']}/{row['modality']}"
                )
            counts[row["modality"]] = counts.get(row["modality"], 0) + 1
    return dict(sorted(counts.items()))


def build_final_evidence(
    source_state: Path,
    evidence_dir: Path,
    final_root: Path,
    archive_name: str,
    archive_report: dict,
) -> dict:
    source_contract = load_json(source_state / "contract.json")
    source_coverage = load_json(source_state / "coverage.json")
    source_contract_hash = source_contract.get("contract_sha256")
    if not isinstance(source_contract_hash, str):
        raise ValueError("Source embedding contract lacks contract_sha256")

    final_contract = dict(source_contract)
    final_contract.pop("contract_sha256", None)
    retired_baseline = final_contract.pop("baseline", None)
    final_contract["finalized_embedding_archive"] = {
        "name": archive_name,
        "path": str((final_root / archive_name).resolve()),
        "sha256": archive_report["archive_sha256"],
        "size_bytes": archive_report["archive_size_bytes"],
        "member_count": archive_report["member_count"],
        "member_content_sha256": archive_report["member_content_sha256"],
    }
    final_contract["source_retry_state"] = {
        "contract_sha256": source_contract_hash,
        "coverage_sha256": sha256_file(source_state / "coverage.json"),
        "pair_status_sha256": sha256_file(source_state / "pair_status.tsv"),
        "baseline_provenance_retired_after_publication": retired_baseline,
    }
    final_contract["finalized_at"] = utc_now()
    final_contract["contract_sha256"] = canonical_sha256(final_contract)

    final_coverage = dict(source_coverage)
    final_coverage["contract_sha256"] = final_contract["contract_sha256"]
    final_coverage["state_root"] = str((final_root / "evidence").resolve())
    final_coverage["finalized_archive_sha256"] = archive_report["archive_sha256"]

    evidence_dir.mkdir(parents=True, exist_ok=False)
    atomic_write_json(evidence_dir / "contract.json", final_contract)
    atomic_write_json(evidence_dir / "coverage.json", final_coverage)
    shutil.copy2(source_state / "targets.tsv", evidence_dir / "targets.tsv")
    shutil.copy2(source_state / "pair_status.tsv", evidence_dir / "pair_status.tsv")
    return {
        "source_contract_sha256": source_contract_hash,
        "final_contract_sha256": final_contract["contract_sha256"],
        "accepted_counts": accepted_counts(evidence_dir / "pair_status.tsv"),
        "embedding_gate_passed": bool(final_coverage.get("embedding_gate_passed")),
    }


@contextmanager
def exclusive_state_lock(state_root: Path) -> Iterator[None]:
    lock_path = state_root / ".state.lock"
    with lock_path.open("a+", encoding="ascii") as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def retire_source_embeddings(
    state_root: Path,
    expected_snapshot: dict[str, str],
    final_root: Path,
) -> dict:
    validated_marker = final_root / "CACHE_ARCHIVE_VALIDATED.json"
    if not validated_marker.is_file():
        raise ValueError("Refusing source retirement without a published validation marker")
    validated = load_json(validated_marker)
    archive_path = final_root / str(validated.get("archive_name", ""))
    if not archive_path.is_file() or sha256_file(archive_path) != validated.get("archive_sha256"):
        raise ValueError("Published final archive no longer matches its validation marker")

    with exclusive_state_lock(state_root):
        if source_snapshot(state_root) != expected_snapshot:
            raise ValueError("Embedding state changed after hydration; refusing source retirement")
        source_contract = load_json(state_root / "contract.json")
        baseline = source_contract.get("baseline", {})
        baseline_path_text = baseline.get("archive", {}).get("path")
        if not isinstance(baseline_path_text, str):
            raise ValueError("Source contract has no baseline archive path")
        baseline_path = Path(baseline_path_text)
        delta_cache = state_root / "cache"
        if not baseline_path.is_file() or baseline_path.is_symlink():
            raise ValueError(f"Contracted baseline archive is unsafe or missing: {baseline_path}")
        if not delta_cache.is_dir() or delta_cache.is_symlink():
            raise ValueError(f"Retry delta cache is unsafe or missing: {delta_cache}")
        if baseline_path.resolve() == archive_path.resolve():
            raise ValueError("Final archive unexpectedly aliases the source baseline")

        token = uuid.uuid4().hex
        baseline_quarantine = baseline_path.with_name(
            f".{baseline_path.name}.retiring-{token}"
        )
        cache_quarantine = state_root / f".cache.retiring-{token}"
        if baseline_quarantine.exists() or cache_quarantine.exists():
            raise ValueError("Retirement quarantine path already exists")

        os.replace(baseline_path, baseline_quarantine)
        try:
            os.replace(delta_cache, cache_quarantine)
        except BaseException:
            os.replace(baseline_quarantine, baseline_path)
            raise

        try:
            shutil.rmtree(cache_quarantine)
            baseline_quarantine.unlink()
        except BaseException:
            # The validated final archive remains authoritative. Any surviving
            # quarantine path is deliberately conspicuous for manual cleanup.
            raise

        result = {
            "schema_version": 1,
            "retired_at": utc_now(),
            "final_root": str(final_root.resolve()),
            "final_archive_sha256": validated["archive_sha256"],
            "removed_baseline_archive": str(baseline_path),
            "removed_retry_delta_cache": str(delta_cache),
            "source_evidence_retained": True,
        }
        atomic_write_json(state_root / "SOURCE_EMBEDDINGS_RETIRED.json", result)
        return result


def validate_cache(
    label: str,
    cache_root: Path,
    data_dir: Path,
    config: Path,
    preparation_report: Path,
    evidence_dir: Path,
    report: Path,
    issues: Path,
    log: Path,
) -> None:
    command = [
        sys.executable,
        str(VALIDATOR),
        "--data-dir",
        str(data_dir),
        "--cache-root",
        str(cache_root),
        "--config",
        str(config),
        "--mode",
        "full",
        "--report",
        str(report),
        "--issues-tsv",
        str(issues),
        "--preparation-report",
        str(preparation_report),
        "--require-embedding-evidence",
    ]
    for name in ("coverage.json", "contract.json", "targets.tsv", "pair_status.tsv"):
        command.extend(("--embedding-evidence", str(evidence_dir / name)))
    run_logged(label, command, log)
    if load_json(report).get("status") != "passed":
        raise ValueError(f"{label} report does not declare status=passed")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--state-root", type=Path, required=True)
    parser.add_argument("--benchmark-dir", type=Path, required=True)
    parser.add_argument("--obo-file", type=Path, required=True)
    parser.add_argument("--pfp-root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--final-root", type=Path, required=True)
    parser.add_argument("--report-dir", type=Path, required=True)
    parser.add_argument("--archive-name", default="contemporary_embedding_cache.tar.gz")
    parser.add_argument("--confirm-retries-finished", action="store_true")
    parser.add_argument("--retire-source-embeddings", action="store_true")
    args = parser.parse_args()

    if not args.confirm_retries_finished:
        raise SystemExit("ERROR: --confirm-retries-finished is required")
    if not args.retire_source_embeddings:
        raise SystemExit("ERROR: --retire-source-embeddings is required")
    if Path(args.archive_name).name != args.archive_name:
        raise SystemExit("ERROR: --archive-name must be a plain filename")
    for path in (args.state_root, args.benchmark_dir, args.pfp_root):
        if not path.is_dir():
            raise SystemExit(f"ERROR: Required directory is missing: {path}")
    for path in (args.obo_file, args.config):
        if not path.is_file():
            raise SystemExit(f"ERROR: Required file is missing: {path}")
    if args.final_root.exists():
        raise SystemExit(f"ERROR: Final root already exists: {args.final_root}")
    if args.work_dir.exists() and any(args.work_dir.iterdir()):
        raise SystemExit(f"ERROR: Work directory is not empty: {args.work_dir}")
    if args.report_dir.exists() and any(args.report_dir.iterdir()):
        raise SystemExit(f"ERROR: Report directory is not empty: {args.report_dir}")

    args.work_dir.mkdir(parents=True, exist_ok=True)
    args.report_dir.mkdir(parents=True, exist_ok=True)
    logs = args.report_dir / "logs"
    reports = args.report_dir / "reports"
    logs.mkdir()
    reports.mkdir()
    staging = args.final_root.parent / (
        f".{args.final_root.name}.staging-{uuid.uuid4().hex}"
    )
    published = False

    try:
        prepared = args.work_dir / "prepared_data"
        preparation_report = reports / "preparation.json"
        run_logged(
            "[1/8] Validate and prepare the benchmark with its frozen ontology",
            [
                sys.executable,
                str(PREPARER),
                "--benchmark-dir",
                str(args.benchmark_dir),
                "--data-dir",
                str(prepared),
                "--obo-file",
                str(args.obo_file),
                "--pfp-root",
                str(args.pfp_root),
                "--config",
                str(args.config),
                "--report",
                str(preparation_report),
                "--log-dir",
                str(logs / "preparation"),
            ],
            logs / "01_prepare.log",
        )

        upgrade_report = reports / "evidence_hash_upgrade.json"
        run_logged(
            "[2/8] Upgrade per-array evidence hashes without changing membership",
            [
                sys.executable,
                str(STATE_MANAGER),
                "upgrade-evidence-hashes",
                "--state-root",
                str(args.state_root),
                "--report",
                str(upgrade_report),
            ],
            logs / "02_evidence_upgrade.log",
        )
        upgrade = load_json(upgrade_report)
        if upgrade.get("accepted_membership_unchanged") is not True:
            raise ValueError("Evidence upgrade did not prove accepted membership unchanged")
        snapshot = source_snapshot(args.state_root)

        hydrated = args.work_dir / "hydrated_cache"
        hydrate_report = reports / "hydrate.json"
        run_logged(
            "[3/8] Hydrate the accepted baseline and retry delta into scratch",
            [
                sys.executable,
                str(STATE_MANAGER),
                "hydrate",
                "--state-root",
                str(args.state_root),
                "--output-cache-root",
                str(hydrated),
                "--report",
                str(hydrate_report),
                "--preserve-evidence",
            ],
            logs / "03_hydrate.log",
        )
        if source_snapshot(args.state_root) != snapshot:
            raise ValueError("Embedding state changed during hydration")

        validate_cache(
            "[4/8] Exhaustively validate the hydrated cache and source evidence",
            hydrated,
            prepared,
            args.config,
            preparation_report,
            args.state_root,
            reports / "hydrated_cache_validation.json",
            reports / "hydrated_cache_issues.tsv",
            logs / "04_validate_hydrated.log",
        )

        scratch_archive = args.work_dir / args.archive_name
        archive_report_path = reports / "archive_creation.json"
        run_logged(
            "[5/8] Create one consolidated embedding archive in scratch",
            [
                sys.executable,
                str(ARCHIVE_MANAGER),
                "create",
                "--cache-root",
                str(hydrated),
                "--archive",
                str(scratch_archive),
                "--config",
                str(args.config),
                "--report",
                str(archive_report_path),
            ],
            logs / "05_create_archive.log",
        )
        archive_report = load_json(archive_report_path)

        args.final_root.parent.mkdir(parents=True, exist_ok=True)
        staging.mkdir()
        shutil.copy2(scratch_archive, staging / args.archive_name)
        if sha256_file(staging / args.archive_name) != archive_report["archive_sha256"]:
            raise ValueError("Archive SHA-256 changed while copying to persistent storage")
        evidence_summary = build_final_evidence(
            args.state_root,
            staging / "evidence",
            args.final_root,
            args.archive_name,
            archive_report,
        )

        roundtrip = args.work_dir / "roundtrip_cache"
        extraction_report_path = reports / "archive_extraction.json"
        run_logged(
            "[6/8] Read the copied SAN archive back into fresh scratch",
            [
                sys.executable,
                str(ARCHIVE_MANAGER),
                "extract",
                "--archive",
                str(staging / args.archive_name),
                "--output-cache-root",
                str(roundtrip),
                "--config",
                str(args.config),
                "--report",
                str(extraction_report_path),
            ],
            logs / "06_extract_roundtrip.log",
        )
        extraction_report = load_json(extraction_report_path)
        for key in ("archive_sha256", "member_count", "members_by_directory", "member_content_sha256"):
            if extraction_report.get(key) != archive_report.get(key):
                raise ValueError(f"Round-trip archive report differs for {key}")

        validate_cache(
            "[7/8] Revalidate every array read back from the copied SAN archive",
            roundtrip,
            prepared,
            args.config,
            preparation_report,
            staging / "evidence",
            reports / "roundtrip_cache_validation.json",
            reports / "roundtrip_cache_issues.tsv",
            logs / "07_validate_roundtrip.log",
        )
        if source_snapshot(args.state_root) != snapshot:
            raise ValueError("Embedding state changed before final publication")

        shutil.copytree(reports, staging / "reports")
        validation_marker = {
            "schema_version": 1,
            "validated": True,
            "validated_at": utc_now(),
            "archive_name": args.archive_name,
            "archive_sha256": archive_report["archive_sha256"],
            "archive_size_bytes": archive_report["archive_size_bytes"],
            "member_count": archive_report["member_count"],
            "accepted_counts": evidence_summary["accepted_counts"],
            "missing_nonsequence_policy": "PFP native zero-vector plus availability mask",
            "embedding_gate_passed": evidence_summary["embedding_gate_passed"],
            "source_contract_sha256": evidence_summary["source_contract_sha256"],
            "final_contract_sha256": evidence_summary["final_contract_sha256"],
        }
        atomic_write_json(staging / "CACHE_ARCHIVE_VALIDATED.json", validation_marker)
        os.replace(staging, args.final_root)
        published = True
        if sha256_file(args.final_root / args.archive_name) != archive_report["archive_sha256"]:
            raise ValueError("Published final archive failed its post-rename checksum")

        print("==> [8/8] Retire only the superseded baseline and retry-delta embedding bytes", flush=True)
        retirement = retire_source_embeddings(args.state_root, snapshot, args.final_root)
        final_marker = dict(validation_marker)
        final_marker.update(
            {
                "complete": True,
                "completed_at": utc_now(),
                "source_embeddings_retired": True,
                "retirement": retirement,
            }
        )
        atomic_write_json(args.final_root / "FINAL_CACHE_COMPLETE.json", final_marker)
        atomic_write_json(args.report_dir / "finalization_report.json", final_marker)
        print(f"Final embedding archive: {args.final_root / args.archive_name}")
        print("The source evidence remains on SAN; duplicate embedding bytes were removed.")
        return 0
    except (OSError, ValueError, RuntimeError, KeyError, json.JSONDecodeError) as error:
        atomic_write_json(
            args.report_dir / "FINALIZATION_FAILED.json",
            {
                "complete": False,
                "failed_at": utc_now(),
                "error": str(error),
                "final_archive_published": published,
                "source_embeddings_retired": False,
            },
        )
        raise SystemExit(f"ERROR: {error}") from error
    finally:
        if staging.exists() and not published:
            shutil.rmtree(staging, ignore_errors=True)


if __name__ == "__main__":
    raise SystemExit(main())
