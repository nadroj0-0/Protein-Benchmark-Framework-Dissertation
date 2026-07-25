#!/usr/bin/env python3
"""Run PFP's temporal text recipe with a configurable historical cutoff."""

from __future__ import annotations

import argparse
from datetime import datetime, timedelta
import functools
import hashlib
import importlib.util
import json
import shutil
import subprocess
import sys
import threading
from pathlib import Path


def load_module(path: Path):
    specification = importlib.util.spec_from_file_location("pfp_extract_uniprot_text", path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"Cannot import PFP text extractor: {path}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pfp-root", type=Path, required=True)
    parser.add_argument("--cafa-assessment-dir", type=Path, required=True)
    parser.add_argument("--cutoff-date", required=True, help="Historical cutoff as YYYY-MM-DD")
    parser.add_argument("--workers", type=int, default=5)
    return parser.parse_args()


def validate_cutoff_date(value: str) -> str:
    try:
        parsed = datetime.strptime(value, "%Y-%m-%d")
    except ValueError as error:
        raise ValueError(f"Invalid historical cutoff date: {value}") from error
    if parsed.strftime("%Y-%m-%d") != value:
        raise ValueError(f"Historical cutoff date is not canonical: {value}")
    return value


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_ids(values: list[str]) -> str:
    payload = "".join(f"{value}\n" for value in sorted(values)).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def sha256_tree(paths: list[Path]) -> dict[str, object]:
    """Hash the names and contents of every regular file below paths."""
    files: list[tuple[str, Path]] = []
    for root in paths:
        if not root.is_dir():
            continue
        files.extend(
            (f"{root.name}/{path.relative_to(root).as_posix()}", path)
            for path in root.rglob("*")
            if path.is_file()
        )
    digest = hashlib.sha256()
    for relative, path in sorted(files):
        digest.update(relative.encode("utf-8"))
        digest.update(b"\0")
        digest.update(sha256_file(path).encode("ascii"))
        digest.update(b"\n")
    return {"file_count": len(files), "sha256": digest.hexdigest()}


def git_revision(path: Path) -> str | None:
    result = subprocess.run(
        ["git", "-C", str(path), "rev-parse", "HEAD"],
        capture_output=True,
        text=True,
        check=False,
        timeout=10,
    )
    if result.returncode != 0:
        return None
    revision = result.stdout.strip()
    return revision if len(revision) == 40 else None


def configure_historical_cutoff(module, cutoff_date: str) -> dict[str, object]:
    """Make PFP's implicit selector use the caller's explicit cutoff."""
    original = module.find_historical_version

    @functools.wraps(original)
    def configured(versions, cutoff_date=cutoff_date):
        if cutoff_date != configured.framework_effective_cutoff:
            raise ValueError(
                "PFP historical selector received a cutoff different from the "
                f"framework contract: {cutoff_date} != "
                f"{configured.framework_effective_cutoff}"
            )
        return original(versions, cutoff_date=cutoff_date)

    configured.framework_effective_cutoff = cutoff_date
    module.find_historical_version = configured
    module.CUTOFF_DATE = cutoff_date

    cutoff = datetime.strptime(cutoff_date, "%Y-%m-%d")
    probe = [
        {
            "firstReleaseDate": (cutoff - timedelta(days=1)).strftime("%d-%b-%Y"),
            "entryVersion": 991,
        },
        {
            "firstReleaseDate": (cutoff + timedelta(days=1)).strftime("%d-%b-%Y"),
            "entryVersion": 992,
        },
    ]
    expected = original(probe, cutoff_date=cutoff_date)
    observed = module.find_historical_version(probe)
    if expected != 991 or observed != expected:
        raise RuntimeError(
            "Historical cutoff probe failed: "
            f"requested={cutoff_date} expected=991 observed={observed}"
        )
    return {
        "strategy": "framework-wrapper-explicit-cutoff",
        "requested_cutoff": cutoff_date,
        "effective_cutoff": configured.framework_effective_cutoff,
        "probe_entry_version": observed,
    }


def atomic_write_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    temporary.replace(path)


def ensure_state_contract(state_dir: Path, payload: dict[str, object]) -> Path:
    state_dir.mkdir(parents=True, exist_ok=True)
    path = state_dir / "historical_state_contract.json"
    if path.is_file():
        observed = json.loads(path.read_text(encoding="utf-8"))
        if observed != payload:
            raise ValueError(
                "Historical text state contract changed; use a fresh cutoff-scoped "
                f"state directory: {state_dir}"
            )
        return path
    unexpected = [candidate for candidate in state_dir.iterdir() if candidate != path]
    if unexpected:
        raise ValueError(
            "Historical text state predates its contract; refusing unsafe resume: "
            f"{state_dir}"
        )
    atomic_write_json(path, payload)
    return path


def configure_version_provenance(module, output_path: Path) -> None:
    """Record the UniSave entry version selected for each historical protein."""
    records: dict[tuple[str, str], int] = {}
    if output_path.is_file():
        for line_number, line in enumerate(
            output_path.read_text(encoding="utf-8").splitlines(), start=1
        ):
            try:
                row = json.loads(line)
                key = (str(row["protein_id"]), str(row["accession"]))
                version = int(row["entry_version"])
            except (KeyError, TypeError, ValueError, json.JSONDecodeError) as error:
                raise ValueError(
                    f"Invalid UniSave selection record at {output_path}:{line_number}"
                ) from error
            existing = records.get(key)
            if existing is not None and existing != version:
                raise ValueError(
                    f"Conflicting UniSave versions for {key}: {existing} != {version}"
                )
            records[key] = version

    original_find = module.find_historical_version
    original_process = module.process_single_historical_protein
    context = threading.local()
    lock = threading.Lock()

    @functools.wraps(original_find)
    def recorded_find(*args, **kwargs):
        version = original_find(*args, **kwargs)
        protein_id = getattr(context, "protein_id", None)
        accession = getattr(context, "accession", None)
        if version is None or protein_id is None or accession is None:
            return version
        key = (protein_id, accession)
        with lock:
            existing = records.get(key)
            if existing is not None and existing != version:
                raise ValueError(
                    f"UniSave version changed within one state for {key}: "
                    f"{existing} != {version}"
                )
            if existing is None:
                output_path.parent.mkdir(parents=True, exist_ok=True)
                with output_path.open("a", encoding="utf-8") as handle:
                    handle.write(
                        json.dumps(
                            {
                                "accession": accession,
                                "entry_version": version,
                                "protein_id": protein_id,
                                "raw_file": str(context.raw_path),
                            },
                            sort_keys=True,
                        )
                        + "\n"
                    )
                    handle.flush()
                records[key] = version
        return version

    @functools.wraps(original_process)
    def recorded_process(protein_id, accession, session, raw_dir):
        raw_path = Path(raw_dir) / f"{protein_id}_{accession}.txt"
        key = (str(protein_id), str(accession))
        if raw_path.is_file() and key not in records:
            raise ValueError(
                "Historical raw record has no selected-version provenance: "
                f"{raw_path}"
            )
        context.protein_id = str(protein_id)
        context.accession = str(accession)
        context.raw_path = raw_path.resolve()
        try:
            return original_process(protein_id, accession, session, raw_dir)
        finally:
            context.protein_id = None
            context.accession = None
            context.raw_path = None

    module.find_historical_version = recorded_find
    module.process_single_historical_protein = recorded_process


def audit_version_provenance(
    historical_tsv: Path, raw_dir: Path, selected_versions: Path
) -> dict[str, object]:
    records = []
    if selected_versions.is_file():
        records = [
            json.loads(line)
            for line in selected_versions.read_text(encoding="utf-8").splitlines()
        ]
    recorded_ids = {str(row["protein_id"]) for row in records}
    recorded_raw = {Path(str(row["raw_file"])).resolve() for row in records}
    observed_raw = (
        {path.resolve() for path in raw_dir.glob("*.txt")} if raw_dir.is_dir() else set()
    )
    orphaned_raw = sorted(str(path) for path in observed_raw - recorded_raw)
    if orphaned_raw:
        raise ValueError(
            "Historical raw records lack selected-version provenance: "
            f"{orphaned_raw[:5]} (total={len(orphaned_raw)})"
        )

    successful_ids = set()
    if historical_tsv.is_file():
        for line_number, line in enumerate(
            historical_tsv.read_text(encoding="utf-8").splitlines(), start=1
        ):
            fields = line.split("\t", 1)
            if len(fields) != 2 or not fields[0]:
                raise ValueError(
                    f"Invalid historical description row at {historical_tsv}:{line_number}"
                )
            successful_ids.add(fields[0])
    missing_success_provenance = sorted(successful_ids - recorded_ids)
    if missing_success_provenance:
        raise ValueError(
            "Successful historical descriptions lack selected-version provenance: "
            f"{missing_success_provenance[:5]} "
            f"(total={len(missing_success_provenance)})"
        )
    return {
        "recorded_selection_count": len(records),
        "recorded_protein_count": len(recorded_ids),
        "raw_record_count": len(observed_raw),
        "successful_description_count": len(successful_ids),
        "all_raw_records_bound": True,
        "all_successful_descriptions_bound": True,
    }


def main() -> int:
    args = parse_args()
    cutoff_date = validate_cutoff_date(args.cutoff_date)
    pfp_root = args.pfp_root.resolve()
    data_dir = pfp_root / "data"
    text_dir = data_dir / "embedding_cache" / "uniprot_text"
    temporal_dir = text_dir / "temporal_recipe"
    current = text_dir / "protein_descriptions.tsv"
    current_checkpoint = text_dir / "processed_checkpoint.txt"
    historical_state = temporal_dir / f"cutoff_{cutoff_date}"
    historical = historical_state / "protein_descriptions_historical.tsv"
    historical_checkpoint = historical_state / "historical_checkpoint.txt"
    historical_raw = historical_state / "historical_raw"
    punct = historical_state / "protein_descriptions_historical_punct_v1_test.tsv"
    mixed = historical_state / "protein_descriptions_mixed.tsv"

    script = pfp_root / "scripts" / "extract_uniprot_text.py"
    if not script.is_file():
        raise SystemExit(f"Missing PFP text extractor: {script}")
    if not args.cafa_assessment_dir.is_dir():
        raise SystemExit(f"Missing CAFA assessment directory: {args.cafa_assessment_dir}")

    module = load_module(script)
    cutoff_contract = configure_historical_cutoff(module, cutoff_date)
    test_ids = sorted(module.get_split_protein_ids(data_dir, splits=["test"]))
    mapping_root = args.cafa_assessment_dir.resolve() / "ID_conversion"
    mapping_contract = sha256_tree(
        [mapping_root / "CAFA_mapping", mapping_root / "uniprot_mapping"]
    )
    state_contract = {
        "schema_version": 1,
        "requested_cutoff": cutoff_date,
        "effective_cutoff": cutoff_contract["effective_cutoff"],
        "test_protein_count": len(test_ids),
        "test_protein_ids_sha256": sha256_ids(test_ids),
        "pfp_text_script": str(script),
        "pfp_text_script_sha256": sha256_file(script),
        "pfp_revision": git_revision(pfp_root),
        "mapping_inputs": mapping_contract,
        "data_dir": str(data_dir),
        "cafa_assessment_dir": str(args.cafa_assessment_dir.resolve()),
    }
    state_contract_path = ensure_state_contract(historical_state, state_contract)
    selected_versions = historical_state / "selected_unisave_versions.jsonl"
    configure_version_provenance(module, selected_versions)
    module.TEXT_BUNDLE_METADATA = historical_state / "metadata.json"

    current_status = module.run_current_extraction(
        data_dir=data_dir,
        cafa_assessment_dir=args.cafa_assessment_dir,
        output_file=current,
        checkpoint_file=current_checkpoint,
    )
    historical_success, historical_failed, historical_status = module.extract_historical_text(
        data_dir=data_dir,
        cafa_assessment_dir=args.cafa_assessment_dir,
        output_file=historical,
        checkpoint_file=historical_checkpoint,
        raw_dir=historical_raw,
        splits=["test"],
        workers=args.workers,
    )
    historical_created_empty = False
    if not historical.exists():
        historical.touch()
        historical_created_empty = True
    version_provenance_audit = audit_version_provenance(
        historical, historical_raw, selected_versions
    )
    punct_metadata = module.build_historical_punct_v1_test_tsv(
        historical_tsv=historical,
        output_tsv=punct,
        data_dir=data_dir,
    )
    mixed_metadata = module.build_mixed_temporal_tsv(
        current_tsv=current,
        hist_test_tsv=punct,
        output_tsv=mixed,
        bundle_dir=historical_state,
        historical_tsv=historical,
        data_dir=data_dir,
    )

    if not mixed.is_file():
        raise SystemExit(f"PFP temporal recipe did not create: {mixed}")
    current_backup = historical_state / "protein_descriptions_current_before_mixed.tsv"
    shutil.copyfile(current, current_backup)
    shutil.copyfile(mixed, current)

    report = {
        "schema_version": 2,
        "pfp_text_script": str(script),
        "pfp_text_script_sha256": sha256_file(script),
        "historical_cutoff": cutoff_date,
        "requested_cutoff": cutoff_date,
        "effective_cutoff": cutoff_contract["effective_cutoff"],
        "cutoff_binding": cutoff_contract,
        "historical_state_dir": str(historical_state),
        "historical_state_contract": str(state_contract_path),
        "historical_state_contract_sha256": sha256_file(state_contract_path),
        "historical_raw_dir": str(historical_raw),
        "historical_checkpoint": str(historical_checkpoint),
        "selected_unisave_versions": str(selected_versions),
        "selected_unisave_versions_count": (
            len(selected_versions.read_text(encoding="utf-8").splitlines())
            if selected_versions.is_file()
            else 0
        ),
        "selected_unisave_versions_sha256": (
            sha256_file(selected_versions) if selected_versions.is_file() else None
        ),
        "selected_unisave_versions_audit": version_provenance_audit,
        "current_backup": str(current_backup),
        "current_status": current_status,
        "historical_success": historical_success,
        "historical_failed": historical_failed,
        "historical_status_counts": historical_status,
        "historical_created_empty": historical_created_empty,
        "punctuation_recipe": punct_metadata,
        "mixed_recipe": mixed_metadata,
        "embedding_input": str(current),
        "embedding_input_source": str(mixed),
    }
    report_path = temporal_dir / "framework_temporal_text_run.json"
    state_report_path = historical_state / "framework_temporal_text_run.json"
    atomic_write_json(state_report_path, report)
    atomic_write_json(report_path, report)
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
