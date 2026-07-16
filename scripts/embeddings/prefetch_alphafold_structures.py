#!/usr/bin/env python3
"""Acquire AlphaFold PDBs with bounded concurrency without editing PFP.

PFP's mapping and API interpretation functions remain the source of truth. This
framework wrapper calls those functions for only PDBs absent from the persistent
cache, downloads files atomically, and copies the requested view into disposable
scratch for ESM-IF1 inference.
"""

from __future__ import annotations

import argparse
import csv
import fcntl
import hashlib
import importlib.util
import json
import os
import shutil
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterator, Optional, Tuple

import requests


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def valid_pdb(path: Path) -> bool:
    if not path.is_file() or path.stat().st_size < 100:
        return False
    with path.open("rb") as handle:
        sample = handle.read(1024 * 1024)
    return b"ATOM  " in sample or b"HETATM" in sample


def atomic_copy(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    descriptor, name = tempfile.mkstemp(
        prefix=f".{destination.name}.", suffix=".tmp", dir=str(destination.parent)
    )
    os.close(descriptor)
    temporary = Path(name)
    try:
        shutil.copyfile(source, temporary)
        if not valid_pdb(temporary):
            raise ValueError(f"Invalid PDB content copied from {source}")
        os.replace(temporary, destination)
    finally:
        if temporary.exists():
            temporary.unlink()


@contextmanager
def cache_lock(cache_dir: Path) -> Iterator[None]:
    cache_dir.mkdir(parents=True, exist_ok=True)
    with (cache_dir / ".alphafold.lock").open("a+", encoding="ascii") as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def load_pfp_module(path: Path):
    specification = importlib.util.spec_from_file_location("pfp_alphafold_coverage", path)
    if specification is None or specification.loader is None:
        raise ValueError(f"Cannot load PFP AlphaFold helper: {path}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


def download_atomic(url: str, destination: Path, attempts: int) -> Tuple[bool, str]:
    destination.parent.mkdir(parents=True, exist_ok=True)
    for attempt in range(1, attempts + 1):
        descriptor, name = tempfile.mkstemp(
            prefix=f".{destination.name}.", suffix=".partial", dir=str(destination.parent)
        )
        os.close(descriptor)
        temporary = Path(name)
        try:
            response = requests.get(url, timeout=(15, 120))
            if response.status_code == 200:
                temporary.write_bytes(response.content)
                if not valid_pdb(temporary):
                    return False, "invalid_pdb_content"
                os.replace(temporary, destination)
                return True, ""
            if response.status_code == 404:
                return False, "pdb_http_404"
            detail = f"pdb_http_{response.status_code}"
        except requests.RequestException as error:
            detail = f"pdb_request_error:{type(error).__name__}:{str(error)[:160]}"
        finally:
            if temporary.exists():
                temporary.unlink()
        if attempt < attempts:
            time.sleep(min(2 ** attempt, 8))
    return False, detail


def load_manifest(path: Path) -> Dict[str, dict]:
    if not path.is_file():
        return {}
    result: Dict[str, dict] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            result[row["protein_id"]] = row
    return result


def write_manifest(path: Path, rows: Dict[str, dict]) -> None:
    columns = [
        "protein_id",
        "sha256",
        "size_bytes",
        "pdb_url",
        "alphafold_version",
        "resolved_accession",
        "acquired_at",
    ]
    lines = ["\t".join(columns)]
    for protein_id in sorted(rows):
        row = rows[protein_id]
        lines.append("\t".join(str(row.get(column, "")) for column in columns))
    content = "\n".join(lines) + "\n"
    descriptor, name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent)
    )
    temporary = Path(name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pfp-root", type=Path, required=True)
    parser.add_argument("--cafa-assessment-dir", type=Path, required=True)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--persistent-cache-dir", type=Path, required=True)
    parser.add_argument("--workspace-pdb-dir", type=Path, required=True)
    parser.add_argument("--coverage-report", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--api-workers", type=int, default=8)
    parser.add_argument("--download-workers", type=int, default=8)
    parser.add_argument("--attempts", type=int, default=3)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.api_workers <= 0 or args.download_workers <= 0 or args.attempts <= 0:
        raise SystemExit("Worker and attempt counts must be positive")
    script = args.pfp_root / "scripts/check_alphafold_coverage.py"
    if not script.is_file():
        raise SystemExit(f"Missing PFP AlphaFold helper: {script}")
    if not args.cafa_assessment_dir.is_dir():
        raise SystemExit(f"Missing CAFA assessment directory: {args.cafa_assessment_dir}")

    module = load_pfp_module(script)
    requested = {str(value) for value in module.get_all_cafa_proteins(args.data_dir)}
    if not requested:
        raise SystemExit("Prepared workspace contains no proteins")
    for path in args.workspace_pdb_dir.glob("*.pdb") if args.workspace_pdb_dir.is_dir() else ():
        if path.stem not in requested:
            path.unlink()

    args.persistent_cache_dir.mkdir(parents=True, exist_ok=True)
    args.workspace_pdb_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = args.persistent_cache_dir / "alphafold_source_manifest.tsv"
    with cache_lock(args.persistent_cache_dir):
        manifest = load_manifest(manifest_path)
        cached = set()
        invalid_cached = []
        for protein_id in sorted(requested):
            path = args.persistent_cache_dir / f"{protein_id}.pdb"
            metadata = manifest.get(protein_id)
            authenticated = bool(
                metadata
                and metadata.get("sha256")
                and valid_pdb(path)
                and sha256_file(path) == metadata["sha256"]
            )
            if authenticated:
                cached.add(protein_id)
            elif path.exists():
                path.unlink()
                manifest.pop(protein_id, None)
                invalid_cached.append(protein_id)

        missing = requested - cached
        results = {"found": [], "not_found": [], "no_mapping": [], "errors": []}
        if missing:
            mapping = module.build_cafa_to_accession_mapping(
                args.cafa_assessment_dir, args.data_dir
            )
            results = module.check_alphafold_coverage(
                missing,
                mapping,
                output_file=args.coverage_report,
                num_workers=args.api_workers,
            )
        else:
            args.coverage_report.parent.mkdir(parents=True, exist_ok=True)
            args.coverage_report.write_text(
                "All requested structures were already present in the persistent cache.\n",
                encoding="utf-8",
            )

        found_metadata = {
            row[0]: {
                "protein_id": row[0],
                "resolved_accession": row[2] or row[3] or "",
                "alphafold_version": row[4],
                "pdb_url": row[5] or "",
            }
            for row in results["found"]
        }
        downloaded = []
        failures: Dict[str, str] = {}

        def fetch(protein_id: str) -> Tuple[str, bool, str]:
            metadata = found_metadata[protein_id]
            url = metadata["pdb_url"]
            if not url:
                return protein_id, False, "missing_pdb_url"
            success, detail = download_atomic(
                url,
                args.persistent_cache_dir / f"{protein_id}.pdb",
                args.attempts,
            )
            return protein_id, success, detail

        with ThreadPoolExecutor(max_workers=args.download_workers) as executor:
            futures = {
                executor.submit(fetch, protein_id): protein_id
                for protein_id in sorted(found_metadata)
            }
            for future in as_completed(futures):
                protein_id, success, detail = future.result()
                if success:
                    downloaded.append(protein_id)
                    path = args.persistent_cache_dir / f"{protein_id}.pdb"
                    metadata = found_metadata[protein_id]
                    manifest[protein_id] = {
                        **metadata,
                        "sha256": sha256_file(path),
                        "size_bytes": path.stat().st_size,
                        "acquired_at": utc_now(),
                    }
                else:
                    failures[protein_id] = detail

        write_manifest(manifest_path, manifest)
        available = sorted(
            protein_id
            for protein_id in requested
            if valid_pdb(args.persistent_cache_dir / f"{protein_id}.pdb")
        )
        for protein_id in available:
            atomic_copy(
                args.persistent_cache_dir / f"{protein_id}.pdb",
                args.workspace_pdb_dir / f"{protein_id}.pdb",
            )

    report = {
        "schema_version": 1,
        "completed_at": utc_now(),
        "pfp_helper": str(script.resolve()),
        "requested": len(requested),
        "cached_before": len(cached),
        "missing_before": len(missing),
        "invalid_cached_removed": invalid_cached,
        "api_workers": args.api_workers,
        "download_workers": args.download_workers,
        "downloaded": len(downloaded),
        "available_for_if1": len(available),
        "not_found": len(results["not_found"]),
        "no_mapping": len(results["no_mapping"]),
        "api_errors": len(results["errors"]),
        "download_failures": failures,
        "persistent_cache_dir": str(args.persistent_cache_dir.resolve()),
        "workspace_pdb_dir": str(args.workspace_pdb_dir.resolve()),
        "source_manifest": str(manifest_path.resolve()),
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
