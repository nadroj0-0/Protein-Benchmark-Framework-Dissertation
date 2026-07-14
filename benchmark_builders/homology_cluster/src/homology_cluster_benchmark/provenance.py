from __future__ import annotations

from contextlib import contextmanager
from datetime import datetime, timezone
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
from typing import Iterator
import uuid

from . import __version__
from .inputs import sha256_file
from .models import ResolvedInput


VARIABLE_FILES = frozenset({
    "disk_preflight.json",
    "input_manifest.json",
    "publication_metadata.json",
    "run_provenance.json",
    "output_manifest.json",
    "RUN_COMPLETE.json",
})


def _is_deterministic_payload(relative_path: Path) -> bool:
    return relative_path.name not in VARIABLE_FILES and "logs" not in relative_path.parts


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _package_version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "not-installed"


def git_state(repository: Path) -> dict[str, object]:
    commit = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=repository, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, check=False,
    )
    status = subprocess.run(
        ["git", "status", "--porcelain"], cwd=repository, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, check=False,
    )
    return {
        "commit": commit.stdout.strip() if commit.returncode == 0 else None,
        "dirty": bool(status.stdout.strip()) if status.returncode == 0 else None,
        "status_porcelain": status.stdout.splitlines() if status.returncode == 0 else [],
    }


def runtime_provenance(
    repository: Path,
    inputs: dict[str, ResolvedInput],
    parameters: dict[str, object],
    mmseqs_runtime: dict[str, object],
) -> dict[str, object]:
    return {
        "schema_version": 1,
        "builder_version": __version__,
        "created_at": utc_now(),
        "command": list(sys.argv),
        "repository": git_state(repository),
        "runtime": {
            "python": platform.python_version(),
            "python_executable": sys.executable,
            "conda_prefix": os.environ.get("CONDA_PREFIX"),
            "numpy": _package_version("numpy"),
            "pandas": _package_version("pandas"),
            "mmseqs2": mmseqs_runtime,
            "platform": platform.platform(),
            "operating_system": os.name,
        },
        "inputs": {key: value.as_dict() for key, value in sorted(inputs.items())},
        "parameters": parameters,
        "determinism": {
            "non_mmseqs_transformations": "deterministic for identical inputs, parameters, and package versions",
            "mmseqs_byte_identity": "not claimed across MMseqs2 versions or environments",
            "explicitly_variable_files": sorted(VARIABLE_FILES),
        },
    }


@contextmanager
def staging_output(final_dir: Path) -> Iterator[Path]:
    final_dir = final_dir.resolve()
    if final_dir.exists():
        raise FileExistsError(
            f"Output directory already exists; resume/overwrite is intentionally unsupported: {final_dir}"
        )
    final_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = final_dir.parent / f".{final_dir.name}.staging-{uuid.uuid4().hex}"
    staging.mkdir()
    try:
        yield staging
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise


def publish(staging: Path, final_dir: Path) -> None:
    if final_dir.exists():
        raise FileExistsError(final_dir)
    os.replace(staging, final_dir)


def output_manifest(directory: Path) -> dict[str, object]:
    files = []
    for path in sorted(directory.rglob("*"), key=lambda item: item.relative_to(directory).as_posix()):
        if not path.is_file() or path.name in {"output_manifest.json", "RUN_COMPLETE.json"}:
            continue
        relative = path.relative_to(directory)
        files.append({
            "path": relative.as_posix(),
            "size_bytes": path.stat().st_size,
            "sha256": sha256_file(path),
            "deterministic_payload": _is_deterministic_payload(relative),
        })
    return {
        "schema_version": 1,
        "created_at": utc_now(),
        "payload_file_count": len(files),
        "files": files,
    }


def write_output_manifest(directory: Path) -> Path:
    path = directory / "output_manifest.json"
    path.write_text(json.dumps(output_manifest(directory), indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def verify_output_manifest(directory: Path) -> None:
    manifest_path = directory / "output_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("schema_version") != 1 or not isinstance(manifest.get("files"), list):
        raise ValueError("Published output manifest has an unsupported schema")
    entries = manifest["files"]
    if manifest.get("payload_file_count") != len(entries):
        raise ValueError("Published output manifest payload_file_count is inconsistent")
    listed: set[str] = set()
    for entry in entries:
        relative = Path(str(entry.get("path", "")))
        if relative.is_absolute() or ".." in relative.parts or not relative.parts:
            raise ValueError(f"Unsafe output manifest path: {relative}")
        relative_text = relative.as_posix()
        if relative_text in listed:
            raise ValueError(f"Duplicate output manifest path: {relative_text}")
        listed.add(relative_text)
        path = directory / entry["path"]
        if not path.is_file():
            raise ValueError(f"Published manifest file is missing: {entry['path']}")
        if path.stat().st_size != entry["size_bytes"] or sha256_file(path) != entry["sha256"]:
            raise ValueError(f"Published manifest hash/size mismatch: {entry['path']}")
    actual = {
        path.relative_to(directory).as_posix()
        for path in directory.rglob("*")
        if path.is_file() and path.name not in {"output_manifest.json", "RUN_COMPLETE.json"}
    }
    if actual != listed:
        raise ValueError(
            "Published manifest does not reconcile with the actual payload file set: "
            f"unlisted={sorted(actual - listed)[:20]}, missing={sorted(listed - actual)[:20]}"
        )


PUBLICATION_MARKER_KEYS = (
    "fixture_mode", "production_eligible", "benchmark_scope", "identity_percent",
    "identities", "split_policy", "training_population", "seed", "min_count",
    "run_input_manifest_sha256", "frozen_input_manifest_sha256",
    "expected_mmseqs_version", "observed_mmseqs_version", "repository_commit",
    "mmseqs_resolved_executable", "mmseqs_executable_sha256", "scientific_fingerprint",
)


def write_publication_metadata(directory: Path, payload: dict[str, object]) -> Path:
    path = directory / "publication_metadata.json"
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def write_completion_marker(directory: Path) -> Path:
    manifest_path = directory / "output_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    publication_path = directory / "publication_metadata.json"
    publication = json.loads(publication_path.read_text(encoding="utf-8"))
    marker = {
        "complete": True,
        "completed_at": utc_now(),
        "manifest": "output_manifest.json",
        "manifest_sha256": sha256_file(manifest_path),
        "publication_metadata": "publication_metadata.json",
        "publication_metadata_sha256": sha256_file(publication_path),
        "payload_file_count": manifest["payload_file_count"],
        "post_publication_hash_verification": True,
        **{key: publication.get(key) for key in PUBLICATION_MARKER_KEYS},
    }
    temporary = directory / f".RUN_COMPLETE.json.part-{uuid.uuid4().hex}"
    temporary.write_text(json.dumps(marker, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    final = directory / "RUN_COMPLETE.json"
    os.replace(temporary, final)
    return final
