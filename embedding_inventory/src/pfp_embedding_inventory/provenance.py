"""Deterministic identities and complete run provenance for inventory plans."""

from __future__ import annotations

import hashlib
import platform
import shlex
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping, Optional, Sequence, Tuple

import numpy as np

from . import __version__
from .benchmark import required_csv_names
from .models import ArtifactVerification, BenchmarkData, CacheCatalog, MODALITIES, PlannerConfig
from .paths import PathSafetyError, resolve_within


CACHE_CATALOG_SCHEMA = "relative-path-tab-size-tab-file-sha256-lf-v1"


class ProvenanceError(ValueError):
    pass


class HashCache:
    """Hash immutable run inputs once, even when paths are used in two roles."""

    def __init__(self) -> None:
        self._values: Dict[Path, str] = {}
        self.compute_count: Dict[Path, int] = {}

    def sha256(self, path: Path, *, remember: bool = True) -> str:
        resolved = path.resolve(strict=True)
        if remember and resolved in self._values:
            return self._values[resolved]
        digest = hashlib.sha256()
        with resolved.open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
        value = digest.hexdigest()
        if remember:
            self.compute_count[resolved] = self.compute_count.get(resolved, 0) + 1
            self._values[resolved] = value
        return value


def compute_cache_catalog(
    embedding_cache: Path,
    config: PlannerConfig,
    hash_cache: Optional[HashCache] = None,
) -> CacheCatalog:
    """Hash every configured ``.npy`` byte stream into a stable cache catalog.

    Catalog lines are globally sorted and encoded as
    ``relative/path<TAB>size<TAB>file_sha256<LF>``. Directory names and mtimes
    cannot establish identity; file bytes and paths can.
    """
    root = embedding_cache.resolve(strict=True)
    hashes = hash_cache or HashCache()
    modality_counts: Dict[str, int] = {modality: 0 for modality in MODALITIES}
    modality_bytes: Dict[str, int] = {modality: 0 for modality in MODALITIES}
    modality_digests = {modality: hashlib.sha256() for modality in MODALITIES}
    global_digest = hashlib.sha256()
    # Sorting directory prefixes, then basenames, is identical to sorting all
    # cache-relative paths and avoids retaining hundreds of thousands of lines.
    for modality in sorted(MODALITIES, key=lambda name: config.modalities[name].directory):
        try:
            source_dir = resolve_within(
                root, Path(config.modalities[modality].directory),
                "%s modality directory" % modality,
            )
        except PathSafetyError as exc:
            raise ProvenanceError(str(exc)) from exc
        paths: Iterable[Path] = (
            sorted(source_dir.glob("*.npy"), key=lambda item: item.name)
            if source_dir.is_dir()
            else ()
        )
        for path in paths:
            if not path.is_file():
                continue
            resolved = path.resolve(strict=True)
            try:
                resolved.relative_to(source_dir)
            except ValueError as exc:
                raise ProvenanceError("cache catalog refuses symlink outside modality directory") from exc
            size = resolved.stat().st_size
            relative = path.relative_to(root).as_posix()
            # Do not retain 265k embedding digests: each one is used once.
            line = "%s\t%d\t%s\n" % (
                relative,
                size,
                hashes.sha256(resolved, remember=False),
            )
            encoded = line.encode("utf-8")
            global_digest.update(encoded)
            modality_digests[modality].update(encoded)
            modality_counts[modality] += 1
            modality_bytes[modality] += size
    modality_fingerprints = {
        modality: modality_digests[modality].hexdigest() for modality in MODALITIES
    }
    return CacheCatalog(
        schema=CACHE_CATALOG_SCHEMA,
        fingerprint=global_digest.hexdigest(),
        modality_fingerprints=modality_fingerprints,
        modality_counts=modality_counts,
        modality_bytes=modality_bytes,
        total_files=sum(modality_counts.values()),
        total_bytes=sum(modality_bytes.values()),
    )


def assert_cache_catalog_unchanged(before: CacheCatalog, after: CacheCatalog) -> None:
    """Abort if cache bytes or membership changed during array validation."""
    if before.as_dict() != after.as_dict():
        raise ProvenanceError(
            "embedding cache changed after its catalog was verified; refusing to write manifests"
        )


def verify_artifact_scope(
    config: PlannerConfig,
    target: BenchmarkData,
    source: BenchmarkData,
    catalog: CacheCatalog,
    embedding_cache: Path,
    artifact_root: Path,
    hash_cache: Optional[HashCache] = None,
) -> ArtifactVerification:
    spec = config.artifact_scope
    if spec.mode == "none":
        return ArtifactVerification(
            configured=False,
            verified=False,
            artifact_id=spec.artifact_id,
            checks={},
            reasons=["no verified artifact scope is configured"],
            expected={},
            observed={
                "target_benchmark_fingerprint": target.fingerprint,
                "source_benchmark_fingerprint": source.fingerprint,
                "cache_catalog_fingerprint": catalog.fingerprint,
                "embedding_cache_root": str(embedding_cache.resolve()),
            },
        )

    hashes = hash_cache or HashCache()
    root = artifact_root.resolve()
    checks: Dict[str, bool] = {
        "source_benchmark_fingerprint": source.fingerprint == spec.expected_benchmark_fingerprint,
        "cache_catalog_fingerprint": catalog.fingerprint == spec.expected_cache_catalog_fingerprint,
        "cache_modality_counts": catalog.modality_counts == spec.expected_modality_counts,
        "cache_total_files": catalog.total_files == spec.expected_total_files,
        "cache_total_bytes": catalog.total_bytes == spec.expected_total_bytes,
    }
    archive_observed: Dict[str, Any] = {}
    for archive in spec.archives:
        try:
            path = resolve_within(root, Path(archive.path), "artifact archive")
        except PathSafetyError as exc:
            raise ProvenanceError(str(exc)) from exc
        exists = path.is_file()
        observed_hash = hashes.sha256(path) if exists else ""
        checks["archive:%s" % archive.path] = exists and observed_hash == archive.sha256
        archive_observed[archive.path] = {
            "path": str(path.resolve()) if exists else str(path),
            "exists": exists,
            "sha256": observed_hash,
            "expected_sha256": archive.sha256,
        }

    commit = _git_value(root, ("rev-parse", "HEAD"))
    checks["reference_commit"] = commit == spec.expected_reference_commit
    reference_observed: Dict[str, Any] = {}
    for reference in spec.reference_files:
        try:
            path = resolve_within(root, Path(reference.path), "artifact reference file")
        except PathSafetyError as exc:
            raise ProvenanceError(str(exc)) from exc
        exists = path.is_file()
        observed_hash = hashes.sha256(path) if exists else ""
        checks["reference_file:%s" % reference.path] = exists and observed_hash == reference.sha256
        reference_observed[reference.path] = {
            "path": str(path.resolve()) if exists else str(path),
            "exists": exists,
            "sha256": observed_hash,
            "expected_sha256": reference.sha256,
        }

    failures = [name for name, passed in checks.items() if not passed]
    return ArtifactVerification(
        configured=True,
        verified=all(checks.values()),
        artifact_id=spec.artifact_id,
        checks=checks,
        reasons=["failed artifact proof: %s" % name for name in failures],
        expected={
            "source_benchmark_fingerprint": spec.expected_benchmark_fingerprint,
            "cache_catalog_fingerprint": spec.expected_cache_catalog_fingerprint,
            "modality_counts": spec.expected_modality_counts,
            "total_files": spec.expected_total_files,
            "total_bytes": spec.expected_total_bytes,
            "reference_commit": spec.expected_reference_commit,
            "metadata_url": spec.metadata_url,
        },
        observed={
            "target_benchmark_fingerprint": target.fingerprint,
            "source_benchmark_fingerprint": source.fingerprint,
            "embedding_cache_root": str(embedding_cache.resolve()),
            "cache_catalog": catalog.as_dict(),
            "archives": archive_observed,
            "reference_commit": commit,
            "reference_files": reference_observed,
        },
    )


def build_run_provenance(
    *,
    command: Sequence[str],
    repository: Path,
    config_path: Path,
    alias_path: Optional[Path],
    target: BenchmarkData,
    source: BenchmarkData,
    embedding_cache: Path,
    artifact_root: Path,
    catalog: CacheCatalog,
    verification: ArtifactVerification,
    policy: str,
    report_level: str,
    runtime_options: Mapping[str, Any],
    hash_cache: Optional[HashCache] = None,
) -> Dict[str, Any]:
    hashes = hash_cache or HashCache()
    repo = repository.resolve()
    git_commit, git_commit_error = _git_probe(repo, ("rev-parse", "HEAD"))
    dirty_worktree, git_status_error = _git_dirty_probe(repo)
    return {
        "schema_version": 1,
        "timestamp_utc": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        "command": shlex.join([str(item) for item in command]),
        "software": {
            "repository": str(repo),
            "git_available": not bool(git_commit_error),
            "git_commit": git_commit,
            "dirty_worktree": None if git_status_error else dirty_worktree,
            "git_commit_error": git_commit_error,
            "git_status_error": git_status_error,
            "python_version": platform.python_version(),
            "python_implementation": platform.python_implementation(),
            "numpy_version": np.__version__,
            "package_version": __version__,
            "executable": sys.executable,
        },
        "inputs": {
            "config": _file_identity(config_path, hashes),
            "alias_map": _file_identity(alias_path, hashes) if alias_path else None,
            "target_csvs": _csv_manifest(target.directory, hashes),
            "source_csvs": _csv_manifest(source.directory, hashes),
            "target_benchmark_fingerprint": target.fingerprint,
            "source_benchmark_fingerprint": source.fingerprint,
            "embedding_cache": str(embedding_cache.resolve()),
            "artifact_root": str(artifact_root.resolve()),
            "cache_catalog": catalog.as_dict(),
            "available_archives": verification.observed.get("archives", {}),
        },
        "artifact_verification": verification.as_dict(),
        "run": {
            "compatibility_policy": policy,
            "report_level": report_level,
            "runtime_options": dict(runtime_options),
        },
    }


def provenance_markdown(provenance: Mapping[str, Any]) -> str:
    software = provenance["software"]
    inputs = provenance["inputs"]
    artifact = provenance["artifact_verification"]
    lines = [
        "# Embedding inventory run provenance",
        "",
        "- UTC timestamp: `%s`" % provenance["timestamp_utc"],
        "- Command: `%s`" % provenance["command"].replace("`", "'"),
        "- Repository commit: `%s`" % (software["git_commit"] or "not available"),
        "- Dirty worktree: `%s`" % str(software["dirty_worktree"]).lower(),
        "- Python / NumPy / package: `%s` / `%s` / `%s`"
        % (software["python_version"], software["numpy_version"], software["package_version"]),
        "- Target benchmark fingerprint: `%s`" % inputs["target_benchmark_fingerprint"],
        "- Source benchmark fingerprint: `%s`" % inputs["source_benchmark_fingerprint"],
        "- Cache catalog fingerprint: `%s`" % inputs["cache_catalog"]["fingerprint"],
        "- Published cache authenticated: `%s`" % str(artifact["verified"]).lower(),
        "- Policy / report level: `%s` / `%s`"
        % (provenance["run"]["compatibility_policy"], provenance["run"]["report_level"]),
        "",
    ]
    if artifact["reasons"]:
        lines.extend(["## Artifact verification exceptions", ""])
        lines.extend("- %s" % reason for reason in artifact["reasons"])
        lines.append("")
    return "\n".join(lines)


def _csv_manifest(directory: Path, hashes: HashCache) -> Dict[str, Dict[str, Any]]:
    return {
        filename: _file_identity(directory / filename, hashes)
        for filename in required_csv_names()
    }


def _file_identity(path: Path, hashes: HashCache) -> Dict[str, Any]:
    resolved = path.resolve(strict=True)
    return {
        "path": str(resolved),
        "size_bytes": resolved.stat().st_size,
        "sha256": hashes.sha256(resolved),
    }


def _git_value(root: Path, args: Tuple[str, ...], strip: bool = True) -> str:
    return _git_probe(root, args, strip=strip)[0]


def _git_probe(
    root: Path, args: Tuple[str, ...], strip: bool = True
) -> Tuple[str, str]:
    try:
        completed = subprocess.run(
            ("git", "-C", str(root), *args),
            check=False,
            capture_output=True,
            text=True,
            timeout=15,
        )
    except subprocess.TimeoutExpired:
        return "", "git command timed out after 15 seconds"
    except OSError as exc:
        return "", "%s: %s" % (type(exc).__name__, exc)
    if completed.returncode != 0:
        message = completed.stderr.strip() or "exit code %d" % completed.returncode
        return "", message
    return (completed.stdout.strip() if strip else completed.stdout), ""


def _git_dirty_probe(root: Path) -> Tuple[bool, str]:
    """Check tracked changes first, then untracked files without full status."""
    try:
        tracked = subprocess.run(
            ("git", "-C", str(root), "diff-index", "--quiet", "HEAD", "--"),
            check=False,
            capture_output=True,
            text=True,
            timeout=15,
        )
        if tracked.returncode == 1:
            return True, ""
        if tracked.returncode != 0:
            return False, tracked.stderr.strip() or "git diff-index failed"
        untracked = subprocess.run(
            (
                "git", "-C", str(root), "ls-files", "--others",
                "--exclude-standard", "--directory",
            ),
            check=False,
            capture_output=True,
            text=True,
            timeout=15,
        )
        if untracked.returncode != 0:
            return False, untracked.stderr.strip() or "git ls-files failed"
        return bool(untracked.stdout), ""
    except subprocess.TimeoutExpired:
        return False, "git dirty-worktree check timed out after 15 seconds"
    except OSError as exc:
        return False, "%s: %s" % (type(exc).__name__, exc)
