#!/usr/bin/env python3
"""Create and safely extract consolidated PFP embedding-cache archives."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import re
import shutil
import tarfile
import tempfile
from collections import Counter
from pathlib import Path, PurePosixPath
from typing import BinaryIO, Iterable


ARCHIVE_PREFIX = PurePosixPath("data/embedding_cache")
SAFE_ID = re.compile(r"^[^\s/\\]+$")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


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


def load_directories(config_path: Path) -> tuple[str, ...]:
    config = json.loads(config_path.read_text(encoding="utf-8"))
    modalities = config.get("modalities")
    if not isinstance(modalities, dict) or set(modalities) != {
        "sequence",
        "text",
        "structure",
        "ppi",
    }:
        raise ValueError("Run config must define exactly the four PFP modalities")
    directories = tuple(str(modalities[name]["directory"]) for name in sorted(modalities))
    if len(set(directories)) != len(directories):
        raise ValueError("Run config repeats an embedding cache directory")
    for directory in directories:
        if Path(directory).name != directory or not SAFE_ID.fullmatch(directory):
            raise ValueError(f"Unsafe embedding cache directory: {directory!r}")
    return directories


def iter_cache_files(cache_root: Path, directories: Iterable[str]) -> list[tuple[str, Path]]:
    expected = set(directories)
    if not cache_root.is_dir() or cache_root.is_symlink():
        raise ValueError(f"Embedding cache root must be a real directory: {cache_root}")
    unexpected = sorted(path.name for path in cache_root.iterdir() if path.name not in expected)
    if unexpected:
        raise ValueError(f"Embedding cache has unexpected top-level entries: {unexpected[:5]}")
    result: list[tuple[str, Path]] = []
    for directory in sorted(expected):
        source_dir = cache_root / directory
        if not source_dir.is_dir() or source_dir.is_symlink():
            raise ValueError(f"Embedding modality directory is missing or unsafe: {source_dir}")
        for path in sorted(source_dir.iterdir(), key=lambda item: item.name):
            if (
                not path.is_file()
                or path.is_symlink()
                or path.suffix != ".npy"
                or not SAFE_ID.fullmatch(path.stem)
            ):
                raise ValueError(f"Unexpected embedding cache entry: {path}")
            result.append((directory, path))
    if not result:
        raise ValueError("Embedding cache contains no arrays")
    return result


def copy_and_hash(source: BinaryIO, destination: BinaryIO) -> tuple[str, int]:
    digest = hashlib.sha256()
    size = 0
    for chunk in iter(lambda: source.read(1024 * 1024), b""):
        digest.update(chunk)
        destination.write(chunk)
        size += len(chunk)
    return digest.hexdigest(), size


def create_archive(
    cache_root: Path,
    archive_path: Path,
    config_path: Path,
) -> dict:
    if archive_path.exists():
        raise ValueError(f"Archive output already exists: {archive_path}")
    directories = load_directories(config_path)
    files = iter_cache_files(cache_root, directories)
    archive_path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{archive_path.name}.", suffix=".tmp", dir=str(archive_path.parent)
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    counts: Counter[str] = Counter()
    content_digest = hashlib.sha256()
    try:
        with temporary.open("wb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
                with tarfile.open(
                    mode="w", fileobj=compressed, format=tarfile.PAX_FORMAT
                ) as archive:
                    for directory, path in files:
                        member_name = str(ARCHIVE_PREFIX / directory / path.name)
                        info = tarfile.TarInfo(member_name)
                        info.size = path.stat().st_size
                        info.mtime = 0
                        info.mode = 0o644
                        info.uid = info.gid = 0
                        info.uname = info.gname = ""
                        with path.open("rb") as source:
                            archive.addfile(info, source)
                        file_sha = sha256_file(path)
                        content_digest.update(member_name.encode("utf-8"))
                        content_digest.update(b"\t")
                        content_digest.update(file_sha.encode("ascii"))
                        content_digest.update(b"\n")
                        counts[directory] += 1
            raw.flush()
            os.fsync(raw.fileno())
        os.replace(temporary, archive_path)
    finally:
        if temporary.exists():
            temporary.unlink()
    return {
        "schema_version": 1,
        "operation": "create",
        "archive": str(archive_path.resolve()),
        "archive_size_bytes": archive_path.stat().st_size,
        "archive_sha256": sha256_file(archive_path),
        "member_count": len(files),
        "members_by_directory": dict(sorted(counts.items())),
        "member_content_sha256": content_digest.hexdigest(),
    }


def normalized_member(
    member: tarfile.TarInfo, directories: set[str]
) -> tuple[str, str] | None:
    name = member.name.removeprefix("./")
    value = PurePosixPath(name)
    allowed_directories = {
        ARCHIVE_PREFIX,
        PurePosixPath("data"),
        *(ARCHIVE_PREFIX / directory for directory in directories),
    }
    if member.isdir() and value in allowed_directories:
        return None
    if not member.isfile():
        raise ValueError(f"Archive member is not a regular file: {member.name}")
    if len(value.parts) != 4 or PurePosixPath(*value.parts[:2]) != ARCHIVE_PREFIX:
        raise ValueError(f"Archive member is outside data/embedding_cache: {member.name}")
    directory, filename = value.parts[2], value.parts[3]
    if directory not in directories:
        raise ValueError(f"Archive member uses an unknown modality directory: {member.name}")
    file_path = PurePosixPath(filename)
    if (
        file_path.name != filename
        or not filename.endswith(".npy")
        or not SAFE_ID.fullmatch(filename[:-4])
    ):
        raise ValueError(f"Archive member has an unsafe embedding filename: {member.name}")
    return directory, filename


def extract_archive(
    archive_path: Path,
    output_cache_root: Path,
    config_path: Path,
) -> dict:
    if not archive_path.is_file() or archive_path.is_symlink():
        raise ValueError(f"Embedding archive is missing or unsafe: {archive_path}")
    if output_cache_root.exists() and any(output_cache_root.iterdir()):
        raise ValueError(f"Extraction directory is not empty: {output_cache_root}")
    directories = set(load_directories(config_path))
    output_cache_root.mkdir(parents=True, exist_ok=True)
    for directory in directories:
        (output_cache_root / directory).mkdir()
    seen: set[tuple[str, str]] = set()
    counts: Counter[str] = Counter()
    content_digest = hashlib.sha256()
    try:
        with tarfile.open(archive_path, mode="r:gz") as archive:
            for member in archive:
                normalized = normalized_member(member, directories)
                if normalized is None:
                    continue
                if normalized in seen:
                    raise ValueError(f"Archive repeats embedding member: {member.name}")
                seen.add(normalized)
                directory, filename = normalized
                destination = output_cache_root / directory / filename
                source = archive.extractfile(member)
                if source is None:
                    raise ValueError(f"Cannot read archive member: {member.name}")
                descriptor, temporary_name = tempfile.mkstemp(
                    prefix=f".{filename}.", suffix=".tmp", dir=str(destination.parent)
                )
                try:
                    with source, os.fdopen(descriptor, "wb") as output:
                        file_sha, size = copy_and_hash(source, output)
                        output.flush()
                        os.fsync(output.fileno())
                    if size != member.size:
                        raise ValueError(f"Archive member size changed while reading: {member.name}")
                    os.replace(temporary_name, destination)
                finally:
                    try:
                        os.unlink(temporary_name)
                    except FileNotFoundError:
                        pass
                canonical_name = str(ARCHIVE_PREFIX / directory / filename)
                content_digest.update(canonical_name.encode("utf-8"))
                content_digest.update(b"\t")
                content_digest.update(file_sha.encode("ascii"))
                content_digest.update(b"\n")
                counts[directory] += 1
    except BaseException:
        shutil.rmtree(output_cache_root, ignore_errors=True)
        raise
    if not seen:
        shutil.rmtree(output_cache_root, ignore_errors=True)
        raise ValueError("Embedding archive contains no arrays")
    return {
        "schema_version": 1,
        "operation": "extract",
        "archive": str(archive_path.resolve()),
        "archive_size_bytes": archive_path.stat().st_size,
        "archive_sha256": sha256_file(archive_path),
        "member_count": len(seen),
        "members_by_directory": dict(sorted(counts.items())),
        "member_content_sha256": content_digest.hexdigest(),
        "output_cache_root": str(output_cache_root.resolve()),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    create = subparsers.add_parser("create")
    create.add_argument("--cache-root", type=Path, required=True)
    create.add_argument("--archive", type=Path, required=True)
    create.add_argument("--config", type=Path, required=True)
    create.add_argument("--report", type=Path)
    extract = subparsers.add_parser("extract")
    extract.add_argument("--archive", type=Path, required=True)
    extract.add_argument("--output-cache-root", type=Path, required=True)
    extract.add_argument("--config", type=Path, required=True)
    extract.add_argument("--report", type=Path)
    args = parser.parse_args()
    try:
        if args.command == "create":
            result = create_archive(args.cache_root, args.archive, args.config)
        else:
            result = extract_archive(args.archive, args.output_cache_root, args.config)
        if args.report:
            atomic_write_json(args.report, result)
    except (OSError, ValueError, json.JSONDecodeError, tarfile.TarError) as error:
        raise SystemExit(f"ERROR: {error}") from error
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
