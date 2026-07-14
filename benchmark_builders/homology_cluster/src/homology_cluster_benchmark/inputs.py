from __future__ import annotations

import gzip
import hashlib
import logging
import os
from pathlib import Path
import re
import shutil
import time
from typing import Iterator, TextIO
from urllib.error import HTTPError, URLError
from urllib.request import urlopen
import uuid

from .models import InputSpec, ResolvedInput


LOGGER = logging.getLogger(__name__)


SHA256_RE = re.compile(r"^[0-9a-fA-F]{64}$")


def open_text(path: str | Path) -> TextIO:
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="strict", newline="")
    return path.open("r", encoding="utf-8", errors="strict", newline="")


def sha256_file(path: str | Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def iter_lines(path: str | Path) -> Iterator[str]:
    with open_text(path) as handle:
        yield from handle


def _download(
    url: str, destination: Path, attempts: int = 3, timeout_seconds: int = 60
) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    for attempt in range(1, attempts + 1):
        temporary = destination.with_name(destination.name + f".part-{uuid.uuid4().hex}")
        try:
            LOGGER.info(
                "Downloading %s to %s (attempt %d/%d, timeout=%ds)",
                url, destination, attempt, attempts, timeout_seconds,
            )
            with urlopen(url, timeout=timeout_seconds) as response, temporary.open("wb") as out:
                expected = response.headers.get("Content-Length", "unknown")
                declared = int(expected) if str(expected).isdigit() else None
                reserve = max(1024 ** 3, int(declared * 0.05)) if declared else 1024 ** 3
                if declared and shutil.disk_usage(destination.parent).free < declared + reserve:
                    raise OSError(
                        f"Insufficient free space to download {url}: declared={declared}, "
                        f"reserve={reserve}, free={shutil.disk_usage(destination.parent).free}"
                    )
                copied = 0
                next_report = 1024 ** 3
                while True:
                    chunk = response.read(1024 * 1024)
                    if not chunk:
                        break
                    out.write(chunk)
                    copied += len(chunk)
                    if copied >= next_report:
                        if shutil.disk_usage(destination.parent).free < reserve:
                            raise OSError(
                                f"Download aborted before filling scratch; free space fell below "
                                f"the {reserve}-byte reserve"
                            )
                        LOGGER.info(
                            "Download progress for %s: %d bytes (declared total %s)",
                            url, copied, expected,
                        )
                        next_report += 1024 ** 3
                out.flush()
                os.fsync(out.fileno())
            os.replace(temporary, destination)
            LOGGER.info("Downloaded %s bytes from %s to %s", copied, url, destination)
            return
        except (HTTPError, URLError, TimeoutError, OSError) as exc:
            temporary.unlink(missing_ok=True)
            if attempt == attempts:
                raise OSError(
                    f"Failed to download {url} after {attempts} attempts: {exc}"
                ) from exc
            delay = attempt
            LOGGER.warning("Download attempt %d failed for %s: %s; retrying", attempt, url, exc)
            time.sleep(delay)


def resolve_input(spec: InputSpec, download_dir: Path, allow_downloads: bool = True) -> ResolvedInput:
    expected = spec.expected_sha256.lower() if spec.expected_sha256 else None
    if expected is not None and SHA256_RE.fullmatch(expected) is None:
        raise ValueError(f"{spec.name}: expected SHA-256 must contain exactly 64 hexadecimal characters")

    acquisition: str
    if spec.path is not None:
        path = spec.path.expanduser()
        if not path.is_file():
            raise FileNotFoundError(f"{spec.name}: configured local input does not exist: {path}")
        acquisition = "local"
    elif spec.url:
        if not allow_downloads:
            raise FileNotFoundError(
                f"{spec.name}: no local path was supplied and downloads are disabled; configured URL={spec.url}"
            )
        basename = Path(spec.url.split("?", 1)[0]).name or spec.name
        path = download_dir / f"{spec.name}-{basename}"
        if path.is_file():
            acquisition = "reused-download"
        else:
            _download(spec.url, path)
            acquisition = "downloaded"
    else:
        raise FileNotFoundError(f"{spec.name}: supply an explicit local path or source URL")

    size = path.stat().st_size
    if size <= 0:
        raise ValueError(f"{spec.name}: input is empty: {path}")
    actual = sha256_file(path)
    if expected is not None and actual != expected:
        raise ValueError(
            f"{spec.name}: SHA-256 mismatch for {path}; expected {expected}, observed {actual}"
        )
    return ResolvedInput(
        name=spec.name,
        resolved_path=path.resolve(),
        source_url=spec.url,
        release=spec.release,
        size_bytes=size,
        sha256=actual,
        expected_sha256=expected,
        acquisition=acquisition,
        source_population=spec.source_population,
    )
