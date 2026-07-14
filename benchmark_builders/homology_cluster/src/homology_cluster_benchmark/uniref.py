from __future__ import annotations

import hashlib
import logging
from pathlib import Path
import re
import sqlite3
import time
from typing import Iterator

from .inputs import open_text


LOGGER = logging.getLogger(__name__)


SEQUENCE_RE = re.compile(r"^[ACDEFGHIKLMNPQRSTVWYBXZJUO]+$", re.IGNORECASE)


def iter_fasta(path: str | Path) -> Iterator[tuple[str, str]]:
    header: str | None = None
    parts: list[str] = []
    with open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(parts)
                header = line[1:]
                if not header:
                    raise ValueError(f"Empty FASTA header at {path}:{line_number}")
                parts = []
            else:
                if header is None:
                    raise ValueError(f"Sequence data appears before a FASTA header at {path}:{line_number}")
                parts.append(line)
        if header is not None:
            yield header, "".join(parts)


def fasta_identifier(header: str) -> str:
    identifier = header.split(None, 1)[0]
    if not identifier:
        raise ValueError("FASTA identifier is empty")
    return identifier


class UniRefIndex:
    """Disk-backed UniRef90 identifier and sequence-hash index."""

    def __init__(self, database: Path):
        self.database = database

    @classmethod
    def build(cls, fasta: Path, database: Path) -> "UniRefIndex":
        started = time.monotonic()
        LOGGER.info("UniRef90 index build started: %s -> %s", fasta, database)
        database.parent.mkdir(parents=True, exist_ok=True)
        database.unlink(missing_ok=True)
        connection = sqlite3.connect(database)
        try:
            # Disposable scratch index: avoid a database-sized WAL; any interrupted load is rebuilt.
            connection.execute("PRAGMA journal_mode=OFF")
            connection.execute("PRAGMA synchronous=OFF")
            connection.execute(
                "CREATE TABLE uniref90 ("
                "uniref90_id TEXT PRIMARY KEY, sequence_sha256 TEXT NOT NULL, "
                "sequence_length INTEGER NOT NULL)"
            )
            batch: list[tuple[str, str, int]] = []
            processed = 0
            for processed, (header, sequence) in enumerate(iter_fasta(fasta), start=1):
                identifier = fasta_identifier(header)
                if not identifier.startswith("UniRef90_"):
                    raise ValueError(
                        f"UniRef90 FASTA identifier must begin with 'UniRef90_': {identifier}"
                    )
                if not sequence or SEQUENCE_RE.fullmatch(sequence) is None:
                    raise ValueError(f"Invalid or empty sequence for {identifier}")
                digest = hashlib.sha256(sequence.upper().encode("ascii")).hexdigest()
                batch.append((identifier, digest, len(sequence)))
                if len(batch) >= 10000:
                    cls._insert_batch(connection, batch)
                    batch.clear()
                if processed % 1_000_000 == 0:
                    LOGGER.info(
                        "UniRef90 index progress: sequences=%d elapsed_seconds=%.1f",
                        processed, time.monotonic() - started,
                    )
            if batch:
                cls._insert_batch(connection, batch)
            count = connection.execute("SELECT COUNT(*) FROM uniref90").fetchone()[0]
            if count == 0:
                raise ValueError(f"UniRef90 FASTA contains no sequences: {fasta}")
            connection.commit()
        finally:
            connection.close()
        LOGGER.info(
            "UniRef90 index completed: sequences=%d elapsed_seconds=%.1f",
            count, time.monotonic() - started,
        )
        return cls(database)

    @staticmethod
    def _insert_batch(connection: sqlite3.Connection, rows: list[tuple[str, str, int]]) -> None:
        try:
            connection.executemany(
                "INSERT INTO uniref90 VALUES (?, ?, ?)", rows
            )
        except sqlite3.IntegrityError as exc:
            raise ValueError("Duplicate UniRef90 FASTA identifier") from exc

    def count(self) -> int:
        with sqlite3.connect(self.database) as connection:
            return int(connection.execute("SELECT COUNT(*) FROM uniref90").fetchone()[0])

    def contains(self, uniref90_id: str) -> bool:
        with sqlite3.connect(self.database) as connection:
            row = connection.execute(
                "SELECT 1 FROM uniref90 WHERE uniref90_id=?", (uniref90_id,)
            ).fetchone()
        return row is not None

    def present_ids(self, uniref90_ids: set[str]) -> set[str]:
        """Return present IDs using one disk-backed join instead of one connection per ID."""
        if not uniref90_ids:
            return set()
        with sqlite3.connect(self.database) as connection:
            connection.execute("CREATE TEMP TABLE requested (uniref90_id TEXT PRIMARY KEY)")
            connection.executemany(
                "INSERT INTO requested VALUES (?)",
                ((identifier,) for identifier in sorted(uniref90_ids)),
            )
            return {
                str(row[0])
                for row in connection.execute(
                    "SELECT r.uniref90_id FROM requested r "
                    "JOIN uniref90 u ON u.uniref90_id=r.uniref90_id"
                )
            }

    def sequence_hash(self, uniref90_id: str) -> str:
        with sqlite3.connect(self.database) as connection:
            row = connection.execute(
                "SELECT sequence_sha256 FROM uniref90 WHERE uniref90_id=?", (uniref90_id,)
            ).fetchone()
        if row is None:
            raise KeyError(uniref90_id)
        return str(row[0])

    def iter_records(self) -> Iterator[tuple[str, str, int]]:
        with sqlite3.connect(self.database) as connection:
            for row in connection.execute(
                "SELECT uniref90_id, sequence_sha256, sequence_length "
                "FROM uniref90 ORDER BY uniref90_id"
            ):
                yield str(row[0]), str(row[1]), int(row[2])
