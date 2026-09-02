from __future__ import annotations

import csv
from dataclasses import dataclass
import logging
from pathlib import Path
import sqlite3
import time

from .inputs import open_text
from .mmseqs import ClusterIndex


LOGGER = logging.getLogger(__name__)


@dataclass(frozen=True)
class MultilinkageStats:
    raw_rows: int
    duplicate_rows: int
    members: int
    clusters: int
    uniparc_members: int

    def as_dict(self) -> dict[str, int]:
        return {
            "raw_rows": self.raw_rows,
            "duplicate_rows": self.duplicate_rows,
            "members": self.members,
            "clusters": self.clusters,
            "uniparc_members": self.uniparc_members,
        }


def build_multilinkage_index(
    memberships_csv: Path, database: Path
) -> tuple[ClusterIndex, MultilinkageStats]:
    """Load Daniel's unordered clique-family CSV into the existing cluster index."""
    source = memberships_csv.expanduser().resolve()
    if not source.is_file() or source.stat().st_size == 0:
        raise ValueError(f"Multi-linkage membership CSV is missing or empty: {source}")
    database.parent.mkdir(parents=True, exist_ok=True)
    database.unlink(missing_ok=True)
    started = time.monotonic()
    connection = sqlite3.connect(database)
    raw_rows = 0
    try:
        connection.execute("PRAGMA journal_mode=OFF")
        connection.execute("PRAGMA synchronous=OFF")
        connection.execute(
            "CREATE TABLE assignments ("
            "member_id TEXT PRIMARY KEY, cluster_id TEXT NOT NULL, conflict INTEGER NOT NULL DEFAULT 0)"
        )
        connection.execute(
            "CREATE TABLE representatives ("
            "cluster_id TEXT PRIMARY KEY, representative_id TEXT NOT NULL, "
            "conflict INTEGER NOT NULL DEFAULT 0)"
        )
        assignment_sql = (
            "INSERT INTO assignments(member_id, cluster_id) VALUES (?, ?) "
            "ON CONFLICT(member_id) DO UPDATE SET conflict=1 "
            "WHERE assignments.cluster_id != excluded.cluster_id"
        )
        representative_sql = (
            "INSERT INTO representatives(cluster_id, representative_id) VALUES (?, ?) "
            "ON CONFLICT(cluster_id) DO UPDATE SET conflict=1 "
            "WHERE representatives.representative_id != excluded.representative_id"
        )
        assignment_batch: list[tuple[str, str]] = []
        representative_batch: list[tuple[str, str]] = []
        with open_text(source) as handle:
            for line_number, columns in enumerate(csv.reader(handle), start=1):
                if not columns or all(not value.strip() for value in columns):
                    continue
                if len(columns) != 3 or any(not value.strip() for value in columns):
                    raise ValueError(
                        f"Malformed multi-linkage row at {source}:{line_number}; expected "
                        "cluster_id,representative_id,member_id"
                    )
                cluster_id, representative_id, member_id = (
                    value.strip() for value in columns
                )
                if not cluster_id.isdecimal():
                    raise ValueError(
                        f"Malformed multi-linkage row at {source}:{line_number}; "
                        "cluster_id must be numeric"
                    )
                raw_rows += 1
                assignment_batch.append((member_id, cluster_id))
                representative_batch.append((cluster_id, representative_id))
                if len(assignment_batch) >= 10_000:
                    connection.executemany(assignment_sql, assignment_batch)
                    connection.executemany(representative_sql, representative_batch)
                    assignment_batch.clear()
                    representative_batch.clear()
                if raw_rows % 1_000_000 == 0:
                    LOGGER.info(
                        "Multi-linkage index progress: rows=%d elapsed_seconds=%.1f",
                        raw_rows,
                        time.monotonic() - started,
                    )
        if assignment_batch:
            connection.executemany(assignment_sql, assignment_batch)
            connection.executemany(representative_sql, representative_batch)
        if raw_rows == 0:
            raise ValueError(f"Multi-linkage membership CSV has no data rows: {source}")
        member_conflicts = int(
            connection.execute(
                "SELECT COUNT(*) FROM assignments WHERE conflict != 0"
            ).fetchone()[0]
        )
        representative_conflicts = int(
            connection.execute(
                "SELECT COUNT(*) FROM representatives WHERE conflict != 0"
            ).fetchone()[0]
        )
        if member_conflicts or representative_conflicts:
            raise ValueError(
                "Invalid multi-linkage membership: "
                f"members_in_multiple_clusters={member_conflicts}, "
                f"clusters_with_multiple_representatives={representative_conflicts}"
            )
        connection.execute(
            "CREATE INDEX assignments_cluster_member_idx "
            "ON assignments(cluster_id, member_id)"
        )
        members = int(connection.execute("SELECT COUNT(*) FROM assignments").fetchone()[0])
        clusters = int(
            connection.execute("SELECT COUNT(*) FROM representatives").fetchone()[0]
        )
        uniparc_members = int(
            connection.execute(
                "SELECT COUNT(*) FROM assignments WHERE member_id LIKE 'UPI%'"
            ).fetchone()[0]
        )
        connection.commit()
    finally:
        connection.close()
    stats = MultilinkageStats(
        raw_rows=raw_rows,
        duplicate_rows=raw_rows - members,
        members=members,
        clusters=clusters,
        uniparc_members=uniparc_members,
    )
    LOGGER.info(
        "Multi-linkage index completed: members=%d clusters=%d duplicates=%d elapsed_seconds=%.1f",
        members,
        clusters,
        stats.duplicate_rows,
        time.monotonic() - started,
    )
    return ClusterIndex(database), stats
