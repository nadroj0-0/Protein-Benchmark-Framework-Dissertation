from __future__ import annotations

import csv
from dataclasses import dataclass
import logging
import os
from pathlib import Path
import re
import shlex
import shutil
import sqlite3
import subprocess
import time
from typing import Iterable, Iterator

from .config import BuildConfig
from .inputs import open_text, sha256_file
from .uniref import UniRefIndex


LOGGER = logging.getLogger(__name__)

GIT_COMMIT_RE = re.compile(r"[0-9a-fA-F]{40}")


@dataclass(frozen=True)
class CommandSpec:
    stage: str
    argv: tuple[str, ...]

    @property
    def display(self) -> str:
        return shlex.join(self.argv)


@dataclass(frozen=True)
class MMseqsRuntime:
    requested_executable: str
    resolved_executable: str | None
    observed_version: str | None
    version_token: str | None
    version_exit_code: int | None
    executable_sha256: str | None

    def as_dict(self, expected_version: str | None) -> dict[str, object]:
        return {
            "expected_version": expected_version,
            "observed_version": self.observed_version,
            "observed_version_token": self.version_token,
            "requested_executable": self.requested_executable,
            "resolved_executable": self.resolved_executable,
            "version_exit_code": self.version_exit_code,
            "executable_sha256": self.executable_sha256,
        }


def build_mmseqs_commands(config: BuildConfig, fasta: Path, work_dir: Path) -> tuple[CommandSpec, ...]:
    database = work_dir / "uniref90_db"
    cluster_database = work_dir / "clusters_db"
    temporary = work_dir / "mmseqs_tmp"
    cluster_tsv = work_dir / "uniref90_clusters.tsv"
    return (
        CommandSpec("createdb", (
            config.mmseqs_bin, "createdb", str(fasta), str(database),
            "--dbtype", "1", "--shuffle", "0",
        )),
        CommandSpec("cluster", (
            config.mmseqs_bin, "cluster", str(database), str(cluster_database), str(temporary),
            "--min-seq-id", f"{config.identity:.2f}",
            "-c", f"{config.coverage:.1f}",
            "--cov-mode", str(config.cov_mode),
            "--cluster-mode", str(config.cluster_mode),
            "--alignment-mode", str(config.alignment_mode),
            "--seq-id-mode", "0",
            "--cluster-reassign", str(config.cluster_reassign),
            "-s", "7.5",
            "-e", "1e-4",
            "--threads", str(config.threads),
        )),
        CommandSpec("createtsv", (
            config.mmseqs_bin, "createtsv", str(database), str(database),
            str(cluster_database), str(cluster_tsv),
        )),
    )


def _version_token(version_text: str) -> str:
    stripped = version_text.strip()
    if GIT_COMMIT_RE.fullmatch(stripped):
        return stripped.lower()
    matches = set(re.findall(r"(?<![A-Za-z0-9])(\d{1,2}(?:[-.][A-Za-z0-9]+)+)(?![A-Za-z0-9])", version_text))
    if len(matches) != 1:
        raise ValueError(
            "Could not parse one exact MMseqs2 release token from version output; "
            f"observed {version_text!r}"
        )
    return next(iter(matches))


def _matches_expected_release(expected: str, observed: str) -> bool:
    if expected == observed:
        return True
    release = re.fullmatch(r"\d{1,2}[-.]([0-9a-fA-F]{5,40})", expected)
    return bool(
        release
        and GIT_COMMIT_RE.fullmatch(observed)
        and observed.lower().startswith(release.group(1).lower())
    )


def resolve_mmseqs_runtime(binary: str) -> MMseqsRuntime:
    candidate = Path(binary).expanduser() if "/" in binary else None
    located = str(candidate) if candidate is not None else shutil.which(binary)
    if located is None:
        return MMseqsRuntime(binary, None, None, None, None, None)
    resolved = Path(located).resolve()
    if not resolved.is_file() or not os.access(resolved, os.X_OK):
        return MMseqsRuntime(binary, str(resolved), None, None, None, None)
    try:
        result = subprocess.run(
            [str(resolved), "version"], text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            check=False,
        )
    except FileNotFoundError:
        return MMseqsRuntime(binary, None, None, None, None, None)
    output = result.stdout.strip()
    token = None
    if result.returncode == 0 and output:
        try:
            token = _version_token(output)
        except ValueError:
            token = None
    digest = sha256_file(resolved) if os.access(resolved, os.R_OK) else None
    return MMseqsRuntime(
        binary, str(resolved), output or None, token, result.returncode, digest
    )


def mmseqs_version(binary: str) -> str:
    runtime = resolve_mmseqs_runtime(binary)
    if runtime.resolved_executable is None:
        return "not-found"
    if runtime.version_exit_code != 0:
        return f"exit-{runtime.version_exit_code}"
    return runtime.observed_version or "empty-output"


def validate_mmseqs_version(version_text: str, minimum_release: int = 12) -> int:
    """Require a release containing the cluster-reassign/cov-mode fix."""
    match = re.search(r"(?i)(?:version\s*[:=]?\s*)?(\d{1,2})(?=[.\-])", version_text)
    if match is None:
        raise ValueError(
            "Could not parse the MMseqs2 release from version output; "
            f"observed {version_text!r}. Pin a modern release >= {minimum_release}."
        )
    release = int(match.group(1))
    if release < minimum_release:
        raise ValueError(
            f"MMseqs2 release {release} is too old; release >= {minimum_release} is required "
            "for correct --cluster-reassign coverage-mode handling"
        )
    return release


def validate_exact_mmseqs_version(expected: str, runtime: MMseqsRuntime) -> str:
    if runtime.resolved_executable is None:
        raise FileNotFoundError(f"MMseqs2 executable is unavailable: {runtime.requested_executable}")
    if runtime.version_exit_code != 0:
        raise ValueError(
            f"MMseqs2 version probe exited {runtime.version_exit_code}; exact production version is unknown"
        )
    if runtime.observed_version is None or runtime.version_token is None:
        raise ValueError("MMseqs2 version output is empty or unparseable")
    expected_token = _version_token(expected.strip())
    validate_mmseqs_version(expected_token)
    if not _matches_expected_release(expected_token, runtime.version_token):
        raise ValueError(
            "MMseqs2 exact version mismatch: "
            f"expected release {expected_token!r}, "
            f"observed binary identity {runtime.version_token!r}"
        )
    return runtime.version_token


def validate_recorded_exact_mmseqs_version(expected: str, observed: str) -> str:
    """Validate immutable recorded tokens without requiring the old executable path to exist."""
    expected_token = _version_token(expected.strip())
    observed_token = _version_token(observed.strip())
    if observed.strip() != observed_token:
        raise ValueError(
            "Recorded observed MMseqs2 version must be exactly one version identity"
        )
    validate_mmseqs_version(expected_token)
    if not _matches_expected_release(expected_token, observed_token):
        raise ValueError(
            "MMseqs2 exact version mismatch in publication metadata: "
            f"expected release {expected_token!r}, "
            f"observed binary identity {observed_token!r}"
        )
    return observed_token


def verify_mmseqs_executable_unchanged(runtime: MMseqsRuntime) -> None:
    if runtime.executable_sha256 is None or runtime.resolved_executable is None:
        return
    path = Path(runtime.resolved_executable)
    if not path.is_file() or sha256_file(path) != runtime.executable_sha256:
        raise ValueError("The resolved MMseqs2 executable changed while the run was in progress")


def execute_commands(commands: Iterable[CommandSpec], log_dir: Path) -> None:
    log_dir.mkdir(parents=True, exist_ok=True)
    for command in commands:
        log_path = log_dir / f"mmseqs_{command.stage}.log"
        started = time.monotonic()
        LOGGER.info("MMseqs2 stage started: %s; log=%s", command.stage, log_path)
        with log_path.open("w", encoding="utf-8") as log:
            result = subprocess.run(
                list(command.argv), stdout=log, stderr=subprocess.STDOUT, text=True, check=False,
            )
        if result.returncode != 0:
            raise RuntimeError(
                f"MMseqs2 stage {command.stage!r} failed with exit code {result.returncode}; "
                f"see {log_path}"
            )
        LOGGER.info(
            "MMseqs2 stage completed: %s elapsed_seconds=%.1f",
            command.stage, time.monotonic() - started,
        )


def write_command_manifest(path: Path, commands: Iterable[CommandSpec]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["stage", "argument_index", "argument", "display_command"])
        for command in commands:
            for index, argument in enumerate(command.argv):
                writer.writerow([command.stage, index, argument, command.display if index == 0 else ""])


class ClusterIndex:
    """Disk-backed, validated representative/member assignment table."""

    def __init__(self, database: Path):
        self.database = database

    @classmethod
    def build(
        cls,
        cluster_tsv: Path,
        uniref: UniRefIndex,
        database: Path,
        *,
        has_header: bool = False,
    ) -> "ClusterIndex":
        started = time.monotonic()
        LOGGER.info("MMseqs2 assignment index started: %s -> %s", cluster_tsv, database)
        if not cluster_tsv.is_file() or cluster_tsv.stat().st_size == 0:
            raise ValueError(f"MMseqs2 cluster TSV is missing or empty: {cluster_tsv}")
        database.parent.mkdir(parents=True, exist_ok=True)
        database.unlink(missing_ok=True)
        connection = sqlite3.connect(database)
        try:
            # Disposable scratch index: avoid a database-sized WAL; any interrupted load is rebuilt.
            connection.execute("PRAGMA journal_mode=OFF")
            connection.execute("PRAGMA synchronous=OFF")
            connection.execute(
                "CREATE TABLE assignments ("
                "member_id TEXT PRIMARY KEY, cluster_id TEXT NOT NULL)"
            )
            batch: list[tuple[str, str]] = []
            first_content = True
            with open_text(cluster_tsv) as handle:
                for line_number, raw_line in enumerate(handle, start=1):
                    if line_number % 1_000_000 == 0:
                        LOGGER.info(
                            "MMseqs2 assignment index progress: rows=%d elapsed_seconds=%.1f",
                            line_number, time.monotonic() - started,
                        )
                    if not raw_line.strip():
                        continue
                    columns = raw_line.rstrip("\n\r").split("\t")
                    if has_header and first_content:
                        first_content = False
                        if columns != ["mmseqs_cluster_id", "uniref90_id"]:
                            raise ValueError(
                                "Cached MMseqs2 membership has an unexpected header"
                            )
                        continue
                    first_content = False
                    if len(columns) != 2 or not columns[0] or not columns[1]:
                        raise ValueError(
                            f"Malformed MMseqs2 cluster row at {cluster_tsv}:{line_number}; "
                            "expected representative<TAB>member"
                        )
                    batch.append((columns[1], columns[0]))
                    if len(batch) >= 10000:
                        cls._insert(connection, batch)
                        batch.clear()
            if batch:
                cls._insert(connection, batch)
            connection.commit()
            connection.execute(
                "CREATE INDEX assignments_cluster_member_idx "
                "ON assignments(cluster_id, member_id)"
            )
            connection.execute(
                "ATTACH DATABASE ? AS u", (str(uniref.database),)
            )
            unknown = connection.execute(
                "SELECT COUNT(*) FROM assignments a LEFT JOIN u.uniref90 u "
                "ON a.member_id=u.uniref90_id WHERE u.uniref90_id IS NULL"
            ).fetchone()[0]
            missing = connection.execute(
                "SELECT COUNT(*) FROM u.uniref90 u LEFT JOIN assignments a "
                "ON a.member_id=u.uniref90_id WHERE a.member_id IS NULL"
            ).fetchone()[0]
            representative_missing = connection.execute(
                "SELECT COUNT(DISTINCT a.cluster_id) FROM assignments a "
                "LEFT JOIN assignments self ON self.member_id=a.cluster_id "
                "AND self.cluster_id=a.cluster_id WHERE self.member_id IS NULL"
            ).fetchone()[0]
            if unknown or missing or representative_missing:
                raise ValueError(
                    "Invalid MMseqs2 mapping: "
                    f"unknown_members={unknown}, missing_members={missing}, "
                    f"clusters_without_representative_self_row={representative_missing}"
                )
            indexed_members = int(
                connection.execute("SELECT COUNT(*) FROM assignments").fetchone()[0]
            )
            indexed_clusters = int(
                connection.execute(
                    "SELECT COUNT(DISTINCT cluster_id) FROM assignments"
                ).fetchone()[0]
            )
        finally:
            connection.close()
        LOGGER.info(
            "MMseqs2 assignment index completed: members=%d clusters=%d elapsed_seconds=%.1f",
            indexed_members, indexed_clusters,
            time.monotonic() - started,
        )
        return cls(database)

    @staticmethod
    def _insert(connection: sqlite3.Connection, rows: list[tuple[str, str]]) -> None:
        try:
            connection.executemany("INSERT INTO assignments VALUES (?, ?)", rows)
        except sqlite3.IntegrityError as exc:
            raise ValueError("A UniRef90 member appears more than once in MMseqs2 output") from exc

    def member_count(self) -> int:
        with sqlite3.connect(self.database) as connection:
            return int(connection.execute("SELECT COUNT(*) FROM assignments").fetchone()[0])

    def cluster_count(self) -> int:
        with sqlite3.connect(self.database) as connection:
            return int(connection.execute("SELECT COUNT(DISTINCT cluster_id) FROM assignments").fetchone()[0])

    def cluster_for(self, member_id: str) -> str | None:
        with sqlite3.connect(self.database) as connection:
            row = connection.execute(
                "SELECT cluster_id FROM assignments WHERE member_id=?", (member_id,)
            ).fetchone()
        return str(row[0]) if row else None

    def clusters_for(self, member_ids: set[str]) -> dict[str, str]:
        """Resolve many members with one SQLite connection and join."""
        if not member_ids:
            return {}
        with sqlite3.connect(self.database) as connection:
            connection.execute("CREATE TEMP TABLE requested (member_id TEXT PRIMARY KEY)")
            connection.executemany(
                "INSERT INTO requested VALUES (?)",
                ((member_id,) for member_id in sorted(member_ids)),
            )
            return {
                str(member): str(cluster)
                for member, cluster in connection.execute(
                    "SELECT r.member_id, a.cluster_id FROM requested r "
                    "JOIN assignments a ON a.member_id=r.member_id"
                )
            }

    def iter_assignments(self) -> Iterator[tuple[str, str]]:
        with sqlite3.connect(self.database) as connection:
            for member, cluster in connection.execute(
                "SELECT member_id, cluster_id FROM assignments ORDER BY cluster_id, member_id"
            ):
                yield str(cluster), str(member)

    def iter_assignments_with_metadata(
        self, uniref: UniRefIndex
    ) -> Iterator[tuple[str, str, str, int]]:
        with sqlite3.connect(self.database) as connection:
            connection.execute("ATTACH DATABASE ? AS u", (str(uniref.database),))
            for cluster, member, digest, length in connection.execute(
                "SELECT a.cluster_id, a.member_id, u.sequence_sha256, u.sequence_length "
                "FROM assignments a JOIN u.uniref90 u ON a.member_id=u.uniref90_id "
                "ORDER BY a.cluster_id, a.member_id"
            ):
                yield str(cluster), str(member), str(digest), int(length)

    def iter_assignments_with_metadata_for_clusters(
        self, uniref: UniRefIndex, cluster_ids: Iterable[str]
    ) -> Iterator[tuple[str, str, str, int]]:
        ordered_cluster_ids = sorted(cluster_ids)
        if not ordered_cluster_ids:
            return
        with sqlite3.connect(self.database) as connection:
            connection.execute("ATTACH DATABASE ? AS u", (str(uniref.database),))
            connection.execute("CREATE TEMP TABLE retained (cluster_id TEXT PRIMARY KEY)")
            connection.executemany(
                "INSERT INTO retained VALUES (?)",
                ((cluster_id,) for cluster_id in ordered_cluster_ids),
            )
            for cluster, member, digest, length in connection.execute(
                "SELECT a.cluster_id, a.member_id, u.sequence_sha256, u.sequence_length "
                "FROM retained r JOIN assignments a ON a.cluster_id=r.cluster_id "
                "JOIN u.uniref90 u ON a.member_id=u.uniref90_id "
                "ORDER BY a.cluster_id, a.member_id"
            ):
                yield str(cluster), str(member), str(digest), int(length)

    def cluster_sizes(self) -> dict[str, int]:
        with sqlite3.connect(self.database) as connection:
            connection.execute("PRAGMA temp_store=FILE")
            return {
                str(cluster): int(count)
                for cluster, count in connection.execute(
                    "SELECT cluster_id, COUNT(*) FROM assignments "
                    "GROUP BY cluster_id ORDER BY cluster_id"
                )
            }

    def iter_cluster_sizes(self) -> Iterator[tuple[str, int]]:
        with sqlite3.connect(self.database) as connection:
            for cluster, count in connection.execute(
                "SELECT cluster_id, COUNT(*) FROM assignments "
                "GROUP BY cluster_id ORDER BY cluster_id"
            ):
                yield str(cluster), int(count)

    def sizes_for(self, cluster_ids: set[str]) -> dict[str, int]:
        if not cluster_ids:
            return {}
        with sqlite3.connect(self.database) as connection:
            connection.execute("CREATE TEMP TABLE retained (cluster_id TEXT PRIMARY KEY)")
            connection.executemany(
                "INSERT INTO retained VALUES (?)",
                ((cluster_id,) for cluster_id in sorted(cluster_ids)),
            )
            return {
                str(cluster): int(count)
                for cluster, count in connection.execute(
                    "SELECT a.cluster_id, COUNT(*) FROM assignments a "
                    "JOIN retained r ON r.cluster_id=a.cluster_id GROUP BY a.cluster_id "
                    "ORDER BY a.cluster_id"
                )
            }

    def exact_sequence_split_conflicts(
        self, uniref: UniRefIndex, cluster_splits: dict[str, str], limit: int = 20
    ) -> list[dict[str, str]]:
        with sqlite3.connect(self.database) as connection:
            connection.execute("ATTACH DATABASE ? AS u", (str(uniref.database),))
            connection.execute(
                "CREATE TEMP TABLE cluster_splits (cluster_id TEXT PRIMARY KEY, split TEXT NOT NULL)"
            )
            connection.executemany(
                "INSERT INTO cluster_splits VALUES (?, ?)",
                sorted(cluster_splits.items()),
            )
            rows = connection.execute(
                "SELECT u.sequence_sha256, GROUP_CONCAT(DISTINCT s.split), "
                "GROUP_CONCAT(a.member_id) "
                "FROM assignments a JOIN u.uniref90 u ON a.member_id=u.uniref90_id "
                "JOIN cluster_splits s ON a.cluster_id=s.cluster_id "
                "GROUP BY u.sequence_sha256 HAVING COUNT(DISTINCT s.split) > 1 "
                "ORDER BY u.sequence_sha256 LIMIT ?",
                (limit,),
            ).fetchall()
        return [
            {"sequence_sha256": str(digest), "splits": str(splits), "uniref90_ids": str(members)}
            for digest, splits, members in rows
        ]

    def global_exact_sequence_split_conflicts(
        self,
        uniref: UniRefIndex,
        cluster_splits: dict[str, str],
        protein_sequences: Iterable[tuple[str, str, str]],
        limit: int = 20,
    ) -> list[dict[str, str]]:
        """Compare retained UniRef scaffold and mapped UniProt sequences in one digest space."""
        with sqlite3.connect(self.database) as connection:
            connection.execute("PRAGMA temp_store=FILE")
            connection.execute("ATTACH DATABASE ? AS u", (str(uniref.database),))
            connection.execute(
                "CREATE TEMP TABLE cluster_splits (cluster_id TEXT PRIMARY KEY, split TEXT NOT NULL)"
            )
            connection.executemany(
                "INSERT INTO cluster_splits VALUES (?, ?)", sorted(cluster_splits.items())
            )
            connection.execute(
                "CREATE TEMP TABLE protein_sequences ("
                "sequence_sha256 TEXT NOT NULL, source_id TEXT NOT NULL, split TEXT NOT NULL)"
            )
            connection.executemany(
                "INSERT INTO protein_sequences VALUES (?, ?, ?)", protein_sequences
            )
            rows = connection.execute(
                "SELECT sequence_sha256, GROUP_CONCAT(DISTINCT split), "
                "GROUP_CONCAT(DISTINCT CASE WHEN source_id LIKE 'UniRef90:%' THEN split END), "
                "GROUP_CONCAT(DISTINCT CASE WHEN source_id LIKE 'UniProtKB:%' THEN split END), "
                "MIN(CASE WHEN source_id LIKE 'UniRef90:%' THEN source_id END), "
                "MIN(CASE WHEN source_id LIKE 'UniProtKB:%' THEN source_id END), "
                "COUNT(*) "
                "FROM ("
                "  SELECT u.sequence_sha256 AS sequence_sha256, "
                "         'UniRef90:' || a.member_id AS source_id, s.split AS split "
                "  FROM assignments a JOIN u.uniref90 u ON a.member_id=u.uniref90_id "
                "  JOIN cluster_splits s ON a.cluster_id=s.cluster_id "
                "  UNION ALL "
                "  SELECT sequence_sha256, 'UniProtKB:' || source_id, split FROM protein_sequences"
                ") GROUP BY sequence_sha256 HAVING COUNT(DISTINCT split) > 1 "
                "ORDER BY sequence_sha256 LIMIT ?",
                (limit,),
            ).fetchall()
        return [
            {
                "sequence_sha256": str(digest),
                "splits": str(splits),
                "uniref90_splits": str(uniref_splits or ""),
                "uniprot_splits": str(uniprot_splits or ""),
                "sample_uniref90_id": str(sample_uniref or ""),
                "sample_uniprot_accession": str(sample_uniprot or ""),
                "source_count": str(source_count),
            }
            for (
                digest, splits, uniref_splits, uniprot_splits,
                sample_uniref, sample_uniprot, source_count,
            ) in rows
        ]
