from __future__ import annotations

import csv
import hashlib
import re
from pathlib import Path
from typing import Dict, Iterable, List, MutableMapping, Sequence, Tuple

from .models import (
    BenchmarkData,
    InputFileIdentity,
    ProteinRecord,
    REQUIRED_CSV_NAMES,
    ReusePlannerError,
)


_BENCHMARK_NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
_SEQUENCE_RE = re.compile(r"^[A-Za-z*.-]+$")


class BenchmarkError(ReusePlannerError):
    pass


def validate_benchmark_name(name: str) -> None:
    if not _BENCHMARK_NAME_RE.fullmatch(name):
        raise BenchmarkError(
            "Benchmark names must match [A-Za-z0-9][A-Za-z0-9._-]*: %r" % name
        )


def parse_benchmark(name: str, directory: Path) -> BenchmarkData:
    validate_benchmark_name(name)
    directory = directory.expanduser().resolve()
    if not directory.is_dir():
        raise BenchmarkError("Benchmark directory does not exist: %s" % directory)

    missing = [filename for filename in REQUIRED_CSV_NAMES if not (directory / filename).is_file()]
    if missing:
        raise BenchmarkError(
            "Benchmark %s is missing required CSVs: %s" % (name, ", ".join(missing))
        )

    mutable_proteins: Dict[str, Tuple[str, set[str]]] = {}
    duplicate_occurrences = 0
    identities: List[InputFileIdentity] = []
    for filename in REQUIRED_CSV_NAMES:
        path = directory / filename
        before = _file_identity(path, filename)
        duplicate_occurrences += _parse_csv(path, filename, mutable_proteins)
        after = _file_identity(path, filename)
        if before != after:
            raise BenchmarkError("Input CSV changed while it was being parsed: %s" % path)
        identities.append(after)

    if not mutable_proteins:
        raise BenchmarkError("Benchmark %s contains no proteins" % name)

    proteins = {
        protein_id: ProteinRecord(
            protein_id=protein_id,
            sequence=sequence,
            sequence_sha256=_sha256_bytes(sequence.encode("utf-8")),
            memberships=tuple(sorted(memberships)),
        )
        for protein_id, (sequence, memberships) in sorted(mutable_proteins.items())
    }
    return BenchmarkData(
        name=name,
        directory=directory,
        proteins=proteins,
        input_files=tuple(identities),
        duplicate_occurrences=duplicate_occurrences,
    )


def verify_input_identities(benchmarks: Sequence[BenchmarkData]) -> None:
    for benchmark in benchmarks:
        for expected in benchmark.input_files:
            path = benchmark.directory / expected.relative_path
            if not path.is_file():
                raise BenchmarkError("Input CSV changed after planning: missing %s" % path)
            observed = _file_identity(path, expected.relative_path)
            if observed != expected:
                raise BenchmarkError("Input CSV changed after planning: %s" % path)


def _parse_csv(
    path: Path,
    filename: str,
    proteins: MutableMapping[str, Tuple[str, set[str]]],
) -> int:
    duplicate_occurrences = 0
    seen_rows: Dict[str, Tuple[str, bytes, int]] = {}
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.reader(handle, strict=True)
            header = _read_header(reader, filename)
            if len(header) < 3 or header[:2] != ["proteins", "sequences"]:
                raise BenchmarkError(
                    "%s must begin with proteins,sequences and contain GO columns" % filename
                )
            go_columns = header[2:]
            if len(go_columns) != len(set(go_columns)):
                raise BenchmarkError("%s contains duplicate GO columns" % filename)
            if any(not column.startswith("GO:") or column == "GO:" for column in go_columns):
                raise BenchmarkError("%s contains malformed GO columns" % filename)

            for line_number, row in _read_rows(reader, filename):
                if len(row) != len(header):
                    raise BenchmarkError(
                        "%s:%d has %d columns; expected %d"
                        % (filename, line_number, len(row), len(header))
                    )
                protein_id, sequence = row[:2]
                _validate_protein_id(protein_id, filename, line_number)
                _validate_sequence(sequence, protein_id, filename, line_number)
                label_digest = hashlib.sha256()
                for column_number, value in enumerate(row[2:], start=3):
                    if value not in {"0", "1"}:
                        raise BenchmarkError(
                            "%s:%d has non-binary GO label in column %d for %s"
                            % (filename, line_number, column_number, protein_id)
                        )
                    label_digest.update(value.encode("ascii"))
                row_labels_digest = label_digest.digest()

                previous_row = seen_rows.get(protein_id)
                if previous_row is not None:
                    if previous_row[:2] != (sequence, row_labels_digest):
                        raise BenchmarkError(
                            "%s has contradictory duplicate rows for %s (lines %d and %d)"
                            % (filename, protein_id, previous_row[2], line_number)
                        )
                    duplicate_occurrences += 1
                    continue
                seen_rows[protein_id] = (sequence, row_labels_digest, line_number)

                previous = proteins.get(protein_id)
                if previous is not None and previous[0] != sequence:
                    raise BenchmarkError(
                        "Protein ID %s has conflicting sequences in benchmark CSVs" % protein_id
                    )
                if previous is None:
                    proteins[protein_id] = (sequence, {filename})
                else:
                    previous[1].add(filename)
    except BenchmarkError:
        raise
    except (OSError, UnicodeError, csv.Error) as exc:
        raise BenchmarkError("Cannot parse %s: %s" % (path, exc)) from exc
    return duplicate_occurrences


def _read_header(reader: Iterable[List[str]], filename: str) -> List[str]:
    try:
        return next(iter(reader))
    except StopIteration as exc:
        raise BenchmarkError("Required CSV is empty: %s" % filename) from exc
    except csv.Error as exc:
        raise BenchmarkError("Malformed CSV %s: %s" % (filename, exc)) from exc


def _read_rows(
    reader: Iterable[List[str]], filename: str
) -> Iterable[Tuple[int, List[str]]]:
    iterator = iter(reader)
    line_number = 1
    while True:
        try:
            row = next(iterator)
        except StopIteration:
            return
        except csv.Error as exc:
            raise BenchmarkError(
                "Malformed CSV %s near record %d: %s" % (filename, line_number + 1, exc)
            ) from exc
        line_number += 1
        yield line_number, row


def _validate_protein_id(protein_id: str, filename: str, line_number: int) -> None:
    unsafe = (
        not protein_id
        or protein_id in {".", ".."}
        or any(character.isspace() for character in protein_id)
        or any(ord(character) < 32 or ord(character) == 127 for character in protein_id)
        or "/" in protein_id
        or "\\" in protein_id
    )
    if unsafe:
        raise BenchmarkError(
            "%s:%d has an empty or unsafe protein ID: %r"
            % (filename, line_number, protein_id)
        )


def _validate_sequence(
    sequence: str, protein_id: str, filename: str, line_number: int
) -> None:
    if not sequence or _SEQUENCE_RE.fullmatch(sequence) is None:
        raise BenchmarkError(
            "%s:%d has an empty or malformed sequence for %s"
            % (filename, line_number, protein_id)
        )


def _file_identity(path: Path, relative_path: str) -> InputFileIdentity:
    try:
        size_bytes = path.stat().st_size
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
    except OSError as exc:
        raise BenchmarkError("Cannot hash input CSV %s: %s" % (path, exc)) from exc
    return InputFileIdentity(
        relative_path=relative_path,
        resolved_path=path.resolve(),
        size_bytes=size_bytes,
        sha256=digest.hexdigest(),
    )


def _sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()
