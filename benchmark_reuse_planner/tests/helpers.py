from __future__ import annotations

import csv
from pathlib import Path
from typing import Dict, Iterable, Mapping, Optional, Sequence, Tuple, Union


CSV_NAMES: Tuple[str, ...] = tuple(
    f"{ontology}-{split}.csv"
    for ontology in ("bp", "cc", "mf")
    for split in ("training", "validation", "test")
)
Labels = Union[str, Sequence[str]]
Row = Union[Tuple[str, str], Tuple[str, str, Labels]]


def write_benchmark(
    directory: Path,
    rows_by_file: Optional[Mapping[str, Sequence[Row]]] = None,
    headers_by_file: Optional[Mapping[str, Sequence[str]]] = None,
) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    if rows_by_file is None:
        rows_by_file = {"bp-training.csv": [("P1", "AAAA")]}
    if headers_by_file is None:
        headers_by_file = {}
    for filename in CSV_NAMES:
        header = list(headers_by_file.get(filename, ("proteins", "sequences", "GO:0000001")))
        with (directory / filename).open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, lineterminator="\n")
            writer.writerow(header)
            for row in rows_by_file.get(filename, ()):
                protein_id, sequence = row[:2]
                if len(row) == 2:
                    labels: Iterable[str] = ("1",)
                else:
                    raw_labels = row[2]
                    labels = (raw_labels,) if isinstance(raw_labels, str) else raw_labels
                writer.writerow((protein_id, sequence, *labels))
    return directory


def rows_in(filename: str, *rows: Row) -> Dict[str, Sequence[Row]]:
    return {filename: rows}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))
