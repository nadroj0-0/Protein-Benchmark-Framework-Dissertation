import csv
import hashlib
import re
from collections import defaultdict
from pathlib import Path
from typing import DefaultDict, Dict, Iterable, List, Set, Tuple

from .models import (
    AliasEntry,
    BenchmarkContract,
    BenchmarkData,
    MODALITIES,
    ONTOLOGIES,
    SPLITS,
    ProteinRecord,
)


class BenchmarkError(ValueError):
    pass


def required_csv_names() -> List[str]:
    return ["%s-%s.csv" % (ontology.lower(), split) for ontology in ONTOLOGIES for split in SPLITS]


def parse_benchmark(directory: Path, contract: BenchmarkContract) -> BenchmarkData:
    directory = directory.resolve()
    missing = [name for name in required_csv_names() if not (directory / name).is_file()]
    if missing:
        raise BenchmarkError(
            "Benchmark is missing required CSVs: %s" % ", ".join(missing)
        )

    id_re = re.compile(contract.protein_id_pattern)
    sequence_re = re.compile(contract.sequence_pattern)
    proteins: Dict[str, ProteinRecord] = {}
    file_members: Dict[Tuple[str, str], Set[str]] = {}
    duplicate_rows = 0

    for ontology in ONTOLOGIES:
        for split in SPLITS:
            filename = "%s-%s.csv" % (ontology.lower(), split)
            path = directory / filename
            members: Set[str] = set()
            seen_rows: Dict[str, Tuple[str, Tuple[str, ...]]] = {}
            with path.open("r", newline="", encoding="utf-8-sig") as handle:
                reader = csv.reader(handle)
                try:
                    header = next(reader)
                except StopIteration as exc:
                    raise BenchmarkError("Required CSV is empty: %s" % path) from exc
                if len(header) < 3 or header[:2] != ["proteins", "sequences"]:
                    raise BenchmarkError(
                        "%s must begin with proteins,sequences and contain GO columns" % filename
                    )
                go_columns = header[2:]
                if len(go_columns) != len(set(go_columns)):
                    raise BenchmarkError("%s contains duplicate GO columns" % filename)
                if any(not value.startswith("GO:") for value in go_columns):
                    raise BenchmarkError("%s contains malformed GO columns" % filename)

                for line_number, row in enumerate(reader, start=2):
                    if len(row) != len(header):
                        raise BenchmarkError(
                            "%s:%d has %d columns; expected %d"
                            % (filename, line_number, len(row), len(header))
                        )
                    protein_id, sequence = row[0], row[1]
                    if (
                        not protein_id
                        or not _safe_protein_id(protein_id)
                        or id_re.fullmatch(protein_id) is None
                    ):
                        raise BenchmarkError(
                            "%s:%d has an empty or malformed protein ID: %r"
                            % (filename, line_number, protein_id)
                        )
                    if not sequence or sequence_re.fullmatch(sequence) is None:
                        raise BenchmarkError(
                            "%s:%d has an empty or malformed sequence for %s"
                            % (filename, line_number, protein_id)
                        )
                    labels = tuple(row[2:])
                    if any(value not in {"0", "1"} for value in labels):
                        raise BenchmarkError(
                            "%s:%d has non-binary GO labels for %s"
                            % (filename, line_number, protein_id)
                        )

                    previous_row = seen_rows.get(protein_id)
                    if previous_row is not None:
                        if previous_row != (sequence, labels):
                            raise BenchmarkError(
                                "%s has contradictory duplicate rows for %s" % (filename, protein_id)
                            )
                        duplicate_rows += 1
                        continue
                    seen_rows[protein_id] = (sequence, labels)

                    existing = proteins.get(protein_id)
                    if existing is not None and existing.sequence != sequence:
                        raise BenchmarkError(
                            "Protein ID %s has conflicting sequences in %s and %s"
                            % (protein_id, ",".join(sorted(existing.source_files)), filename)
                        )
                    if existing is None:
                        existing = ProteinRecord(
                            protein_id=protein_id,
                            sequence=sequence,
                            sequence_sha256=hashlib.sha256(sequence.encode("utf-8")).hexdigest(),
                            sequence_length=len(sequence),
                        )
                        proteins[protein_id] = existing
                    existing.ontologies.add(ontology)
                    existing.splits.add(split)
                    existing.source_files.add(filename)
                    existing.memberships.add((ontology, split))
                    members.add(protein_id)
            file_members[(ontology, split)] = members

    if not proteins:
        raise BenchmarkError("Benchmark contains no proteins")
    _validate_overlap(proteins, file_members, contract)
    return BenchmarkData(
        directory=directory,
        proteins=proteins,
        file_members=file_members,
        duplicate_rows=duplicate_rows,
    )


def load_aliases(
    path: Path, protein_id_pattern: str = r"^[^\s/\\]+$"
) -> Dict[Tuple[str, str], List[AliasEntry]]:
    aliases: DefaultDict[Tuple[str, str], List[AliasEntry]] = defaultdict(list)
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "protein_id",
            "source_protein_id",
            "modality",
            "mapping_route",
            "source_identity",
            "mapping_evidence",
        }
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise BenchmarkError(
                "Alias TSV must contain: %s" % ", ".join(sorted(required))
            )
        for line_number, row in enumerate(reader, start=2):
            protein_id = (row.get("protein_id") or "").strip()
            source_id = (row.get("source_protein_id") or "").strip()
            modality = (row.get("modality") or "").strip().lower()
            route = (row.get("mapping_route") or "").strip()
            source_identity = (row.get("source_identity") or "").strip()
            mapping_evidence = (row.get("mapping_evidence") or "").strip()
            if not protein_id or not source_id or not route or not source_identity or not mapping_evidence:
                raise BenchmarkError("Alias TSV line %d contains an empty required value" % line_number)
            id_re = re.compile(protein_id_pattern)
            if (
                not _safe_protein_id(protein_id)
                or not _safe_protein_id(source_id)
                or id_re.fullmatch(protein_id) is None
                or id_re.fullmatch(source_id) is None
            ):
                raise BenchmarkError(
                    "Alias TSV line %d contains an unsafe or malformed protein ID" % line_number
                )
            target_modalities: Iterable[str]
            if modality in {"*", "all"}:
                target_modalities = MODALITIES
            elif modality in MODALITIES:
                target_modalities = (modality,)
            else:
                raise BenchmarkError("Alias TSV line %d has unknown modality %r" % (line_number, modality))
            for target_modality in target_modalities:
                entry = AliasEntry(
                    protein_id,
                    source_id,
                    target_modality,
                    route,
                    source_identity,
                    mapping_evidence,
                )
                if entry not in aliases[(protein_id, target_modality)]:
                    aliases[(protein_id, target_modality)].append(entry)
    return dict(aliases)


def _safe_protein_id(protein_id: str) -> bool:
    return (
        protein_id not in {".", ".."}
        and not any(character.isspace() for character in protein_id)
        and "/" not in protein_id
        and "\\" not in protein_id
    )


def sequence_index(benchmark: BenchmarkData) -> Dict[str, List[str]]:
    index: DefaultDict[str, List[str]] = defaultdict(list)
    for protein_id, protein in benchmark.proteins.items():
        index[protein.sequence_sha256].append(protein_id)
    return {digest: sorted(ids) for digest, ids in index.items()}


def temporal_text_role(protein: ProteinRecord) -> str:
    current = bool(protein.splits & {"training", "validation"})
    test = "test" in protein.splits
    if current and test:
        return "mixed-current-and-test"
    if current:
        return "current-train-validation"
    if test:
        return "historical-test"
    return "unknown"


def _validate_overlap(
    proteins: Dict[str, ProteinRecord],
    file_members: Dict[Tuple[str, str], Set[str]],
    contract: BenchmarkContract,
) -> None:
    _check_overlap_kind(proteins, file_members, contract.id_overlap, use_sequences=False)
    _check_overlap_kind(proteins, file_members, contract.sequence_overlap, use_sequences=True)


def _check_overlap_kind(
    proteins: Dict[str, ProteinRecord],
    file_members: Dict[Tuple[str, str], Set[str]],
    policy: str,
    use_sequences: bool,
) -> None:
    if policy == "allow":
        return
    label = "exact sequences" if use_sequences else "protein IDs"
    pairs = (("training", "validation"), ("training", "test"), ("validation", "test"))
    groups: List[Tuple[str, Dict[str, Set[str]]]] = []
    if policy == "global-disjoint":
        groups.append(
            (
                "global",
                {
                    split: set().union(*(file_members[(ontology, split)] for ontology in ONTOLOGIES))
                    for split in SPLITS
                },
            )
        )
    else:
        for ontology in ONTOLOGIES:
            groups.append(
                (ontology, {split: set(file_members[(ontology, split)]) for split in SPLITS})
            )

    for group_name, split_ids in groups:
        values: Dict[str, Set[str]] = {}
        for split, ids in split_ids.items():
            if use_sequences:
                values[split] = {proteins[protein_id].sequence_sha256 for protein_id in ids}
            else:
                values[split] = ids
        for left, right in pairs:
            overlap = values[left] & values[right]
            if overlap:
                example = sorted(overlap)[0]
                raise BenchmarkError(
                    "Benchmark contract %s violated: %s %s overlap between %s and %s "
                    "(%d values; example %s)"
                    % (policy, group_name, label, left, right, len(overlap), example)
                )
