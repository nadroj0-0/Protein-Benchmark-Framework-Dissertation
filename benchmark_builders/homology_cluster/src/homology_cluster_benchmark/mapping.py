from __future__ import annotations

from dataclasses import replace
import logging
from pathlib import Path
import re
import time
from typing import Iterator

from .inputs import open_text
from .models import AnnotationRecord, GoaLoadResult, ProteinCatalog, ProteinRecord
from .uniref import SEQUENCE_RE, fasta_identifier, iter_fasta


LOGGER = logging.getLogger(__name__)


OX_RE = re.compile(r"\bOX=(\d+)\b")
DAT_TAXON_RE = re.compile(r"NCBI_TaxID=(\d+)")


def _fasta_accession(header: str) -> str:
    first = fasta_identifier(header)
    fields = first.split("|")
    if len(fields) >= 3 and fields[0] in {"sp", "tr"}:
        return fields[1]
    return first


def _iter_uniprot_fasta(path: Path) -> Iterator[ProteinRecord]:
    for header, sequence in iter_fasta(path):
        protein_id = _fasta_accession(header)
        if not sequence or SEQUENCE_RE.fullmatch(sequence) is None:
            raise ValueError(f"Invalid or empty UniProt sequence for {protein_id}")
        taxon = OX_RE.search(header)
        yield ProteinRecord(
            protein_id=protein_id,
            sequence=sequence.upper(),
            taxon_id=taxon.group(1) if taxon else "",
            aliases=(protein_id,),
            source_header=header,
        )


def _iter_uniprot_dat(path: Path) -> Iterator[ProteinRecord]:
    entry_name = ""
    accessions: list[str] = []
    taxon = ""
    in_sequence = False
    sequence_parts: list[str] = []
    with open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\n\r")
            tag = line[:2]
            if tag == "ID":
                if entry_name or accessions or in_sequence:
                    raise ValueError(
                        f"Nested UniProt DAT ID before // terminator at {path}:{line_number}"
                    )
                fields = line.split()
                entry_name = fields[1] if len(fields) > 1 else ""
            elif tag == "AC":
                accessions.extend(item.strip() for item in line[5:].split(";") if item.strip())
            elif tag == "OX":
                match = DAT_TAXON_RE.search(line)
                if match:
                    taxon = match.group(1)
            elif tag == "SQ":
                in_sequence = True
            elif line == "//":
                sequence = "".join(sequence_parts).upper()
                if not accessions:
                    raise ValueError(f"UniProt DAT record ending at {path}:{line_number} has no accession")
                if not sequence or SEQUENCE_RE.fullmatch(sequence) is None:
                    raise ValueError(f"Invalid or empty UniProt DAT sequence for {accessions[0]}")
                yield ProteinRecord(
                    protein_id=accessions[0], sequence=sequence, taxon_id=taxon,
                    aliases=tuple(accessions), source_header=entry_name,
                )
                entry_name = ""
                accessions = []
                taxon = ""
                in_sequence = False
                sequence_parts = []
            elif in_sequence:
                sequence_parts.append("".join(character for character in line if character.isalpha()))
    if entry_name or accessions or in_sequence or sequence_parts:
        raise ValueError(f"Unterminated UniProt DAT record at end of {path}")


def iter_uniprot(path: Path) -> Iterator[ProteinRecord]:
    suffixes = "".join(path.suffixes)
    if suffixes.endswith(".dat") or suffixes.endswith(".dat.gz"):
        yield from _iter_uniprot_dat(path)
    else:
        yield from _iter_uniprot_fasta(path)


def load_requested_proteins(path: Path, requested_accessions: set[str]) -> ProteinCatalog:
    """Stream UniProt and retain only records touching a qualifying GOA ID."""
    started = time.monotonic()
    LOGGER.info(
        "UniProt sequence scan started: %s requested_accessions=%d",
        path, len(requested_accessions),
    )
    catalog = ProteinCatalog()
    seen_requested: dict[str, str] = {}
    processed = 0
    for processed, record in enumerate(iter_uniprot(path), start=1):
        if processed % 1_000_000 == 0:
            LOGGER.info(
                "UniProt sequence progress: records=%d retained=%d elapsed_seconds=%.1f",
                processed, len(catalog.records), time.monotonic() - started,
            )
        aliases = set(record.aliases) | {record.protein_id}
        matches = aliases & requested_accessions
        if not matches:
            continue
        existing = catalog.records.get(record.protein_id)
        if existing is not None and existing.sequence != record.sequence:
            raise ValueError(f"Conflicting sequences for UniProt accession {record.protein_id}")
        catalog.records.setdefault(record.protein_id, record)
        for alias in sorted(aliases):
            previous = catalog.alias_to_primary.get(alias)
            if previous is None or previous == record.protein_id:
                catalog.alias_to_primary[alias] = record.protein_id
            else:
                previous_record = catalog.records.get(previous)
                if previous_record is not None and previous_record.sequence != record.sequence:
                    raise ValueError(
                        f"Conflicting sequences for duplicate UniProt alias {alias}: "
                        f"{previous} versus {record.protein_id}"
                    )
                catalog.ambiguous_aliases.add(alias)
                catalog.alias_to_primary.pop(alias, None)
        for requested in matches:
            previous = seen_requested.get(requested)
            if previous is not None and previous != record.protein_id:
                previous_record = catalog.records.get(previous)
                if previous_record is not None and previous_record.sequence != record.sequence:
                    raise ValueError(
                        f"Conflicting sequences for requested UniProt accession {requested}: "
                        f"{previous} versus {record.protein_id}"
                    )
                catalog.ambiguous_aliases.add(requested)
                catalog.alias_to_primary.pop(requested, None)
            else:
                seen_requested[requested] = record.protein_id
    LOGGER.info(
        "UniProt sequence scan completed: records=%d retained=%d ambiguous_aliases=%d "
        "elapsed_seconds=%.1f",
        processed, len(catalog.records), len(catalog.ambiguous_aliases),
        time.monotonic() - started,
    )
    return catalog


def canonicalize_goa_accessions(result: GoaLoadResult, catalog: ProteinCatalog) -> GoaLoadResult:
    records: list[AnnotationRecord] = []
    annotations: dict[str, set[str]] = {}
    canonical_accessions: set[str] = set()
    accession_map: dict[str, tuple[str, str]] = {}
    raw_accessions = result.qualifying_accessions or set(result.annotations)
    for raw in sorted(raw_accessions):
        if raw in catalog.ambiguous_aliases:
            protein_id = raw
            action = "ambiguous-secondary"
        else:
            protein_id = catalog.alias_to_primary.get(raw, raw)
            action = "exact" if protein_id == raw else "secondary-to-primary"
        accession_map[raw] = (protein_id, action)
        canonical_accessions.add(protein_id)
        if result.record_spool is None:
            annotations.setdefault(protein_id, set()).update(result.annotations[raw])

    if result.record_spool is not None:
        result.records = []
        result.annotations = annotations
        result.qualifying_accessions = canonical_accessions
        result.accession_map = accession_map
        return result

    for record in result.records:
        raw = record.raw_accession
        protein_id, action = accession_map[raw]
        updated = replace(record, protein_id=protein_id, accession_action=action)
        records.append(updated)
    result.records = records
    result.annotations = annotations
    result.qualifying_accessions = canonical_accessions
    result.accession_map = accession_map
    return result
