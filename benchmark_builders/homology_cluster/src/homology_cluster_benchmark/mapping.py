from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import replace
import hashlib
import logging
from pathlib import Path
import re
import sqlite3
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


def _iter_uniprot_fasta(path: Path, source_population: str) -> Iterator[ProteinRecord]:
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
            source_population=source_population,
        )


def _iter_uniprot_dat(path: Path, source_population: str) -> Iterator[ProteinRecord]:
    entry_name = ""
    accessions: list[str] = []
    taxon = ""
    in_sequence = False
    sequence_parts: list[str] = []
    expected_review_status = {
        "sprot": "Reviewed",
        "trembl": "Unreviewed",
    }.get(source_population)
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
                review_status = fields[2].rstrip(";") if len(fields) > 2 else ""
                if not entry_name or not review_status:
                    raise ValueError(
                        f"Malformed UniProt DAT ID line at {path}:{line_number}"
                    )
                if expected_review_status and review_status != expected_review_status:
                    raise ValueError(
                        "UniProt DAT source-role mismatch at "
                        f"{path}:{line_number}: source population {source_population!r} "
                        f"requires ID status {expected_review_status!r}, observed "
                        f"{review_status or 'missing'!r}"
                    )
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
                if not entry_name:
                    raise ValueError(
                        f"UniProt DAT record ending at {path}:{line_number} has no ID line"
                    )
                if not accessions:
                    raise ValueError(f"UniProt DAT record ending at {path}:{line_number} has no accession")
                if not sequence or SEQUENCE_RE.fullmatch(sequence) is None:
                    raise ValueError(f"Invalid or empty UniProt DAT sequence for {accessions[0]}")
                yield ProteinRecord(
                    protein_id=accessions[0], sequence=sequence, taxon_id=taxon,
                    aliases=tuple(accessions), source_header=entry_name,
                    source_population=source_population,
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


def iter_uniprot(path: Path, source_population: str = "unknown") -> Iterator[ProteinRecord]:
    suffixes = "".join(path.suffixes)
    if suffixes.endswith(".dat") or suffixes.endswith(".dat.gz"):
        yield from _iter_uniprot_dat(path, source_population)
    else:
        yield from _iter_uniprot_fasta(path, source_population)


def load_requested_proteins(path: Path, requested_accessions: set[str]) -> ProteinCatalog:
    """Compatibility helper for one explicitly synthetic/diagnostic source."""
    return load_requested_proteins_from_sources(
        {"fixture": path}, requested_accessions, strict_collisions=False
    )


def load_requested_proteins_from_sources(
    sources: dict[str, Path],
    requested_accessions: set[str],
    *,
    strict_collisions: bool,
    collision_database: Path | None = None,
) -> ProteinCatalog:
    """Stream selected UniProt populations in fixed order with a disk-backed accession audit."""
    started = time.monotonic()
    LOGGER.info(
        "UniProt sequence scan started: sources=%s requested_accessions=%d",
        sorted(sources), len(requested_accessions),
    )
    if not sources:
        raise ValueError("At least one selected UniProt source is required")
    catalog = ProteinCatalog()
    source_collision_counts: dict[str, Counter[str]] = defaultdict(Counter)
    seen_requested: dict[str, str] = {}
    processed = 0
    database = collision_database or Path(":memory:")
    if collision_database is not None:
        collision_database.parent.mkdir(parents=True, exist_ok=True)
        collision_database.unlink(missing_ok=True)
    connection = sqlite3.connect(str(database))
    try:
        connection.execute("PRAGMA journal_mode=OFF")
        connection.execute("PRAGMA synchronous=OFF")
        connection.execute(
            "CREATE TABLE identifiers (accession TEXT PRIMARY KEY, primary_accession TEXT NOT NULL, "
            "sequence_sha256 TEXT NOT NULL, source_population TEXT NOT NULL, kind TEXT NOT NULL)"
        )
        source_order = [name for name in ("sprot", "trembl", "fixture") if name in sources]
        source_order.extend(sorted(set(sources) - set(source_order)))
        for source_population in source_order:
            path = sources[source_population]
            format_name = "dat" if "".join(path.suffixes).endswith((".dat", ".dat.gz")) else "fasta-diagnostic"
            source_primary = 0
            source_secondary = 0
            source_retained = 0
            for record in iter_uniprot(path, source_population):
                processed += 1
                source_primary += 1
                source_secondary += max(0, len(set(record.aliases)) - 1)
                if processed % 1_000_000 == 0:
                    LOGGER.info(
                        "UniProt sequence progress: records=%d retained=%d elapsed_seconds=%.1f",
                        processed, len(catalog.records), time.monotonic() - started,
                    )
                aliases = set(record.aliases) | {record.protein_id}
                sequence_digest = hashlib.sha256(record.sequence.encode("ascii")).hexdigest()
                for accession in sorted(aliases):
                    kind = "primary" if accession == record.protein_id else "secondary"
                    inserted = connection.execute(
                        "INSERT OR IGNORE INTO identifiers VALUES (?, ?, ?, ?, ?)",
                        (accession, record.protein_id, sequence_digest, source_population, kind),
                    ).rowcount
                    if inserted:
                        continue
                    previous = connection.execute(
                        "SELECT primary_accession, sequence_sha256, source_population, kind "
                        "FROM identifiers WHERE accession=?", (accession,),
                    ).fetchone()
                    assert previous is not None
                    previous_primary, previous_digest, previous_source, previous_kind = previous
                    if (
                        previous_primary == record.protein_id
                        and previous_digest == sequence_digest
                        and previous_source == source_population
                        and previous_kind == kind
                    ):
                        collision_kind = "duplicate-record-identical"
                    elif previous_digest != sequence_digest:
                        collision_kind = "conflicting-sequence"
                    elif kind == "primary" or previous_kind == "primary":
                        collision_kind = "duplicate-primary-identical"
                    else:
                        collision_kind = "ambiguous-secondary-identical"
                    catalog.collision_counts[collision_kind] += 1
                    for involved_source in {previous_source, source_population}:
                        source_collision_counts[involved_source][collision_kind] += 1
                        if kind == "primary" and previous_kind == "primary":
                            source_collision_counts[involved_source][
                                "duplicate-primary-accession"
                            ] += 1
                    if collision_kind == "ambiguous-secondary-identical":
                        # A retired secondary accession can legitimately remain attached to
                        # multiple isolate-specific records with the same sequence. Choosing one
                        # primary would invent provenance, so keep the alias explicitly ambiguous
                        # and exclude it from supervision if GOA requests it.
                        catalog.ambiguous_aliases.add(accession)
                    elif strict_collisions:
                        raise ValueError(
                            f"UniProt accession collision ({collision_kind}) for {accession}: "
                            f"{previous_source}/{previous_primary} versus "
                            f"{source_population}/{record.protein_id}"
                        )

                matches = aliases & requested_accessions
                if not matches:
                    continue
                source_retained += 1
                existing = catalog.records.get(record.protein_id)
                if existing is not None and existing.sequence != record.sequence:
                    raise ValueError(f"Conflicting sequences for UniProt accession {record.protein_id}")
                catalog.records.setdefault(record.protein_id, record)
                catalog.primary_source.setdefault(record.protein_id, source_population)
                for alias in sorted(aliases):
                    previous = catalog.alias_to_primary.get(alias)
                    if previous is None or previous == record.protein_id:
                        catalog.alias_to_primary[alias] = record.protein_id
                        catalog.alias_source[alias] = source_population
                    else:
                        previous_record = catalog.records.get(previous)
                        if previous_record is not None and previous_record.sequence != record.sequence:
                            raise ValueError(
                                f"Conflicting sequences for duplicate UniProt alias {alias}: "
                                f"{previous} versus {record.protein_id}"
                            )
                        catalog.ambiguous_aliases.add(alias)
                        catalog.alias_to_primary.pop(alias, None)
                        catalog.alias_source.pop(alias, None)
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
                        catalog.alias_source.pop(requested, None)
                    else:
                        seen_requested[requested] = record.protein_id
            catalog.source_counts[source_population] = {
                "primary_accessions_read": source_primary,
                "secondary_aliases_read": source_secondary,
                "records_matching_qualifying_goa": source_retained,
                "retained_primary_accessions": sum(
                    source == source_population for source in catalog.primary_source.values()
                ),
                "sequence_format": format_name,
            }
        for source_population in source_order:
            collisions = source_collision_counts[source_population]
            catalog.source_counts[source_population].update({
                "conflicting_sequences": int(collisions["conflicting-sequence"]),
                "duplicate_primary_accessions": int(
                    collisions["duplicate-primary-accession"]
                ),
                "ambiguous_secondary_aliases": int(
                    collisions["ambiguous-secondary-identical"]
                ),
                "accession_collision_counts": {
                    key: int(value) for key, value in sorted(collisions.items())
                },
            })
        connection.commit()
    finally:
        connection.close()
    ambiguous_secondary_collisions = catalog.collision_counts[
        "ambiguous-secondary-identical"
    ]
    if ambiguous_secondary_collisions:
        LOGGER.warning(
            "Observed %d identical-sequence ambiguous secondary-accession collisions; "
            "qualifying ambiguous aliases are excluded rather than assigned to an arbitrary "
            "primary accession",
            ambiguous_secondary_collisions,
        )
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
