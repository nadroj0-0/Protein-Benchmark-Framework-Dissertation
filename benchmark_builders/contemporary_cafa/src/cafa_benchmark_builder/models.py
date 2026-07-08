from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class ProteinRecord:
    protein_id: str
    sequence: str
    taxon_id: str | None = None
    reviewed: bool | None = None
    entry_name: str | None = None
    accessions: tuple[str, ...] = ()


@dataclass(frozen=True)
class GafRecord:
    db: str
    db_object_id: str
    db_object_symbol: str
    qualifier: tuple[str, ...]
    go_id: str
    db_reference: str
    evidence: str
    with_from: str
    aspect: str
    db_object_name: str
    synonym: str
    db_object_type: str
    taxon_id: str
    date: str
    assigned_by: str
    annotation_extension: str = ""
    gene_product_form_id: str = ""
