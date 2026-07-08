from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import Iterator

from .config import CAFA3_FINAL_EXP_CODES
from .io_utils import open_text
from .models import GafRecord


GAF_FIELDS = (
    "DB",
    "DB_Object_ID",
    "DB_Object_Symbol",
    "Qualifier",
    "GO_ID",
    "DB:Reference",
    "Evidence",
    "With",
    "Aspect",
    "DB_Object_Name",
    "Synonym",
    "DB_Object_Type",
    "Taxon_ID",
    "Date",
    "Assigned_By",
    "Annotation_Extension",
    "Gene_Product_Form_ID",
)


def split_multi(value: str) -> tuple[str, ...]:
    if value:
        return tuple(x for x in value.split("|") if x)
    return ()


def parse_taxon(value: str) -> str:
    first = value.split("|", 1)[0]
    if first.startswith("taxon:"):
        return first.split(":", 1)[1]
    return first


def iter_gaf(path: str | Path, max_records: int | None = None) -> Iterator[GafRecord]:
    """Stream GAF 2.x records.

    This mirrors the CAFA_benchmark use of GAF rows but avoids loading the
    complete GOA file into memory.
    """
    emitted = 0
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("!"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 15:
                raise ValueError(f"Invalid GAF row at {path}:{line_number}: expected at least 15 columns")
            if len(cols) < 17:
                cols = cols + [""] * (17 - len(cols))
            yield GafRecord(
                db=cols[0],
                db_object_id=cols[1],
                db_object_symbol=cols[2],
                qualifier=split_multi(cols[3]),
                go_id=cols[4],
                db_reference=cols[5],
                evidence=cols[6],
                with_from=cols[7],
                aspect=cols[8],
                db_object_name=cols[9],
                synonym=cols[10],
                db_object_type=cols[11],
                taxon_id=parse_taxon(cols[12]),
                date=cols[13],
                assigned_by=cols[14],
                annotation_extension=cols[15],
                gene_product_form_id=cols[16],
            )
            emitted += 1
            if max_records is not None and emitted >= max_records:
                return


def keep_cafa_annotation(
    rec: GafRecord,
    evidence_codes: frozenset[str] = CAFA3_FINAL_EXP_CODES,
    target_taxa: frozenset[str] = frozenset(),
) -> bool:
    if rec.db != "UniProtKB":
        return False
    if rec.evidence not in evidence_codes:
        return False
    if "NOT" in rec.qualifier:
        return False
    if rec.aspect not in {"P", "C", "F"}:
        return False
    if target_taxa and rec.taxon_id not in target_taxa:
        return False
    return True


def load_annotation_map(
    path: str | Path,
    evidence_codes: frozenset[str] = CAFA3_FINAL_EXP_CODES,
    target_taxa: frozenset[str] = frozenset(),
    max_records: int | None = None,
) -> dict[str, set[str]]:
    annots: dict[str, set[str]] = defaultdict(set)
    for rec in iter_gaf(path, max_records=max_records):
        if keep_cafa_annotation(rec, evidence_codes=evidence_codes, target_taxa=target_taxa):
            annots[rec.db_object_id].add(rec.go_id)
    return dict(annots)
