from __future__ import annotations

from contextlib import contextmanager
from collections import defaultdict
import logging
import os
from pathlib import Path
import shutil
import subprocess
from typing import Iterator

from .config import CAFA3_FINAL_EXP_CODES
from .io_utils import open_text
from .models import GafRecord

LOGGER = logging.getLogger(__name__)


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


@contextmanager
def open_gaf_text(path: str | Path):
    """Open GAF text, optionally using pigz for faster gzip decompression."""
    path = Path(path)
    use_pigz = (
        path.suffix == ".gz"
        and os.environ.get("CAFA_BUILDER_USE_PIGZ", "1") != "0"
        and shutil.which("pigz") is not None
    )
    if not use_pigz:
        with open_text(path) as handle:
            yield handle
        return

    LOGGER.info("Streaming gzip with pigz: %s", path)
    proc = subprocess.Popen(["pigz", "-dc", str(path)], stdout=subprocess.PIPE, text=True)
    if proc.stdout is None:
        raise RuntimeError("pigz did not provide stdout")
    try:
        yield proc.stdout
    finally:
        proc.stdout.close()
        return_code = proc.wait()
        # Early smoke-test exits can close the pipe before pigz reaches EOF.
        if return_code not in (0, -13, 141):
            raise RuntimeError(f"pigz failed for {path} with exit code {return_code}")


def progress_interval_from_env() -> int:
    value = os.environ.get("CAFA_BUILDER_GOA_PROGRESS_INTERVAL", "1000000")
    try:
        return max(0, int(value))
    except ValueError:
        LOGGER.warning("Invalid CAFA_BUILDER_GOA_PROGRESS_INTERVAL=%r; using 1000000", value)
        return 1_000_000


def iter_gaf(path: str | Path, max_records: int | None = None) -> Iterator[GafRecord]:
    """Stream GAF 2.x records.

    This mirrors the CAFA_benchmark use of GAF rows but avoids loading the
    complete GOA file into memory.
    """
    emitted = 0
    with open_gaf_text(path) as handle:
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
    allowed_proteins: set[str] | frozenset[str] | None = None,
    progress_interval: int | None = None,
) -> dict[str, set[str]]:
    annots: dict[str, set[str]] = defaultdict(set)
    processed = 0
    kept = 0
    skipped_not_allowed = 0
    skipped_evidence = 0
    skipped_not = 0
    skipped_aspect = 0
    skipped_taxon = 0
    skipped_db = 0
    progress_every = progress_interval_from_env() if progress_interval is None else progress_interval

    if allowed_proteins is not None:
        LOGGER.info("Filtering GOA %s to %d loaded UniProt sequence IDs", path, len(allowed_proteins))

    with open_gaf_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("!"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 15:
                raise ValueError(f"Invalid GAF row at {path}:{line_number}: expected at least 15 columns")

            processed += 1

            if cols[0] != "UniProtKB":
                skipped_db += 1
            else:
                protein_id = cols[1]
                if allowed_proteins is not None and protein_id not in allowed_proteins:
                    skipped_not_allowed += 1
                elif cols[6] not in evidence_codes:
                    skipped_evidence += 1
                elif cols[3] and "NOT" in cols[3].split("|"):
                    skipped_not += 1
                elif cols[8] not in {"P", "C", "F"}:
                    skipped_aspect += 1
                elif target_taxa and parse_taxon(cols[12]) not in target_taxa:
                    skipped_taxon += 1
                else:
                    annots[protein_id].add(cols[4])
                    kept += 1

            if progress_every and processed % progress_every == 0:
                LOGGER.info(
                    "GOA progress for %s: processed=%d kept_rows=%d proteins=%d skipped_outside_sequences=%d",
                    path,
                    processed,
                    kept,
                    len(annots),
                    skipped_not_allowed,
                )
            if max_records is not None and processed >= max_records:
                break

    LOGGER.info(
        "Loaded GOA annotations from %s: processed=%d kept_rows=%d proteins=%d "
        "skipped_db=%d skipped_outside_sequences=%d skipped_evidence=%d skipped_not=%d "
        "skipped_aspect=%d skipped_taxon=%d",
        path,
        processed,
        kept,
        len(annots),
        skipped_db,
        skipped_not_allowed,
        skipped_evidence,
        skipped_not,
        skipped_aspect,
        skipped_taxon,
    )
    return dict(annots)
