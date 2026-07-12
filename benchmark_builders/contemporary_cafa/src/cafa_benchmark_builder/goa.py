from __future__ import annotations

from contextlib import contextmanager
from collections import Counter, defaultdict
import logging
import os
from pathlib import Path
import shutil
import subprocess
from typing import Iterator

from .config import CAFA3_FINAL_EXP_CODES
from .io_utils import open_text
from .models import AnnotationLoadResult, GafRecord
from .ontology import Ontology

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


def load_normalized_annotation_map(
    path: str | Path,
    *,
    alias_to_primary: dict[str, str],
    source_ontology: Ontology,
    benchmark_ontology: Ontology,
    other_ontology: Ontology | None = None,
    snapshot: str = "",
    allow_frozen_source_fallback: bool = True,
    evidence_codes: frozenset[str] = CAFA3_FINAL_EXP_CODES,
    target_taxa: frozenset[str] = frozenset(),
    exclude_on_or_before: str | None = None,
    include_on_or_before: str | None = None,
    require_valid_dates: bool = True,
    max_records: int | None = None,
    progress_interval: int | None = None,
) -> AnnotationLoadResult:
    """Stream, filter and normalise a GOA snapshot.

    UniProt secondary accessions are collapsed onto the snapshot's primary
    accession. GO IDs are first resolved in the source snapshot and then mapped
    into the frozen benchmark ontology. Terms introduced after the frozen
    ontology are counted separately from malformed/unresolvable GO IDs.
    """
    annots: dict[str, set[str]] = defaultdict(set)
    counters = Counter()
    evidence_counts = Counter()
    taxon_counts = Counter()
    unmapped_terms = Counter()
    out_of_benchmark_terms = Counter()
    source_diagnostics: list[dict[str, object]] = []
    outside_frozen_diagnostics: list[dict[str, object]] = []
    date_filter_counts: dict[tuple[str, str], Counter] = defaultdict(Counter)
    progress_every = progress_interval_from_env() if progress_interval is None else progress_interval

    with open_gaf_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("!"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 15:
                raise ValueError(f"Invalid GAF row at {path}:{line_number}: expected at least 15 columns")
            counters["processed"] += 1

            if cols[0] != "UniProtKB":
                counters["skipped_db"] += 1
            else:
                protein_id = alias_to_primary.get(cols[1])
                taxon_id = parse_taxon(cols[12])
                if protein_id is None:
                    counters["skipped_outside_sequences"] += 1
                elif cols[6] not in evidence_codes:
                    counters["skipped_evidence"] += 1
                elif cols[3] and "NOT" in cols[3].split("|"):
                    counters["skipped_not"] += 1
                elif cols[8] not in {"P", "C", "F"}:
                    counters["skipped_aspect"] += 1
                elif target_taxa and taxon_id not in target_taxa:
                    counters["skipped_taxon"] += 1
                elif (exclude_on_or_before is not None or include_on_or_before is not None) and (
                    not cols[13].isdigit() or len(cols[13]) != 8
                ):
                    counters["invalid_date"] += 1
                    date_filter_counts[(protein_id, cols[8])]["invalid_date"] += 1
                elif exclude_on_or_before is not None and cols[13] <= exclude_on_or_before:
                    counters["skipped_backfill"] += 1
                    date_filter_counts[(protein_id, cols[8])]["skipped_backfill"] += 1
                elif include_on_or_before is not None and cols[13] > include_on_or_before:
                    counters["skipped_after_cutoff"] += 1
                    date_filter_counts[(protein_id, cols[8])]["skipped_after_cutoff"] += 1
                else:
                    date_filter_counts[(protein_id, cols[8])]["passed_date_filter"] += 1
                    source_term = source_ontology.resolve_term(cols[4])
                    if source_term is None:
                        frozen_term = benchmark_ontology.resolve_term(cols[4])
                        source_description = source_ontology.describe_id(cols[4])
                        frozen_description = benchmark_ontology.describe_id(cols[4])
                        other_description = (
                            other_ontology.describe_id(cols[4]) if other_ontology else {}
                        )
                        use_frozen = bool(allow_frozen_source_fallback and frozen_term)
                        classification = (
                            "source_snapshot_mismatch_valid_in_frozen"
                            if use_frozen else
                            "obsolete_without_unique_replacement"
                            if source_description["is_obsolete"] else
                            "absent_from_source_ontology"
                        )
                        source_diagnostics.append({
                            "snapshot": snapshot,
                            "raw_protein_identifier": cols[1],
                            "canonical_protein_identifier": protein_id,
                            "raw_go_id": cols[4],
                            "evidence_code": cols[6],
                            "qualifier": cols[3],
                            "assigned_annotation_date": cols[13],
                            "source_database": cols[0],
                            "taxon": taxon_id,
                            "gaf_line_number": line_number,
                            "source_ontology_file": str(source_ontology.filename),
                            "source_ontology_date": source_ontology.data_version or "",
                            "frozen_benchmark_ontology_file": str(benchmark_ontology.filename),
                            "frozen_benchmark_ontology_date": benchmark_ontology.data_version or "",
                            "exists_in_other_ontology": int(bool(other_description.get("exists"))),
                            "other_ontology_canonical_id": other_description.get("canonical_id", ""),
                            "alt_id_mapping": source_description["canonical_id"] if source_description["is_alt_id"] else "",
                            "obsolete_status": int(bool(source_description["is_obsolete"])),
                            "replaced_by": "|".join(source_description["replaced_by"]),
                            "consider": "|".join(source_description["consider"]),
                            "exists_in_frozen_ontology": int(bool(frozen_description["exists"])),
                            "frozen_canonical_id": frozen_term or "",
                            "final_classification": classification,
                            "final_action": "use_frozen_term" if use_frozen else "fail_strict_qc",
                        })
                        if use_frozen:
                            source_term = frozen_term
                            counters["resolved_in_frozen_fallback"] += 1
                        else:
                            counters["unmapped_source_go"] += 1
                            unmapped_terms[cols[4]] += 1
                            source_term = None
                    if source_term is not None:
                        benchmark_term = benchmark_ontology.resolve_term(source_term)
                        if benchmark_term is None:
                            # If the source snapshot has replaced a t0 term,
                            # the original GAF ID may still be the resolvable
                            # identifier in the frozen benchmark graph.
                            benchmark_term = benchmark_ontology.resolve_term(cols[4])
                        if benchmark_term is None:
                            counters["outside_frozen_ontology"] += 1
                            out_of_benchmark_terms[source_term] += 1
                            source_description = source_ontology.describe_id(cols[4])
                            other_description = (
                                other_ontology.describe_id(cols[4]) if other_ontology else {}
                            )
                            outside_frozen_diagnostics.append({
                                "snapshot": snapshot,
                                "raw_protein_identifier": cols[1],
                                "canonical_protein_identifier": protein_id,
                                "raw_go_id": cols[4],
                                "source_canonical_go_id": source_term,
                                "evidence_code": cols[6],
                                "qualifier": cols[3],
                                "assigned_annotation_date": cols[13],
                                "source_database": cols[0],
                                "taxon": taxon_id,
                                "gaf_line_number": line_number,
                                "source_ontology_file": str(source_ontology.filename),
                                "source_ontology_date": source_ontology.data_version or "",
                                "frozen_benchmark_ontology_file": str(benchmark_ontology.filename),
                                "frozen_benchmark_ontology_date": benchmark_ontology.data_version or "",
                                "exists_in_other_ontology": int(bool(other_description.get("exists"))),
                                "alt_id_mapping": source_description["canonical_id"] if source_description["is_alt_id"] else "",
                                "obsolete_status": int(bool(source_description["is_obsolete"])),
                                "replaced_by": "|".join(source_description["replaced_by"]),
                                "consider": "|".join(source_description["consider"]),
                                "final_classification": "valid_source_term_outside_frozen_graph",
                                "final_action": "exclude_from_frozen_label_space",
                            })
                        else:
                            annots[protein_id].add(benchmark_term)
                            counters["kept_rows"] += 1
                            evidence_counts[cols[6]] += 1
                            taxon_counts[taxon_id] += 1

            if progress_every and counters["processed"] % progress_every == 0:
                LOGGER.info(
                    "GOA progress for %s: processed=%d kept_rows=%d proteins=%d "
                    "skipped_outside_sequences=%d skipped_backfill=%d skipped_after_cutoff=%d",
                    path,
                    counters["processed"],
                    counters["kept_rows"],
                    len(annots),
                    counters["skipped_outside_sequences"],
                    counters["skipped_backfill"],
                    counters["skipped_after_cutoff"],
                )
            if max_records is not None and counters["processed"] >= max_records:
                break

    if require_valid_dates and counters["invalid_date"]:
        raise ValueError(
            f"{path} contains {counters['invalid_date']} otherwise-eligible rows with invalid GAF dates; "
            "backfill filtering cannot be applied safely"
        )

    LOGGER.info(
        "Loaded normalized GOA annotations from %s: processed=%d kept_rows=%d proteins=%d "
        "skipped_outside_sequences=%d skipped_evidence=%d skipped_not=%d "
        "skipped_backfill=%d skipped_after_cutoff=%d outside_frozen_ontology=%d "
        "unmapped_source_go=%d",
        path,
        counters["processed"],
        counters["kept_rows"],
        len(annots),
        counters["skipped_outside_sequences"],
        counters["skipped_evidence"],
        counters["skipped_not"],
        counters["skipped_backfill"],
        counters["skipped_after_cutoff"],
        counters["outside_frozen_ontology"],
        counters["unmapped_source_go"],
    )
    return AnnotationLoadResult(
        annotations={protein_id: set(terms) for protein_id, terms in annots.items()},
        counters=counters,
        evidence_counts=evidence_counts,
        taxon_counts=taxon_counts,
        unmapped_terms=unmapped_terms,
        out_of_benchmark_terms=out_of_benchmark_terms,
        source_diagnostics=source_diagnostics,
        outside_frozen_diagnostics=outside_frozen_diagnostics,
        date_filter_counts={key: Counter(value) for key, value in date_filter_counts.items()},
    )
