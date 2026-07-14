from __future__ import annotations

from collections import Counter
from contextlib import ExitStack
from dataclasses import asdict
import gzip
import json
import logging
from pathlib import Path
import re
import time

from .config import ASPECT_TO_NAMESPACE, ROOT_TERMS, SUPERVISOR_EVIDENCE_CODES
from .inputs import open_text
from .models import AnnotationRecord, ExcludedAnnotation, GoaLoadResult
from .ontology import Ontology


LOGGER = logging.getLogger(__name__)


GO_ID_RE = re.compile(r"^GO:\d{7}$")
DATE_RE = re.compile(r"^\d{8}$")
ISOFORM_ACCESSION_RE = re.compile(r"^.+-\d+$")
REQUIRED_GAF_FIELDS = {
    0: "DB",
    1: "DB_Object_ID",
    2: "DB_Object_Symbol",
    3: "Relation_or_Qualifier",
    4: "GO_ID",
    5: "DB_Reference",
    6: "Evidence_Code",
    8: "Aspect",
    11: "DB_Object_Type",
    12: "Taxon",
    13: "Date",
    14: "Assigned_By",
}


def _taxon(value: str) -> str:
    first = value.split("|", 1)[0]
    return first.split(":", 1)[1] if first.startswith("taxon:") else first


def _header_value(line: str) -> tuple[str, str] | None:
    text = line[1:].strip()
    if ":" not in text:
        return None
    key, value = text.split(":", 1)
    return key.strip().lower().replace("-", "_"), value.strip()


def load_goa(
    path: Path,
    ontology: Ontology,
    evidence_codes: frozenset[str] = SUPERVISOR_EVIDENCE_CODES,
    strict_malformed: bool = True,
    spool_dir: Path | None = None,
    excluded_sample_per_reason: int = 1000,
) -> GoaLoadResult:
    if excluded_sample_per_reason < 0:
        raise ValueError("excluded_sample_per_reason cannot be negative")
    started = time.monotonic()
    LOGGER.info("GOA scan started: %s", path)
    records: list[AnnotationRecord] = []
    excluded: list[ExcludedAnnotation] = []
    annotations: dict[str, set[str]] = {}
    qualifying_accessions: set[str] = set()
    counters = Counter()
    evidence_counts = Counter()
    taxonomy_counts = Counter()
    direct_go_counts = Counter()
    rejected_evidence_reason_counts = Counter()
    rejection_reason_counts = Counter()
    rejection_evidence_counts = Counter()
    rejection_database_counts = Counter()
    rejection_object_type_counts = Counter()
    rejection_aspect_counts = Counter()
    accepted_database_counts = Counter()
    accepted_object_type_counts = Counter()
    accepted_aspect_counts = Counter()
    annotation_decision_counts = Counter()
    excluded_sample_counts = Counter()
    candidate_accessions: set[str] = set()
    headers: dict[str, str] = {}
    record_spool = spool_dir / "qualifying_annotations.raw.jsonl.gz" if spool_dir else None
    excluded_spool = spool_dir / "excluded_annotations.sample.jsonl.gz" if spool_dir else None
    if spool_dir:
        spool_dir.mkdir(parents=True, exist_ok=True)

    with ExitStack() as stack:
        # These are private scratch spools, so favor throughput over archival compression ratio.
        record_handle = stack.enter_context(
            gzip.open(record_spool, "wt", encoding="utf-8", compresslevel=1)
        ) if record_spool else None
        excluded_handle = stack.enter_context(
            gzip.open(excluded_spool, "wt", encoding="utf-8", compresslevel=1)
        ) if excluded_spool else None

        def emit_excluded(record: ExcludedAnnotation) -> None:
            reason = record.rejection_reason
            evidence = record.evidence or "missing"
            database = record.database or "missing"
            object_type = record.object_type or "missing"
            aspect = record.aspect or "missing"
            rejected_evidence_reason_counts[(evidence, reason)] += 1
            rejection_reason_counts[reason] += 1
            rejection_evidence_counts[evidence] += 1
            rejection_database_counts[database] += 1
            rejection_object_type_counts[object_type] += 1
            rejection_aspect_counts[aspect] += 1
            annotation_decision_counts[
                ("rejected", reason, evidence, database, object_type, aspect)
            ] += 1
            if excluded_sample_counts[reason] >= excluded_sample_per_reason:
                return
            excluded_sample_counts[reason] += 1
            if excluded_handle:
                excluded_handle.write(json.dumps(asdict(record), separators=(",", ":")) + "\n")
            else:
                excluded.append(record)

        handle = stack.enter_context(open_text(path))
        for line_number, raw_line in enumerate(handle, start=1):
            if line_number % 1_000_000 == 0:
                LOGGER.info(
                    "GOA progress: lines=%d accepted=%d rejected=%d malformed=%d",
                    line_number,
                    counters["kept_rows"],
                    sum(value for key, value in counters.items() if key.startswith("rejected_")),
                    counters["malformed"],
                )
            if raw_line.startswith("!"):
                parsed = _header_value(raw_line)
                if parsed:
                    headers[parsed[0]] = parsed[1]
                continue
            if not raw_line.strip():
                continue
            columns = raw_line.rstrip("\n\r").split("\t")
            if len(columns) != 17:
                counters["malformed"] += 1
                excluded_record = ExcludedAnnotation(
                    raw_accession=columns[1] if len(columns) > 1 else "",
                    raw_go_id=columns[4] if len(columns) > 4 else "",
                    evidence=columns[6] if len(columns) > 6 else "",
                    qualifier=columns[3] if len(columns) > 3 else "",
                    aspect=columns[8] if len(columns) > 8 else "",
                    object_type=columns[11] if len(columns) > 11 else "",
                    taxon_id="", line_number=line_number,
                    rejection_reason="malformed",
                    detail=f"GAF 2.2 requires exactly 17 columns; observed {len(columns)}",
                    database=columns[0] if columns else "",
                    symbol=columns[2] if len(columns) > 2 else "",
                    reference=columns[5] if len(columns) > 5 else "",
                    with_from=columns[7] if len(columns) > 7 else "",
                    assigned_date=columns[13] if len(columns) > 13 else "",
                    assigned_by=columns[14] if len(columns) > 14 else "",
                    annotation_extension=columns[15] if len(columns) > 15 else "",
                    gene_product_form=columns[16] if len(columns) > 16 else "",
                )
                emit_excluded(excluded_record)
                continue
            counters["processed"] += 1
            accession = columns[1].strip()
            symbol = columns[2].strip()
            qualifier_tokens = tuple(item for item in columns[3].split("|") if item)
            raw_go_id = columns[4].strip()
            reference = columns[5].strip()
            evidence = columns[6].strip()
            with_from = columns[7].strip()
            aspect = columns[8].strip()
            object_type = columns[11].strip()
            taxon_id = _taxon(columns[12].strip())
            assigned_date = columns[13].strip()
            assigned_by = columns[14].strip()
            annotation_extension = columns[15].strip()
            gene_product_form = columns[16].strip()

            reason = ""
            detail = ""
            missing_required = [
                name for index, name in REQUIRED_GAF_FIELDS.items() if not columns[index].strip()
            ]
            if missing_required or not GO_ID_RE.fullmatch(raw_go_id) or not DATE_RE.fullmatch(assigned_date):
                reason = "malformed_required_fields"
                details = []
                if missing_required:
                    details.append("blank=" + "|".join(missing_required))
                if raw_go_id and not GO_ID_RE.fullmatch(raw_go_id):
                    details.append("invalid_GO_ID=" + raw_go_id)
                if assigned_date and not DATE_RE.fullmatch(assigned_date):
                    details.append("invalid_date=" + assigned_date)
                detail = ";".join(details)
                counters["malformed"] += 1
            elif columns[0] != "UniProtKB":
                reason = "database"
            elif object_type != "protein":
                reason = "object_type"
            elif not accession:
                reason = "missing_accession"
            elif ISOFORM_ACCESSION_RE.fullmatch(accession):
                reason = "isoform_accession"
                detail = "isoform suffix is preserved and excluded; no canonicalization policy is assumed"
            elif evidence not in evidence_codes:
                reason = "evidence_code"
            elif "NOT" in qualifier_tokens:
                reason = "not_qualifier"
            elif aspect not in ASPECT_TO_NAMESPACE:
                reason = "aspect"
            elif gene_product_form:
                reason = "isoform_specific"
                detail = gene_product_form

            term_description: dict[str, object] = {}
            canonical = None
            if not reason:
                candidate_accessions.add(accession)
                term_description = ontology.describe(raw_go_id)
                canonical = str(term_description["canonical_id"] or "")
                if not canonical:
                    reason = str(term_description["status"])
                    detail = "consider=" + "|".join(term_description["consider"])  # type: ignore[arg-type]
                elif ontology.namespace(canonical) != ASPECT_TO_NAMESPACE[aspect]:
                    reason = "namespace_mismatch"
                    detail = (
                        f"GAF aspect={aspect}; ontology namespace={ontology.namespace(canonical)}"
                    )

            if reason:
                counters[f"rejected_{reason}"] += 1
                excluded_record = ExcludedAnnotation(
                    raw_accession=accession,
                    raw_go_id=raw_go_id,
                    evidence=evidence,
                    qualifier="|".join(qualifier_tokens),
                    aspect=aspect,
                    object_type=object_type,
                    taxon_id=taxon_id,
                    line_number=line_number,
                    rejection_reason=reason,
                    detail=detail,
                    database=columns[0].strip(),
                    symbol=symbol,
                    reference=reference,
                    with_from=with_from,
                    assigned_date=assigned_date,
                    assigned_by=assigned_by,
                    annotation_extension=annotation_extension,
                    gene_product_form=gene_product_form,
                )
                emit_excluded(excluded_record)
                continue

            assert canonical is not None
            term_action = str(term_description["status"])
            if canonical in ROOT_TERMS:
                counters["accepted_direct_root"] += 1
            record = AnnotationRecord(
                database=columns[0].strip(),
                raw_accession=accession,
                protein_id=accession,
                symbol=symbol,
                raw_go_id=raw_go_id,
                go_id=canonical,
                namespace=ontology.namespace(canonical),
                aspect=aspect,
                evidence=evidence,
                qualifier="|".join(qualifier_tokens),
                reference=reference,
                with_from=with_from,
                taxon_id=taxon_id,
                assigned_date=assigned_date,
                assigned_by=assigned_by,
                annotation_extension=annotation_extension,
                gene_product_form=gene_product_form,
                line_number=line_number,
                term_action=term_action,
            )
            if record_handle:
                record_handle.write(json.dumps(asdict(record), separators=(",", ":")) + "\n")
            else:
                records.append(record)
            # In production-spool mode only accession membership is needed in memory; direct terms
            # remain in the disk spool and are consumed from there. Tests without a spool retain the
            # convenient complete mapping.
            qualifying_accessions.add(accession)
            if record_spool is None:
                annotations.setdefault(accession, set()).add(canonical)
            counters["kept_rows"] += 1
            evidence_counts[evidence] += 1
            accepted_database_counts[record.database] += 1
            accepted_object_type_counts[object_type] += 1
            accepted_aspect_counts[aspect] += 1
            annotation_decision_counts[
                ("accepted", "qualifying", evidence, record.database, object_type, aspect)
            ] += 1
            taxonomy_counts[taxon_id or "unknown"] += 1
            direct_go_counts[canonical] += 1

    if headers.get("gaf_version") != "2.2":
        counters["malformed_header"] += 1
    counters["qualifying_proteins"] = len(qualifying_accessions)
    counters["rejected_malformed"] = counters["malformed"]
    if strict_malformed and counters["malformed"]:
        raise ValueError(
            f"GOA contains {counters['malformed']} malformed data row(s); "
            "rerun a parser-only audit with strict_malformed=False to inspect exclusions"
        )
    if strict_malformed and counters["malformed_header"]:
        raise ValueError(
            "GOA does not declare the required !gaf-version: 2.2 header"
        )
    LOGGER.info(
        "GOA scan completed: lines=%d accepted_rows=%d qualifying_accessions=%d "
        "malformed=%d elapsed_seconds=%.1f",
        line_number if 'line_number' in locals() else 0,
        counters["kept_rows"],
        len(qualifying_accessions),
        counters["malformed"],
        time.monotonic() - started,
    )
    return GoaLoadResult(
        records=records,
        excluded=excluded,
        annotations=annotations,
        qualifying_accessions=qualifying_accessions,
        counters=counters,
        evidence_counts=evidence_counts,
        taxonomy_counts=taxonomy_counts,
        direct_go_counts=direct_go_counts,
        rejected_evidence_reason_counts=rejected_evidence_reason_counts,
        rejection_reason_counts=rejection_reason_counts,
        rejection_evidence_counts=rejection_evidence_counts,
        rejection_database_counts=rejection_database_counts,
        rejection_object_type_counts=rejection_object_type_counts,
        rejection_aspect_counts=rejection_aspect_counts,
        accepted_database_counts=accepted_database_counts,
        accepted_object_type_counts=accepted_object_type_counts,
        accepted_aspect_counts=accepted_aspect_counts,
        annotation_decision_counts=annotation_decision_counts,
        excluded_sample_counts=excluded_sample_counts,
        candidate_accessions=candidate_accessions,
        headers=headers,
        record_spool=record_spool,
        excluded_spool=excluded_spool,
    )


def iter_annotation_records(result: GoaLoadResult):
    if result.record_spool is None:
        yield from result.records
        return
    with gzip.open(result.record_spool, "rt", encoding="utf-8") as handle:
        for line in handle:
            values = json.loads(line)
            protein_id, action = result.accession_map.get(
                values["raw_accession"], (values["protein_id"], values["accession_action"])
            )
            values["protein_id"] = protein_id
            values["accession_action"] = action
            yield AnnotationRecord(**values)


def iter_excluded_annotations(result: GoaLoadResult):
    if result.excluded_spool is None:
        yield from result.excluded
        return
    with gzip.open(result.excluded_spool, "rt", encoding="utf-8") as handle:
        for line in handle:
            yield ExcludedAnnotation(**json.loads(line))
