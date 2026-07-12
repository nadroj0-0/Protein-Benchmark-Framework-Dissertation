from __future__ import annotations

from collections import Counter, defaultdict
import csv
import hashlib
import importlib.metadata
import json
import logging
import platform
from pathlib import Path
import subprocess

import pandas as pd

from .builder import (
    build_deepgoplus_dataframe_from_fasta,
    build_training_dataframe,
    export_pfp_csvs,
    load_deepgoplus_annotation_file,
    make_terms_dataframe,
    propagate_annotations,
    split_train_valid,
)
from .config import BuildConfig, PREFIX_TO_NAMESPACE
from .goa import load_normalized_annotation_map
from .models import AnnotationLoadResult, IdentityMatch, ProteinCatalog
from .official_targets import load_official_target_catalog
from .ontology import Ontology
from .parsers import load_protein_catalog
from .training_annotations import load_released_training_annotations


LOGGER = logging.getLogger(__name__)
PROTEIN_BINDING = "GO:0005515"
SOURCE_DIAGNOSTIC_FIELDS = [
    "snapshot", "raw_protein_identifier", "canonical_protein_identifier", "raw_go_id",
    "evidence_code", "qualifier", "assigned_annotation_date", "source_database", "taxon",
    "gaf_line_number", "source_ontology_file", "source_ontology_date",
    "frozen_benchmark_ontology_file", "frozen_benchmark_ontology_date",
    "exists_in_other_ontology", "other_ontology_canonical_id", "alt_id_mapping",
    "obsolete_status", "replaced_by", "consider", "exists_in_frozen_ontology",
    "frozen_canonical_id", "final_classification", "final_action",
]
OUTSIDE_FROZEN_FIELDS = [
    "snapshot", "raw_protein_identifier", "canonical_protein_identifier", "raw_go_id",
    "source_canonical_go_id", "evidence_code", "qualifier", "assigned_annotation_date",
    "source_database", "taxon", "gaf_line_number", "source_ontology_file",
    "source_ontology_date", "frozen_benchmark_ontology_file",
    "frozen_benchmark_ontology_date", "exists_in_other_ontology", "alt_id_mapping",
    "obsolete_status", "replaced_by", "consider", "final_classification", "final_action",
]
OFFICIAL_TARGET_FIELDS = [
    "snapshot", "target_id", "source_identifiers", "mapping_files", "taxon_id",
    "sequence_length", "special_or_custom_source", "status", "reason", "mapping_method",
    "uniprot_accession", "uniprot_entry_name", "reviewed", "source_candidate_count",
    "exact_sequence_candidate_count", "present_in_snapshot",
]


def _catalog_from_records(catalog: ProteinCatalog, protein_ids: set[str]) -> ProteinCatalog:
    subset = ProteinCatalog()
    for protein_id in sorted(protein_ids):
        record = catalog.records[protein_id]
        subset.records[protein_id] = record
        for alias in set(record.accessions) | {protein_id}:
            if alias in catalog.ambiguous_aliases:
                subset.ambiguous_aliases.add(alias)
            elif catalog.alias_to_primary.get(alias) == protein_id:
                subset.alias_to_primary[alias] = protein_id
    return subset


def _merge_catalogs(*catalogs: ProteinCatalog) -> ProteinCatalog:
    merged = ProteinCatalog()
    for catalog in catalogs:
        for protein_id, record in catalog.records.items():
            existing = merged.records.get(protein_id)
            if existing is not None and existing.sequence != record.sequence:
                raise ValueError(f"Conflicting sequences for UniProt accession {protein_id}")
            merged.records.setdefault(protein_id, record)
        for alias in catalog.ambiguous_aliases:
            merged.ambiguous_aliases.add(alias)
            merged.alias_to_primary.pop(alias, None)
        for alias, protein_id in catalog.alias_to_primary.items():
            if alias in merged.ambiguous_aliases:
                continue
            previous = merged.alias_to_primary.get(alias)
            if previous is None or previous == protein_id:
                merged.alias_to_primary[alias] = protein_id
            else:
                merged.ambiguous_aliases.add(alias)
                del merged.alias_to_primary[alias]
    return merged


def _target_catalog(
    paths: tuple[Path, ...],
    training_catalog: ProteinCatalog,
    training_taxa: frozenset[str],
    target_taxa: frozenset[str],
    training_reviewed_only: bool,
    target_reviewed_only: bool,
) -> ProteinCatalog:
    matching = {
        protein_id for protein_id, record in training_catalog.records.items()
        if not target_taxa or record.taxon_id in target_taxa
    }
    # Reuse the training pass only when it is guaranteed to contain the full
    # target universe. Otherwise perform a separate streaming pass.
    training_covers_targets = (
        (not training_taxa or target_taxa <= training_taxa)
        and (not training_reviewed_only or target_reviewed_only)
    )
    if training_covers_targets:
        return _catalog_from_records(training_catalog, matching)
    return load_protein_catalog(paths, target_taxa, target_reviewed_only)


def _build_identity_crosswalk(
    t0_catalog: ProteinCatalog,
    t1_catalog: ProteinCatalog,
    sequence_change_policy: str,
) -> tuple[list[IdentityMatch], dict[str, str]]:
    preliminary: dict[str, tuple[str | None, str, bool]] = {}
    reverse: dict[str, list[str]] = defaultdict(list)

    for t0_id in sorted(t0_catalog.records):
        record = t0_catalog.records[t0_id]
        candidates = {
            t1_catalog.alias_to_primary[alias]
            for alias in set(record.accessions) | {t0_id}
            if alias in t1_catalog.alias_to_primary
        }
        if not candidates:
            preliminary[t0_id] = (None, "missing_at_t1", False)
        elif len(candidates) > 1:
            preliminary[t0_id] = (None, "ambiguous_t1_identity", False)
        else:
            t1_id = next(iter(candidates))
            changed = record.sequence != t1_catalog.records[t1_id].sequence
            preliminary[t0_id] = (t1_id, "matched", changed)
            reverse[t1_id].append(t0_id)

    matches: list[IdentityMatch] = []
    t1_to_t0: dict[str, str] = {}
    for t0_id in sorted(preliminary):
        t1_id, reason, changed = preliminary[t0_id]
        if t1_id is None:
            matches.append(IdentityMatch(t0_id, None, "excluded", reason, changed))
            continue
        if len(reverse[t1_id]) > 1:
            matches.append(IdentityMatch(t0_id, t1_id, "excluded", "many_t0_ids_map_to_one_t1", changed))
            continue
        if changed and sequence_change_policy == "error":
            raise ValueError(f"Sequence changed between t0 and t1 for {t0_id} -> {t1_id}")
        if changed and sequence_change_policy == "exclude":
            matches.append(IdentityMatch(t0_id, t1_id, "excluded", "sequence_changed", True))
            continue
        matches.append(IdentityMatch(t0_id, t1_id, "matched", "matched", changed))
        t1_to_t0[t1_id] = t0_id

    mapped_t1 = set(t1_to_t0)
    for t1_id in sorted(set(t1_catalog.records) - mapped_t1):
        if not any(match.t1_id == t1_id for match in matches):
            matches.append(IdentityMatch("", t1_id, "excluded", "not_present_at_t0", False))
    return matches, t1_to_t0


def _drop_protein_binding_only(terms: set[str], go: Ontology, policy: str) -> set[str]:
    if policy == "keep":
        return set(terms)
    if policy != "drop-mf-protein-binding-only":
        raise ValueError(f"Unknown protein-binding policy: {policy}")
    mf_terms = {term for term in terms if go.get_namespace(term) == "molecular_function"}
    if mf_terms == {PROTEIN_BINDING}:
        return set(terms) - mf_terms
    return set(terms)


def _test_eligible_terms(
    t0_terms: set[str],
    t1_terms: set[str],
    go: Ontology,
    policy: str,
) -> tuple[set[str], set[str]]:
    """Return newly assigned terms and t0-blocked namespaces for a test target."""
    if policy == "global-no-knowledge":
        blocked = {go.get_namespace(term) for term in t0_terms if go.has_term(term)}
        return (set() if t0_terms else set(t1_terms)), blocked
    if policy != "ontology-no-knowledge":
        raise ValueError(f"Unknown test eligibility policy: {policy}")

    blocked = {go.get_namespace(term) for term in t0_terms if go.has_term(term)}
    eligible = {
        term for term in t1_terms
        if go.has_term(term) and go.get_namespace(term) not in blocked
    }
    return eligible, blocked


def _build_test_dataframe(
    go: Ontology,
    t0_catalog: ProteinCatalog,
    t1_catalog: ProteinCatalog,
    matches: list[IdentityMatch],
    t1_to_t0: dict[str, str],
    t0_annots: dict[str, set[str]],
    t1_annots: dict[str, set[str]],
    protein_binding_policy: str,
    test_eligibility_policy: str,
) -> tuple[pd.DataFrame, list[dict[str, object]]]:
    translated_t1 = {
        t1_to_t0[t1_id]: set(terms)
        for t1_id, terms in t1_annots.items()
        if t1_id in t1_to_t0
    }
    rows = []
    flow: list[dict[str, object]] = []

    for match in sorted(matches, key=lambda item: (item.t0_id, item.t1_id or "")):
        t0_record = t0_catalog.records.get(match.t0_id)
        t1_record = t1_catalog.records.get(match.t1_id or "")
        entry = {
            "t0_id": match.t0_id,
            "t1_id": match.t1_id or "",
            "taxon_id": (t0_record or t1_record).taxon_id if (t0_record or t1_record) else "",
            "sequence_changed": int(match.sequence_changed),
            "t0_annotation_count": len(t0_annots.get(match.t0_id, set())),
            "t1_annotation_count": len(translated_t1.get(match.t0_id, set())),
            "blocked_ontologies": "",
            "eligible_ontologies": "",
            "status": match.status,
            "reason": match.reason,
        }
        if match.status != "matched":
            flow.append(entry)
            continue
        t0_terms = t0_annots.get(match.t0_id, set())
        if test_eligibility_policy == "global-no-knowledge" and t0_terms:
            entry.update(status="excluded", reason="qualifying_annotation_at_t0")
            flow.append(entry)
            continue

        direct_terms = translated_t1.get(match.t0_id, set())
        if not direct_terms:
            entry.update(status="excluded", reason="no_new_qualifying_annotation_at_t1")
            flow.append(entry)
            continue
        direct_terms, blocked = _test_eligible_terms(
            t0_terms, direct_terms, go, test_eligibility_policy
        )
        entry["blocked_ontologies"] = ",".join(sorted(blocked))
        if not direct_terms:
            entry.update(
                status="excluded",
                reason="qualifying_annotation_at_t0_in_gained_ontology",
            )
            flow.append(entry)
            continue
        direct_terms = _drop_protein_binding_only(direct_terms, go, protein_binding_policy)
        if not direct_terms:
            entry.update(status="excluded", reason="protein_binding_only")
            flow.append(entry)
            continue
        entry["eligible_ontologies"] = ",".join(sorted({
            go.get_namespace(term) for term in direct_terms if go.has_term(term)
        }))
        propagated = propagate_annotations(go, direct_terms)
        if not propagated:
            entry.update(status="excluded", reason="empty_annotation_after_propagation")
            flow.append(entry)
            continue

        rows.append({
            "proteins": match.t0_id,
            "sequences": t0_record.sequence,
            "annotations": tuple(sorted(propagated)),
        })
        entry.update(status="selected", reason="new_annotation_after_t0")
        flow.append(entry)

    frame = pd.DataFrame(rows, columns=["proteins", "sequences", "annotations"])
    if not frame.empty:
        frame = frame.sort_values("proteins", kind="stable").reset_index(drop=True)
    return frame, flow


def _build_released_test_dataframe(
    go: Ontology,
    fasta_paths: tuple[Path, ...],
    annotations_path: Path,
    target_catalog: ProteinCatalog,
) -> tuple[pd.DataFrame, list[dict[str, object]], AnnotationLoadResult]:
    """Recreate DeepGOPlus test_data.pkl from released CAFA ground truth."""
    direct_annotations = load_deepgoplus_annotation_file(annotations_path)
    frames = []
    for fasta_path in fasta_paths:
        frame, _ = build_deepgoplus_dataframe_from_fasta(
            go=go,
            sequences_file=fasta_path,
            annots=direct_annotations,
            count_terms=False,
        )
        frames.append(frame)
    test_df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(
        columns=["proteins", "sequences", "annotations"]
    )
    if test_df["proteins"].duplicated().any():
        duplicates = sorted(test_df.loc[test_df["proteins"].duplicated(), "proteins"].unique())
        raise ValueError(f"Duplicate released CAFA test IDs: {duplicates[:10]}")

    selected_ids = set(test_df["proteins"])
    unknown_ids = selected_ids - set(target_catalog.records)
    if unknown_ids:
        raise ValueError(f"Released test labels reference unknown CAFA targets: {sorted(unknown_ids)[:10]}")
    flow = []
    for protein_id in test_df["proteins"]:
        record = target_catalog.records[protein_id]
        flow.append({
            "t0_id": protein_id,
            "t1_id": protein_id,
            "taxon_id": record.taxon_id or "",
            "sequence_changed": 0,
            "t0_annotation_count": 0,
            "t1_annotation_count": len(direct_annotations[protein_id]),
            "blocked_ontologies": "",
            "eligible_ontologies": "released-groundtruth",
            "status": "selected",
            "reason": "released_official_groundtruth",
        })

    kept_rows = sum(len(direct_annotations[protein_id]) for protein_id in selected_ids)
    skipped_rows = sum(
        len(terms) for protein_id, terms in direct_annotations.items()
        if protein_id not in selected_ids
    )
    result = AnnotationLoadResult(
        annotations={protein_id: set(direct_annotations[protein_id]) for protein_id in selected_ids},
        counters=Counter({
            "processed": kept_rows + skipped_rows,
            "kept_rows": kept_rows,
            "proteins": len(selected_ids),
            "skipped_outside_sequences": skipped_rows,
            "released_groundtruth": kept_rows,
        }),
    )
    return test_df, flow, result


def _write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _diagnostic_markdown(t0_result: AnnotationLoadResult, t1_result: AnnotationLoadResult) -> str:
    rows = t0_result.source_diagnostics + t1_result.source_diagnostics
    grouped: dict[tuple[object, object], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[(row["snapshot"], row["raw_go_id"])].append(row)
    lines = [
        "## GO source-resolution diagnostics",
        "",
        f"- Source-unresolvable rows observed before policy handling: {len(rows)}",
        f"- Rows resolved through an exact raw-ID match in the frozen graph: "
        f"{sum(row['final_action'] == 'use_frozen_term' for row in rows)}",
        f"- Rows remaining as strict-QC failures: "
        f"{sum(row['final_action'] == 'fail_strict_qc' for row in rows)}",
        f"- Valid source rows excluded outside the frozen graph: "
        f"{len(t0_result.outside_frozen_diagnostics) + len(t1_result.outside_frozen_diagnostics)}",
        "",
    ]
    if grouped:
        lines.extend([
            "|snapshot|GO ID|annotations|proteins|classification|action|",
            "|---|---|---:|---:|---|---|",
        ])
        for (snapshot, go_id), group in sorted(grouped.items()):
            proteins = {row["canonical_protein_identifier"] for row in group}
            classifications = "<br>".join(sorted({str(row["final_classification"]) for row in group}))
            actions = "<br>".join(sorted({str(row["final_action"]) for row in group}))
            lines.append(
                f"|{snapshot}|{go_id}|{len(group)}|{len(proteins)}|{classifications}|{actions}|"
            )
        lines.append("")
    lines.append(
        "No `consider` term is selected automatically; only primary/alt IDs, a unique "
        "`replaced_by`, or an exact raw-ID match in the frozen benchmark graph can resolve a row."
    )
    return "\n".join(lines) + "\n"


def _official_target_markdown(rows: list[dict[str, object]]) -> str:
    if not rows:
        return ""
    status = Counter(str(row["status"]) for row in rows)
    reasons = Counter(str(row["reason"]) for row in rows if row["status"] != "mapped")
    reviewed = Counter(
        "unknown" if row["reviewed"] == "" else "reviewed" if row["reviewed"] == 1 else "unreviewed"
        for row in rows
    )
    special = sum(int(row["special_or_custom_source"]) for row in rows)
    absent = sum(not int(row["present_in_snapshot"]) for row in rows)
    lines = [
        "## Official CAFA3 target mapping",
        "",
        f"- Mapping rows (t0 and t1): {len(rows)}",
        f"- Mapped: {status['mapped']}",
        f"- Unmapped: {status['unmapped']}",
        f"- Ambiguous: {status['ambiguous']}",
        f"- Special/custom-source rows: {special}",
        f"- Targets absent from the available UniProt mapping snapshot: {absent}",
        f"- Reviewed mappings: {reviewed['reviewed']}",
        f"- Unreviewed mappings: {reviewed['unreviewed']}",
        f"- Mapping status unknown: {reviewed['unknown']}",
        "",
    ]
    if reasons:
        lines.extend(["|unmapped/ambiguous reason|rows|", "|---|---:|"])
        lines.extend(f"|{reason}|{count}|" for reason, count in sorted(reasons.items()))
        lines.append("")
    lines.append(
        "The released CAFA IDs and exact FASTA sequences remain authoritative even when "
        "a UniProt mapping is absent or ambiguous; no target is silently discarded."
    )
    return "\n".join(lines) + "\n"


def _write_preflight_diagnostics(
    config: BuildConfig,
    t0_result: AnnotationLoadResult,
    t1_result: AnnotationLoadResult,
    official_target_rows: list[dict[str, object]],
) -> dict[str, Path]:
    config.reports.mkdir(parents=True, exist_ok=True)
    unresolved = config.reports / "unresolved_source_go_annotations.tsv"
    outside = config.reports / "outside_frozen_go_annotations.tsv"
    target_report = config.reports / "official_target_mapping.tsv"
    _write_tsv(
        unresolved,
        SOURCE_DIAGNOSTIC_FIELDS,
        t0_result.source_diagnostics + t1_result.source_diagnostics,
    )
    _write_tsv(
        outside,
        OUTSIDE_FROZEN_FIELDS,
        t0_result.outside_frozen_diagnostics + t1_result.outside_frozen_diagnostics,
    )
    if official_target_rows:
        _write_tsv(target_report, OFFICIAL_TARGET_FIELDS, official_target_rows)
    report = config.reports / "benchmark_build_report.md"
    report.write_text(
        "# CAFA benchmark build report\n\n"
        "The annotation inputs were parsed. This preflight report is written before strict QC "
        "so diagnostics survive a failed build.\n\n"
        + _diagnostic_markdown(t0_result, t1_result)
        + _official_target_markdown(official_target_rows)
    )
    written = {
        "unresolved_source_go_annotations": unresolved,
        "outside_frozen_go_annotations": outside,
        "build_report": report,
    }
    if official_target_rows:
        written["official_target_mapping"] = target_report
    return written


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _write_checksum_file(path: Path, files: list[Path]) -> None:
    with path.open("w") as handle:
        for file_path in sorted(set(files), key=lambda item: str(item)):
            handle.write(f"{_sha256(file_path)}  {file_path}\n")


def _validate_csv_outputs(written: dict[str, Path], strict: bool) -> dict[str, dict[str, int]]:
    csv.field_size_limit(1_000_000_000)
    stats: dict[str, dict[str, int]] = {}
    split_data: dict[tuple[str, str], tuple[set[str], set[str]]] = {}
    for prefix in sorted(PREFIX_TO_NAMESPACE):
        for split in ("training", "validation", "test"):
            key = f"{prefix}-{split}"
            path = written[key]
            proteins: set[str] = set()
            sequences: set[str] = set()
            rows = 0
            with path.open(newline="") as handle:
                reader = csv.reader(handle)
                header = next(reader, None)
                if header is None or header[:2] != ["proteins", "sequences"] or len(header) < 3:
                    raise ValueError(f"Invalid PFP CSV schema: {path}")
                if len(header[2:]) != len(set(header[2:])):
                    raise ValueError(f"Duplicate GO columns in {path}")
                for row_number, row in enumerate(reader, start=2):
                    if len(row) != len(header):
                        raise ValueError(f"Malformed CSV row at {path}:{row_number}")
                    if row[0] in proteins:
                        raise ValueError(f"Duplicate protein ID {row[0]} in {path}")
                    if any(value not in {"0", "1"} for value in row[2:]):
                        raise ValueError(f"Non-binary GO label at {path}:{row_number}")
                    proteins.add(row[0])
                    sequences.add(row[1])
                    rows += 1
            if strict and rows == 0:
                raise ValueError(f"Strict QC rejected empty ontology/split output: {path}")
            stats[key] = {"rows": rows, "go_terms": len(header) - 2}
            split_data[(prefix, split)] = (proteins, sequences)

        train_ids, train_sequences = split_data[(prefix, "training")]
        valid_ids, valid_sequences = split_data[(prefix, "validation")]
        test_ids, test_sequences = split_data[(prefix, "test")]
        if train_ids & test_ids or valid_ids & test_ids:
            raise ValueError(f"Protein overlap across train/validation and test for {prefix}")
        if train_sequences & test_sequences or valid_sequences & test_sequences:
            raise ValueError(f"Exact sequence overlap across train/validation and test for {prefix}")
        if train_sequences & valid_sequences:
            raise ValueError(f"Exact sequence overlap across train and validation for {prefix}")
    return stats


def _git_commit() -> str | None:
    repo_root = Path(__file__).resolve().parents[4]
    result = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=repo_root, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, check=False,
    )
    return result.stdout.strip() or None if result.returncode == 0 else None


def _package_versions() -> dict[str, str]:
    versions = {}
    for package in ("numpy", "pandas"):
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = "not-installed"
    return versions


def _counter_dict(counter: Counter) -> dict[str, int]:
    return {str(key): int(counter[key]) for key in sorted(counter)}


def _write_reports(
    config: BuildConfig,
    go: Ontology,
    written: dict[str, Path],
    flow: list[dict[str, object]],
    train_df: pd.DataFrame,
    test_df: pd.DataFrame,
    terms_df: pd.DataFrame,
    training_catalog: ProteinCatalog,
    target_catalog: ProteinCatalog,
    t0_result: AnnotationLoadResult,
    t1_result: AnnotationLoadResult,
    csv_stats: dict[str, dict[str, int]],
    training_annotation_result: AnnotationLoadResult | None = None,
    official_target_rows: list[dict[str, object]] | None = None,
) -> dict[str, Path]:
    report_dir = config.reports
    report_dir.mkdir(parents=True, exist_ok=True)
    reports: dict[str, Path] = {}

    flow_path = report_dir / "protein_flow.tsv"
    _write_tsv(flow_path, [
        "t0_id", "t1_id", "taxon_id", "sequence_changed", "t0_annotation_count",
        "t1_annotation_count", "blocked_ontologies", "eligible_ontologies", "status", "reason",
    ], flow)
    reports["protein_flow_report"] = flow_path

    exclusions = Counter(row["reason"] for row in flow if row["status"] != "selected")
    exclusion_path = report_dir / "exclusion_reasons.tsv"
    _write_tsv(exclusion_path, ["reason", "count"], [
        {"reason": reason, "count": exclusions[reason]} for reason in sorted(exclusions)
    ])
    reports["exclusion_report"] = exclusion_path

    evidence_rows = []
    for timepoint, result in (("t0", t0_result), ("t1", t1_result)):
        for evidence in sorted(result.evidence_counts):
            evidence_rows.append({
                "timepoint": timepoint,
                "evidence_code": evidence,
                "kept_rows": result.evidence_counts[evidence],
            })
    evidence_path = report_dir / "evidence_summary.tsv"
    _write_tsv(evidence_path, ["timepoint", "evidence_code", "kept_rows"], evidence_rows)
    reports["evidence_report"] = evidence_path

    taxon_rows = []
    train_ids = set(train_df["proteins"].tolist())
    test_ids = set(test_df["proteins"].tolist())
    for stage, protein_ids, catalog in (
        ("training_before_split", train_ids, training_catalog),
        ("test_before_ontology_export", test_ids, target_catalog),
    ):
        counts = Counter(catalog.records[protein_id].taxon_id or "unknown" for protein_id in protein_ids)
        for taxon_id in sorted(counts):
            taxon_rows.append({"stage": stage, "taxon_id": taxon_id, "proteins": counts[taxon_id]})
    taxon_path = report_dir / "taxon_summary.tsv"
    _write_tsv(taxon_path, ["stage", "taxon_id", "proteins"], taxon_rows)
    reports["taxon_report"] = taxon_path

    selected_flow = [row for row in flow if row["status"] == "selected"]
    gain_rows = []
    for prefix, namespace in sorted(PREFIX_TO_NAMESPACE.items()):
        namespace_terms = {
            term for term in terms_df["terms"].tolist()
            if go.has_term(term) and go.get_namespace(term) == namespace
        }
        proteins = sum(
            1 for annotations in test_df["annotations"].tolist()
            if set(annotations) & namespace_terms
        )
        gain_rows.append({
            "ontology": prefix,
            "selected_test_proteins": proteins,
            "training_defined_terms": len(namespace_terms),
        })
    gain_path = report_dir / "annotation_gain_summary.tsv"
    _write_tsv(gain_path, ["ontology", "selected_test_proteins", "training_defined_terms"], gain_rows)
    reports["annotation_gain_report"] = gain_path

    statistics = {
        "profile": config.profile_name,
        "training_proteins_before_split": len(train_df),
        "selected_test_proteins_before_ontology_export": len(test_df),
        "selected_target_flow_rows": len(selected_flow),
        "training_defined_terms": len(terms_df),
        "csv_outputs": csv_stats,
        "t0_goa": _counter_dict(t0_result.counters),
        "t1_goa": _counter_dict(t1_result.counters),
        "excluded_target_proteins": _counter_dict(exclusions),
        "t0_unmapped_terms": _counter_dict(t0_result.unmapped_terms),
        "t1_unmapped_terms": _counter_dict(t1_result.unmapped_terms),
        "t0_terms_outside_frozen_ontology": _counter_dict(t0_result.out_of_benchmark_terms),
        "t1_terms_outside_frozen_t0_ontology": _counter_dict(t1_result.out_of_benchmark_terms),
        "training_annotation_source": (
            _counter_dict(training_annotation_result.counters)
            if training_annotation_result else "t0_goa"
        ),
        "test_annotation_source": (
            "released_official_groundtruth" if config.test_annotations_file else "temporal_goa"
        ),
        "official_target_mapping": _counter_dict(Counter(
            row["status"] for row in (official_target_rows or [])
        )),
    }
    statistics_path = report_dir / "benchmark_statistics.json"
    statistics_path.write_text(json.dumps(statistics, indent=2, sort_keys=True) + "\n")
    reports["statistics"] = statistics_path

    manifest = {
        "builder_version": "0.2.0",
        "git_commit": _git_commit(),
        "python": platform.python_version(),
        "packages": _package_versions(),
        "profile": config.profile_name,
        "target_universe_policy": config.target_universe_policy,
        "training_snapshot_id": config.training_snapshot_id,
        "training_snapshot_date": config.training_snapshot_date,
        "evidence_codes": sorted(config.evidence_codes),
        "training_taxa": sorted(config.training_taxa),
        "target_taxa": sorted(config.target_taxa),
        "t0_cutoff": config.t0_cutoff,
        "t1_cutoff": config.t1_cutoff,
        "test_eligibility_policy": config.test_eligibility_policy,
        "exclude_t1_backfill": config.exclude_t1_backfill,
        "require_t0_presence": config.require_t0_presence,
        "sequence_change_policy": config.sequence_change_policy,
        "protein_binding_policy": config.protein_binding_policy,
        "min_count": config.min_count,
        "split": config.split,
        "seed": config.seed,
        "training_reviewed_only": config.reviewed_only,
        "target_reviewed_only": config.target_reviewed_only,
        "allow_frozen_source_fallback": config.allow_frozen_source_fallback,
        "include_relationships": config.include_rels,
        "ontology_policy": (
            "Resolve each GAF against its source product, then map labels into the frozen "
            "t0 benchmark ontology; an exact raw-ID match in the frozen graph may resolve a "
            "nearest-source snapshot mismatch when enabled; never select consider terms; "
            "report and exclude valid terms outside the frozen graph."
        ),
        "inputs": {
            "uniprot_t0": [str(path) for path in config.uniprot_t0],
            "uniprot_t1": [str(path) for path in config.uniprot_t1],
            "target_uniprot_t0": [str(path) for path in config.target_uniprot_t0],
            "target_uniprot_t1": [str(path) for path in config.target_uniprot_t1],
            "official_target_fastas": [str(path) for path in config.official_target_fastas],
            "official_target_mapping_dir": (
                str(config.official_target_mapping_dir)
                if config.official_target_mapping_dir else None
            ),
            "training_annotations_file": (
                str(config.training_annotations_file)
                if config.training_annotations_file else None
            ),
            "test_annotations_file": (
                str(config.test_annotations_file)
                if config.test_annotations_file else None
            ),
            "goa_t0": str(config.goa_t0) if config.goa_t0 else None,
            "goa_t1": str(config.goa_t1) if config.goa_t1 else None,
            "benchmark_ontology": str(config.go_obo),
            "t0_ontology": str(config.ontology_t0),
            "t1_ontology": str(config.ontology_t1),
        },
    }
    manifest_path = report_dir / "build_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    reports["manifest"] = manifest_path

    if config.write_checksums:
        input_paths = (
            list(config.uniprot_t0) + list(config.uniprot_t1)
            + list(config.target_uniprot_t0) + list(config.target_uniprot_t1)
            + list(config.official_target_fastas) + [
            config.go_obo, config.ontology_t0, config.ontology_t1,
        ])
        input_paths.extend(path for path in (config.goa_t0, config.goa_t1) if path)
        if config.training_annotations_file:
            input_paths.append(config.training_annotations_file)
        if config.test_annotations_file:
            input_paths.append(config.test_annotations_file)
        if config.official_target_mapping_dir:
            input_paths.extend(sorted(config.official_target_mapping_dir.glob("*")))
        input_checksums = report_dir / "input_checksums.sha256"
        _write_checksum_file(input_checksums, input_paths)
        reports["input_checksums"] = input_checksums

    output_checksums = report_dir / "output_checksums.sha256"
    _write_checksum_file(output_checksums, list(written.values()))
    reports["output_checksums"] = output_checksums

    report_path = report_dir / "benchmark_build_report.md"
    report_path.write_text(
        "# CAFA benchmark build report\n\n"
        f"- Profile: `{config.profile_name}`\n"
        f"- Training proteins before split: {len(train_df)}\n"
        f"- Test proteins before ontology export: {len(test_df)}\n"
        f"- Training-defined GO terms: {len(terms_df)}\n"
        f"- Test annotation source: `{'released-official-groundtruth' if config.test_annotations_file else 'temporal-goa'}`\n"
        f"- t1 rows excluded as backfill: {t1_result.counters['skipped_backfill']}\n"
        f"- t1 rows after the endpoint cutoff: {t1_result.counters['skipped_after_cutoff']}\n"
        f"- Test eligibility policy: `{config.test_eligibility_policy}`\n"
        f"- Target universe policy: `{config.target_universe_policy}`\n"
        f"- t0 terms outside the frozen benchmark ontology: {t0_result.counters['outside_frozen_ontology']}\n"
        f"- t1 terms outside the frozen t0 ontology: {t1_result.counters['outside_frozen_ontology']}\n"
        "- All nine PFP CSVs passed schema, duplicate, binary-label and overlap checks.\n\n"
        + _diagnostic_markdown(t0_result, t1_result)
        + _official_target_markdown(official_target_rows or [])
    )
    reports["build_report"] = report_path
    return reports


def _validate_config(config: BuildConfig) -> None:
    if not 0 < config.split < 1:
        raise ValueError("split must be between 0 and 1")
    if config.min_count < 1:
        raise ValueError("min_count must be positive")
    if config.sequence_change_policy not in {"exclude", "use-t0", "error"}:
        raise ValueError("sequence_change_policy must be exclude, use-t0 or error")
    if config.test_eligibility_policy not in {
        "global-no-knowledge", "ontology-no-knowledge",
    }:
        raise ValueError(
            "test_eligibility_policy must be global-no-knowledge or ontology-no-knowledge"
        )
    if config.t0_cutoff and config.t1_cutoff and config.t1_cutoff <= config.t0_cutoff:
        raise ValueError("t1_cutoff must be later than t0_cutoff")
    if not config.require_t0_presence:
        raise ValueError("This builder requires t0 presence; disabling it would violate the benchmark contract")
    if config.target_universe_policy not in {
        "reconstructed-all-qualifying", "official-cafa3-targets",
    }:
        raise ValueError("Unknown target universe policy")
    if config.target_universe_policy == "official-cafa3-targets":
        if config.profile_name != "cafa3-reconstructed":
            raise ValueError("Official CAFA3 target mode is historical-validation only")
        if not config.official_target_fastas:
            raise ValueError("Official CAFA3 target mode requires FASTA files")
        if config.test_annotations_file is None and config.official_target_mapping_dir is None:
            raise ValueError("Raw-GOA official target mode requires a mapping directory")
    if config.test_annotations_file:
        if config.profile_name != "cafa3-reconstructed":
            raise ValueError("Released CAFA3 test ground truth is historical-validation only")
        if config.target_universe_policy != "official-cafa3-targets":
            raise ValueError("Released CAFA3 test ground truth requires the official target universe")
        if not config.test_annotations_file.is_file():
            raise FileNotFoundError(config.test_annotations_file)
    if config.training_annotations_file is None and config.goa_t0 is None:
        raise ValueError("A t0 GOA is required when released training annotations are not supplied")
    if config.test_annotations_file is None and (config.goa_t0 is None or config.goa_t1 is None):
        raise ValueError("Both GOA snapshots are required for temporal test construction")
    paths = (
        list(config.uniprot_t0) + list(config.uniprot_t1)
        + list(config.target_uniprot_t0) + list(config.target_uniprot_t1)
        + list(config.official_target_fastas)
        + [config.go_obo, config.ontology_t0, config.ontology_t1]
        + [path for path in (config.goa_t0, config.goa_t1) if path]
    )
    for path in paths:
        if not Path(path).is_file():
            raise FileNotFoundError(path)
    if config.training_annotations_file and not config.training_annotations_file.is_file():
        raise FileNotFoundError(config.training_annotations_file)
    if config.official_target_mapping_dir and not config.official_target_mapping_dir.is_dir():
        raise FileNotFoundError(config.official_target_mapping_dir)


def build_snapshot_benchmark(config: BuildConfig) -> dict[str, Path]:
    _validate_config(config)
    config.output_dir.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s:%(name)s:%(message)s")

    LOGGER.info("Loading frozen benchmark ontology from %s", config.go_obo)
    benchmark_go = Ontology(config.go_obo, with_rels=config.include_rels)
    t0_go = benchmark_go if config.ontology_t0 == config.go_obo else Ontology(
        config.ontology_t0, with_rels=config.include_rels
    )
    t1_go = benchmark_go if config.ontology_t1 == config.go_obo else Ontology(
        config.ontology_t1, with_rels=config.include_rels
    )
    LOGGER.info("Loading t0 training UniProt catalogue")
    training_paths = tuple(
        path for path in config.uniprot_t0
        if not (config.reviewed_only and "trembl" in path.name.lower())
    )
    if not training_paths:
        raise ValueError("No t0 UniProt input remains after applying the reviewed training policy")
    training_catalog = load_protein_catalog(
        training_paths, config.training_taxa, config.reviewed_only
    )
    LOGGER.info("Loaded %d canonical t0 training proteins", len(training_catalog.records))
    official_target_rows: list[dict[str, object]] = []
    if config.target_universe_policy == "official-cafa3-targets":
        released_groundtruth = config.test_annotations_file is not None
        t0_reference_paths = config.target_uniprot_t0 or config.uniprot_t0
        t1_reference_paths = config.target_uniprot_t1 or config.uniprot_t1
        t0_reference = (
            ProteinCatalog() if released_groundtruth
            else load_protein_catalog(t0_reference_paths, frozenset(), False)
        )
        t1_reference = (
            ProteinCatalog() if released_groundtruth
            else load_protein_catalog(t1_reference_paths, frozenset(), False)
        )
        t0_official = load_official_target_catalog(
            config.official_target_fastas,
            config.official_target_mapping_dir,
            t0_reference,
            config.target_taxa,
            "t0",
        )
        t1_official = load_official_target_catalog(
            config.official_target_fastas,
            config.official_target_mapping_dir,
            t1_reference,
            config.target_taxa,
            "t1",
        )
        t0_target_catalog = t0_official.catalog
        t1_target_catalog = t1_official.catalog
        if not released_groundtruth:
            official_target_rows = t0_official.rows + t1_official.rows
    else:
        reconstructed_t0_paths = config.target_uniprot_t0 or config.uniprot_t0
        reconstructed_t1_paths = config.target_uniprot_t1 or config.uniprot_t1
        t0_target_catalog = _target_catalog(
            reconstructed_t0_paths,
            training_catalog,
            config.training_taxa,
            config.target_taxa,
            config.reviewed_only,
            config.target_reviewed_only,
        )
        t1_target_catalog = load_protein_catalog(
            reconstructed_t1_paths, config.target_taxa, config.target_reviewed_only
        )
    LOGGER.info("Loaded %d canonical t0 target proteins", len(t0_target_catalog.records))
    LOGGER.info("Loaded %d canonical t1 target proteins", len(t1_target_catalog.records))

    annotation_catalog = (
        t0_target_catalog
        if config.training_annotations_file
        else _merge_catalogs(training_catalog, t0_target_catalog)
    )
    if config.goa_t0:
        LOGGER.info("Loading and normalising t0 GOA annotations from %s", config.goa_t0)
        t0_result = load_normalized_annotation_map(
            config.goa_t0,
            alias_to_primary=annotation_catalog.alias_to_primary,
            source_ontology=t0_go,
            benchmark_ontology=benchmark_go,
            other_ontology=t1_go,
            snapshot="t0",
            allow_frozen_source_fallback=config.allow_frozen_source_fallback,
            evidence_codes=config.evidence_codes,
            max_records=config.max_gaf_records,
        )
    else:
        LOGGER.info("Bypassing t0 GOA because released training and test annotations were supplied")
        t0_result = AnnotationLoadResult(annotations={})

    released_test_df: pd.DataFrame | None = None
    released_test_flow: list[dict[str, object]] | None = None
    if config.test_annotations_file:
        LOGGER.info("Loading released CAFA3 test ground truth from %s", config.test_annotations_file)
        released_test_df, released_test_flow, t1_result = _build_released_test_dataframe(
            benchmark_go,
            config.official_target_fastas,
            config.test_annotations_file,
            t0_target_catalog,
        )
    else:
        LOGGER.info("Loading and normalising t1 GOA annotations from %s", config.goa_t1)
        t1_result = load_normalized_annotation_map(
            config.goa_t1,
            alias_to_primary=t1_target_catalog.alias_to_primary,
            source_ontology=t1_go,
            benchmark_ontology=benchmark_go,
            other_ontology=t0_go,
            snapshot="t1",
            allow_frozen_source_fallback=config.allow_frozen_source_fallback,
            evidence_codes=config.evidence_codes,
            target_taxa=config.target_taxa,
            exclude_on_or_before=config.t0_cutoff if config.exclude_t1_backfill else None,
            include_on_or_before=config.t1_cutoff,
            max_records=config.max_gaf_records,
        )
    preflight_reports = _write_preflight_diagnostics(
        config, t0_result, t1_result, official_target_rows
    )
    if config.strict_qc and (t0_result.unmapped_terms or t1_result.unmapped_terms):
        raise ValueError("Strict QC found GO IDs that cannot be resolved in their source ontology")

    training_ids = set(training_catalog.records)
    target_t0_ids = set(t0_target_catalog.records)
    training_annotation_result: AnnotationLoadResult | None = None
    if config.training_annotations_file:
        training_annotation_result = load_released_training_annotations(
            config.training_annotations_file,
            training_catalog.alias_to_primary,
            benchmark_go,
        )
        if config.strict_qc and training_annotation_result.unmapped_terms:
            raise ValueError("Strict QC found released training GO IDs outside the benchmark ontology")
        train_annots = training_annotation_result.annotations
    else:
        train_annots = {
            protein_id: terms for protein_id, terms in t0_result.annotations.items()
            if protein_id in training_ids
        }
    target_t0_annots = {
        protein_id: terms for protein_id, terms in t0_result.annotations.items()
        if protein_id in target_t0_ids
    }

    LOGGER.info("Building deterministic DeepGOPlus-style training dataframe")
    train_all_df, counts = build_training_dataframe(
        benchmark_go, training_catalog.sequences, train_annots
    )
    terms_df = make_terms_dataframe(counts, config.min_count)
    train_df, valid_df = split_train_valid(train_all_df, split=config.split, seed=config.seed)
    LOGGER.info(
        "Training proteins=%d terms=%d split=%d/%d",
        len(train_all_df), len(terms_df), len(train_df), len(valid_df),
    )

    if released_test_df is not None and released_test_flow is not None:
        test_df, flow = released_test_df, released_test_flow
    else:
        matches, t1_to_t0 = _build_identity_crosswalk(
            t0_target_catalog, t1_target_catalog, config.sequence_change_policy
        )
        test_df, flow = _build_test_dataframe(
            benchmark_go,
            t0_target_catalog,
            t1_target_catalog,
            matches,
            t1_to_t0,
            target_t0_annots,
            t1_result.annotations,
            config.protein_binding_policy,
            config.test_eligibility_policy,
        )

    test_ids = set(test_df["proteins"].tolist())
    if not test_ids <= target_t0_ids:
        raise ValueError("Test contains a protein absent from the t0 target snapshot")
    if (
        config.test_annotations_file is None
        and
        config.test_eligibility_policy == "global-no-knowledge"
        and test_ids & set(target_t0_annots)
    ):
        raise ValueError("Test contains a protein with a qualifying t0 annotation")
    if (
        config.test_annotations_file is None
        and
        config.test_eligibility_policy == "global-no-knowledge"
        and test_ids & set(train_all_df["proteins"].tolist())
    ):
        raise ValueError("Protein ID overlap between training and test")
    LOGGER.info("Selected %d temporal test proteins", len(test_df))

    written: dict[str, Path] = {}
    if config.write_intermediates:
        intermediates = {
            "train_data": config.output_dir / "train_data.pkl",
            "train_data_train": config.output_dir / "train_data_train.pkl",
            "train_data_valid": config.output_dir / "train_data_valid.pkl",
            "test_data": config.output_dir / "test_data.pkl",
            "terms": config.output_dir / "terms.pkl",
        }
        train_all_df.to_pickle(intermediates["train_data"])
        train_df.to_pickle(intermediates["train_data_train"])
        valid_df.to_pickle(intermediates["train_data_valid"])
        test_df.to_pickle(intermediates["test_data"])
        terms_df.to_pickle(intermediates["terms"])
        written.update(intermediates)

    written.update(export_pfp_csvs(
        benchmark_go, train_df, valid_df, test_df, terms_df, config.output_dir
    ))
    csv_stats = _validate_csv_outputs(written, strict=config.strict_qc)
    reports = _write_reports(
        config,
        benchmark_go,
        written,
        flow,
        train_all_df,
        test_df,
        terms_df,
        training_catalog,
        t0_target_catalog,
        t0_result,
        t1_result,
        csv_stats,
        training_annotation_result,
        official_target_rows,
    )
    reports.update(preflight_reports)
    written.update(reports)
    return written
