"""Storage-safe reports for embedding inventory and reuse plans."""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import shutil
import tempfile
from collections import Counter, defaultdict
from contextlib import contextmanager
from dataclasses import asdict
from pathlib import Path
from typing import Any, Dict, Iterable, Iterator, List, Mapping, Optional, Sequence, Set, TextIO

from .models import InventoryRecord, InventoryResult, MODALITIES, ONTOLOGIES, SPLITS
from .provenance import provenance_markdown


class ReportError(ValueError):
    pass


INVENTORY_FIELDS = [
    "protein_id", "sequence_sha256", "modality", "source_directory", "source_file",
    "exists", "observed_shape", "expected_shape", "dtype", "finite", "valid",
    "scientifically_eligible", "source_protein_id", "match_route", "sequence_match",
    "provenance", "factual_status", "requested_action", "reason",
]
ACTION_FIELDS = [
    "protein_id", "sequence_sha256", "modality", "source_protein_id", "source_file",
    "match_route", "factual_status", "requested_action", "provenance", "reason",
]
def write_reports(
    result: InventoryResult,
    output_dir: Path,
    embedding_cache: Path,
    *,
    report_level: str = "compact",
    provenance: Optional[Mapping[str, Any]] = None,
    protected_roots: Sequence[Path] = (),
) -> Dict[str, Any]:
    if report_level not in {"compact", "full"}:
        raise ReportError("report_level must be compact or full")
    _prepare_output(output_dir, result, embedding_cache, protected_roots)
    staging = Path(
        tempfile.mkdtemp(prefix=".%s.staging-" % output_dir.name, dir=str(output_dir.parent))
    )
    published = False
    try:
        summary = _write_reports_staged(
            result,
            staging,
            embedding_cache,
            report_level=report_level,
            provenance=provenance,
        )
        _write_output_integrity(staging)
        staging.replace(output_dir)
        published = True
    finally:
        if not published:
            shutil.rmtree(staging, ignore_errors=True)
    return summary


def _write_reports_staged(
    result: InventoryResult,
    output_dir: Path,
    embedding_cache: Path,
    *,
    report_level: str,
    provenance: Optional[Mapping[str, Any]],
) -> Dict[str, Any]:
    records_by_protein: Dict[str, Dict[str, InventoryRecord]] = defaultdict(dict)
    for record in result.records:
        records_by_protein[record.protein_id][record.modality] = record

    _write_benchmark_proteins(result, output_dir / "benchmark_proteins.tsv.gz", False)
    _write_inventory(result.records, output_dir / "embedding_inventory.tsv.gz")
    _write_protein_summary(
        result, records_by_protein, output_dir / "protein_embedding_summary.tsv.gz"
    )
    _write_action_manifest(result.records, output_dir / "reuse_manifest.tsv.gz", "reuse")
    _write_action_manifest(
        result.records, output_dir / "generation_manifest.tsv", "generate"
    )
    _write_action_manifest(
        result.records, output_dir / "manual_review.tsv.gz", "manual-review"
    )
    extras = _write_cache_extras(result, output_dir / "cache_extras.tsv.gz")
    _write_modality_lists(result, output_dir)
    _write_exact_sequence_reuse(result, output_dir / "exact_sequence_reuse.tsv")
    _write_reason_counts(result, output_dir / "manual_review_reasons.tsv")
    _write_errors(result, output_dir / "errors.tsv")
    if report_level == "full":
        _write_benchmark_proteins(
            result, output_dir / "benchmark_proteins_full.tsv.gz", True
        )

    summary = _build_summary(result, records_by_protein, extras, embedding_cache, report_level)
    (output_dir / "embedding_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (output_dir / "embedding_summary.md").write_text(
        _summary_markdown(summary), encoding="utf-8"
    )
    provenance_payload: Mapping[str, Any] = provenance or {
        "schema_version": 1,
        "timestamp_utc": "not-recorded (direct library call)",
        "command": "not-recorded (direct library call)",
        "software": {
            "git_commit": "", "dirty_worktree": False, "python_version": "unknown",
            "numpy_version": "unknown", "package_version": "unknown",
        },
        "inputs": {
            "target_benchmark_fingerprint": result.benchmark.fingerprint,
            "source_benchmark_fingerprint": result.source_benchmark.fingerprint,
            "cache_catalog": {"fingerprint": "not-recorded"},
        },
        "artifact_verification": result.artifact_verification.as_dict(),
        "run": {"compatibility_policy": result.policy, "report_level": report_level},
    }
    (output_dir / "run_provenance.json").write_text(
        json.dumps(provenance_payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (output_dir / "run_provenance.md").write_text(
        provenance_markdown(provenance_payload), encoding="utf-8"
    )
    return summary


def _write_output_integrity(staging: Path) -> None:
    excluded = {
        "run_provenance.json", "run_provenance.md",
        "output_manifest.json", "RUN_COMPLETE.json",
    }
    payload_files = _file_manifest(staging, excluded)
    manifest = {
        "schema_version": 1,
        "catalog_schema": "relative-path-tab-size-tab-sha256-lf-v1",
        "files": payload_files,
        "file_count": len(payload_files),
        "total_bytes": sum(item["size_bytes"] for item in payload_files),
        "catalog_sha256": _manifest_digest(payload_files),
    }
    manifest_path = staging / "output_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    provenance_path = staging / "run_provenance.json"
    provenance_payload = json.loads(provenance_path.read_text(encoding="utf-8"))
    provenance_payload["outputs"] = {
        "payload_manifest": _path_identity(staging, manifest_path),
        "payload_file_count": manifest["file_count"],
        "payload_total_bytes": manifest["total_bytes"],
        "payload_catalog_sha256": manifest["catalog_sha256"],
    }
    provenance_path.write_text(
        json.dumps(provenance_payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    (staging / "run_provenance.md").write_text(
        provenance_markdown(provenance_payload), encoding="utf-8"
    )
    completion = {
        "schema_version": 1,
        "complete": True,
        "payload_manifest": _path_identity(staging, manifest_path),
        "run_provenance_json": _path_identity(staging, provenance_path),
        "run_provenance_markdown": _path_identity(staging, staging / "run_provenance.md"),
    }
    (staging / "RUN_COMPLETE.json").write_text(
        json.dumps(completion, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def _file_manifest(root: Path, excluded: Set[str]) -> List[Dict[str, Any]]:
    return [
        _path_identity(root, path)
        for path in sorted(root.iterdir(), key=lambda item: item.name)
        if path.is_file() and path.name not in excluded
    ]


def _path_identity(root: Path, path: Path) -> Dict[str, Any]:
    return {
        "path": path.relative_to(root).as_posix(),
        "size_bytes": path.stat().st_size,
        "sha256": _sha256_path(path),
    }


def _sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _manifest_digest(files: Sequence[Mapping[str, Any]]) -> str:
    digest = hashlib.sha256()
    for item in files:
        digest.update(
            ("%s\t%d\t%s\n" % (item["path"], item["size_bytes"], item["sha256"])).encode(
                "utf-8"
            )
        )
    return digest.hexdigest()


def _prepare_output(
    output_dir: Path,
    result: InventoryResult,
    embedding_cache: Path,
    protected_roots: Sequence[Path],
) -> None:
    resolved = output_dir.resolve()
    protected = [
        result.benchmark.directory.resolve(),
        result.source_benchmark.directory.resolve(),
        embedding_cache.resolve(),
        *(root.resolve() for root in protected_roots),
    ]
    for root in protected:
        if resolved == root or root in resolved.parents:
            raise ReportError("output directory cannot be inside a benchmark or embedding cache")
    if output_dir.exists():
        raise ReportError("output directory already exists; choose a new path: %s" % output_dir)
    output_dir.parent.mkdir(parents=True, exist_ok=True)


@contextmanager
def _text_writer(path: Path) -> Iterator[TextIO]:
    if path.suffix == ".gz":
        with path.open("wb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as zipped:
                with io.TextIOWrapper(zipped, encoding="utf-8", newline="") as text:
                    yield text
    else:
        with path.open("w", newline="", encoding="utf-8") as text:
            yield text


def _write_benchmark_proteins(
    result: InventoryResult, path: Path, include_sequence: bool
) -> None:
    fields = ["protein_id"]
    if include_sequence:
        fields.append("sequence")
    fields.extend(
        ["sequence_sha256", "sequence_length", "ontologies", "splits", "source_csv_files"]
    )
    with _text_writer(path) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for protein_id in sorted(result.benchmark.proteins):
            protein = result.benchmark.proteins[protein_id]
            row = {
                "protein_id": protein.protein_id,
                "sequence_sha256": protein.sequence_sha256,
                "sequence_length": protein.sequence_length,
                "ontologies": ";".join(sorted(protein.ontologies)),
                "splits": ";".join(sorted(protein.splits)),
                "source_csv_files": ";".join(sorted(protein.source_files)),
            }
            if include_sequence:
                row["sequence"] = protein.sequence
            writer.writerow(row)


def _write_inventory(records: Iterable[InventoryRecord], path: Path) -> None:
    with _text_writer(path) as handle:
        writer = csv.DictWriter(handle, fieldnames=INVENTORY_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for record in records:
            writer.writerow(record.as_dict())


def _write_action_manifest(
    records: Iterable[InventoryRecord], path: Path, action: str
) -> None:
    with _text_writer(path) as handle:
        writer = csv.DictWriter(handle, fieldnames=ACTION_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for record in records:
            if record.requested_action == action:
                row = record.as_dict()
                writer.writerow({field: row[field] for field in ACTION_FIELDS})


def _write_protein_summary(
    result: InventoryResult,
    records_by_protein: Mapping[str, Mapping[str, InventoryRecord]],
    path: Path,
) -> None:
    fields = [
        "protein_id", "sequence_sha256", "sequence_length", "has_any_embedding",
        "has_prott5", "has_text", "has_structure", "has_ppi", "valid_modalities",
        "reusable_modalities", "missing_modalities", "generation_modalities",
        "masked_modalities", "unavailable_modalities", "manual_review_modalities",
    ]
    with _text_writer(path) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for protein_id in sorted(records_by_protein):
            records = records_by_protein[protein_id]
            protein = result.benchmark.proteins[protein_id]
            present = [m for m in MODALITIES if records[m].exists]
            writer.writerow(
                {
                    "protein_id": protein_id,
                    "sequence_sha256": protein.sequence_sha256,
                    "sequence_length": protein.sequence_length,
                    "has_any_embedding": _bool(bool(present)),
                    "has_prott5": _bool("prott5" in present),
                    "has_text": _bool("text" in present),
                    "has_structure": _bool("structure" in present),
                    "has_ppi": _bool("ppi" in present),
                    "valid_modalities": ";".join(m for m in MODALITIES if records[m].valid),
                    "reusable_modalities": ";".join(_modalities_with_action(records, "reuse")),
                    "missing_modalities": ";".join(m for m in MODALITIES if records[m].factual_status == "missing"),
                    "generation_modalities": ";".join(_modalities_with_action(records, "generate")),
                    "masked_modalities": ";".join(_modalities_with_action(records, "leave-masked")),
                    "unavailable_modalities": ";".join(_modalities_with_action(records, "unavailable")),
                    "manual_review_modalities": ";".join(_modalities_with_action(records, "manual-review")),
                }
            )


def _write_cache_extras(result: InventoryResult, path: Path) -> Dict[str, int]:
    target_ids = set(result.benchmark.proteins)
    counts: Dict[str, int] = {}
    fields = ["modality", "protein_id", "source_file", "used_as_reuse_source"]
    with _text_writer(path) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for modality in MODALITIES:
            extra_ids = sorted(result.cache_ids[modality] - target_ids)
            counts[modality] = len(extra_ids)
            for protein_id in extra_ids:
                writer.writerow(
                    {
                        "modality": modality,
                        "protein_id": protein_id,
                        "source_file": str(Path(result.config.modalities[modality].directory) / (protein_id + ".npy")),
                        "used_as_reuse_source": _bool(protein_id in result.used_source_ids[modality]),
                    }
                )
    return counts


def _write_modality_lists(result: InventoryResult, output_dir: Path) -> None:
    by_modality: Dict[str, List[InventoryRecord]] = defaultdict(list)
    for record in result.records:
        by_modality[record.modality].append(record)
    for modality in MODALITIES:
        _write_id_list(output_dir / ("reuse_%s.txt" % modality), by_modality[modality], "reuse")
        _write_status_list(output_dir / ("missing_%s.txt" % modality), by_modality[modality], "missing")
        _write_id_list(output_dir / ("masked_%s.txt" % modality), by_modality[modality], "leave-masked")
        _write_id_list(output_dir / ("%s_manual_review.txt" % modality), by_modality[modality], "manual-review")
    _write_fasta(output_dir / "generate_prott5.fasta", by_modality["prott5"], result)
    _write_id_list(output_dir / "generate_text.txt", by_modality["text"], "generate")
    _write_id_list(output_dir / "generate_structure.txt", by_modality["structure"], "generate")
    _write_id_list(output_dir / "structure_unavailable.txt", by_modality["structure"], "unavailable")
    _write_id_list(output_dir / "extract_ppi.txt", by_modality["ppi"], "generate")
    _write_id_list(output_dir / "ppi_unavailable.txt", by_modality["ppi"], "unavailable")


def _write_id_list(path: Path, records: Sequence[InventoryRecord], action: str) -> None:
    ids = sorted(record.protein_id for record in records if record.requested_action == action)
    path.write_text("".join(protein_id + "\n" for protein_id in ids), encoding="utf-8")


def _write_status_list(path: Path, records: Sequence[InventoryRecord], status: str) -> None:
    ids = sorted(record.protein_id for record in records if record.factual_status == status)
    path.write_text("".join(protein_id + "\n" for protein_id in ids), encoding="utf-8")


def _write_fasta(path: Path, records: Sequence[InventoryRecord], result: InventoryResult) -> None:
    ids = sorted(record.protein_id for record in records if record.requested_action == "generate")
    with path.open("w", encoding="utf-8") as handle:
        for protein_id in ids:
            sequence = result.benchmark.proteins[protein_id].sequence
            handle.write(">%s\n" % protein_id)
            for start in range(0, len(sequence), 60):
                handle.write(sequence[start : start + 60] + "\n")


def _write_exact_sequence_reuse(result: InventoryResult, path: Path) -> None:
    fields = ["protein_id", "sequence_sha256", "modality", "source_protein_id", "match_route"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for record in result.records:
            if record.requested_action == "reuse" and record.match_route == "sequence-sha256":
                writer.writerow({field: record.as_dict()[field] for field in fields})


def _write_reason_counts(result: InventoryResult, path: Path) -> None:
    counts = Counter(
        (record.modality, record.factual_status, record.reason)
        for record in result.records if record.requested_action == "manual-review"
    )
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["modality", "factual_status", "reason", "count"])
        for (modality, status, reason), count in sorted(counts.items()):
            writer.writerow([modality, status, reason, count])


def _write_errors(result: InventoryResult, path: Path) -> None:
    error_statuses = {"unreadable", "wrong-dimension", "non-finite", "unsupported-dtype", "sequence-mismatch"}
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=ACTION_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for record in result.records:
            if record.factual_status in error_statuses:
                row = record.as_dict()
                writer.writerow({field: row[field] for field in ACTION_FIELDS})


def _build_summary(
    result: InventoryResult,
    records_by_protein: Mapping[str, Mapping[str, InventoryRecord]],
    extras: Dict[str, int],
    embedding_cache: Path,
    report_level: str,
) -> Dict[str, Any]:
    all_ids = set(result.benchmark.proteins)
    by_split = {split: {pid for pid, p in result.benchmark.proteins.items() if split in p.splits} for split in SPLITS}
    by_ontology = {ontology: {pid for pid, p in result.benchmark.proteins.items() if ontology in p.ontologies} for ontology in ONTOLOGIES}
    by_ontology_split = {
        "%s/%s" % (ontology, split): set(result.benchmark.file_members[(ontology, split)])
        for ontology in ONTOLOGIES for split in SPLITS
    }
    return {
        "schema_version": 2,
        "config_name": result.config.name,
        "policy": result.policy,
        "report_level": report_level,
        "benchmark_directory": str(result.benchmark.directory),
        "source_benchmark_directory": str(result.source_benchmark.directory),
        "embedding_cache": str(embedding_cache.resolve()),
        "target_benchmark_fingerprint": result.benchmark.fingerprint,
        "source_benchmark_fingerprint": result.source_benchmark.fingerprint,
        "artifact_verification": result.artifact_verification.as_dict(),
        "population": len(all_ids),
        "unique_sequences": len({p.sequence_sha256 for p in result.benchmark.proteins.values()}),
        "duplicate_identical_rows": result.benchmark.duplicate_rows,
        "pfp_missing_behavior": "Absent or unreadable arrays become zero vectors with mask 0.0 in PFP mmfp.dataset.MultiModalDataset._load_embedding_raw().",
        "configuration": {
            "target_benchmark_contract": asdict(result.config.target_benchmark_contract),
            "source_benchmark_contract": asdict(result.config.source_benchmark_contract),
            "modalities": {m: asdict(result.config.modalities[m]) for m in MODALITIES},
        },
        "coverage": {
            "global": _group_coverage(all_ids, records_by_protein),
            "by_split": {name: _group_coverage(ids, records_by_protein) for name, ids in by_split.items()},
            "by_ontology": {name: _group_coverage(ids, records_by_protein) for name, ids in by_ontology.items()},
            "by_ontology_split": {name: _group_coverage(ids, records_by_protein) for name, ids in by_ontology_split.items()},
        },
        "cache_files": {m: len(result.cache_ids[m]) for m in MODALITIES},
        "cache_extras": extras,
    }


def _group_coverage(
    protein_ids: Set[str], records_by_protein: Mapping[str, Mapping[str, InventoryRecord]]
) -> Dict[str, Any]:
    total = len(protein_ids)
    by_modality: Dict[str, Any] = {}
    for modality in MODALITIES:
        records = [records_by_protein[pid][modality] for pid in protein_ids]
        statuses = Counter(record.factual_status for record in records)
        actions = Counter(record.requested_action for record in records)
        by_modality[modality] = {
            "present": _metric(sum(record.exists for record in records), total),
            "valid": _metric(sum(record.valid for record in records), total),
            "reusable": _metric(actions["reuse"], total),
            "missing": _metric(statuses["missing"], total),
            "generation": _metric(actions["generate"], total),
            "masked": _metric(actions["leave-masked"], total),
            "unavailable": _metric(actions["unavailable"], total),
            "manual_review": _metric(actions["manual-review"], total),
            "factual_status_counts": dict(sorted(statuses.items())),
            "requested_action_counts": dict(sorted(actions.items())),
        }
    at_least_one = sum(any(records_by_protein[p][m].exists for m in MODALITIES) for p in protein_ids)
    complete_present = sum(all(records_by_protein[p][m].exists for m in MODALITIES) for p in protein_ids)
    complete_valid = sum(all(records_by_protein[p][m].valid for m in MODALITIES) for p in protein_ids)
    complete_reusable = sum(all(records_by_protein[p][m].requested_action == "reuse" for m in MODALITIES) for p in protein_ids)
    actions = Counter(records_by_protein[p][m].requested_action for p in protein_ids for m in MODALITIES)
    return {
        "population": total,
        "at_least_one_modality": _metric(at_least_one, total),
        "complete_four_modalities_present": _metric(complete_present, total),
        "complete_four_modalities_valid": _metric(complete_valid, total),
        "complete_four_modalities_reusable": _metric(complete_reusable, total),
        "by_modality": by_modality,
        "requested_action_counts": dict(sorted(actions.items())),
    }


def _summary_markdown(summary: Mapping[str, Any]) -> str:
    global_coverage = summary["coverage"]["global"]
    lines = [
        "# Embedding inventory summary", "",
        "- Configuration: `%s`" % summary["config_name"],
        "- Policy / report level: `%s` / `%s`" % (summary["policy"], summary["report_level"]),
        "- Benchmark proteins / unique sequences: %s / %s" % (_integer(summary["population"]), _integer(summary["unique_sequences"])),
        "- Verified exact artifact scope: `%s`" % str(summary["artifact_verification"]["verified"]).lower(),
        "- At least one physical modality: %s" % _format_metric(global_coverage["at_least_one_modality"]),
        "- Complete four-modality physical coverage: %s" % _format_metric(global_coverage["complete_four_modalities_present"]),
        "- Complete four-modality reusable coverage: %s" % _format_metric(global_coverage["complete_four_modalities_reusable"]),
        "", "## Global modality coverage", "",
        "| Modality | Present | Valid | Reuse | Missing | Generate | Masked | Unavailable | Manual review |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for modality in MODALITIES:
        values = global_coverage["by_modality"][modality]
        lines.append("| %s | %s | %s | %s | %s | %s | %s | %s | %s |" % (
            modality, _format_metric(values["present"]), _format_metric(values["valid"]),
            _format_metric(values["reusable"]), _format_metric(values["missing"]),
            _format_metric(values["generation"]), _format_metric(values["masked"]),
            _format_metric(values["unavailable"]), _format_metric(values["manual_review"]),
        ))
    lines.extend(["", "## Interpretation", "", summary["pfp_missing_behavior"], "",
        "Physical validity and scientific reuse eligibility are separate. Exact published-artifact reuse requires the recorded cryptographic proof. Cross-benchmark temporal text, structure, and PPI remain conservative when their required context or mapping evidence is absent.", ""])
    return "\n".join(lines)


def _modalities_with_action(records: Mapping[str, InventoryRecord], action: str) -> List[str]:
    return [m for m in MODALITIES if records[m].requested_action == action]


def _metric(count: int, total: int) -> Dict[str, Any]:
    return {"count": int(count), "total": int(total), "fraction": count / total if total else 0.0}


def _format_metric(metric: Mapping[str, Any]) -> str:
    return "%s/%s (%.2f%%)" % (_integer(metric["count"]), _integer(metric["total"]), 100.0 * metric["fraction"])


def _integer(value: int) -> str:
    return format(int(value), ",")


def _bool(value: bool) -> str:
    return "true" if value else "false"
