import csv
import json
from collections import Counter, defaultdict
from dataclasses import asdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Sequence, Set

from .models import InventoryRecord, InventoryResult, MODALITIES, ONTOLOGIES, SPLITS


INVENTORY_FIELDS = [
    "protein_id",
    "sequence_sha256",
    "modality",
    "source_directory",
    "source_file",
    "exists",
    "observed_shape",
    "expected_shape",
    "dtype",
    "finite",
    "valid",
    "scientifically_eligible",
    "source_protein_id",
    "match_route",
    "sequence_match",
    "provenance",
    "factual_status",
    "requested_action",
    "reason",
]


def write_reports(result: InventoryResult, output_dir: Path, embedding_cache: Path) -> Dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    records_by_protein: Dict[str, Dict[str, InventoryRecord]] = defaultdict(dict)
    for record in result.records:
        records_by_protein[record.protein_id][record.modality] = record

    _write_benchmark_proteins(result, output_dir / "benchmark_proteins.tsv")
    _write_inventory(result.records, output_dir / "embedding_inventory.tsv")
    _write_protein_summary(result, records_by_protein, output_dir / "protein_embedding_summary.tsv")
    _write_inventory(
        (record for record in result.records if record.requested_action == "reuse"),
        output_dir / "reuse_manifest.tsv",
    )
    _write_inventory(
        (record for record in result.records if record.requested_action == "generate"),
        output_dir / "generation_manifest.tsv",
    )
    _write_inventory(
        (record for record in result.records if record.requested_action == "manual-review"),
        output_dir / "manual_review.tsv",
    )
    extras = _write_cache_extras(result, output_dir / "cache_extras.tsv")
    _write_modality_lists(result, output_dir)

    summary = _build_summary(result, records_by_protein, extras, embedding_cache)
    (output_dir / "embedding_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (output_dir / "embedding_summary.md").write_text(
        _summary_markdown(summary), encoding="utf-8"
    )
    return summary


def _write_benchmark_proteins(result: InventoryResult, path: Path) -> None:
    fields = [
        "protein_id",
        "sequence",
        "sequence_sha256",
        "sequence_length",
        "ontologies",
        "splits",
        "source_csv_files",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for protein_id in sorted(result.benchmark.proteins):
            protein = result.benchmark.proteins[protein_id]
            writer.writerow(
                {
                    "protein_id": protein.protein_id,
                    "sequence": protein.sequence,
                    "sequence_sha256": protein.sequence_sha256,
                    "sequence_length": protein.sequence_length,
                    "ontologies": ";".join(sorted(protein.ontologies)),
                    "splits": ";".join(sorted(protein.splits)),
                    "source_csv_files": ";".join(sorted(protein.source_files)),
                }
            )


def _write_inventory(records: Iterable[InventoryRecord], path: Path) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=INVENTORY_FIELDS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for record in records:
            writer.writerow(record.as_dict())


def _write_protein_summary(
    result: InventoryResult,
    records_by_protein: Mapping[str, Mapping[str, InventoryRecord]],
    path: Path,
) -> None:
    fields = [
        "protein_id",
        "has_any_embedding",
        "has_prott5",
        "has_text",
        "has_structure",
        "has_ppi",
        "valid_modalities",
        "reusable_modalities",
        "missing_modalities",
        "generation_modalities",
        "masked_modalities",
        "unavailable_modalities",
        "manual_review_modalities",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for protein_id in sorted(records_by_protein):
            records = records_by_protein[protein_id]
            present = [modality for modality in MODALITIES if records[modality].exists]
            valid = [modality for modality in MODALITIES if records[modality].valid]
            reusable = _modalities_with_action(records, "reuse")
            generation = _modalities_with_action(records, "generate")
            masked = _modalities_with_action(records, "leave-masked")
            unavailable = _modalities_with_action(records, "unavailable")
            manual = _modalities_with_action(records, "manual-review")
            missing = [
                modality for modality in MODALITIES if records[modality].factual_status == "missing"
            ]
            writer.writerow(
                {
                    "protein_id": protein_id,
                    "has_any_embedding": _bool(bool(present)),
                    "has_prott5": _bool("prott5" in present),
                    "has_text": _bool("text" in present),
                    "has_structure": _bool("structure" in present),
                    "has_ppi": _bool("ppi" in present),
                    "valid_modalities": ";".join(valid),
                    "reusable_modalities": ";".join(reusable),
                    "missing_modalities": ";".join(missing),
                    "generation_modalities": ";".join(generation),
                    "masked_modalities": ";".join(masked),
                    "unavailable_modalities": ";".join(unavailable),
                    "manual_review_modalities": ";".join(manual),
                }
            )


def _write_cache_extras(result: InventoryResult, path: Path) -> Dict[str, int]:
    target_ids = set(result.benchmark.proteins)
    counts: Dict[str, int] = {}
    fields = ["modality", "protein_id", "source_file", "used_as_reuse_source"]
    with path.open("w", newline="", encoding="utf-8") as handle:
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
                        "source_file": str(
                            Path(result.config.modalities[modality].directory)
                            / (protein_id + ".npy")
                        ),
                        "used_as_reuse_source": _bool(
                            protein_id in result.used_source_ids[modality]
                        ),
                    }
                )
    return counts


def _write_modality_lists(result: InventoryResult, output_dir: Path) -> None:
    by_modality: Dict[str, List[InventoryRecord]] = defaultdict(list)
    for record in result.records:
        by_modality[record.modality].append(record)

    _write_id_list(output_dir / "reuse_prott5.txt", by_modality["prott5"], "reuse")
    _write_fasta(
        output_dir / "generate_prott5.fasta",
        by_modality["prott5"],
        result,
    )
    _write_id_list(output_dir / "reuse_text.txt", by_modality["text"], "reuse")
    _write_id_list(output_dir / "generate_text.txt", by_modality["text"], "generate")
    _write_id_list(
        output_dir / "text_manual_review.txt", by_modality["text"], "manual-review"
    )
    _write_id_list(output_dir / "reuse_structure.txt", by_modality["structure"], "reuse")
    _write_id_list(
        output_dir / "generate_structure.txt", by_modality["structure"], "generate"
    )
    _write_id_list(
        output_dir / "structure_unavailable.txt", by_modality["structure"], "unavailable"
    )
    _write_id_list(output_dir / "reuse_ppi.txt", by_modality["ppi"], "reuse")
    _write_id_list(output_dir / "extract_ppi.txt", by_modality["ppi"], "generate")
    _write_id_list(output_dir / "ppi_unavailable.txt", by_modality["ppi"], "unavailable")


def _write_id_list(path: Path, records: Sequence[InventoryRecord], action: str) -> None:
    ids = sorted(record.protein_id for record in records if record.requested_action == action)
    path.write_text("".join(protein_id + "\n" for protein_id in ids), encoding="utf-8")


def _write_fasta(path: Path, records: Sequence[InventoryRecord], result: InventoryResult) -> None:
    ids = sorted(record.protein_id for record in records if record.requested_action == "generate")
    with path.open("w", encoding="utf-8") as handle:
        for protein_id in ids:
            sequence = result.benchmark.proteins[protein_id].sequence
            handle.write(">%s\n" % protein_id)
            for start in range(0, len(sequence), 60):
                handle.write(sequence[start : start + 60] + "\n")


def _build_summary(
    result: InventoryResult,
    records_by_protein: Mapping[str, Mapping[str, InventoryRecord]],
    extras: Dict[str, int],
    embedding_cache: Path,
) -> Dict[str, Any]:
    all_ids = set(result.benchmark.proteins)
    by_split = {
        split: {
            protein_id
            for protein_id, protein in result.benchmark.proteins.items()
            if split in protein.splits
        }
        for split in SPLITS
    }
    by_ontology = {
        ontology: {
            protein_id
            for protein_id, protein in result.benchmark.proteins.items()
            if ontology in protein.ontologies
        }
        for ontology in ONTOLOGIES
    }
    by_ontology_split = {
        "%s/%s" % (ontology, split): set(result.benchmark.file_members[(ontology, split)])
        for ontology in ONTOLOGIES
        for split in SPLITS
    }
    return {
        "schema_version": 1,
        "config_name": result.config.name,
        "policy": result.policy,
        "benchmark_directory": str(result.benchmark.directory),
        "source_benchmark_directory": str(result.source_benchmark.directory),
        "embedding_cache": str(embedding_cache.resolve()),
        "population": len(all_ids),
        "unique_sequences": len(
            {protein.sequence_sha256 for protein in result.benchmark.proteins.values()}
        ),
        "duplicate_identical_rows": result.benchmark.duplicate_rows,
        "pfp_missing_behavior": (
            "Absent or unreadable arrays become zero vectors with mask 0.0 in "
            "PFP mmfp.dataset.MultiModalDataset._load_embedding_raw()."
        ),
        "configuration": {
            "benchmark_contract": asdict(result.config.benchmark_contract),
            "modalities": {
                modality: asdict(result.config.modalities[modality]) for modality in MODALITIES
            },
        },
        "coverage": {
            "global": _group_coverage(all_ids, records_by_protein),
            "by_split": {
                split: _group_coverage(ids, records_by_protein) for split, ids in by_split.items()
            },
            "by_ontology": {
                ontology: _group_coverage(ids, records_by_protein)
                for ontology, ids in by_ontology.items()
            },
            "by_ontology_split": {
                name: _group_coverage(ids, records_by_protein)
                for name, ids in by_ontology_split.items()
            },
        },
        "cache_files": {modality: len(result.cache_ids[modality]) for modality in MODALITIES},
        "cache_extras": extras,
    }


def _group_coverage(
    protein_ids: Set[str],
    records_by_protein: Mapping[str, Mapping[str, InventoryRecord]],
) -> Dict[str, Any]:
    total = len(protein_ids)
    by_modality: Dict[str, Any] = {}
    for modality in MODALITIES:
        records = [records_by_protein[protein_id][modality] for protein_id in protein_ids]
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
    at_least_one = sum(
        any(records_by_protein[protein_id][modality].exists for modality in MODALITIES)
        for protein_id in protein_ids
    )
    complete_present = sum(
        all(records_by_protein[protein_id][modality].exists for modality in MODALITIES)
        for protein_id in protein_ids
    )
    complete_valid = sum(
        all(records_by_protein[protein_id][modality].valid for modality in MODALITIES)
        for protein_id in protein_ids
    )
    complete_reusable = sum(
        all(
            records_by_protein[protein_id][modality].requested_action == "reuse"
            for modality in MODALITIES
        )
        for protein_id in protein_ids
    )
    actions = Counter(
        records_by_protein[protein_id][modality].requested_action
        for protein_id in protein_ids
        for modality in MODALITIES
    )
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
        "# Embedding inventory summary",
        "",
        "- Configuration: `%s`" % summary["config_name"],
        "- Policy: `%s`" % summary["policy"],
        "- Benchmark proteins: %s" % _integer(summary["population"]),
        "- Unique complete sequences: %s" % _integer(summary["unique_sequences"]),
        "- At least one physical modality: %s" % _format_metric(global_coverage["at_least_one_modality"]),
        "- Complete four-modality physical coverage: %s"
        % _format_metric(global_coverage["complete_four_modalities_present"]),
        "- Complete four-modality reusable coverage: %s"
        % _format_metric(global_coverage["complete_four_modalities_reusable"]),
        "",
        "## Global modality coverage",
        "",
        "| Modality | Present | Valid | Reuse | Missing | Generate | Masked | Unavailable | Manual review |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for modality in MODALITIES:
        values = global_coverage["by_modality"][modality]
        lines.append(
            "| %s | %s | %s | %s | %s | %s | %s | %s | %s |"
            % (
                modality,
                _format_metric(values["present"]),
                _format_metric(values["valid"]),
                _format_metric(values["reusable"]),
                _format_metric(values["missing"]),
                _format_metric(values["generation"]),
                _format_metric(values["masked"]),
                _format_metric(values["unavailable"]),
                _format_metric(values["manual_review"]),
            )
        )

    lines.extend(["", "## Population slices", "", "| Slice | Proteins | Complete four present |", "|---|---:|---:|"])
    for split in SPLITS:
        values = summary["coverage"]["by_split"][split]
        lines.append(
            "| %s | %s | %s |"
            % (split, _integer(values["population"]), _format_metric(values["complete_four_modalities_present"]))
        )
    for ontology in ONTOLOGIES:
        values = summary["coverage"]["by_ontology"][ontology]
        lines.append(
            "| %s | %s | %s |"
            % (ontology, _integer(values["population"]), _format_metric(values["complete_four_modalities_present"]))
        )

    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            summary["pfp_missing_behavior"],
            "",
            "Physical validity and scientific reuse eligibility are separate. `paper-faithful` leaves "
            "missing or unreadable published modalities masked; `maximize-coverage` applies each modality's "
            "configured generation or unavailability rule. Provenance uncertainty always routes to manual review.",
            "",
            "Text source selection is singular and explicit in the configuration; no current/temporal fallback is used.",
            "",
        ]
    )
    return "\n".join(lines)


def _modalities_with_action(records: Mapping[str, InventoryRecord], action: str) -> List[str]:
    return [modality for modality in MODALITIES if records[modality].requested_action == action]


def _metric(count: int, total: int) -> Dict[str, Any]:
    return {"count": int(count), "total": int(total), "fraction": (count / total if total else 0.0)}


def _format_metric(metric: Mapping[str, Any]) -> str:
    return "%s/%s (%.2f%%)" % (
        _integer(metric["count"]),
        _integer(metric["total"]),
        100.0 * metric["fraction"],
    )


def _integer(value: int) -> str:
    return format(int(value), ",")


def _bool(value: bool) -> str:
    return "true" if value else "false"
