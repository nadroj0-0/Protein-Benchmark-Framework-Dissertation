#!/usr/bin/env python3
"""Combine validated PFP specificity analyses without erasing benchmark scope."""

from __future__ import annotations

import argparse
import csv
import io
import json
import os
import shutil
import tempfile
from pathlib import Path
from typing import Any, Iterable, Mapping

from label_space_common import (
    atomic_write_json,
    atomic_write_text,
    output_manifest,
    sha256_file,
)


ASPECTS = ("BPO", "CCO", "MFO")
MEASURES = ("ia", "xu_totipotency_raw")
BOOTSTRAP_METRICS = {
    "f": "fixed_unweighted_f",
    "weighted_f": "fixed_weighted_f",
    "jaccard_set_agreement": "fixed_unweighted_jaccard",
    "weighted_jaccard_set_agreement": "fixed_weighted_jaccard",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare complete IA/Xu specificity analyses descriptively."
    )
    parser.add_argument(
        "--source",
        action="append",
        required=True,
        metavar="LABEL=PATH",
        help="Repeat for each completed specificity analysis.",
    )
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def parse_sources(values: Iterable[str]) -> list[tuple[str, Path]]:
    result: list[tuple[str, Path]] = []
    for value in values:
        label, separator, raw_path = value.partition("=")
        if not separator or not label or not raw_path:
            raise ValueError(f"--source must be LABEL=PATH, received: {value}")
        if not all(character.isalnum() or character in "._-" for character in label):
            raise ValueError(f"Unsafe source label: {label}")
        path = Path(raw_path).expanduser().resolve()
        if (path / "analysis" / "specificity_bins.tsv").is_file():
            path = path / "analysis"
        result.append((label, path))
    labels = [label for label, _path in result]
    if len(result) < 2:
        raise ValueError("At least two specificity sources are required")
    if len(labels) != len(set(labels)):
        raise ValueError("Specificity source labels must be unique")
    return result


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def tsv(rows: list[Mapping[str, Any]], fields: tuple[str, ...]) -> str:
    output = io.StringIO()
    writer = csv.DictWriter(
        output, fieldnames=fields, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows({field: row.get(field, "") for field in fields} for row in rows)
    return output.getvalue()


def verify_manifest(root: Path) -> dict[str, Any]:
    manifest_path = root / "output_manifest.json"
    marker_path = root / "RUN_COMPLETE.json"
    if not manifest_path.is_file() or not marker_path.is_file():
        raise FileNotFoundError(f"Incomplete specificity source: {root}")
    marker = json.loads(marker_path.read_text(encoding="utf-8"))
    if marker.get("complete") is not True:
        raise ValueError(
            f"Specificity completion marker is not complete: {marker_path}"
        )
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    for item in manifest.get("files", []):
        relative = Path(item["path"])
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError(f"Unsafe manifest path in {manifest_path}: {relative}")
        path = root / relative
        if not path.is_file():
            raise FileNotFoundError(f"Manifest payload is missing: {path}")
        if path.stat().st_size != item["bytes"] or sha256_file(path) != item["sha256"]:
            raise ValueError(f"Manifest payload changed: {path}")
    return {
        "path": str(root),
        "run_complete_sha256": sha256_file(marker_path),
        "output_manifest_sha256": sha256_file(manifest_path),
    }


def load_source(label: str, root: Path) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    provenance = verify_manifest(root)
    bins = read_tsv(root / "specificity_bins.tsv")
    bootstraps = read_tsv(root / "bootstrap_intervals.tsv")
    report = json.loads(
        (root / "specificity_analysis.json").read_text(encoding="utf-8")
    )
    if report.get("status") != "complete":
        raise ValueError(f"Specificity report is not complete: {root}")
    if report.get("specificity_measures") != list(MEASURES):
        raise ValueError(
            f"{label} does not contain the required separate IA and raw Xu panels"
        )
    benchmark_ids = {row["benchmark_id"] for row in bins}
    modes = {row["mode"] for row in bins}
    if len(benchmark_ids) != 1 or modes != {"full"}:
        raise ValueError(f"{label} must contain one full-model benchmark")

    bootstrap_by_key = {
        (
            row["aspect"],
            row["specificity_measure"],
            row["specificity_bin"],
            row["metric_name"],
        ): row
        for row in bootstraps
    }
    combined: list[dict[str, Any]] = []
    for row in bins:
        value: dict[str, Any] = {"source_label": label, **row}
        for metric_name, output_name in BOOTSTRAP_METRICS.items():
            bootstrap = bootstrap_by_key.get(
                (
                    row["aspect"],
                    row["specificity_measure"],
                    row["specificity_bin"],
                    metric_name,
                )
            )
            value[f"{output_name}_ci_low"] = bootstrap["ci_low"] if bootstrap else ""
            value[f"{output_name}_ci_high"] = bootstrap["ci_high"] if bootstrap else ""
        combined.append(value)
    provenance.update(
        {
            "label": label,
            "benchmark_id": next(iter(benchmark_ids)),
            "mode": "full",
            "specificity_analysis_sha256": sha256_file(
                root / "specificity_analysis.json"
            ),
            "specificity_bins_sha256": sha256_file(root / "specificity_bins.tsv"),
            "bootstrap_intervals_sha256": sha256_file(root / "bootstrap_intervals.tsv"),
        }
    )
    return provenance, combined


def broad_to_specific(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    indexed = {
        (
            row["source_label"],
            row["aspect"],
            row["specificity_measure"],
            row["specificity_bin"],
        ): row
        for row in rows
    }
    output: list[dict[str, Any]] = []
    for label in sorted({row["source_label"] for row in rows}):
        for aspect in ASPECTS:
            for measure in MEASURES:
                q1 = indexed.get((label, aspect, measure, "specificity_q1"))
                q4 = indexed.get((label, aspect, measure, "specificity_q4"))
                if q1 is None or q4 is None:
                    raise ValueError(
                        f"Missing Q1/Q4 bins for {label}/{aspect}/{measure}"
                    )
                q1_f = float(q1["fixed_unweighted_f"])
                q4_f = float(q4["fixed_unweighted_f"])
                q1_wf = float(q1["fixed_weighted_f"])
                q4_wf = float(q4["fixed_weighted_f"])
                output.append(
                    {
                        "source_label": label,
                        "benchmark_id": q1["benchmark_id"],
                        "aspect": aspect,
                        "specificity_measure": measure,
                        "fixed_threshold": q1["fixed_threshold"],
                        "q1_fixed_f": q1_f,
                        "q1_fixed_f_ci_low": q1["fixed_unweighted_f_ci_low"],
                        "q1_fixed_f_ci_high": q1["fixed_unweighted_f_ci_high"],
                        "q4_fixed_f": q4_f,
                        "q4_fixed_f_ci_low": q4["fixed_unweighted_f_ci_low"],
                        "q4_fixed_f_ci_high": q4["fixed_unweighted_f_ci_high"],
                        "q4_minus_q1_fixed_f": q4_f - q1_f,
                        "q4_over_q1_fixed_f": q4_f / q1_f if q1_f else "",
                        "q1_fixed_weighted_f": q1_wf,
                        "q4_fixed_weighted_f": q4_wf,
                        "q4_minus_q1_fixed_weighted_f": q4_wf - q1_wf,
                        "q1_positive_proteins": q1["target_positive_proteins"],
                        "q4_positive_proteins": q4["target_positive_proteins"],
                    }
                )
    return output


def markdown(rows: list[dict[str, Any]]) -> str:
    lines = [
        "# Three-Way PFP Specificity Comparison",
        "",
        "All values below are fixed-threshold flat diagnostics from independently ",
        "validated full-model prediction arrays. Q1 is broad and Q4 is specific for ",
        "both IA and raw Xu totipotency.",
        "",
        "| Benchmark | Aspect | Measure | Q1 F | Q4 F | Q4 - Q1 | Q1 positives | Q4 positives |",
        "|---|---|---|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            "| {source_label} | {aspect} | {specificity_measure} | {q1_fixed_f:.3f} | "
            "{q4_fixed_f:.3f} | {q4_minus_q1_fixed_f:+.3f} | "
            "{q1_positive_proteins} | {q4_positive_proteins} |".format(**row)
        )
    lines.extend(
        [
            "",
            "## Interpretation limits",
            "",
            "- These are robustness measurements, not a controlled trend line.",
            "- Each benchmark has its own training-derived IA values, term universe and fixed threshold.",
            "- Quartiles are benchmark-local rank groups; they do not contain identical GO terms.",
            "- CAFA3 also uses an older ontology snapshot. Xu panels are therefore comparable in direction, not term-for-term identity.",
            "- Bootstrap intervals quantify protein-sampling uncertainty, not training-seed uncertainty.",
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> int:
    args = parse_args()
    sources = parse_sources(args.source)
    output_dir = args.output_dir.expanduser().resolve()
    if output_dir.exists():
        raise FileExistsError(f"Output already exists: {output_dir}")
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(
        tempfile.mkdtemp(prefix=f".{output_dir.name}.staging-", dir=output_dir.parent)
    )
    try:
        provenance: list[dict[str, Any]] = []
        rows: list[dict[str, Any]] = []
        for label, path in sources:
            source_provenance, source_rows = load_source(label, path)
            provenance.append(source_provenance)
            rows.extend(source_rows)
        summary = broad_to_specific(rows)
        long_fields = tuple(rows[0])
        summary_fields = tuple(summary[0])
        atomic_write_text(
            stage / "specificity_comparison_long.tsv", tsv(rows, long_fields)
        )
        atomic_write_text(
            stage / "broad_to_specific_summary.tsv", tsv(summary, summary_fields)
        )
        atomic_write_text(
            stage / "specificity_three_way_comparison.md", markdown(summary)
        )
        atomic_write_json(
            stage / "comparison.json",
            {
                "schema_version": 1,
                "complete": True,
                "analysis_kind": "cross_benchmark_specificity_comparison",
                "sources": provenance,
                "source_count": len(provenance),
                "comparison_policy": "descriptive benchmark-local quartiles; no controlled causal ordering",
            },
        )
        atomic_write_json(
            stage / "output_manifest.json",
            output_manifest(
                stage, exclude={"output_manifest.json", "RUN_COMPLETE.json"}
            ),
        )
        atomic_write_json(
            stage / "RUN_COMPLETE.json",
            {
                "schema_version": 1,
                "complete": True,
                "analysis_kind": "cross_benchmark_specificity_comparison",
                "source_labels": [label for label, _path in sources],
                "output_manifest_sha256": sha256_file(stage / "output_manifest.json"),
            },
        )
        os.replace(stage, output_dir)
    except BaseException:
        shutil.rmtree(stage, ignore_errors=True)
        raise
    print(json.dumps({"complete": True, "output_dir": str(output_dir)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
