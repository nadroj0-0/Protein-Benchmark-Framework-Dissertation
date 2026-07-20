#!/usr/bin/env python3
"""Compare any number of benchmark-agnostic PFP label-space audit reports."""

from __future__ import annotations

import argparse
import csv
import io
import json
import os
import shutil
import tempfile
from pathlib import Path
from typing import Any

from label_space_common import (
    ASPECTS,
    SCHEMA_VERSION,
    atomic_write_json,
    atomic_write_text,
    file_snapshot,
    output_manifest,
    require_unchanged,
    sha256_file,
)


def load_report(path: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    completion_path = path.parent / "RUN_COMPLETE.json"
    output_path = path.parent / "output_manifest.json"
    report_snapshot = file_snapshot(path)
    completion_snapshot = file_snapshot(completion_path)
    output_snapshot = file_snapshot(output_path)
    completion = json.loads(completion_path.read_text(encoding="utf-8"))
    if not completion.get("complete"):
        raise ValueError(f"Label-space completion marker is false: {path.parent}")
    if completion.get("output_manifest_sha256") != sha256_file(output_path):
        raise ValueError(f"Label-space completion marker does not bind its manifest: {path}")
    output = json.loads(output_path.read_text(encoding="utf-8"))
    listed = {item.get("path"): item for item in output.get("files", [])}
    if path.name not in listed:
        raise ValueError(f"Label-space output manifest does not bind its report: {path}")
    if (
        listed[path.name].get("bytes") != report_snapshot["bytes"]
        or listed[path.name].get("sha256") != report_snapshot["sha256"]
    ):
        raise ValueError(f"Label-space manifest binds different report bytes: {path}")
    for item in output.get("files", []):
        relative = Path(item["path"])
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError(f"Label-space manifest contains unsafe path: {relative}")
        artifact = path.parent / relative
        if not artifact.is_file():
            raise FileNotFoundError(f"Label-space artifact is missing: {artifact}")
        if (
            artifact.stat().st_size != item["bytes"]
            or sha256_file(artifact) != item["sha256"]
        ):
            raise ValueError(f"Label-space artifact changed after publication: {artifact}")
    value = json.loads(path.read_text(encoding="utf-8"))
    if (
        not isinstance(value, dict)
        or value.get("schema_version") != SCHEMA_VERSION
        or value.get("status") != "passed"
    ):
        raise ValueError(f"Unsupported or incomplete label-space audit: {path}")
    require_unchanged(path, report_snapshot, "Label-space report")
    require_unchanged(completion_path, completion_snapshot, "Completion marker")
    require_unchanged(output_path, output_snapshot, "Output manifest")
    return value, report_snapshot


def markdown_report(comparison: dict[str, Any]) -> str:
    lines = [
        "# PFP Label-Space Comparison",
        "",
        "| Benchmark | Aspect | Test proteins | Root-only | Root-only % | Mean non-root labels | Mean max depth | Root-only baseline F |",
        "|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in comparison["rows"]:
        mean_depth = row["mean_max_shortest_depth"]
        lines.append(
            "| {benchmark_id} | {aspect} | {test_proteins:,} | {root_only_rows:,} | "
            "{root_only_fraction:.2%} | {mean_non_root_labels:.2f} | "
            "{mean_max_depth} | {root_only_baseline_f:.3f} |".format(
                **row,
                mean_max_depth=(f"{mean_depth:.2f}" if mean_depth is not None else "n/a"),
            )
        )
    lines.extend(
        [
            "",
            "Root-only prevalence and the arithmetic baseline describe benchmark composition; they do not rank model quality.",
            "All comparisons use the immutable hashes recorded by each source audit.",
            "",
        ]
    )
    return "\n".join(lines)


def tsv_text(rows: list[dict[str, Any]]) -> str:
    fields = [
        "benchmark_id",
        "benchmark_fingerprint",
        "aspect",
        "test_proteins",
        "terms",
        "root_only_rows",
        "root_only_fraction",
        "rows_with_non_root_labels",
        "mean_labels",
        "mean_non_root_labels",
        "mean_max_shortest_depth",
        "median_max_shortest_depth",
        "maximum_shortest_depth",
        "root_only_baseline_f",
    ]
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(
        buffer, fieldnames=fields, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--report", action="append", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    if len(args.report) < 2:
        raise ValueError("At least two --report inputs are required")
    output_dir = args.output_dir.resolve()
    if output_dir.exists():
        raise ValueError(f"Output directory already exists: {output_dir}")
    loaded = [
        (path.resolve(), *load_report(path.resolve())) for path in args.report
    ]
    ids = [report["benchmark_id"] for _, report, _ in loaded]
    if len(ids) != len(set(ids)):
        raise ValueError("Each compared benchmark_id must be unique")

    rows = []
    sources = []
    for path, report, snapshot in loaded:
        sources.append(
            {
                "path": str(path),
                "sha256": snapshot["sha256"],
                "benchmark_id": report["benchmark_id"],
                "benchmark_fingerprint": report["benchmark_fingerprint"],
                "metadata": report.get("metadata", {}),
            }
        )
        for aspect in ASPECTS:
            item = report["files"][f"{aspect}:test"]
            rows.append(
                {
                    "benchmark_id": report["benchmark_id"],
                    "benchmark_fingerprint": report["benchmark_fingerprint"],
                    "aspect": aspect,
                    "test_proteins": item["rows"],
                    "terms": item["terms"],
                    "root_only_rows": item["root_only_rows"],
                    "root_only_fraction": item["root_only_fraction"],
                    "rows_with_non_root_labels": item["rows_with_non_root_labels"],
                    "mean_labels": item["labels_per_protein"]["mean"],
                    "mean_non_root_labels": item["non_root_labels_per_protein"]["mean"],
                    "mean_max_shortest_depth": item[
                        "max_shortest_depth_per_non_root_row"
                    ]["mean"],
                    "median_max_shortest_depth": item[
                        "max_shortest_depth_per_non_root_row"
                    ]["median"],
                    "maximum_shortest_depth": item[
                        "max_shortest_depth_per_non_root_row"
                    ]["maximum"],
                    "root_only_baseline_f": item["root_only_diagnostic_baseline"]["macro_f"],
                }
            )
    comparison = {
        "schema_version": SCHEMA_VERSION,
        "status": "passed",
        "sources": sources,
        "rows": rows,
        "interpretation_policy": {
            "root_only_is_failure_criterion": False,
            "model_ranking_permitted": False,
            "note": "This report compares benchmark label spaces, not trained models.",
        },
    }

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(
        tempfile.mkdtemp(prefix=f".{output_dir.name}.stage-", dir=str(output_dir.parent))
    )
    try:
        atomic_write_json(stage / "label_space_comparison.json", comparison)
        atomic_write_text(stage / "label_space_comparison.md", markdown_report(comparison))
        atomic_write_text(stage / "label_space_comparison.tsv", tsv_text(rows))
        manifest = output_manifest(
            stage, exclude={"output_manifest.json", "RUN_COMPLETE.json"}
        )
        atomic_write_json(stage / "output_manifest.json", manifest)
        atomic_write_json(
            stage / "RUN_COMPLETE.json",
            {
                "schema_version": 1,
                "complete": True,
                "compared_benchmarks": ids,
                "output_manifest_sha256": sha256_file(stage / "output_manifest.json"),
            },
        )
        os.replace(stage, output_dir)
    except BaseException:
        shutil.rmtree(stage, ignore_errors=True)
        raise
    print(json.dumps({"status": "passed", "output_dir": str(output_dir)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
