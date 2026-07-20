#!/usr/bin/env python3
"""Audit the label space of any benchmark implementing PFP's nine-CSV contract."""

from __future__ import annotations

import argparse
import csv
import io
import json
import os
import shutil
import tempfile
import time
from pathlib import Path
from typing import Any

from label_space_common import (
    ASPECTS,
    ASPECT_TO_PREFIX,
    ASPECT_TO_ROOT,
    CSV_SPLITS,
    SCHEMA_VERSION,
    atomic_write_json,
    atomic_write_text,
    audit_csv,
    file_snapshot,
    output_manifest,
    peak_rss_bytes,
    read_ia_file,
    read_obo,
    require_unchanged,
    required_csv_names,
    sha256_file,
    sha256_json,
    verify_prepared_data,
)


def parse_assignment(value: str, option: str) -> tuple[str, str]:
    if "=" not in value:
        raise ValueError(f"{option} must use NAME=VALUE syntax: {value!r}")
    name, raw = value.split("=", 1)
    if not name or not raw:
        raise ValueError(f"{option} contains an empty name or value: {value!r}")
    return name, raw


def load_policy(config_path: Path | None) -> tuple[dict[str, Any], dict[str, Any] | None]:
    if config_path is None:
        return {
            "allow_legacy_singular_protein_header": False,
            "allow_all_zero_rows": False,
        }, None
    config = json.loads(config_path.read_text(encoding="utf-8"))
    if not isinstance(config, dict) or config.get("schema_version") != 1:
        raise ValueError(f"Unsupported run config: {config_path}")
    contract = config.get("benchmark_contract")
    if not isinstance(contract, dict):
        raise ValueError(f"Run config has no benchmark_contract: {config_path}")
    for key in (
        "allow_legacy_singular_protein_header",
        "allow_all_zero_rows",
    ):
        if not isinstance(contract.get(key, False), bool):
            raise ValueError(f"Run config {key} must be boolean: {config_path}")
    return {
        "allow_legacy_singular_protein_header": bool(
            contract.get("allow_legacy_singular_protein_header", False)
        ),
        "allow_all_zero_rows": bool(contract.get("allow_all_zero_rows", False)),
    }, config


def tsv_text(fieldnames: list[str], rows: list[dict[str, Any]]) -> str:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(
        buffer,
        fieldnames=fieldnames,
        delimiter="\t",
        lineterminator="\n",
        extrasaction="ignore",
    )
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue()


def markdown_report(report: dict[str, Any]) -> str:
    lines = [
        f"# PFP Label-Space Audit: {report['benchmark_id']}",
        "",
        f"- Status: `{report['status']}`",
        f"- Benchmark fingerprint: `{report['benchmark_fingerprint']}`",
        f"- OBO SHA-256: `{report['inputs']['obo']['sha256']}`",
        f"- Prepared-data references verified: {len(report['prepared_data_verification'])}",
        "",
        "## Test Label Profile",
        "",
        "| Aspect | Proteins | Terms | Root-only | Root-only % | Mean labels | Mean non-root | Root-only baseline F |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for aspect in ASPECTS:
        item = report["files"][f"{aspect}:test"]
        lines.append(
            "| {aspect} | {rows:,} | {terms:,} | {root_only_rows:,} | "
            "{root_only_fraction:.2%} | {labels:.2f} | {non_root:.2f} | {baseline:.3f} |".format(
                aspect=aspect,
                rows=item["rows"],
                terms=item["terms"],
                root_only_rows=item["root_only_rows"],
                root_only_fraction=item["root_only_fraction"],
                labels=item["labels_per_protein"]["mean"],
                non_root=item["non_root_labels_per_protein"]["mean"],
                baseline=item["root_only_diagnostic_baseline"]["macro_f"],
            )
        )
    lines.extend(
        [
            "",
            "The root-only baseline is an arithmetic diagnostic over the retained label matrix, not a cafaeval result.",
            "Root-only prevalence is reported as quality-control evidence and is not a failure criterion.",
            "",
            "## Files",
            "",
            "- `label_space_audit.json`: complete machine-readable report.",
            "- `root_only_targets.tsv`: exact root-only target membership.",
            "- `term_support.tsv`: support, depth and optional IA for every retained term.",
            "- `label_depth_histogram.tsv`: per-file maximum shortest-depth histogram.",
            "- `input_manifest.json`: immutable input hashes and policy.",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--benchmark-id", required=True)
    parser.add_argument("--benchmark-dir", type=Path, required=True)
    parser.add_argument("--obo-file", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--config", type=Path)
    parser.add_argument(
        "--prepared-data",
        action="append",
        default=[],
        metavar="LABEL=DIR",
        help="Repeatable prepared-data reference to verify against the CSVs.",
    )
    parser.add_argument(
        "--ia-file-dir",
        type=Path,
        help="Optional directory containing BPO_ia.txt, CCO_ia.txt and MFO_ia.txt.",
    )
    parser.add_argument(
        "--metadata", action="append", default=[], metavar="KEY=VALUE"
    )
    parser.add_argument(
        "--source-evidence", action="append", default=[], type=Path
    )
    args = parser.parse_args()
    started = time.perf_counter()

    benchmark_dir = args.benchmark_dir.resolve()
    obo_file = args.obo_file.resolve()
    output_dir = args.output_dir.resolve()
    ia_file_dir = args.ia_file_dir.resolve() if args.ia_file_dir else None
    if output_dir.exists():
        raise ValueError(f"Output directory already exists: {output_dir}")
    missing = [name for name in required_csv_names() if not (benchmark_dir / name).is_file()]
    if missing:
        raise FileNotFoundError(f"Benchmark is missing required CSVs: {missing}")
    extras = sorted(
        path.name
        for path in benchmark_dir.glob("*.csv")
        if path.name not in required_csv_names()
    )
    if extras:
        raise ValueError(f"Benchmark directory contains unexpected CSVs: {extras}")

    config_path = args.config.resolve() if args.config else None
    config_sha256 = sha256_file(config_path) if config_path else None
    policy, config = load_policy(config_path)
    if config_path and sha256_file(config_path) != config_sha256:
        raise ValueError(f"Run config changed while it was being loaded: {config_path}")
    obo_bytes = obo_file.stat().st_size
    obo_sha256 = sha256_file(obo_file)
    graph = read_obo(obo_file)
    metadata_pairs = [
        parse_assignment(value, "--metadata") for value in args.metadata
    ]
    if len({key for key, _ in metadata_pairs}) != len(metadata_pairs):
        raise ValueError("Each --metadata key must be unique")
    metadata = dict(metadata_pairs)
    prepared = [
        (label, Path(raw).expanduser().resolve())
        for label, raw in (
            parse_assignment(value, "--prepared-data") for value in args.prepared_data
        )
    ]
    if len({label for label, _ in prepared}) != len(prepared):
        raise ValueError("Each --prepared-data label must be unique")
    evidence = []
    for path in args.source_evidence:
        resolved = path.resolve()
        if not resolved.is_file():
            raise FileNotFoundError(f"Source-evidence file is missing: {resolved}")
        evidence.append(
            {"path": str(resolved), "bytes": resolved.stat().st_size, "sha256": sha256_file(resolved)}
        )

    csv_reports: dict[str, dict[str, Any]] = {}
    root_only_rows: list[dict[str, Any]] = []
    term_rows: list[dict[str, Any]] = []
    histogram_rows: list[dict[str, Any]] = []
    aliases = []
    ia_snapshots: dict[str, dict[str, Any]] = {}
    for aspect in ASPECTS:
        ia_path = ia_file_dir / f"{aspect}_ia.txt" if ia_file_dir else None
        if ia_path is not None and not ia_path.is_file():
            raise FileNotFoundError(f"IA file is missing: {ia_path}")
        if ia_path is not None:
            ia_snapshots[aspect] = file_snapshot(ia_path)
        ia_values = read_ia_file(ia_path)
        aspect_terms = None
        for split in CSV_SPLITS:
            path = benchmark_dir / f"{ASPECT_TO_PREFIX[aspect]}-{split}.csv"
            item, roots, supports, histogram = audit_csv(
                path=path,
                aspect=aspect,
                split=split,
                graph=graph,
                ia_values=ia_values,
                allow_singular_header=policy["allow_legacy_singular_protein_header"],
                allow_all_zero_rows=policy["allow_all_zero_rows"],
            )
            if aspect_terms is None:
                aspect_terms = item["ordered_terms"]
            elif item["ordered_terms"] != aspect_terms:
                raise ValueError(
                    f"{aspect} GO columns differ between training, validation and test"
                )
            csv_reports[f"{aspect}:{split}"] = item
            root_only_rows.extend(roots)
            term_rows.extend(supports)
            aliases.extend([item["header_alias"]] if item["header_alias"] else [])
            for depth, count in sorted(histogram.items()):
                histogram_rows.append(
                    {"aspect": aspect, "split": split, "depth": depth, "proteins": count}
                )

    prepared_results = [
        verify_prepared_data(path, label, csv_reports) for label, path in prepared
    ]
    if obo_file.stat().st_size != obo_bytes or sha256_file(obo_file) != obo_sha256:
        raise ValueError(f"GO OBO changed while the benchmark was being audited: {obo_file}")
    for aspect, snapshot in ia_snapshots.items():
        require_unchanged(Path(snapshot["path"]), snapshot, f"{aspect} IA file")
    input_manifest = {
        "schema_version": SCHEMA_VERSION,
        "benchmark_id": args.benchmark_id,
        "benchmark_dir": str(benchmark_dir),
        "obo": {
            "path": str(obo_file),
            "bytes": obo_bytes,
            "sha256": obo_sha256,
        },
        "csvs": {
            item["file"]: {"bytes": item["bytes"], "sha256": item["sha256"]}
            for item in csv_reports.values()
        },
        "config": (
            {
                "path": str(args.config.resolve()),
                "sha256": config_sha256,
                "name": config.get("name"),
            }
            if args.config
            else None
        ),
        "policy": policy,
        "metadata": metadata,
        "source_evidence": evidence,
        "ia_files": ia_snapshots,
        "prepared_data": [
            {"label": label, "path": str(path)} for label, path in prepared
        ],
    }
    benchmark_fingerprint = sha256_json(
        {
            "csvs": input_manifest["csvs"],
            "obo": input_manifest["obo"]["sha256"],
            "ia_files": {
                aspect: value["sha256"] for aspect, value in ia_snapshots.items()
            },
            "policy": policy,
        }
    )
    report = {
        "schema_version": SCHEMA_VERSION,
        "status": "passed",
        "benchmark_id": args.benchmark_id,
        "benchmark_fingerprint": benchmark_fingerprint,
        "metadata": metadata,
        "inputs": input_manifest,
        "roots": ASPECT_TO_ROOT,
        "header_compatibility_aliases": aliases,
        "files": csv_reports,
        "prepared_data_verification": prepared_results,
        "resource_usage": {
            "wall_seconds": time.perf_counter() - started,
            "peak_rss_bytes": peak_rss_bytes(),
        },
        "root_only_policy": {
            "definition": "A labelled row with the ontology root positive and zero positive non-root terms.",
            "is_failure_criterion": False,
        },
    }

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(
        tempfile.mkdtemp(prefix=f".{output_dir.name}.stage-", dir=str(output_dir.parent))
    )
    try:
        atomic_write_json(stage / "input_manifest.json", input_manifest)
        atomic_write_json(stage / "label_space_audit.json", report)
        atomic_write_text(stage / "label_space_audit.md", markdown_report(report))
        atomic_write_text(
            stage / "root_only_targets.tsv",
            tsv_text(
                ["aspect", "split", "protein_id", "sequence_sha256"],
                sorted(root_only_rows, key=lambda row: (row["aspect"], row["split"], row["protein_id"])),
            ),
        )
        atomic_write_text(
            stage / "term_support.tsv",
            tsv_text(
                [
                    "aspect",
                    "split",
                    "term",
                    "is_root",
                    "support",
                    "shortest_depth",
                    "longest_depth",
                    "ia",
                ],
                term_rows,
            ),
        )
        atomic_write_text(
            stage / "label_depth_histogram.tsv",
            tsv_text(["aspect", "split", "depth", "proteins"], histogram_rows),
        )
        manifest = output_manifest(
            stage, exclude={"output_manifest.json", "RUN_COMPLETE.json"}
        )
        atomic_write_json(stage / "output_manifest.json", manifest)
        atomic_write_json(
            stage / "RUN_COMPLETE.json",
            {
                "schema_version": 1,
                "complete": True,
                "benchmark_id": args.benchmark_id,
                "benchmark_fingerprint": benchmark_fingerprint,
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
