#!/usr/bin/env python3
"""Audit the direct evidence-code composition of homology benchmarks."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import re


ASPECTS = ("bp", "cc", "mf")
SPLITS = ("training", "validation", "test")
ASPECT_MAP = {"P": "bp", "C": "cc", "F": "mf"}
EVIDENCE_GROUPS = {
    "experimental": frozenset({"EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "HTP", "HDA", "HMP", "HGI", "HEP"}),
    "author_statement": frozenset({"TAS", "NAS"}),
    "curator_inference": frozenset({"IC"}),
    "computational_or_phylogenetic": frozenset({"IGC", "RCA"}),
    "no_biological_data": frozenset({"ND"}),
}
EXPECTED_CODES = frozenset().union(*EVIDENCE_GROUPS.values())
SCHEMA_NAME = "homology-direct-evidence-policy-audit"
SCHEMA_VERSION = 1


def _utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _spec(raw: str) -> tuple[str, Path]:
    if "=" not in raw:
        raise argparse.ArgumentTypeError("benchmark must be LABEL=PATH")
    label, path = raw.split("=", 1)
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", label):
        raise argparse.ArgumentTypeError(f"Invalid label: {label!r}")
    return label, Path(path).expanduser().resolve()


def _group(code: str) -> str:
    matches = [name for name, codes in EVIDENCE_GROUPS.items() if code in codes]
    if len(matches) != 1:
        raise ValueError(f"Evidence code has no unique local audit group: {code}")
    return matches[0]


def _csv_population(root: Path, aspect: str, split_name: str) -> int:
    path = root / f"{aspect}-{split_name}.csv"
    with path.open("rb") as handle:
        count = sum(1 for _ in handle) - 1
    return count


def _audit_one(label: str, root: Path) -> dict[str, object]:
    annotation_path = root / "qualifying_annotations.tsv.gz"
    if not annotation_path.is_file():
        raise FileNotFoundError(annotation_path)
    row_counts: Counter[tuple[str, str, str, str]] = Counter()
    protein_sets: defaultdict[tuple[str, str, str, str], set[str]] = defaultdict(set)
    protein_masks: defaultdict[tuple[str, str, str], set[str]] = defaultdict(set)
    observed_codes: set[str] = set()
    with gzip.open(annotation_path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"uniprot_accession", "aspect", "evidence_code", "split"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"{annotation_path} lacks required columns")
        for row in reader:
            aspect = ASPECT_MAP.get(row["aspect"])
            split_name = row["split"]
            code = row["evidence_code"]
            protein = row["uniprot_accession"]
            if aspect is None or split_name not in SPLITS or not protein:
                raise ValueError(f"Invalid annotation row in {annotation_path}")
            group = _group(code)
            observed_codes.add(code)
            row_counts[aspect, split_name, code, group] += 1
            protein_sets[aspect, split_name, code, group].add(protein)
            protein_masks[aspect, split_name, protein].add(group)
    if observed_codes - EXPECTED_CODES:
        raise ValueError(
            "Observed evidence codes outside the expected 17-code contract: "
            f"{sorted(observed_codes - EXPECTED_CODES)}"
        )

    populations = {
        aspect: {split_name: _csv_population(root, aspect, split_name) for split_name in SPLITS}
        for aspect in ASPECTS
    }
    category_rows: list[dict[str, object]] = []
    code_rows: list[dict[str, object]] = []
    population_rows: list[dict[str, object]] = []
    for aspect in ASPECTS:
        for split_name in SPLITS:
            denominator = populations[aspect][split_name]
            masks = {
                protein: groups for (a, s, protein), groups in protein_masks.items()
                if a == aspect and s == split_name
            }
            experimental = {
                protein for protein, groups in masks.items() if "experimental" in groups
            }
            non_experimental_only = {
                protein for protein, groups in masks.items()
                if groups and "experimental" not in groups
            }
            population_rows.append({
                "aspect": aspect,
                "split": split_name,
                "benchmark_proteins": denominator,
                "directly_supported_proteins": len(masks),
                "experimental_proteins": len(experimental),
                "non_experimental_only_proteins": len(non_experimental_only),
            })
            for category in EVIDENCE_GROUPS:
                proteins = {
                    protein for (a, s, protein), groups in protein_masks.items()
                    if a == aspect and s == split_name and category in groups
                }
                exclusive = {
                    protein for (a, s, protein), groups in protein_masks.items()
                    if a == aspect and s == split_name and groups == {category}
                }
                rows = sum(
                    count for (a, s, _, group), count in row_counts.items()
                    if a == aspect and s == split_name and group == category
                )
                category_rows.append({
                    "aspect": aspect, "split": split_name, "category": category,
                    "annotation_rows": rows, "proteins": len(proteins),
                    "exclusive_proteins": len(exclusive), "benchmark_proteins": denominator,
                    "protein_fraction": len(proteins) / denominator if denominator else None,
                })
            for code in sorted(EXPECTED_CODES):
                category = _group(code)
                code_rows.append({
                    "aspect": aspect, "split": split_name, "evidence_code": code,
                    "category": category,
                    "annotation_rows": row_counts[aspect, split_name, code, category],
                    "proteins": len(protein_sets[aspect, split_name, code, category]),
                })
    return {
        "label": label,
        "path": str(root),
        "qualifying_annotations": {
            "path": str(annotation_path), "sha256": _sha256(annotation_path),
            "size_bytes": annotation_path.stat().st_size,
        },
        "observed_evidence_codes": sorted(observed_codes),
        "contract_codes_not_observed": sorted(EXPECTED_CODES - observed_codes),
        "benchmark_populations": populations,
        "population_rows": population_rows,
        "category_rows": category_rows,
        "evidence_code_rows": code_rows,
    }


def run(args: argparse.Namespace) -> dict[str, object]:
    specs = [_spec(raw) for raw in args.benchmark]
    if not specs:
        raise ValueError("At least one benchmark is required")
    output = args.output_dir.resolve()
    if output.exists():
        raise FileExistsError(output)
    audits = [_audit_one(label, root) for label, root in specs]
    payload = {
        "schema_name": SCHEMA_NAME, "schema_version": SCHEMA_VERSION,
        "created_at_utc": _utc_now(),
        "classification_boundary": {
            "description": "Local analytical grouping of the benchmark's exact 17-code contract; not a claim that all qualifying codes are experimental.",
            "groups": {name: sorted(codes) for name, codes in EVIDENCE_GROUPS.items()},
            "propagated_truth_caveat": "Counts describe direct qualifying annotation sources. Propagated GO labels cannot be uniquely assigned to one evidence row.",
        },
        "benchmarks": audits,
    }
    output.mkdir(parents=True)
    audit_path = output / "audit.json"
    audit_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    with (output / "evidence_code_counts.tsv").open("w", encoding="utf-8") as handle:
        handle.write("benchmark\taspect\tsplit\tevidence_code\tcategory\tannotation_rows\tproteins\n")
        for audit in audits:
            for row in audit["evidence_code_rows"]:
                handle.write(f"{audit['label']}\t{row['aspect']}\t{row['split']}\t{row['evidence_code']}\t{row['category']}\t{row['annotation_rows']}\t{row['proteins']}\n")
    with (output / "evidence_category_counts.tsv").open("w", encoding="utf-8") as handle:
        handle.write("benchmark\taspect\tsplit\tcategory\tannotation_rows\tproteins\texclusive_proteins\tbenchmark_proteins\tprotein_fraction\n")
        for audit in audits:
            for row in audit["category_rows"]:
                fraction = "" if row["protein_fraction"] is None else f"{row['protein_fraction']:.9f}"
                handle.write(f"{audit['label']}\t{row['aspect']}\t{row['split']}\t{row['category']}\t{row['annotation_rows']}\t{row['proteins']}\t{row['exclusive_proteins']}\t{row['benchmark_proteins']}\t{fraction}\n")
    lines = [
        "# Homology evidence-policy audit", "", "## Scope", "",
        "This read-only audit classifies the direct annotations admitted by the exact 17-code supervisor policy. `Experimental` means the eleven explicitly experimental/high-throughput codes listed below; the remaining six qualifying codes are reported separately rather than described as experimental.",
    ]
    for audit in audits:
        lines.extend(["", f"## {audit['label']}", "", "| Aspect | Split | Experimental proteins | Benchmark proteins | Experimental coverage | Non-experimental-only proteins |", "|---|---|---:|---:|---:|---:|"])
        lookup = {(row["aspect"], row["split"], row["category"]): row for row in audit["category_rows"]}
        population_lookup = {
            (row["aspect"], row["split"]): row for row in audit["population_rows"]
        }
        for aspect in ASPECTS:
            for split_name in SPLITS:
                exp = lookup[aspect, split_name, "experimental"]
                population = population_lookup[aspect, split_name]
                lines.append(
                    f"| {aspect.upper()} | {split_name} | {exp['proteins']:,} | "
                    f"{exp['benchmark_proteins']:,} | {100 * exp['protein_fraction']:.3f}% | "
                    f"{population['non_experimental_only_proteins']:,} |"
                )
    lines.extend(["", "## Interpretation boundary", "", "- Protein counts can overlap categories when a protein has annotations supported by more than one evidence class.", "- Exclusive counts identify proteins supported only by one class.", "- These counts audit label provenance; they do not re-evaluate model performance."])
    summary_path = output / "summary.md"
    summary_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    manifest_files = (audit_path, output / "evidence_code_counts.tsv", output / "evidence_category_counts.tsv", summary_path)
    manifest = {path.name: {"sha256": _sha256(path), "size_bytes": path.stat().st_size} for path in manifest_files}
    (output / "output_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (output / "RUN_COMPLETE.json").write_text(json.dumps({"complete": True, "schema_name": SCHEMA_NAME, "schema_version": SCHEMA_VERSION, "audit_sha256": _sha256(audit_path)}, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return payload


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark", action="append", required=True, metavar="LABEL=PATH")
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
