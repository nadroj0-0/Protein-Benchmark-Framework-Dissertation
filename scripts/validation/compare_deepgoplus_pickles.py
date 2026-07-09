#!/usr/bin/env python3
"""Compare regenerated DeepGOPlus CAFA pickles with released references."""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

from compare_cafa3_outputs import compare_pickles, markdown_table, write_tsv


REQUIRED_PICKLES = [
    "train_data.pkl",
    "test_data.pkl",
    "terms.pkl",
]


def status_from_pickle_rows(rows: list[dict[str, Any]]) -> str:
    required = {row["file"]: row for row in rows if row.get("file") in REQUIRED_PICKLES}
    if set(required) != set(REQUIRED_PICKLES):
        return "FAIL"
    for row in required.values():
        if not row.get("generated_exists") or not row.get("reference_exists"):
            return "FAIL"
        if row.get("status") != "compared":
            return "FAIL"
        if row.get("generated_shape") != row.get("reference_shape"):
            return "FAIL"
        if row.get("columns_equal") is False:
            return "FAIL"

    comparable_rates = []
    for name in ("train_data.pkl", "test_data.pkl"):
        row = required[name]
        comparable_rates.extend([
            row.get("protein_jaccard"),
            row.get("sequence_agreement_rate"),
            row.get("annotation_row_agreement_rate"),
            row.get("annotation_or_term_jaccard"),
        ])
    comparable_rates.append(required["terms.pkl"].get("annotation_or_term_jaccard"))
    rates = [rate for rate in comparable_rates if rate is not None]
    if rates and min(rates) >= 0.9999:
        return "PASS"
    if rates and min(rates) >= 0.95:
        return "PASS WITH MINOR DIFFERENCES"
    return "FAIL"


def write_report(
    path: Path,
    status: str,
    manifest_md: Path | None,
    pickle_rows: list[dict[str, Any]],
) -> None:
    lines = ["# CAFA3 DeepGOPlus Pickle Generation Validation Report", ""]
    lines.append(f"Final status: **{status}**")
    lines.append("")
    lines.append("This validates the historical DeepGOPlus `cafa3_data.py` layer:")
    lines.append("")
    lines.append("```text")
    lines.append("CAFA3 training FASTA/annotations + target FASTA/leaf-only ground truth + GO OBO")
    lines.append("    -> train_data.pkl")
    lines.append("    -> test_data.pkl")
    lines.append("    -> terms.pkl")
    lines.append("```")
    lines.append("")
    lines.append("Status thresholds:")
    lines.append("- PASS: required pickles exist, shapes/columns match, and protein/sequence/annotation/term agreement >= 0.9999.")
    lines.append("- PASS WITH MINOR DIFFERENCES: required pickles exist and minimum comparable agreement >= 0.95.")
    lines.append("- FAIL: missing required pickles, mismatched shapes/columns, or larger discrepancies.")
    lines.append("")
    if manifest_md and manifest_md.exists():
        lines.append("## Run Manifest")
        lines.append("")
        lines.append(manifest_md.read_text())
        lines.append("")
    lines.append("## Pickle Summary")
    lines.append("")
    lines.append(markdown_table(pickle_rows, [
        "file", "status", "generated_exists", "reference_exists", "generated_shape", "reference_shape",
        "columns_equal", "protein_jaccard", "sequence_agreement_rate",
        "annotation_row_agreement_rate", "annotation_or_term_jaccard",
    ]))
    lines.append("## Interpretation Notes")
    lines.append("")
    lines.append("- This comparator checks semantic equality of dataframe contents, not byte-for-byte pickle identity.")
    lines.append("- DeepGOPlus stores annotations as Python sets, so pickle byte order and some term ordering can depend on Python/hash/runtime details.")
    lines.append("- The downstream CSV validation remains the authoritative check for the PFP-facing benchmark interface.")
    lines.append("")
    path.write_text("\n".join(lines))


def main() -> None:
    parser = argparse.ArgumentParser(description="Compare regenerated DeepGOPlus CAFA pickles.")
    parser.add_argument("--generated-dir", required=True, type=Path)
    parser.add_argument("--reference-pickle-dir", required=True, type=Path)
    parser.add_argument("--reports-dir", required=True, type=Path)
    parser.add_argument("--manifest-md", type=Path)
    args = parser.parse_args()

    args.reports_dir.mkdir(parents=True, exist_ok=True)
    pickle_rows, protein_rows, term_rows = compare_pickles(args.generated_dir, args.reference_pickle_dir)
    pickle_rows = [row for row in pickle_rows if row.get("file") in REQUIRED_PICKLES]
    protein_rows = [row for row in protein_rows if row.get("file") in REQUIRED_PICKLES]
    term_rows = [row for row in term_rows if row.get("file") in REQUIRED_PICKLES]

    write_tsv(args.reports_dir / "pickle_generation_comparison.tsv", pickle_rows, [
        "file", "status", "generated_exists", "reference_exists", "generated_type", "reference_type",
        "generated_shape", "reference_shape", "generated_columns", "reference_columns",
        "column_overlap", "columns_equal", "generated_proteins", "reference_proteins",
        "protein_overlap", "protein_jaccard", "sequence_agreement_count", "sequence_agreement_rate",
        "annotation_row_agreement_count", "annotation_row_agreement_rate",
        "generated_annotation_or_term_count", "reference_annotation_or_term_count",
        "annotation_or_term_overlap", "annotation_or_term_jaccard",
    ])
    write_tsv(args.reports_dir / "pickle_generation_protein_overlap.tsv", protein_rows, [
        "source", "file", "generated_count", "reference_count", "overlap_count",
        "generated_only_count", "reference_only_count", "jaccard",
        "generated_only_preview", "reference_only_preview",
    ])
    write_tsv(args.reports_dir / "pickle_generation_go_term_overlap.tsv", term_rows, [
        "source", "file", "generated_count", "reference_count", "overlap_count",
        "generated_only_count", "reference_only_count", "jaccard",
        "generated_only_preview", "reference_only_preview",
    ])

    status = status_from_pickle_rows(pickle_rows)
    write_report(args.reports_dir / "cafa3_deepgoplus_pickle_generation_report.md", status, args.manifest_md, pickle_rows)
    print(f"CAFA3 DeepGOPlus pickle generation validation status: {status}")
    print(f"Reports written to: {args.reports_dir}")
    if status == "FAIL":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
