#!/usr/bin/env python3
"""Compare generated CAFA3 historical validation outputs with references."""
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


CSV_FILES = [
    "bp-training.csv",
    "bp-validation.csv",
    "bp-test.csv",
    "cc-training.csv",
    "cc-validation.csv",
    "cc-test.csv",
    "mf-training.csv",
    "mf-validation.csv",
    "mf-test.csv",
]

PICKLE_FILES = [
    "train_data.pkl",
    "test_data.pkl",
    "terms.pkl",
    "train_data_train.pkl",
    "train_data_valid.pkl",
]


def fmt(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        if math.isnan(value):
            return ""
        return f"{value:.6g}"
    return str(value)


def write_tsv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: fmt(row.get(field)) for field in fields})


def jaccard(left: set[str], right: set[str]) -> float:
    union = left | right
    if not union:
        return 1.0
    return len(left & right) / len(union)


def summary(values: pd.Series | np.ndarray | list[float]) -> dict[str, Any]:
    arr = np.asarray(values, dtype=float)
    if arr.size == 0:
        return {
            "min": None,
            "p25": None,
            "median": None,
            "mean": None,
            "p75": None,
            "max": None,
        }
    return {
        "min": float(np.min(arr)),
        "p25": float(np.percentile(arr, 25)),
        "median": float(np.median(arr)),
        "mean": float(np.mean(arr)),
        "p75": float(np.percentile(arr, 75)),
        "max": float(np.max(arr)),
    }


def locate_file(root: Path | None, name: str) -> Path | None:
    if root is None or not root.exists():
        return None
    direct = root / name
    if direct.exists():
        return direct
    matches = sorted(root.rglob(name))
    return matches[0] if matches else None


def load_csv(path: Path) -> tuple[pd.DataFrame, list[str], str]:
    df = pd.read_csv(path)
    protein_col = "proteins" if "proteins" in df.columns else "protein"
    if protein_col not in df.columns:
        raise ValueError(f"{path} has no proteins/protein column")
    if protein_col != "proteins":
        df = df.rename(columns={protein_col: "proteins"})
    terms = [col for col in df.columns if col.startswith("GO:")]
    for term in terms:
        df[term] = pd.to_numeric(df[term], errors="coerce").fillna(0).astype(int).clip(0, 1)
    return df, terms, protein_col


def grouped_labels(df: pd.DataFrame, terms: list[str]) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame(columns=terms)
    return df[["proteins", *terms]].groupby("proteins", sort=False)[terms].max()


def compare_csvs(generated_dir: Path, reference_dir: Path) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    csv_rows: list[dict[str, Any]] = []
    protein_rows: list[dict[str, Any]] = []
    term_rows: list[dict[str, Any]] = []

    for name in CSV_FILES:
        gen_path = locate_file(generated_dir, name)
        ref_path = locate_file(reference_dir, name)
        base = {"file": name, "generated_exists": bool(gen_path), "reference_exists": bool(ref_path)}

        if gen_path is None or ref_path is None:
            csv_rows.append(base)
            protein_rows.append({"source": "csv", "file": name, **base})
            term_rows.append({"source": "csv", "file": name, **base})
            continue

        gen_df, gen_terms, gen_protein_col = load_csv(gen_path)
        ref_df, ref_terms, ref_protein_col = load_csv(ref_path)

        gen_proteins = set(gen_df["proteins"].astype(str))
        ref_proteins = set(ref_df["proteins"].astype(str))
        gen_term_set = set(gen_terms)
        ref_term_set = set(ref_terms)
        common_terms = sorted(gen_term_set & ref_term_set)
        common_proteins = sorted(gen_proteins & ref_proteins)

        gen_matrix = grouped_labels(gen_df, common_terms)
        ref_matrix = grouped_labels(ref_df, common_terms)
        gen_common = gen_matrix.reindex(common_proteins).fillna(0).astype(int)
        ref_common = ref_matrix.reindex(common_proteins).fillna(0).astype(int)

        if common_terms and common_proteins:
            g = gen_common.to_numpy(dtype=int)
            r = ref_common.to_numpy(dtype=int)
            total_cells = int(g.size)
            matches = int((g == r).sum())
            tp = int(((g == 1) & (r == 1)).sum())
            fp = int(((g == 1) & (r == 0)).sum())
            fn = int(((g == 0) & (r == 1)).sum())
            exact_rows = int((g == r).all(axis=1).sum())
        else:
            total_cells = matches = tp = fp = fn = exact_rows = 0

        precision = tp / (tp + fp) if (tp + fp) else None
        recall = tp / (tp + fn) if (tp + fn) else None
        f1 = 2 * precision * recall / (precision + recall) if precision is not None and recall is not None and (precision + recall) else None
        matrix_agreement = matches / total_cells if total_cells else None

        gen_grouped = grouped_labels(gen_df, gen_terms)
        ref_grouped = grouped_labels(ref_df, ref_terms)
        gen_lpp = gen_grouped.sum(axis=1) if gen_terms else pd.Series(dtype=float)
        ref_lpp = ref_grouped.sum(axis=1) if ref_terms else pd.Series(dtype=float)
        gen_ppt = gen_grouped.sum(axis=0) if gen_terms else pd.Series(dtype=float)
        ref_ppt = ref_grouped.sum(axis=0) if ref_terms else pd.Series(dtype=float)

        gen_density = float(gen_grouped.to_numpy().sum() / gen_grouped.size) if gen_grouped.size else 0.0
        ref_density = float(ref_grouped.to_numpy().sum() / ref_grouped.size) if ref_grouped.size else 0.0

        lpp_gen = summary(gen_lpp)
        lpp_ref = summary(ref_lpp)
        ppt_gen = summary(gen_ppt)
        ppt_ref = summary(ref_ppt)

        csv_rows.append({
            **base,
            "generated_rows": len(gen_df),
            "reference_rows": len(ref_df),
            "row_diff": len(gen_df) - len(ref_df),
            "generated_columns": len(gen_df.columns),
            "reference_columns": len(ref_df.columns),
            "generated_protein_column": gen_protein_col,
            "reference_protein_column": ref_protein_col,
            "generated_proteins": len(gen_proteins),
            "reference_proteins": len(ref_proteins),
            "protein_overlap": len(gen_proteins & ref_proteins),
            "protein_jaccard": jaccard(gen_proteins, ref_proteins),
            "generated_go_terms": len(gen_term_set),
            "reference_go_terms": len(ref_term_set),
            "go_term_overlap": len(gen_term_set & ref_term_set),
            "go_term_jaccard": jaccard(gen_term_set, ref_term_set),
            "generated_label_density": gen_density,
            "reference_label_density": ref_density,
            "matrix_cells_compared": total_cells,
            "matrix_agreement": matrix_agreement,
            "exact_row_agreement_count": exact_rows,
            "exact_row_agreement_rate": exact_rows / len(common_proteins) if common_proteins else None,
            "true_positives": tp,
            "false_positives": fp,
            "false_negatives": fn,
            "precision": precision,
            "recall": recall,
            "f1": f1,
            "generated_labels_per_protein_mean": lpp_gen["mean"],
            "reference_labels_per_protein_mean": lpp_ref["mean"],
            "generated_proteins_per_term_mean": ppt_gen["mean"],
            "reference_proteins_per_term_mean": ppt_ref["mean"],
            "generated_labels_per_protein_median": lpp_gen["median"],
            "reference_labels_per_protein_median": lpp_ref["median"],
            "generated_proteins_per_term_median": ppt_gen["median"],
            "reference_proteins_per_term_median": ppt_ref["median"],
        })
        protein_rows.append({
            "source": "csv",
            "file": name,
            "generated_count": len(gen_proteins),
            "reference_count": len(ref_proteins),
            "overlap_count": len(gen_proteins & ref_proteins),
            "generated_only_count": len(gen_proteins - ref_proteins),
            "reference_only_count": len(ref_proteins - gen_proteins),
            "jaccard": jaccard(gen_proteins, ref_proteins),
            "generated_only_preview": ",".join(sorted(gen_proteins - ref_proteins)[:20]),
            "reference_only_preview": ",".join(sorted(ref_proteins - gen_proteins)[:20]),
        })
        term_rows.append({
            "source": "csv",
            "file": name,
            "generated_count": len(gen_term_set),
            "reference_count": len(ref_term_set),
            "overlap_count": len(gen_term_set & ref_term_set),
            "generated_only_count": len(gen_term_set - ref_term_set),
            "reference_only_count": len(ref_term_set - gen_term_set),
            "jaccard": jaccard(gen_term_set, ref_term_set),
            "generated_only_preview": ",".join(sorted(gen_term_set - ref_term_set)[:20]),
            "reference_only_preview": ",".join(sorted(ref_term_set - gen_term_set)[:20]),
        })

    return csv_rows, protein_rows, term_rows


def normalise_annotation(value: Any) -> set[str]:
    if value is None:
        return set()
    if isinstance(value, float) and math.isnan(value):
        return set()
    if isinstance(value, (set, list, tuple, frozenset)):
        return {str(x) for x in value}
    text = str(value).strip()
    if not text or text == "nan":
        return set()
    for sep in ["|", ";", ","]:
        if sep in text:
            return {part.strip().strip("'\"{}[]()") for part in text.split(sep) if part.strip()}
    return {text.strip("'\"{}[]()")}


def pickle_terms(obj: Any) -> set[str]:
    if isinstance(obj, pd.DataFrame):
        if "terms" in obj.columns:
            return set(obj["terms"].astype(str))
        if "annotations" in obj.columns:
            out: set[str] = set()
            for value in obj["annotations"]:
                out |= normalise_annotation(value)
            return out
    if isinstance(obj, pd.Series):
        return set(obj.astype(str))
    if isinstance(obj, (set, list, tuple, frozenset)):
        return {str(x) for x in obj}
    return set()


def compare_pickles(generated_dir: Path, reference_dir: Path | None) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    pickle_rows: list[dict[str, Any]] = []
    protein_rows: list[dict[str, Any]] = []
    term_rows: list[dict[str, Any]] = []

    for name in PICKLE_FILES:
        gen_path = locate_file(generated_dir, name)
        ref_path = locate_file(reference_dir, name) if reference_dir else None
        base = {"file": name, "generated_exists": bool(gen_path), "reference_exists": bool(ref_path)}
        if gen_path is None or ref_path is None:
            pickle_rows.append({**base, "status": "skipped_missing_reference_or_generated"})
            continue

        gen_obj = pd.read_pickle(gen_path)
        ref_obj = pd.read_pickle(ref_path)
        row = {**base, "status": "compared"}
        row["generated_type"] = type(gen_obj).__name__
        row["reference_type"] = type(ref_obj).__name__

        if isinstance(gen_obj, pd.DataFrame) and isinstance(ref_obj, pd.DataFrame):
            gen_cols = list(gen_obj.columns)
            ref_cols = list(ref_obj.columns)
            row.update({
                "generated_shape": str(gen_obj.shape),
                "reference_shape": str(ref_obj.shape),
                "generated_columns": ",".join(gen_cols),
                "reference_columns": ",".join(ref_cols),
                "column_overlap": len(set(gen_cols) & set(ref_cols)),
                "columns_equal": gen_cols == ref_cols,
            })

            if "proteins" in gen_obj.columns and "proteins" in ref_obj.columns:
                gen_proteins = set(gen_obj["proteins"].astype(str))
                ref_proteins = set(ref_obj["proteins"].astype(str))
                common = sorted(gen_proteins & ref_proteins)
                row["generated_proteins"] = len(gen_proteins)
                row["reference_proteins"] = len(ref_proteins)
                row["protein_overlap"] = len(common)
                row["protein_jaccard"] = jaccard(gen_proteins, ref_proteins)
                protein_rows.append({
                    "source": "pickle",
                    "file": name,
                    "generated_count": len(gen_proteins),
                    "reference_count": len(ref_proteins),
                    "overlap_count": len(common),
                    "generated_only_count": len(gen_proteins - ref_proteins),
                    "reference_only_count": len(ref_proteins - gen_proteins),
                    "jaccard": jaccard(gen_proteins, ref_proteins),
                    "generated_only_preview": ",".join(sorted(gen_proteins - ref_proteins)[:20]),
                    "reference_only_preview": ",".join(sorted(ref_proteins - gen_proteins)[:20]),
                })

                gen_idx = gen_obj.drop_duplicates("proteins").set_index("proteins")
                ref_idx = ref_obj.drop_duplicates("proteins").set_index("proteins")
                if "sequences" in gen_idx.columns and "sequences" in ref_idx.columns and common:
                    seq_equal = sum(str(gen_idx.loc[p, "sequences"]) == str(ref_idx.loc[p, "sequences"]) for p in common)
                    row["sequence_agreement_count"] = seq_equal
                    row["sequence_agreement_rate"] = seq_equal / len(common)
                if "annotations" in gen_idx.columns and "annotations" in ref_idx.columns and common:
                    ann_equal = sum(
                        normalise_annotation(gen_idx.loc[p, "annotations"]) == normalise_annotation(ref_idx.loc[p, "annotations"])
                        for p in common
                    )
                    row["annotation_row_agreement_count"] = ann_equal
                    row["annotation_row_agreement_rate"] = ann_equal / len(common)

            gen_terms = pickle_terms(gen_obj)
            ref_terms = pickle_terms(ref_obj)
            row["generated_annotation_or_term_count"] = len(gen_terms)
            row["reference_annotation_or_term_count"] = len(ref_terms)
            row["annotation_or_term_overlap"] = len(gen_terms & ref_terms)
            row["annotation_or_term_jaccard"] = jaccard(gen_terms, ref_terms)
            term_rows.append({
                "source": "pickle",
                "file": name,
                "generated_count": len(gen_terms),
                "reference_count": len(ref_terms),
                "overlap_count": len(gen_terms & ref_terms),
                "generated_only_count": len(gen_terms - ref_terms),
                "reference_only_count": len(ref_terms - gen_terms),
                "jaccard": jaccard(gen_terms, ref_terms),
                "generated_only_preview": ",".join(sorted(gen_terms - ref_terms)[:20]),
                "reference_only_preview": ",".join(sorted(ref_terms - gen_terms)[:20]),
            })

        pickle_rows.append(row)
    return pickle_rows, protein_rows, term_rows


def status_from_csv_rows(rows: list[dict[str, Any]]) -> str:
    if any(not row.get("generated_exists") or not row.get("reference_exists") for row in rows):
        return "FAIL"
    f1s = [row.get("f1") for row in rows if row.get("f1") is not None]
    protein_js = [row.get("protein_jaccard") for row in rows if row.get("protein_jaccard") is not None]
    term_js = [row.get("go_term_jaccard") for row in rows if row.get("go_term_jaccard") is not None]
    if f1s and min(f1s) >= 0.9999 and min(protein_js or [1.0]) >= 0.9999 and min(term_js or [1.0]) >= 0.9999:
        return "PASS"
    if f1s and min(f1s) >= 0.95:
        return "PASS WITH MINOR DIFFERENCES"
    return "FAIL"


def markdown_table(rows: list[dict[str, Any]], fields: list[str]) -> str:
    if not rows:
        return "_No rows._\n"
    out = ["|" + "|".join(fields) + "|", "|" + "|".join(["---"] * len(fields)) + "|"]
    for row in rows:
        out.append("|" + "|".join(fmt(row.get(field)) for field in fields) + "|")
    return "\n".join(out) + "\n"


def write_report(path: Path, status: str, manifest_md: Path | None, csv_rows: list[dict[str, Any]], pickle_rows: list[dict[str, Any]]) -> None:
    lines = ["# CAFA3 Historical Validation Report", ""]
    lines.append(f"Final status: **{status}**")
    lines.append("")
    lines.append("Status thresholds:")
    lines.append("- PASS: every CSV exists, minimum label F1/protein Jaccard/GO-term Jaccard >= 0.9999.")
    lines.append("- PASS WITH MINOR DIFFERENCES: every CSV exists and minimum label F1 >= 0.95.")
    lines.append("- FAIL: missing required CSVs or larger discrepancies.")
    lines.append("")
    if manifest_md and manifest_md.exists():
        lines.append("## Run Manifest")
        lines.append("")
        lines.append(manifest_md.read_text())
        lines.append("")
    lines.append("## CSV Summary")
    lines.append("")
    lines.append(markdown_table(csv_rows, [
        "file", "generated_rows", "reference_rows", "row_diff",
        "protein_jaccard", "go_term_jaccard", "matrix_agreement", "precision", "recall", "f1",
    ]))
    lines.append("## Pickle Summary")
    lines.append("")
    lines.append(markdown_table(pickle_rows, [
        "file", "status", "generated_exists", "reference_exists", "generated_shape", "reference_shape",
        "protein_jaccard", "annotation_or_term_jaccard",
    ]))
    lines.append("## Known Policy Gaps To Consider")
    lines.append("")
    lines.append("- CAFA3 README includes TAS and IC while some public benchmark code keeps only EXP, IDA, IPI, IMP, IGI, IEP.")
    lines.append("- CGD Candida backfill removal may be undocumented outside official materials.")
    lines.append("- MF protein-binding-only removal for GO:0005515 may not be implemented by the public Python path.")
    lines.append("- GO ontology release choice and obsolete/alternate GO handling may affect propagated labels.")
    lines.append("- Reference CSVs may include preprocessing decisions made outside the public benchmark-construction code.")
    lines.append("")
    path.write_text("\n".join(lines))


def main() -> None:
    parser = argparse.ArgumentParser(description="Compare CAFA3 historical validation outputs.")
    parser.add_argument("--generated-dir", required=True, type=Path)
    parser.add_argument("--reference-csv-dir", required=True, type=Path)
    parser.add_argument("--reference-pickle-dir", type=Path)
    parser.add_argument("--reports-dir", required=True, type=Path)
    parser.add_argument("--manifest-md", type=Path)
    args = parser.parse_args()

    args.reports_dir.mkdir(parents=True, exist_ok=True)
    csv_rows, protein_rows, term_rows = compare_csvs(args.generated_dir, args.reference_csv_dir)
    pickle_rows, pickle_protein_rows, pickle_term_rows = compare_pickles(args.generated_dir, args.reference_pickle_dir)
    protein_rows.extend(pickle_protein_rows)
    term_rows.extend(pickle_term_rows)

    write_tsv(args.reports_dir / "csv_comparison.tsv", csv_rows, [
        "file", "generated_exists", "reference_exists", "generated_rows", "reference_rows", "row_diff",
        "generated_columns", "reference_columns", "generated_proteins", "reference_proteins",
        "protein_overlap", "protein_jaccard", "generated_go_terms", "reference_go_terms",
        "go_term_overlap", "go_term_jaccard", "generated_label_density", "reference_label_density",
        "matrix_cells_compared", "matrix_agreement", "exact_row_agreement_count", "exact_row_agreement_rate",
        "true_positives", "false_positives", "false_negatives", "precision", "recall", "f1",
        "generated_labels_per_protein_mean", "reference_labels_per_protein_mean",
        "generated_proteins_per_term_mean", "reference_proteins_per_term_mean",
        "generated_labels_per_protein_median", "reference_labels_per_protein_median",
        "generated_proteins_per_term_median", "reference_proteins_per_term_median",
    ])
    write_tsv(args.reports_dir / "pickle_comparison.tsv", pickle_rows, [
        "file", "status", "generated_exists", "reference_exists", "generated_type", "reference_type",
        "generated_shape", "reference_shape", "generated_columns", "reference_columns",
        "column_overlap", "columns_equal", "generated_proteins", "reference_proteins",
        "protein_overlap", "protein_jaccard", "sequence_agreement_count", "sequence_agreement_rate",
        "annotation_row_agreement_count", "annotation_row_agreement_rate",
        "generated_annotation_or_term_count", "reference_annotation_or_term_count",
        "annotation_or_term_overlap", "annotation_or_term_jaccard",
    ])
    write_tsv(args.reports_dir / "protein_overlap.tsv", protein_rows, [
        "source", "file", "generated_count", "reference_count", "overlap_count",
        "generated_only_count", "reference_only_count", "jaccard",
        "generated_only_preview", "reference_only_preview",
    ])
    write_tsv(args.reports_dir / "go_term_overlap.tsv", term_rows, [
        "source", "file", "generated_count", "reference_count", "overlap_count",
        "generated_only_count", "reference_only_count", "jaccard",
        "generated_only_preview", "reference_only_preview",
    ])
    status = status_from_csv_rows(csv_rows)
    write_report(args.reports_dir / "cafa3_historical_validation_report.md", status, args.manifest_md, csv_rows, pickle_rows)
    print(f"CAFA3 historical validation status: {status}")
    print(f"Reports written to: {args.reports_dir}")


if __name__ == "__main__":
    main()
