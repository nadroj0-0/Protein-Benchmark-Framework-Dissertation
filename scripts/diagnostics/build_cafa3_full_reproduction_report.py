#!/usr/bin/env python3
"""Build compact JSON and Markdown reports for the full CAFA3 reproduction."""

from __future__ import annotations

import argparse
import csv
import json
import platform
import subprocess
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any


ASPECTS = ("BPO", "CCO", "MFO")
EXPECTED = {
    "BPO": {"fmax": 0.601, "wfmax": 0.515},
    "CCO": {"fmax": 0.706, "wfmax": 0.566},
    "MFO": {"fmax": 0.702, "wfmax": 0.605},
}
PACKAGES = (
    "torch",
    "numpy",
    "pandas",
    "scipy",
    "tqdm",
    "scikit-learn",
    "cafaeval",
    "obonet",
    "networkx",
    "transformers",
    "sentencepiece",
    "biopython",
    "h5py",
    "requests",
    "biotite",
    "fair-esm",
    "torch-geometric",
    "torch-scatter",
    "torch-sparse",
)


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def git_commit(path: Path) -> str:
    return subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=path, text=True
    ).strip()


def package_versions() -> dict[str, str]:
    observed = {}
    for package in PACKAGES:
        try:
            observed[package] = version(package)
        except PackageNotFoundError:
            observed[package] = "missing"
    return observed


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def training_results(pfp_root: Path) -> dict[str, Any]:
    results = {}
    for aspect in ASPECTS:
        path = (
            pfp_root
            / "results/full_model/fusion_comparison/prott5"
            / aspect
            / "gated_bilinear/results.json"
        )
        if not path.is_file():
            raise ValueError(f"Missing training result: {path}")
        results[aspect] = read_json(path)
    return results


def fmt(value: Any, digits: int = 6) -> str:
    if value is None:
        return "n/a"
    if isinstance(value, float):
        return f"{value:.{digits}f}"
    return str(value)


def build_markdown(report: dict[str, Any]) -> str:
    lines = [
        "# Full CAFA3 From-Scratch Reproduction Report",
        "",
        "## Outcome",
        "",
        "The workflow regenerated all four embedding modalities from the canonical "
        "Zenodo 7409660 CSVs, compared them with Zijian's authenticated published "
        "embedding cache, trained fresh PFP models, and evaluated those checkpoints.",
        "",
        f"- Workflow complete: `{str(report['complete']).lower()}`",
        f"- PFP commit: `{report['provenance']['pfp_commit']}`",
        f"- Framework commit: `{report['provenance']['framework_commit']}`",
        f"- Text temporal cutoff: `{report['provenance']['text_cutoff_date']}`",
        f"- Published cache discarded after comparison: "
        f"`{str(report['published_cache_discarded']).lower()}`",
        "",
        "## Embedding Comparison",
        "",
        "| Modality | Generated | Published | Common | Byte-exact | Numeric match | Different | Missing generated | Missing published |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for modality, summary in report["embedding_comparison"].items():
        statuses = summary.get("statuses", {})
        lines.append(
            "| {modality} | {generated} | {published} | {common} | {exact} | "
            "{numeric} | {different} | {missing_generated} | {missing_published} |".format(
                modality=modality,
                generated=summary.get("generated_count", "n/a"),
                published=summary.get("published_count", "n/a"),
                common=summary.get("common_count", "n/a"),
                exact=statuses.get("exact_match", 0),
                numeric=statuses.get("numeric_match", 0),
                different=statuses.get("different", 0),
                missing_generated=statuses.get("missing_generated", 0),
                missing_published=statuses.get("missing_published", 0),
            )
        )

    lines.extend(
        [
            "",
            "Array differences are observational results, not workflow failures. "
            "The row-level compressed comparison contains shape, dtype, SHA-256, "
            "cosine, L2, maximum-absolute and mean-absolute comparisons.",
            "",
            "## Evaluation Comparison",
            "",
            "| Ontology | Published Fmax | Reproduced Fmax | Delta | Published wFmax | Reproduced wFmax | Delta | Within 0.0005 |",
            "|---|---:|---:|---:|---:|---:|---:|---|",
        ]
    )
    for row in report["evaluation"]["results"]:
        lines.append(
            f"| {row['aspect']} | {fmt(row['expected_fmax'], 3)} | "
            f"{fmt(row['actual_fmax'], 6)} | {fmt(row['delta_fmax'], 6)} | "
            f"{fmt(row['expected_wfmax'], 3)} | {fmt(row['actual_wfmax'], 6)} | "
            f"{fmt(row['delta_wfmax'], 6)} | {row['passed']} |"
        )

    lines.extend(
        [
            "",
            f"PFP evaluation exit status: `{report['evaluation_exit_status']}`. A "
            "non-zero status is retained as a reference-metric mismatch when the "
            "complete evaluation summaries exist; missing summaries remain fatal.",
            "",
            "## Training",
            "",
            "| Ontology | Best epoch | Total epochs | Best validation Fmax | Test Fmax | CAFA Fmax | CAFA wFmax |",
            "|---|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for aspect, result in report["training"].items():
        lines.append(
            f"| {aspect} | {result.get('best_epoch', 'n/a')} | "
            f"{result.get('total_epochs', 'n/a')} | "
            f"{fmt(result.get('best_val_fmax'))} | {fmt(result.get('test_fmax'))} | "
            f"{fmt(result.get('cafa_fmax'))} | {fmt(result.get('cafa_wfmax'))} |"
        )

    lines.extend(
        [
            "",
            "## Interpretation Boundaries",
            "",
            "- The benchmark CSVs and CAFA3 ontology are authenticated historical inputs.",
            "- STRING PPI uses the fixed v12 files downloaded by the PFP workflow.",
            "- ProtT5, PubMedBERT, UniProt text retrieval, AlphaFold retrieval, CUDA, "
            "and unpinned fair-esm/PyG dependencies can make regenerated bytes "
            "environment- or service-dependent.",
            "- PFP training seeds are set, but the repository does not request all "
            "PyTorch deterministic algorithms; cross-environment bit identity is not expected.",
            "- This report distinguishes input/output contract failures from legitimate "
            "scientific differences. It does not hide either class.",
            "",
            "## Artifacts",
            "",
            "- `embedding_comparison.csv.gz`: row-level embedding comparison",
            "- `embedding_comparison_summary.json`: modality summaries",
            "- `evaluation/reproduction_summary.csv` and `.json`: metric comparison",
            "- `training/*_results.json`: fresh training summaries",
            "- `logs/`: complete stage logs",
            "- `input_acquisition.tsv`: source and checksum trail",
            "- `modality_status.tsv`: preflight and full modality exit statuses",
            "",
        ]
    )
    return "\n".join(lines)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pfp-root", type=Path, required=True)
    parser.add_argument("--framework-root", type=Path, required=True)
    parser.add_argument("--embedding-summary", type=Path, required=True)
    parser.add_argument("--evaluation-summary", type=Path, required=True)
    parser.add_argument("--modality-status", type=Path, required=True)
    parser.add_argument("--input-acquisition", type=Path, required=True)
    parser.add_argument("--evaluation-exit-status", type=int, required=True)
    parser.add_argument("--text-cutoff-date", required=True)
    parser.add_argument("--published-cache-discarded", action="store_true")
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-md", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    evaluation = read_json(args.evaluation_summary)
    if set(row["aspect"] for row in evaluation.get("results", [])) != set(ASPECTS):
        raise ValueError("Evaluation summary does not contain all three ontologies")
    report = {
        "schema_version": 1,
        "complete": True,
        "completed_at": datetime.now(timezone.utc).isoformat(),
        "provenance": {
            "pfp_commit": git_commit(args.pfp_root),
            "framework_commit": git_commit(args.framework_root),
            "text_cutoff_date": args.text_cutoff_date,
        },
        "environment": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "packages": package_versions(),
        },
        "embedding_comparison": read_json(args.embedding_summary),
        "evaluation": evaluation,
        "evaluation_exit_status": args.evaluation_exit_status,
        "training": training_results(args.pfp_root),
        "modality_status": read_tsv(args.modality_status),
        "input_acquisition": read_tsv(args.input_acquisition),
        "published_reference": EXPECTED,
        "published_cache_discarded": args.published_cache_discarded,
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_md.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    args.output_md.write_text(build_markdown(report), encoding="utf-8")
    print(args.output_md)
    print(args.output_json)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
