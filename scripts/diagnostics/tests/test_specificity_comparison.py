from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


FRAMEWORK = Path(__file__).parents[3]
SCRIPT = FRAMEWORK / "scripts" / "diagnostics" / "compare_pfp_specificity_runs.py"
sys.path.insert(0, str(FRAMEWORK / "scripts" / "diagnostics"))

from label_space_common import output_manifest, sha256_file  # noqa: E402


BIN_FIELDS = (
    "benchmark_id",
    "mode",
    "aspect",
    "specificity_measure",
    "specificity_bin",
    "term_count",
    "target_positive_proteins",
    "value_min",
    "value_max",
    "best_unweighted_f_threshold",
    "best_unweighted_f",
    "best_weighted_f_threshold",
    "best_weighted_f",
    "fixed_threshold_type",
    "fixed_threshold",
    "fixed_weighted_f",
    "fixed_unweighted_f",
    "fixed_weighted_jaccard",
    "fixed_unweighted_jaccard",
)


def write_tsv(
    path: Path, fields: tuple[str, ...], rows: list[dict[str, object]]
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_source(root: Path, benchmark_id: str, offset: float) -> None:
    root.mkdir()
    bins = []
    bootstrap = []
    bootstrap_fields = (
        "benchmark_id",
        "mode",
        "aspect",
        "specificity_measure",
        "specificity_bin",
        "threshold_type",
        "threshold",
        "metric_name",
        "estimate",
        "ci_low",
        "ci_high",
        "bootstrap_replicates",
        "bootstrap_seed",
        "resampling_unit",
        "method",
    )
    for aspect in ("BPO", "CCO", "MFO"):
        for measure in ("ia", "xu_totipotency_raw"):
            for number in range(1, 5):
                value = 0.8 - number * 0.1 + offset
                row = {
                    "benchmark_id": benchmark_id,
                    "mode": "full",
                    "aspect": aspect,
                    "specificity_measure": measure,
                    "specificity_bin": f"specificity_q{number}",
                    "term_count": 10,
                    "target_positive_proteins": 20,
                    "value_min": 0.1,
                    "value_max": 0.2,
                    "best_unweighted_f_threshold": 0.5,
                    "best_unweighted_f": value,
                    "best_weighted_f_threshold": 0.5,
                    "best_weighted_f": value,
                    "fixed_threshold_type": "descriptive_test_oracle_fixed",
                    "fixed_threshold": 0.5,
                    "fixed_weighted_f": value,
                    "fixed_unweighted_f": value,
                    "fixed_weighted_jaccard": value / 2,
                    "fixed_unweighted_jaccard": value / 2,
                }
                bins.append(row)
                for metric in (
                    "f",
                    "weighted_f",
                    "jaccard_set_agreement",
                    "weighted_jaccard_set_agreement",
                ):
                    bootstrap.append(
                        {
                            "benchmark_id": benchmark_id,
                            "mode": "full",
                            "aspect": aspect,
                            "specificity_measure": measure,
                            "specificity_bin": f"specificity_q{number}",
                            "threshold_type": "descriptive_test_oracle_fixed",
                            "threshold": 0.5,
                            "metric_name": metric,
                            "estimate": value,
                            "ci_low": value - 0.01,
                            "ci_high": value + 0.01,
                            "bootstrap_replicates": 20,
                            "bootstrap_seed": 42,
                            "resampling_unit": "protein",
                            "method": "percentile",
                        }
                    )
    write_tsv(root / "specificity_bins.tsv", BIN_FIELDS, bins)
    write_tsv(root / "bootstrap_intervals.tsv", bootstrap_fields, bootstrap)
    (root / "specificity_analysis.json").write_text(
        json.dumps(
            {
                "status": "complete",
                "specificity_measures": ["ia", "xu_totipotency_raw"],
            }
        )
    )
    (root / "specificity_analysis.md").write_text("fixture\n")
    manifest = output_manifest(
        root, exclude={"output_manifest.json", "RUN_COMPLETE.json"}
    )
    (root / "output_manifest.json").write_text(json.dumps(manifest))
    (root / "RUN_COMPLETE.json").write_text(
        json.dumps(
            {
                "complete": True,
                "output_manifest_sha256": sha256_file(root / "output_manifest.json"),
            }
        )
    )


class SpecificityComparisonTests(unittest.TestCase):
    def test_combines_three_complete_sources_and_reports_q1_q4_change(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            sources = []
            for number, label in enumerate(("cafa3", "global-nk", "nk-lk")):
                path = root / label
                write_source(path, label, number * 0.01)
                sources.extend(["--source", f"{label}={path}"])
            output = root / "output"
            subprocess.run(
                [sys.executable, str(SCRIPT), *sources, "--output-dir", str(output)],
                check=True,
                capture_output=True,
                text=True,
            )
            marker = json.loads((output / "RUN_COMPLETE.json").read_text())
            self.assertTrue(marker["complete"])
            self.assertEqual(marker["source_labels"], ["cafa3", "global-nk", "nk-lk"])
            with (output / "broad_to_specific_summary.tsv").open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 18)
            self.assertAlmostEqual(float(rows[0]["q4_minus_q1_fixed_f"]), -0.3)
            report = (output / "specificity_three_way_comparison.md").read_text()
            self.assertIn("benchmark-local rank groups", report)

    def test_duplicate_labels_fail(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "source"
            write_source(source, "fixture", 0.0)
            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--source",
                    f"same={source}",
                    "--source",
                    f"same={source}",
                    "--output-dir",
                    str(root / "output"),
                ],
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("labels must be unique", result.stderr)


if __name__ == "__main__":
    unittest.main()
