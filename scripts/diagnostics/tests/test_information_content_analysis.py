from __future__ import annotations

import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np


FRAMEWORK = Path(__file__).parents[3]
DIAGNOSTICS = FRAMEWORK / "scripts" / "diagnostics"
MODEL_EXECUTION = FRAMEWORK / "scripts" / "model_execution"
sys.path.insert(0, str(MODEL_EXECUTION))
sys.path.insert(0, str(DIAGNOSTICS))

import test_pfp_label_sensitivity as sensitivity_tests  # noqa: E402


ANALYSIS = DIAGNOSTICS / "evaluate_pfp_information_content.py"


def load_analysis():
    specification = importlib.util.spec_from_file_location(
        "test_evaluate_pfp_information_content", ANALYSIS
    )
    if specification is None or specification.loader is None:
        raise RuntimeError(f"Cannot import {ANALYSIS}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


class PfpInformationContentAnalysisTests(unittest.TestCase):
    def test_generic_binning_supports_lower_values_as_more_specific(self) -> None:
        from specificity_common import SpecificityMeasure, assign_specificity_bins

        terms = [f"GO:{index:07d}" for index in range(4)]
        measure = SpecificityMeasure(
            name="fixture_totipotency",
            values=np.asarray([1.0, 0.8, 0.2, 0.1]),
            higher_is_more_specific=False,
            zero_bin_label=None,
            source={},
        )
        _, assignments = assign_specificity_bins(
            terms, measure, 2, bin_prefix="specificity_q"
        )
        self.assertEqual(
            assignments,
            ["specificity_q1", "specificity_q1", "specificity_q2", "specificity_q2"],
        )

    def test_bins_are_deterministic_and_keep_zero_ia_separate(self) -> None:
        module = load_analysis()
        terms = [f"GO:{index:07d}" for index in range(6)]
        values = np.asarray([0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
        bins, assignments = module.assign_information_bins(terms, values, 2)
        self.assertEqual(assignments[0], "zero_ia")
        self.assertEqual(assignments[1:4], ["positive_q1"] * 3)
        self.assertEqual(assignments[4:], ["positive_q2"] * 2)
        self.assertEqual([value["term_count"] for value in bins], [1, 3, 2])

    def test_weighted_metrics_penalize_high_information_errors(self) -> None:
        module = load_analysis()
        truth = np.asarray([[1, 0], [0, 1]], dtype=np.uint8)
        scores = np.asarray([[0.9, 0.8], [0.1, 0.4]], dtype=np.float64)
        weights = np.asarray([1.0, 4.0], dtype=np.float64)
        metrics = module.threshold_metrics(truth, scores, weights, 0.5)
        self.assertAlmostEqual(metrics["micro_weighted_recall"], 0.2)
        self.assertLess(metrics["micro_weighted_f"], metrics["micro_f"])
        self.assertLess(
            metrics["weighted_jaccard_set_agreement"],
            metrics["jaccard_set_agreement"],
        )

    def test_xu_totipotency_chain_and_diamond(self) -> None:
        from specificity_common import compute_xu_totipotency, read_xu_ontology

        with tempfile.TemporaryDirectory() as name:
            path = Path(name) / "go.obo"
            path.write_text(
                "\n".join(
                    [
                        "format-version: 1.2",
                        "data-version: fixture",
                        "",
                        "[Term]",
                        "id: GO:0008150",
                        "namespace: biological_process",
                        "",
                        "[Term]",
                        "id: GO:0000001",
                        "namespace: biological_process",
                        "is_a: GO:0008150 ! root",
                        "",
                        "[Term]",
                        "id: GO:0000002",
                        "namespace: biological_process",
                        "is_a: GO:0000001 ! a",
                        "",
                        "[Term]",
                        "id: GO:0000003",
                        "namespace: biological_process",
                        "is_a: GO:0008150 ! root",
                        "",
                        "[Term]",
                        "id: GO:0000004",
                        "namespace: biological_process",
                        "is_a: GO:0000001 ! a",
                        "is_a: GO:0000003 ! b",
                        "",
                        "[Term]",
                        "id: GO:0005575",
                        "namespace: cellular_component",
                        "",
                        "[Term]",
                        "id: GO:0003674",
                        "namespace: molecular_function",
                        "",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            ontology = read_xu_ontology(path)
            terms = [
                "GO:0008150",
                "GO:0000001",
                "GO:0000002",
                "GO:0000003",
                "GO:0000004",
            ]
            raw, neglog, rows = compute_xu_totipotency(ontology, terms, "BPO")
            self.assertTrue(np.allclose(raw.values, [1.0, 0.6, 0.2, 0.4, 0.2]))
            self.assertTrue(np.allclose(neglog.values, -np.log2(raw.values)))
            self.assertEqual(rows[4]["descendant_count"], 1)
            self.assertEqual(rows[0]["aspect_root_descendant_count"], 5)

    def test_bin_without_positive_targets_is_not_reported_as_model_failure(
        self,
    ) -> None:
        module = load_analysis()
        truth = np.asarray([[1, 0], [1, 0]], dtype=np.uint8)
        scores = np.asarray([[0.9, 0.8], [0.8, 0.7]], dtype=np.float64)
        metrics, rows = module.evaluate_bin(
            truth,
            scores,
            np.asarray([1.0, 2.0], dtype=np.float64),
            [1],
            0.5,
        )
        self.assertEqual(metrics["status"], "not_evaluable_no_positive_targets")
        self.assertEqual(rows, [])

    def test_cli_publishes_auditable_bin_reports(self) -> None:
        helper = sensitivity_tests.PfpLabelSensitivityTests(methodName="runTest")
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            obo = root / "go.obo"
            helper.make_obo(obo)
            manifest = helper.make_prediction_artifact(root, obo)
            output = root / "information-content"
            result = subprocess.run(
                [
                    sys.executable,
                    str(ANALYSIS),
                    "--prediction-manifest",
                    str(manifest),
                    "--output-dir",
                    str(output),
                    "--positive-bins",
                    "2",
                    "--bootstrap-replicates",
                    "10",
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            report = json.loads(
                (output / "specificity_analysis.json").read_text(encoding="utf-8")
            )
            self.assertEqual(report["status"], "complete")
            bins = report["aspects"]["BPO"]["measures"]["ia"]["bins"]
            self.assertEqual(bins["zero_ia"]["term_count"], 0)
            self.assertEqual(bins["specificity_q1"]["term_count"], 1)
            self.assertTrue(
                bins["specificity_q1"]["metrics"]["weighted_metrics_available"]
            )
            for filename in (
                "specificity_analysis.md",
                "term_specificity.tsv",
                "specificity_bins.tsv",
                "specificity_thresholds.tsv",
                "bootstrap_intervals.tsv",
                "output_manifest.json",
                "RUN_COMPLETE.json",
            ):
                self.assertTrue((output / filename).is_file(), filename)


if __name__ == "__main__":
    unittest.main()
