from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import numpy as np


FRAMEWORK = Path(__file__).parents[3]
DIAGNOSTICS = FRAMEWORK / "scripts" / "diagnostics"
MODEL_EXECUTION = FRAMEWORK / "scripts" / "model_execution"
sys.path.insert(0, str(MODEL_EXECUTION))
sys.path.insert(0, str(DIAGNOSTICS))

from calibration_common import (  # noqa: E402
    CalibrationPolicy,
    FFPredPlattPolicy,
    apply_calibrator,
    apply_ffpred_platt_calibrator,
    ffpred_term_reliability,
    fit_ffpred_platt_calibrator,
    fit_monotone_hierarchical_calibrator,
)
import calibration_common  # noqa: E402
import test_pfp_label_sensitivity as sensitivity_tests  # noqa: E402


CALIBRATE = DIAGNOSTICS / "calibrate_pfp_predictions.py"
CONFIDENCE = DIAGNOSTICS / "generate_pfp_confidence.py"


class PfpCalibrationTests(unittest.TestCase):
    def test_positive_slope_and_fallback_are_deterministic(self) -> None:
        scores = np.asarray(
            [
                [0.05, 0.20, 0.10],
                [0.20, 0.40, 0.25],
                [0.70, 0.60, 0.80],
                [0.90, 0.85, 0.95],
            ],
            dtype=np.float64,
        )
        truth = np.asarray(
            [
                [0, 0, 0],
                [0, 0, 0],
                [1, 1, 1],
                [1, 1, 1],
            ],
            dtype=np.uint8,
        )
        terms = ["GO:0000001", "GO:0000002", "GO:0000003"]
        bins = ["low", "low", "high"]
        policy = CalibrationPolicy(
            minimum_bin_positives=1,
            minimum_bin_negatives=1,
            minimum_term_positives=1,
            minimum_term_negatives=1,
            maximum_iterations=100,
            protein_chunk_size=2,
        )
        model = fit_monotone_hierarchical_calibrator(scores, truth, terms, bins, policy)
        calibrated, fallback = apply_calibrator(scores, terms, bins, model)
        self.assertGreater(model["positive_slope"], 0)
        self.assertEqual(fallback, ["term_shrinkage"] * 3)
        self.assertTrue(np.all(np.diff(calibrated[:, 0]) >= 0))
        self.assertLess(calibrated[0, 0], calibrated[-1, 0])
        second = fit_monotone_hierarchical_calibrator(
            scores, truth, terms, bins, policy
        )
        self.assertEqual(model["model_sha256"], second["model_sha256"])

    def test_one_class_population_is_explicitly_uncalibrated(self) -> None:
        scores = np.asarray([[0.1], [0.9]], dtype=np.float64)
        truth = np.zeros((2, 1), dtype=np.uint8)
        policy = CalibrationPolicy(
            minimum_bin_positives=1,
            minimum_bin_negatives=1,
            minimum_term_positives=1,
            minimum_term_negatives=1,
        )
        model = fit_monotone_hierarchical_calibrator(
            scores, truth, ["GO:0000001"], ["low"], policy
        )
        calibrated, fallback = apply_calibrator(scores, ["GO:0000001"], ["low"], model)
        self.assertEqual(model["status"], "uncalibrated_insufficient_support")
        self.assertEqual(fallback, ["uncalibrated"])
        self.assertTrue(np.isnan(calibrated).all())

    def test_ffpred_style_is_per_term_monotone_and_reports_reliability(self) -> None:
        scores = np.asarray(
            [
                [0.05, 0.10],
                [0.20, 0.20],
                [0.70, 0.80],
                [0.95, 0.90],
            ],
            dtype=np.float64,
        )
        truth = np.asarray([[0, 0], [0, 0], [1, 1], [1, 1]], dtype=np.uint8)
        terms = ["GO:0000001", "GO:0000002"]
        policy = FFPredPlattPolicy(maximum_iterations=100, protein_chunk_size=2)
        model = fit_ffpred_platt_calibrator(scores, truth, terms, policy)
        posterior, statuses = apply_ffpred_platt_calibrator(
            scores, terms, model, protein_chunk_size=2
        )
        self.assertEqual(statuses, ["term_platt", "term_platt"])
        self.assertTrue(np.all(np.diff(posterior[:, 0]) >= 0))
        self.assertEqual(
            [row["reliability_level"] for row in ffpred_term_reliability(
                posterior, truth, terms, statuses
            )],
            ["high", "high"],
        )
        repeated = fit_ffpred_platt_calibrator(scores, truth, terms, policy)
        self.assertEqual(model["model_sha256"], repeated["model_sha256"])

    def test_ffpred_style_continues_after_iteration_limit(self) -> None:
        scores = np.asarray(
            [[0.05], [0.20], [0.70], [0.95]], dtype=np.float64
        )
        truth = np.asarray([[0], [0], [1], [1]], dtype=np.uint8)
        policy = FFPredPlattPolicy(
            maximum_iterations=100,
            maximum_restarts=1,
            protein_chunk_size=2,
        )
        real_minimize = calibration_common.minimize
        call_count = 0

        def iteration_limited_then_real(*args, **kwargs):
            nonlocal call_count
            call_count += 1
            if call_count == 1:
                initial = np.asarray(args[1], dtype=np.float64)
                return SimpleNamespace(
                    success=False,
                    status=1,
                    message="TOTAL NO. of ITERATIONS REACHED LIMIT",
                    nit=100,
                    nfev=101,
                    fun=1.0,
                    x=initial,
                )
            return real_minimize(*args, **kwargs)

        with mock.patch.object(
            calibration_common,
            "minimize",
            side_effect=iteration_limited_then_real,
        ):
            model = fit_ffpred_platt_calibrator(
                scores, truth, ["GO:0000001"], policy
            )
        self.assertEqual(model["status"], "complete")
        self.assertEqual(len(model["optimizer"]["attempts"]), 2)
        self.assertEqual(
            model["optimizer"]["parameterization"], "bounded_positive_slope"
        )

    def test_cli_fits_validation_only_and_publishes_no_p_values(self) -> None:
        helper = sensitivity_tests.PfpLabelSensitivityTests(methodName="runTest")
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            obo = root / "go.obo"
            helper.make_obo(obo)
            valid = helper.make_prediction_artifact(
                root,
                obo,
                evaluation_split="valid",
                protein_ids_override=["V1", "V2", "V3", "V4"],
                scores_override=np.asarray(
                    [[0.9, 0.1], [0.9, 0.3], [0.9, 0.7], [0.9, 0.9]]
                ),
                truth_override=np.asarray(
                    [[1, 0], [1, 0], [1, 1], [1, 1]], dtype=np.uint8
                ),
            )
            test = helper.make_prediction_artifact(
                root,
                obo,
                evaluation_split="test",
                protein_ids_override=["T1", "T2", "T3"],
                scores_override=np.asarray([[0.9, 0.2], [0.9, 0.6], [0.9, 0.8]]),
                truth_override=np.asarray([[1, 0], [1, 1], [1, 1]], dtype=np.uint8),
            )
            output = root / "calibration"
            result = subprocess.run(
                [
                    sys.executable,
                    str(CALIBRATE),
                    "--validation-prediction-manifest",
                    str(valid),
                    "--test-prediction-manifest",
                    str(test),
                    "--obo",
                    str(obo),
                    "--output-dir",
                    str(output),
                    "--positive-ia-bins",
                    "2",
                    "--reliability-bins",
                    "2",
                    "--minimum-bin-positives",
                    "1",
                    "--minimum-bin-negatives",
                    "1",
                    "--minimum-term-positives",
                    "1",
                    "--minimum-term-negatives",
                    "1",
                    "--maximum-iterations",
                    "100",
                    "--protein-chunk-size",
                    "2",
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            model = json.loads(
                (output / "calibration_model.json").read_text(encoding="utf-8")
            )
            self.assertEqual(
                model["analysis_label"], "post_selection_validation_calibration"
            )
            self.assertIn("not a p-value", model["probability_interpretation"])
            self.assertEqual(model["policy"]["p_values"], "prohibited")
            with np.load(
                output / "BPO_calibration_predictions.npz", allow_pickle=False
            ) as archive:
                self.assertEqual(
                    set(archive.files),
                    {
                        "raw_scores",
                        "postprop_scores",
                        "calibrated_q",
                        "truth",
                        "protein_ids",
                        "go_terms",
                        "information_accretion",
                        "ia_bins",
                        "fallback_level",
                    },
                )
                self.assertEqual(archive["calibrated_q"].shape, (3, 1))
            for filename in (
                "calibration_reliability.tsv",
                "calibration_shift.tsv",
                "calibration_hierarchy_audit.tsv",
                "calibration_analysis.json",
                "output_manifest.json",
                "RUN_COMPLETE.json",
            ):
                self.assertTrue((output / filename).is_file(), filename)

    def test_confidence_cli_publishes_three_explicit_methods(self) -> None:
        helper = sensitivity_tests.PfpLabelSensitivityTests(methodName="runTest")
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            obo = root / "go.obo"
            helper.make_obo(obo)
            valid = helper.make_prediction_artifact(
                root,
                obo,
                evaluation_split="valid",
                protein_ids_override=["V1", "V2", "V3", "V4"],
                scores_override=np.asarray(
                    [[0.9, 0.1], [0.9, 0.3], [0.9, 0.7], [0.9, 0.9]]
                ),
                truth_override=np.asarray(
                    [[1, 0], [1, 0], [1, 1], [1, 1]], dtype=np.uint8
                ),
            )
            test = helper.make_prediction_artifact(
                root,
                obo,
                evaluation_split="test",
                protein_ids_override=["T1", "T2", "T3"],
                scores_override=np.asarray([[0.9, 0.2], [0.9, 0.6], [0.9, 0.8]]),
                truth_override=np.asarray([[1, 0], [1, 1], [1, 1]], dtype=np.uint8),
            )
            output = root / "confidence"
            result = subprocess.run(
                [
                    sys.executable,
                    str(CONFIDENCE),
                    "--validation-prediction-manifest",
                    str(valid),
                    "--test-prediction-manifest",
                    str(test),
                    "--obo",
                    str(obo),
                    "--output-dir",
                    str(output),
                    "--positive-ia-bins",
                    "2",
                    "--reliability-bins",
                    "2",
                    "--minimum-bin-positives",
                    "1",
                    "--minimum-bin-negatives",
                    "1",
                    "--minimum-term-positives",
                    "1",
                    "--minimum-term-negatives",
                    "1",
                    "--maximum-iterations",
                    "100",
                    "--protein-chunk-size",
                    "2",
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            model = json.loads(
                (output / "confidence_model.json").read_text(encoding="utf-8")
            )
            self.assertEqual(
                model["selected_methods"],
                ["raw-sigmoid", "ffpred-style", "hierarchy-calibrated"],
            )
            self.assertEqual(model["policy"]["p_values"], "prohibited")
            with np.load(
                output / "BPO_confidence_outputs.npz", allow_pickle=False
            ) as archive:
                self.assertIn("raw_sigmoid_score", archive.files)
                self.assertIn("ffpred_style_posterior", archive.files)
                self.assertIn("hierarchy_calibrated_probability", archive.files)
                self.assertEqual(
                    archive["hierarchy_calibrated_probability"].shape, (3, 1)
                )
            audit = (output / "confidence_hierarchy_audit.tsv").read_text(
                encoding="utf-8"
            )
            projected = [
                line.split("\t")
                for line in audit.splitlines()[1:]
                if "hierarchy-calibrated-projected" in line
            ]
            self.assertTrue(projected)
            header = audit.splitlines()[0].split("\t")
            violation_index = header.index("violation_events")
            self.assertTrue(all(int(row[violation_index]) == 0 for row in projected))


if __name__ == "__main__":
    unittest.main()
