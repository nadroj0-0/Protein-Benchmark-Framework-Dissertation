from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
import types
import unittest
from pathlib import Path

import numpy as np


FRAMEWORK = Path(__file__).parents[3]
DIAGNOSTICS = FRAMEWORK / "scripts" / "diagnostics"
MODEL_EXECUTION = FRAMEWORK / "scripts" / "model_execution"
sys.path.insert(0, str(MODEL_EXECUTION))
sys.path.insert(0, str(DIAGNOSTICS))

from prediction_artifacts import (  # noqa: E402
    EvaluationArrayCapture,
    publish_prediction_artifacts,
)
from pfp_sensitivity_common import (  # noqa: E402
    cohort_masks,
    flat_non_root_metrics,
    load_aspect_bundle,
    verify_artifact_manifest,
)


SENSITIVITY = DIAGNOSTICS / "evaluate_pfp_label_sensitivity.py"
COMPARE = DIAGNOSTICS / "compare_pfp_label_sensitivity.py"


class PfpLabelSensitivityTests(unittest.TestCase):
    def make_obo(self, path: Path) -> None:
        values = [
            ("GO:0008150", "biological_process"),
            ("GO:0009987", "biological_process"),
            ("GO:0005575", "cellular_component"),
            ("GO:0003674", "molecular_function"),
        ]
        lines = ["format-version: 1.2", ""]
        for term, namespace in values:
            lines.extend(["[Term]", f"id: {term}", f"namespace: {namespace}"])
            if term == "GO:0009987":
                lines.append("is_a: GO:0008150 ! root")
            lines.append("")
        path.write_text("\n".join(lines), encoding="utf-8")

    def make_prediction_artifact(
        self,
        root: Path,
        obo: Path,
        mode: str = "full",
        seed: int = 42,
        has_non_root: bool = True,
        mutate_ia_before_persist: bool = False,
    ) -> Path:
        destination = root / f"prediction-artifacts-{mode}"
        stage = root / f"prediction-stage-{mode}"
        stage.mkdir()
        module = types.SimpleNamespace()
        writer_calls = {"predictions": 0, "truth": 0, "ia": 0}

        def save_predictions(predictions, protein_ids, go_terms, output_file):
            writer_calls["predictions"] += 1
            path = Path(output_file)
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("predictions\n", encoding="utf-8")

        def save_truth(labels, protein_ids, go_terms, output_file):
            writer_calls["truth"] += 1
            path = Path(output_file)
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("truth\n", encoding="utf-8")

        def save_ia(values, output_file):
            writer_calls["ia"] += 1
            Path(output_file).write_text("GO:0008150\t0\n", encoding="utf-8")

        module.save_predictions_cafa_format = save_predictions
        module.save_ground_truth_cafa_format = save_truth
        module.save_ia_file = save_ia
        original_prediction_writer = module.save_predictions_cafa_format
        scores = np.asarray(
            [[0.9, 0.1], [0.9, 0.8], [0.2, 0.1]], dtype=np.float32
        )
        truth = np.asarray(
            [[1, 0], [1, 1 if has_non_root else 0], [0, 0]], dtype=np.uint8
        )
        protein_ids = ["ROOT_ONLY", "DEEP", "ALL_ZERO"]
        go_terms = ["GO:0008150", "GO:0009987"]
        with EvaluationArrayCapture(module, "BPO", stage) as capture:
            module.save_predictions_cafa_format(
                scores, protein_ids, go_terms, root / "predictions.tsv"
            )
            module.save_ground_truth_cafa_format(
                truth, protein_ids, go_terms, root / "truth.tsv"
            )
        self.assertIs(module.save_predictions_cafa_format, original_prediction_writer)
        self.assertEqual(writer_calls, {"predictions": 1, "truth": 1, "ia": 0})
        checkpoint = root / "best_model.pt"
        checkpoint.write_bytes(b"checkpoint")
        ia = root / "BPO_ia.txt"
        ia.write_text("GO:0008150\t0\nGO:0009987\t1\n", encoding="utf-8")
        expected_ia_sha256 = self.sha256(ia)
        if mutate_ia_before_persist:
            ia.write_text("GO:0008150\t0\nGO:0009987\t2\n", encoding="utf-8")
        aspect = capture.persist(
            expected_protein_ids=protein_ids,
            expected_go_terms=go_terms,
            checkpoint=checkpoint,
            expected_checkpoint_sha256=self.sha256(checkpoint),
            cafa_metrics={"fmax": 0.6, "wfmax": 0.5, "smin": 1.2, "threshold": 0.5},
            ia_file=ia,
            expected_ia_sha256=expected_ia_sha256,
        )
        preparation = stage / "preparation_report.json"
        preparation.write_text(
            json.dumps(
                {
                    "status": "passed",
                    "benchmark_fingerprint": "benchmark-fixture",
                    "source_csv_sha256": {"bp-test.csv": "1" * 64},
                },
                sort_keys=True,
            )
            + "\n",
            encoding="utf-8",
        )
        embedding = stage / "embedding_validation_report.json"
        embedding.write_text(
            json.dumps({"status": "passed", "mode": mode}, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        manifest = {
            "schema_version": 2,
            "status": "complete",
            "benchmark_id": "fixture",
            "mode": mode,
            "seed": seed,
            "selected_aspects": ["BPO"],
            "config": {"path": "fixture.json", "sha256": "2" * 64},
            "obo": {"path": str(obo), "sha256": self.sha256(obo)},
            "provenance": {
                "framework_commit": "3" * 40,
                "pfp_commit": "4" * 40,
                "benchmark_fingerprint": "benchmark-fixture",
                "source_csv_sha256": {"bp-test.csv": "1" * 64},
                "preparation_report": {
                    "artifact_file": preparation.name,
                    "source_path": str(preparation),
                    "bytes": preparation.stat().st_size,
                    "sha256": self.sha256(preparation),
                },
                "embedding_validation_report": {
                    "artifact_file": embedding.name,
                    "source_path": str(embedding),
                    "bytes": embedding.stat().st_size,
                    "sha256": self.sha256(embedding),
                },
            },
            "aspects": {"BPO": aspect},
        }
        publish_prediction_artifacts(stage, destination, manifest)
        return destination / "prediction_artifact_manifest.json"

    @staticmethod
    def sha256(path: Path) -> str:
        import hashlib

        return hashlib.sha256(path.read_bytes()).hexdigest()

    def make_fake_cafaeval(self, root: Path) -> Path:
        package = root / "cafaeval"
        package.mkdir(parents=True)
        (package / "__init__.py").write_text("", encoding="utf-8")
        (package / "evaluation.py").write_text(
            "import json\n"
            "from pathlib import Path\n"
            "LAST={}\n"
            "def cafa_eval(obo, predictions, truth, **kwargs):\n"
            " global LAST\n"
            " truth_lines=Path(truth).read_text().splitlines()\n"
            " pred_lines=[]\n"
            " for p in sorted(Path(predictions).glob('*.tsv')): pred_lines.extend(p.read_text().splitlines())\n"
            " LAST={'truth_ids':sorted({x.split('\\t')[0] for x in truth_lines if x}), 'prediction_ids':sorted({x.split('\\t')[0] for x in pred_lines if x}), 'prediction_terms':sorted({x.split('\\t')[1] for x in pred_lines if x}), 'kwargs':kwargs}\n"
            " return ()\n"
            "def write_results(*args, out_dir):\n"
            " p=Path(out_dir); p.mkdir(parents=True, exist_ok=True)\n"
            " (p/'observed_inputs.json').write_text(json.dumps(LAST, sort_keys=True))\n"
            " (p/'evaluation_all.tsv').write_text('tau\\tf\\tpr\\trc\\tcov\\tf_w\\tpr_w\\trc_w\\ts\\n0.500000\\t0.600000\\t0.600000\\t0.600000\\t1.000000\\t0.500000\\t0.500000\\t0.500000\\t1.200000\\n')\n"
            " (p/'evaluation_best_s.tsv').write_text('s\\n1.200000\\n')\n"
            " (p/'evaluation_best_f_w.tsv').write_text('tau\\tf_w\\tpr_w\\trc_w\\n0.500000\\t0.500000\\t0.500000\\t0.500000\\n')\n",
            encoding="utf-8",
        )
        return root

    def test_bundle_validation_and_flat_metric_cohort(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            obo = root / "go.obo"
            self.make_obo(obo)
            manifest_path = self.make_prediction_artifact(root, obo)
            manifest, artifact_root = verify_artifact_manifest(manifest_path)
            bundle = load_aspect_bundle(manifest, artifact_root, "BPO")
            masks = cohort_masks(bundle["truth"], bundle["root_index"])
            self.assertEqual(masks["root_only"].tolist(), [True, False, False])
            self.assertEqual(masks["eligible_non_root"].tolist(), [False, True, False])
            self.assertEqual(masks["all_zero"].tolist(), [False, False, True])
            flat = flat_non_root_metrics(
                bundle["truth"], bundle["scores"], bundle["root_index"], 0.5
            )
            self.assertEqual(flat["status"], "complete")
            self.assertAlmostEqual(
                flat["fixed_at_canonical_threshold"]["macro_f"], 1.0
            )

    def test_sensitivity_cli_is_separate_and_reproduces_canonical_metrics(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            obo = root / "go.obo"
            self.make_obo(obo)
            manifest = self.make_prediction_artifact(root, obo)
            fake = self.make_fake_cafaeval(root / "fake")
            output = root / "sensitivity"
            environment = os.environ.copy()
            environment["PYTHONDONTWRITEBYTECODE"] = "1"
            environment["PYTHONPATH"] = os.pathsep.join(
                [str(fake), environment.get("PYTHONPATH", "")]
            )
            result = subprocess.run(
                [
                    sys.executable,
                    str(SENSITIVITY),
                    "--prediction-manifest",
                    str(manifest),
                    "--obo-file",
                    str(obo),
                    "--output-dir",
                    str(output),
                ],
                check=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                env=environment,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            report = json.loads(
                (output / "root_exclusion_sensitivity.json").read_text()
            )
            bpo = report["aspects"]["BPO"]
            self.assertEqual(bpo["cohorts"]["root_only"], 1)
            self.assertEqual(bpo["cohorts"]["eligible_non_root"], 1)
            self.assertEqual(bpo["cohorts"]["captured_rows"], 3)
            self.assertEqual(bpo["cohorts"]["cafaeval_evaluable_targets"], 2)
            self.assertEqual(bpo["cohorts"]["all_zero"], 1)
            self.assertEqual(bpo["canonical_recheck_absolute_deltas"]["fmax"], 0.0)
            self.assertEqual(bpo["root_only_excluded"]["status"], "complete")
            self.assertTrue((output / "RUN_COMPLETE.json").is_file())
            canonical_inputs = json.loads(
                (output / "BPO" / "canonical_all_targets" / "observed_inputs.json").read_text()
            )
            excluded_inputs = json.loads(
                (output / "BPO" / "root_only_excluded" / "observed_inputs.json").read_text()
            )
            baseline_inputs = json.loads(
                (output / "BPO" / "root_only_prediction_baseline" / "observed_inputs.json").read_text()
            )
            self.assertEqual(canonical_inputs["truth_ids"], ["DEEP", "ROOT_ONLY"])
            self.assertEqual(canonical_inputs["prediction_ids"], ["ALL_ZERO", "DEEP", "ROOT_ONLY"])
            self.assertEqual(excluded_inputs["truth_ids"], ["DEEP"])
            self.assertEqual(excluded_inputs["prediction_ids"], ["DEEP"])
            self.assertEqual(baseline_inputs["truth_ids"], ["DEEP", "ROOT_ONLY"])
            self.assertEqual(baseline_inputs["prediction_ids"], ["DEEP", "ROOT_ONLY"])
            self.assertEqual(baseline_inputs["prediction_terms"], ["GO:0008150"])
            self.assertFalse(baseline_inputs["kwargs"]["no_orphans"])

            sequence_manifest = self.make_prediction_artifact(
                root, obo, mode="sequence-only"
            )
            sequence_output = root / "sensitivity-sequence"
            result = subprocess.run(
                [
                    sys.executable,
                    str(SENSITIVITY),
                    "--prediction-manifest",
                    str(sequence_manifest),
                    "--obo-file",
                    str(obo),
                    "--output-dir",
                    str(sequence_output),
                ],
                check=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                env=environment,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            comparison = root / "comparison"
            result = subprocess.run(
                [
                    sys.executable,
                    str(COMPARE),
                    "--report",
                    str(output / "root_exclusion_sensitivity.json"),
                    "--report",
                    str(sequence_output / "root_exclusion_sensitivity.json"),
                    "--output-dir",
                    str(comparison),
                ],
                check=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                env=environment,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            compared = json.loads(
                (comparison / "label_sensitivity_comparison.json").read_text()
            )
            self.assertEqual(len(compared["mode_delta_rows"]), 1)
            self.assertEqual(
                compared["mode_delta_rows"][0][
                    "root_excluded_at_mode_canonical_threshold_delta"
                ],
                0.0,
            )

    def test_comparator_rejects_different_scientific_provenance(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            obo = root / "go.obo"
            self.make_obo(obo)
            full = self.make_prediction_artifact(root, obo, mode="full", seed=42)
            sequence = self.make_prediction_artifact(
                root, obo, mode="sequence-only", seed=7
            )
            fake = self.make_fake_cafaeval(root / "fake")
            environment = {
                **os.environ,
                "PYTHONDONTWRITEBYTECODE": "1",
                "PYTHONPATH": os.pathsep.join(
                    [str(fake), os.environ.get("PYTHONPATH", "")]
                ),
            }
            reports = []
            for manifest, label in ((full, "full"), (sequence, "sequence")):
                output = root / f"sensitivity-{label}"
                result = subprocess.run(
                    [
                        sys.executable,
                        str(SENSITIVITY),
                        "--prediction-manifest",
                        str(manifest),
                        "--obo-file",
                        str(obo),
                        "--output-dir",
                        str(output),
                    ],
                    check=False,
                    capture_output=True,
                    text=True,
                    env=environment,
                )
                self.assertEqual(result.returncode, 0, result.stderr)
                reports.append(output / "root_exclusion_sensitivity.json")
            result = subprocess.run(
                [
                    sys.executable,
                    str(COMPARE),
                    "--report",
                    str(reports[0]),
                    "--report",
                    str(reports[1]),
                    "--output-dir",
                    str(root / "comparison"),
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("provenance differs", result.stderr)

    def test_non_evaluable_aspect_remains_visible_in_comparison(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            obo = root / "go.obo"
            self.make_obo(obo)
            manifests = [
                self.make_prediction_artifact(
                    root, obo, mode=mode, has_non_root=False
                )
                for mode in ("full", "sequence-only")
            ]
            fake = self.make_fake_cafaeval(root / "fake")
            environment = {
                **os.environ,
                "PYTHONDONTWRITEBYTECODE": "1",
                "PYTHONPATH": os.pathsep.join(
                    [str(fake), os.environ.get("PYTHONPATH", "")]
                ),
            }
            reports = []
            for manifest, label in zip(manifests, ("full", "sequence")):
                output = root / f"sensitivity-{label}"
                result = subprocess.run(
                    [
                        sys.executable,
                        str(SENSITIVITY),
                        "--prediction-manifest",
                        str(manifest),
                        "--obo-file",
                        str(obo),
                        "--output-dir",
                        str(output),
                    ],
                    check=False,
                    capture_output=True,
                    text=True,
                    env=environment,
                )
                self.assertEqual(result.returncode, 0, result.stderr)
                reports.append(output / "root_exclusion_sensitivity.json")
            comparison = root / "comparison"
            result = subprocess.run(
                [
                    sys.executable,
                    str(COMPARE),
                    "--report",
                    str(reports[0]),
                    "--report",
                    str(reports[1]),
                    "--output-dir",
                    str(comparison),
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            report = json.loads(
                (comparison / "label_sensitivity_comparison.json").read_text()
            )
            self.assertEqual(len(report["run_rows"]), 2)
            self.assertEqual(len(report["mode_delta_rows"]), 1)
            self.assertEqual(
                report["mode_delta_rows"][0]["status"],
                "not_evaluable_no_non_root_targets",
            )
            self.assertIsNone(
                report["mode_delta_rows"][0]["root_excluded_fmax_delta"]
            )

    def test_capture_restores_writers_after_exception(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            module = types.SimpleNamespace()
            module.save_predictions_cafa_format = lambda *args, **kwargs: None
            module.save_ground_truth_cafa_format = lambda *args, **kwargs: None
            module.save_ia_file = lambda *args, **kwargs: None
            originals = {
                key: getattr(module, key)
                for key in (
                    "save_predictions_cafa_format",
                    "save_ground_truth_cafa_format",
                    "save_ia_file",
                )
            }
            with self.assertRaisesRegex(RuntimeError, "fixture failure"):
                with EvaluationArrayCapture(module, "BPO", Path(name)):
                    raise RuntimeError("fixture failure")
            for key, value in originals.items():
                self.assertIs(getattr(module, key), value)

    def test_capture_rejects_precomputed_ia_mutation(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            obo = root / "go.obo"
            self.make_obo(obo)
            with self.assertRaisesRegex(ValueError, "Precomputed IA changed"):
                self.make_prediction_artifact(
                    root, obo, mutate_ia_before_persist=True
                )


if __name__ == "__main__":
    unittest.main()
