from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


FRAMEWORK = Path(__file__).parents[3]
COMPARE = FRAMEWORK / "scripts" / "diagnostics" / "compare_pfp_modality_runs.py"
MODALITIES = {
    "full": ("sequence", "text", "structure", "ppi"),
    "sequence-only": ("sequence",),
    "sequence-text": ("sequence", "text"),
    "sequence-structure": ("sequence", "structure"),
    "sequence-ppi": ("sequence", "ppi"),
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class PfpModalityComparisonTests(unittest.TestCase):
    def embedding_report(self, mode: str) -> dict[str, object]:
        return {
            "status": "passed",
            "mode": mode,
            "modalities": {
                modality: {"valid_content_sha256": (modality[0] * 64)[:64]}
                for modality in MODALITIES[mode]
            },
            "information_accretion": {},
            "embedding_evidence_binding": {
                "contract_sha256": "b" * 64,
                "target_manifest_sha256": "c" * 64,
                "pair_status_sha256": "d" * 64,
            },
        }

    def make_run(
        self,
        root: Path,
        mode: str,
        fmax: float,
        wfmax: float,
        smin: float,
        seed: int = 42,
        execution_mode: str = "train-eval",
        framework_commit: str = "1" * 40,
    ) -> Path:
        run_root = root / mode
        reports = run_root / "reports"
        reports.mkdir(parents=True)
        preparation = {
            "status": "passed",
            "benchmark_fingerprint": "benchmark-fingerprint",
            "source_csv_sha256": {"bp-test.csv": "e" * 64},
        }
        embedding = self.embedding_report(mode)
        files = {
            "reports/preparation.json": preparation,
            "reports/embedding_cache.json": embedding,
            "reports/embedding_cache_post.json": embedding,
            "run_config.json": {"name": "fixture"},
        }
        for relative, value in files.items():
            path = run_root / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(json.dumps(value, sort_keys=True) + "\n", encoding="utf-8")
        report = {
            "schema_version": 1,
            "status": "passed",
            "benchmark_id": "fixture",
            "benchmark_fingerprint": "benchmark-fingerprint",
            "execution_mode": execution_mode,
            "modality_mode": mode,
            "seed": seed,
            "aspects": ["BPO"],
            "framework_commit": framework_commit,
            "pfp_commit": "2" * 40,
            "environment": {"python": "fixture", "packages": {"cafaeval": "fixture"}},
            "preparation_report_sha256": sha256(run_root / "reports/preparation.json"),
            "embedding_report_sha256": sha256(run_root / "reports/embedding_cache.json"),
            "embedding_post_report_sha256": sha256(
                run_root / "reports/embedding_cache_post.json"
            ),
            "results": {
                "BPO": {
                    "checkpoint_sha256": "3" * 64,
                    "metrics": {
                        "cafa_fmax": fmax,
                        "cafa_wfmax": wfmax,
                        "cafa_smin": smin,
                        "cafa_threshold": 0.5,
                        "cafa_evaluator_policy": (
                            "strict-ia-norm-cafa-prop-max-no-fallback"
                        ),
                    },
                }
            },
        }
        report_path = reports / "run_report.json"
        report_path.write_text(json.dumps(report, sort_keys=True) + "\n", encoding="utf-8")
        manifest_path = run_root / "output_manifest.json"
        manifest_files = []
        for path in sorted(run_root.rglob("*")):
            if path.is_file() and path.name not in {"output_manifest.json", "WORKFLOW_COMPLETE.json"}:
                manifest_files.append(
                    {
                        "path": path.relative_to(run_root).as_posix(),
                        "bytes": path.stat().st_size,
                        "sha256": sha256(path),
                    }
                )
        manifest_path.write_text(
            json.dumps({"schema_version": 1, "files": manifest_files}, sort_keys=True)
            + "\n",
            encoding="utf-8",
        )
        (run_root / "WORKFLOW_COMPLETE.json").write_text(
            json.dumps(
                {
                    "complete": True,
                    "manifest": "output_manifest.json",
                    "manifest_sha256": sha256(manifest_path),
                    "run_report": "reports/run_report.json",
                },
                sort_keys=True,
            )
            + "\n",
            encoding="utf-8",
        )
        return report_path

    def make_prediction(
        self,
        root: Path,
        mode: str,
        checkpoint: str = "3" * 64,
        framework_commit: str = "1" * 40,
    ) -> Path:
        artifact_root = root / f"{mode}-capture" / "evaluation" / "prediction_artifacts"
        artifact_root.mkdir(parents=True)
        preparation = artifact_root / "preparation_report.json"
        preparation.write_text(
            json.dumps(
                {
                    "status": "passed",
                    "benchmark_fingerprint": "benchmark-fingerprint",
                    "source_csv_sha256": {"bp-test.csv": "e" * 64},
                },
                sort_keys=True,
            )
            + "\n",
            encoding="utf-8",
        )
        embedding = artifact_root / "embedding_validation_report.json"
        embedding.write_text(
            json.dumps(self.embedding_report(mode), sort_keys=True) + "\n",
            encoding="utf-8",
        )
        manifest_path = artifact_root / "prediction_artifact_manifest.json"
        manifest = {
            "schema_version": 2,
            "status": "complete",
            "benchmark_id": "fixture",
            "mode": mode,
            "seed": 42,
            "selected_aspects": ["BPO"],
            "config": {"path": "fixture", "sha256": sha256(root / mode / "run_config.json")},
            "obo": {"path": "fixture.obo", "sha256": "f" * 64},
            "provenance": {
                "framework_commit": framework_commit,
                "pfp_commit": "2" * 40,
                "benchmark_fingerprint": "benchmark-fingerprint",
                "source_csv_sha256": {"bp-test.csv": "e" * 64},
                "preparation_report": {
                    "artifact_file": preparation.name,
                    "bytes": preparation.stat().st_size,
                    "sha256": sha256(preparation),
                },
                "embedding_validation_report": {
                    "artifact_file": embedding.name,
                    "bytes": embedding.stat().st_size,
                    "sha256": sha256(embedding),
                },
            },
            "aspects": {
                "BPO": {
                    "checkpoint_sha256": checkpoint,
                    "ia_file_sha256": "a" * 64,
                    "canonical_cafa_metrics": {
                        "fmax": 0.6 if mode == "sequence-only" else 0.63,
                        "wfmax": 0.5 if mode == "sequence-only" else 0.52,
                        "smin": 1.2 if mode == "sequence-only" else 1.1,
                        "threshold": 0.5,
                    },
                }
            },
        }
        manifest_path.write_text(json.dumps(manifest, sort_keys=True) + "\n", encoding="utf-8")
        output_path = artifact_root / "output_manifest.json"
        output_files = []
        for path in sorted(artifact_root.iterdir()):
            if path.is_file() and path.name not in {"output_manifest.json", "RUN_COMPLETE.json"}:
                output_files.append(
                    {"path": path.name, "bytes": path.stat().st_size, "sha256": sha256(path)}
                )
        output_path.write_text(
            json.dumps({"schema_version": 1, "files": output_files}, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        (artifact_root / "RUN_COMPLETE.json").write_text(
            json.dumps(
                {"complete": True, "output_manifest_sha256": sha256(output_path)},
                sort_keys=True,
            )
            + "\n",
            encoding="utf-8",
        )
        return manifest_path

    def run_compare(self, root: Path, reports: list[Path], predictions: list[tuple[str, Path]] = ()) -> subprocess.CompletedProcess[str]:
        command = [sys.executable, str(COMPARE)]
        for report in reports:
            command.extend(("--run-report", str(report)))
        for mode, path in predictions:
            command.extend(("--prediction-manifest", f"{mode}={path}"))
        command.extend(("--output-dir", str(root / "comparison")))
        return subprocess.run(command, capture_output=True, text=True, check=False)

    def test_compares_retrained_modes_and_binds_predictions(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            reports = [
                self.make_run(root, "sequence-only", 0.60, 0.50, 1.20),
                self.make_run(root, "sequence-text", 0.63, 0.52, 1.10),
            ]
            predictions = [
                ("sequence-only", self.make_prediction(root, "sequence-only")),
                ("sequence-text", self.make_prediction(root, "sequence-text")),
            ]
            result = self.run_compare(root, reports, predictions)
            self.assertEqual(result.returncode, 0, result.stderr)
            comparison = json.loads(
                (root / "comparison/modality_comparison.json").read_text(encoding="utf-8")
            )
            self.assertEqual(comparison["schema_version"], 2)
            self.assertEqual(len(comparison["prediction_sources"]), 2)
            self.assertEqual(
                comparison["prediction_sources"][0]["ia_binding"]["BPO"],
                "computed_ia_bound_by_checkpoint_and_exact_canonical_metrics",
            )
            self.assertAlmostEqual(
                comparison["delta_rows"][0]["cafa_fmax_improvement"], 0.03
            )

    def test_rejects_mismatched_seed(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            baseline = self.make_run(root, "sequence-only", 0.60, 0.50, 1.20)
            candidate = self.make_run(root, "sequence-text", 0.63, 0.52, 1.10, seed=7)
            result = self.run_compare(root, [baseline, candidate])
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("seed", result.stderr)

    def test_rejects_eval_only_as_canonical(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            baseline = self.make_run(root, "sequence-only", 0.60, 0.50, 1.20)
            candidate = self.make_run(
                root, "sequence-text", 0.63, 0.52, 1.10, execution_mode="eval-only"
            )
            result = self.run_compare(root, [baseline, candidate])
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("train-eval", result.stderr)

    def test_rejects_prediction_checkpoint_mismatch(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            reports = [
                self.make_run(root, "sequence-only", 0.60, 0.50, 1.20),
                self.make_run(root, "sequence-text", 0.63, 0.52, 1.10),
            ]
            predictions = [
                ("sequence-only", self.make_prediction(root, "sequence-only")),
                ("sequence-text", self.make_prediction(root, "sequence-text", "9" * 64)),
            ]
            result = self.run_compare(root, reports, predictions)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("checkpoint differs", result.stderr)

    def test_rejects_prediction_framework_drift_without_override(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            reports = [
                self.make_run(root, "sequence-only", 0.60, 0.50, 1.20),
                self.make_run(root, "sequence-text", 0.63, 0.52, 1.10),
            ]
            predictions = [
                ("sequence-only", self.make_prediction(root, "sequence-only")),
                (
                    "sequence-text",
                    self.make_prediction(
                        root, "sequence-text", framework_commit="8" * 40
                    ),
                ),
            ]
            result = self.run_compare(root, reports, predictions)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("framework_commit", result.stderr)


if __name__ == "__main__":
    unittest.main()
