from __future__ import annotations

import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[3]
ANALYZER = REPO_ROOT / "scripts/embeddings/analyze_embedding_reproducibility.py"
RUNTIME = REPO_ROOT / "scripts/embeddings/record_embedding_runtime.py"
WORKFLOW = REPO_ROOT / "scripts/embeddings/run_contemporary_embedding_reproducibility.sh"
WRAPPER = REPO_ROOT / "hpc_jobs/active/hpc_contemporary_embedding_reproducibility.sh"
TEXT_RECIPE = REPO_ROOT / "scripts/embeddings/run_pfp_temporal_text.py"


class ReproducibilityAnalysisTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.contract = self.root / "contract.json"
        self.controls = self.root / "controls.tsv"
        self.input_file = self.root / "input.tsv"
        self.output = self.root / "output"
        self.contract.write_text(
            json.dumps(
                {
                    "policy": {
                        "modalities": {
                            "text": {"cache_directory": "text", "dimension": 3}
                        }
                    }
                }
            ),
            encoding="utf-8",
        )
        self.controls.write_text(
            "protein_id\tmodality\tsequence_sha256\n"
            "P1\ttext\tsha1\n"
            "P2\ttext\tsha2\n",
            encoding="utf-8",
        )
        self.input_file.write_text("P1\tone\nP2\ttwo\n", encoding="utf-8")
        self.roots = {
            name: self.root / name
            for name in ("baseline", "repeat_one", "repeat_two")
        }
        for root in self.roots.values():
            (root / "text").mkdir(parents=True)

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def command(self) -> list[str]:
        return [
            sys.executable,
            str(ANALYZER),
            "--contract",
            str(self.contract),
            "--controls",
            str(self.controls),
            "--modality",
            "text",
            "--baseline-root",
            str(self.roots["baseline"]),
            "--repeat-one-root",
            str(self.roots["repeat_one"]),
            "--repeat-two-root",
            str(self.roots["repeat_two"]),
            "--input-file",
            str(self.input_file),
            "--output-dir",
            str(self.output),
            "--minimum-compared",
            "2",
        ]

    def write_arrays(self) -> None:
        baseline = {
            "P1": np.asarray([1.0, 2.0, 3.0], dtype=np.float32),
            "P2": np.asarray([0.0, 1.0, 2.0], dtype=np.float32),
        }
        repeat_one = {
            "P1": baseline["P1"].copy(),
            "P2": baseline["P2"] + np.asarray([2e-6, 0.0, 0.0], dtype=np.float32),
        }
        repeat_two = {
            "P1": baseline["P1"].copy(),
            "P2": baseline["P2"] + np.asarray([3e-6, 0.0, 0.0], dtype=np.float32),
        }
        for protein_id in baseline:
            np.save(self.roots["baseline"] / "text" / f"{protein_id}.npy", baseline[protein_id])
            np.save(
                self.roots["repeat_one"] / "text" / f"{protein_id}.npy",
                repeat_one[protein_id],
            )
            np.save(
                self.roots["repeat_two"] / "text" / f"{protein_id}.npy",
                repeat_two[protein_id],
            )

    def test_numeric_differences_are_reported_without_becoming_integrity_failures(self) -> None:
        self.write_arrays()
        result = subprocess.run(self.command(), capture_output=True, text=True, check=False)
        self.assertEqual(result.returncode, 0, result.stderr)
        report = json.loads(
            (self.output / "embedding_reproducibility.json").read_text(encoding="utf-8")
        )
        self.assertTrue(report["integrity_passed"])
        self.assertEqual(report["summaries"]["repeat_1_vs_repeat_2"]["compared"], 2)
        self.assertGreater(
            report["summaries"]["baseline_vs_repeat_2"]["max_abs_difference_max"],
            0.0,
        )
        self.assertTrue((self.output / "embedding_reproducibility.tsv").is_file())
        self.assertTrue((self.output / "embedding_reproducibility.md").is_file())
        with (self.output / "input_manifest.tsv").open(encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 1)
        self.assertEqual(len(rows[0]["sha256"]), 64)

    def test_missing_repeat_array_is_an_integrity_failure(self) -> None:
        self.write_arrays()
        (self.roots["repeat_two"] / "text/P2.npy").unlink()
        result = subprocess.run(self.command(), capture_output=True, text=True, check=False)
        self.assertEqual(result.returncode, 1)
        report = json.loads(
            (self.output / "embedding_reproducibility.json").read_text(encoding="utf-8")
        )
        self.assertFalse(report["integrity_passed"])
        self.assertEqual(
            report["summaries"]["repeat_1_vs_repeat_2"]["integrity_failures"], 1
        )

    def test_runtime_report_survives_a_host_without_nvidia_smi(self) -> None:
        output = self.root / "runtime.json"
        result = subprocess.run(
            [
                sys.executable,
                str(RUNTIME),
                "--output",
                str(output),
                "--source-file",
                f"analysis={ANALYZER}",
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        report = json.loads(output.read_text(encoding="utf-8"))
        self.assertIn("hostname", report)
        self.assertEqual(report["sources"][0]["label"], "analysis")


class ReproducibilityWorkflowContractTest(unittest.TestCase):
    def test_wrapper_pins_animal_and_always_cleans_owned_scratch(self) -> None:
        source = WRAPPER.read_text(encoding="utf-8")
        self.assertIn("#$ -l hostname=animal-206-2.local", source)
        self.assertIn("#$ -pe gpu 1", source)
        self.assertIn('EXPECTED_HOST="animal-206-2.local"', source)
        self.assertIn('WORK="/scratch0/contemporary_embedding_reproducibility_${JOB_TOKEN}"', source)
        self.assertIn('rm -rf "$WORK"', source)

    def test_workflow_runs_two_repeats_and_has_no_merge_path(self) -> None:
        source = WORKFLOW.read_text(encoding="utf-8")
        self.assertIn("run_text_repeat repeat_1", source)
        self.assertIn("run_text_repeat repeat_2", source)
        self.assertIn("run_structure_repeat repeat_1", source)
        self.assertIn("run_structure_repeat repeat_2", source)
        self.assertIn('"accepted_embedding_state_modified": False', source)
        self.assertIn('"source_cache_writes_allowed": True', source)
        self.assertIn("--balance-global-splits", source)
        self.assertNotIn("manage_resumable_embedding_state.py\" merge", source)

    def test_text_recipe_materializes_an_empty_historical_file(self) -> None:
        source = TEXT_RECIPE.read_text(encoding="utf-8")
        self.assertIn("if not historical.exists():", source)
        self.assertIn("historical.touch()", source)


if __name__ == "__main__":
    unittest.main()
