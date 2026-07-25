from __future__ import annotations

import csv
import importlib.util
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


def load_text_recipe():
    specification = importlib.util.spec_from_file_location(
        "test_run_pfp_temporal_text", TEXT_RECIPE
    )
    if specification is None or specification.loader is None:
        raise RuntimeError(f"Cannot import {TEXT_RECIPE}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


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

    def test_text_recipe_passes_requested_cutoff_to_pfp_selector(self) -> None:
        module = load_text_recipe()

        class StubPfp:
            CUTOFF_DATE = "2016-02-17"

            @staticmethod
            def find_historical_version(versions, cutoff_date=CUTOFF_DATE):
                cutoff = module.datetime.strptime(cutoff_date, "%Y-%m-%d")
                eligible = [
                    row
                    for row in versions
                    if module.datetime.strptime(
                        row["firstReleaseDate"], "%d-%b-%Y"
                    )
                    <= cutoff
                ]
                return max(
                    (row["entryVersion"] for row in eligible), default=None
                )

        stub = StubPfp()
        binding = module.configure_historical_cutoff(stub, "2025-03-08")
        versions = [
            {"firstReleaseDate": "01-Jan-2020", "entryVersion": 1},
            {"firstReleaseDate": "01-Jan-2026", "entryVersion": 2},
        ]
        self.assertEqual(stub.find_historical_version(versions), 1)
        self.assertEqual(stub.CUTOFF_DATE, "2025-03-08")
        self.assertEqual(binding["effective_cutoff"], "2025-03-08")
        with self.assertRaisesRegex(ValueError, "different from the framework contract"):
            stub.find_historical_version(versions, cutoff_date="2016-02-17")

    def test_historical_state_is_cutoff_scoped_and_contract_bound(self) -> None:
        module = load_text_recipe()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            first = root / "cutoff_2016-02-17"
            second = root / "cutoff_2025-03-08"
            first_payload = {
                "schema_version": 1,
                "requested_cutoff": "2016-02-17",
            }
            second_payload = {
                "schema_version": 1,
                "requested_cutoff": "2025-03-08",
            }
            first_contract = module.ensure_state_contract(first, first_payload)
            second_contract = module.ensure_state_contract(second, second_payload)
            self.assertNotEqual(first_contract.parent, second_contract.parent)
            self.assertEqual(
                json.loads(first_contract.read_text(encoding="utf-8")), first_payload
            )
            module.ensure_state_contract(first, first_payload)
            with self.assertRaisesRegex(ValueError, "contract changed"):
                module.ensure_state_contract(first, second_payload)

    def test_historical_state_rejects_uncontracted_legacy_content(self) -> None:
        module = load_text_recipe()
        with tempfile.TemporaryDirectory() as temporary:
            state = Path(temporary) / "cutoff_2025-03-08"
            state.mkdir()
            (state / "historical_checkpoint.txt").write_text("P1\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "unsafe resume"):
                module.ensure_state_contract(
                    state,
                    {"schema_version": 1, "requested_cutoff": "2025-03-08"},
                )

    def test_text_recipe_end_to_end_uses_separate_cutoff_state(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            pfp = root / "pfp"
            data = pfp / "data"
            script_dir = pfp / "scripts"
            assessment = root / "assessment"
            script_dir.mkdir(parents=True)
            data.mkdir()
            assessment.mkdir()
            extractor = script_dir / "extract_uniprot_text.py"
            extractor.write_text(
                "from datetime import datetime\n"
                "import json\n"
                "from pathlib import Path\n"
                "CUTOFF_DATE = '2016-02-17'\n"
                "TEXT_BUNDLE_METADATA = Path('unused.json')\n"
                "def get_split_protein_ids(data_dir, splits): return ['TEST']\n"
                "def find_historical_version(versions, cutoff_date=CUTOFF_DATE):\n"
                "    cutoff = datetime.strptime(cutoff_date, '%Y-%m-%d')\n"
                "    valid = [v for v in versions if datetime.strptime(v['firstReleaseDate'], '%d-%b-%Y') <= cutoff]\n"
                "    return max((v['entryVersion'] for v in valid), default=None)\n"
                "def run_current_extraction(data_dir, cafa_assessment_dir, output_file, checkpoint_file):\n"
                "    output_file.parent.mkdir(parents=True, exist_ok=True)\n"
                "    output_file.write_text('TRAIN\\tcurrent\\nTEST\\tcurrent\\n')\n"
                "    checkpoint_file.write_text('TRAIN\\nTEST\\n')\n"
                "    return {'records': 2}\n"
                "def process_single_historical_protein(cafa_id, accession, session, raw_dir):\n"
                "    versions = [{'firstReleaseDate': '01-Jan-2010', 'entryVersion': 1}, {'firstReleaseDate': '01-Jan-2024', 'entryVersion': 2}, {'firstReleaseDate': '01-Jan-2026', 'entryVersion': 3}]\n"
                "    selected = find_historical_version(versions)\n"
                "    raw_dir.mkdir(parents=True, exist_ok=True)\n"
                "    (raw_dir / f'{cafa_id}_{accession}.txt').write_text(str(selected))\n"
                "    return cafa_id, f'historical-{selected}', 'ok'\n"
                "def extract_historical_text(data_dir, cafa_assessment_dir, output_file, checkpoint_file, raw_dir, splits, workers):\n"
                "    _, description, _ = process_single_historical_protein('TEST', 'ACC', None, raw_dir)\n"
                "    selected = int(description.rsplit('-', 1)[-1])\n"
                "    output_file.parent.mkdir(parents=True, exist_ok=True)\n"
                "    output_file.write_text(f'TEST\\thistorical-{selected}\\n')\n"
                "    checkpoint_file.write_text('TEST\\n')\n"
                "    return 1, 0, {'selected': selected}\n"
                "def build_historical_punct_v1_test_tsv(historical_tsv, output_tsv, data_dir):\n"
                "    output_tsv.write_text(historical_tsv.read_text())\n"
                "    return {'records': 1}\n"
                "def build_mixed_temporal_tsv(current_tsv, hist_test_tsv, output_tsv, bundle_dir, historical_tsv, data_dir):\n"
                "    output_tsv.write_text('TRAIN\\tcurrent\\n' + hist_test_tsv.read_text())\n"
                "    TEXT_BUNDLE_METADATA.parent.mkdir(parents=True, exist_ok=True)\n"
                "    TEXT_BUNDLE_METADATA.write_text(json.dumps({'mixed': str(output_tsv)}))\n"
                "    return {'records': 2}\n",
                encoding="utf-8",
            )

            def run(cutoff: str) -> dict:
                result = subprocess.run(
                    [
                        sys.executable,
                        str(TEXT_RECIPE),
                        "--pfp-root",
                        str(pfp),
                        "--cafa-assessment-dir",
                        str(assessment),
                        "--cutoff-date",
                        cutoff,
                        "--workers",
                        "1",
                    ],
                    capture_output=True,
                    text=True,
                    check=False,
                )
                self.assertEqual(result.returncode, 0, result.stderr)
                return json.loads(result.stdout)

            recent = run("2025-03-08")
            older = run("2020-03-08")
            self.assertEqual(recent["historical_status_counts"], {"selected": 2})
            self.assertEqual(older["historical_status_counts"], {"selected": 1})
            self.assertEqual(recent["effective_cutoff"], "2025-03-08")
            self.assertEqual(older["effective_cutoff"], "2020-03-08")
            recent_state = Path(recent["historical_state_dir"])
            older_state = Path(older["historical_state_dir"])
            self.assertNotEqual(recent_state, older_state)
            self.assertEqual(
                (recent_state / "historical_raw/TEST_ACC.txt").read_text(), "2"
            )
            self.assertEqual(
                (older_state / "historical_raw/TEST_ACC.txt").read_text(), "1"
            )
            recent_versions = [
                json.loads(line)
                for line in (
                    recent_state / "selected_unisave_versions.jsonl"
                ).read_text(encoding="utf-8").splitlines()
            ]
            older_versions = [
                json.loads(line)
                for line in (
                    older_state / "selected_unisave_versions.jsonl"
                ).read_text(encoding="utf-8").splitlines()
            ]
            self.assertEqual(recent_versions[0]["entry_version"], 2)
            self.assertEqual(older_versions[0]["entry_version"], 1)
            self.assertEqual(recent["selected_unisave_versions_count"], 1)
            self.assertEqual(older["selected_unisave_versions_count"], 1)
            self.assertTrue(
                recent["selected_unisave_versions_audit"][
                    "all_successful_descriptions_bound"
                ]
            )

    def test_version_provenance_rejects_an_unbound_raw_record(self) -> None:
        module = load_text_recipe()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            historical = root / "historical.tsv"
            raw = root / "raw"
            selected = root / "selected.jsonl"
            historical.write_text("P1\tdescription\n", encoding="utf-8")
            raw.mkdir()
            (raw / "P1_A1.txt").write_text("record", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "lack selected-version provenance"):
                module.audit_version_provenance(historical, raw, selected)


if __name__ == "__main__":
    unittest.main()
