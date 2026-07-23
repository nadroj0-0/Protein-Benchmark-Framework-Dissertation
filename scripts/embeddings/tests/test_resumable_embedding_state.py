from __future__ import annotations

import csv
import gzip
import hashlib
import json
import subprocess
import sys
import tarfile
import tempfile
import unittest
from pathlib import Path

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[3]
STATE_SCRIPT = REPO_ROOT / "scripts/embeddings/manage_resumable_embedding_state.py"
RETRY_WORKSPACE = REPO_ROOT / "scripts/embeddings/prepare_embedding_retry_workspace.py"
EQUIVALENCE = REPO_ROOT / "scripts/embeddings/verify_embedding_subset_equivalence.py"
PREFETCH = REPO_ROOT / "scripts/embeddings/prefetch_alphafold_structures.py"


class ResumableEmbeddingStateTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.benchmark = self.root / "benchmark"
        self.data = self.root / "data"
        self.state = self.root / "state"
        self.benchmark.mkdir()
        self.data.mkdir()
        for aspect in ("bp", "cc", "mf"):
            for split in ("training", "validation", "test"):
                (self.benchmark / f"{aspect}-{split}.csv").write_text(
                    "Entry,Sequence,GO:0000001\nP1,ACDE,1\n", encoding="utf-8"
                )
        sequences = {"P1": "ACDE", "P2": "FGHI", "P3": "KLMN"}
        (self.data / "BPO_train_sequences.json").write_text(
            json.dumps(sequences), encoding="utf-8"
        )
        self.policy = self.root / "policy.json"
        self.policy.write_text(
            json.dumps(
                {
                    "schema_version": 1,
                    "modalities": {
                        "sequence": {
                            "cache_directory": "prott5",
                            "dimension": 4,
                            "min_accepted_count": 3,
                        },
                        "text": {
                            "cache_directory": "exp_text_embeddings_temporal",
                            "dimension": 3,
                            "min_accepted_count": 2,
                        },
                        "structure": {
                            "cache_directory": "IF1",
                            "dimension": 2,
                            "min_accepted_count": 2,
                        },
                        "ppi": {
                            "cache_directory": "ppi",
                            "dimension": 2,
                            "min_accepted_count": 1,
                        },
                    },
                }
            ),
            encoding="utf-8",
        )
        self.environment = self.root / "environment.txt"
        self.environment.write_text("python=3.9.23\nnumpy=2.0.2\n", encoding="utf-8")
        self.source = self.root / "source.dat"
        self.source.write_text("frozen-source\n", encoding="utf-8")

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def run_state(self, *arguments: str, check: bool = True) -> subprocess.CompletedProcess:
        result = subprocess.run(
            [sys.executable, str(STATE_SCRIPT), *arguments],
            check=False,
            capture_output=True,
            text=True,
        )
        if check and result.returncode != 0:
            self.fail(
                f"State command failed ({result.returncode}):\n"
                f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
            )
        return result

    def initialize(self, benchmark_id: str = "fixture") -> None:
        self.run_state(
            "initialize",
            "--state-root",
            str(self.state),
            "--benchmark-id",
            benchmark_id,
            "--benchmark-dir",
            str(self.benchmark),
            "--data-dir",
            str(self.data),
            "--policy",
            str(self.policy),
            "--pfp-commit",
            "1" * 40,
            "--framework-commit",
            "2" * 40,
            "--environment-report",
            str(self.environment),
            "--source-file",
            f"fixture={self.source}",
            "--runtime-value",
            "text_cutoff_date=2016-02-17",
        )

    @staticmethod
    def save(cache: Path, directory: str, protein_id: str, dimension: int, value: float) -> None:
        path = cache / directory
        path.mkdir(parents=True, exist_ok=True)
        np.save(path / f"{protein_id}.npy", np.full(dimension, value, dtype=np.float32))

    @staticmethod
    def write_pairs(path: Path, pairs: list[tuple[str, str]]) -> None:
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["protein_id", "modality"])
            writer.writerows(pairs)

    def test_partial_merge_retry_and_gate_preserve_one_cumulative_cache(self) -> None:
        self.initialize()
        first = self.root / "generated-first"
        self.save(first, "prott5", "P1", 4, 1)
        self.save(first, "prott5", "P2", 4, 2)
        self.save(first, "exp_text_embeddings_temporal", "P1", 3, 1)
        self.save(first, "IF1", "P1", 2, 1)
        self.save(first, "ppi", "P1", 2, 1)
        result = self.run_state(
            "merge",
            "--state-root",
            str(self.state),
            "--generated-cache-root",
            str(first),
            "--attempt-id",
            "attempt-1",
        )
        summary = json.loads(result.stdout)
        self.assertFalse(summary["embedding_gate_passed"])
        self.assertEqual(summary["newly_accepted"], 5)
        self.assertTrue((self.state / "GENERATION_INCOMPLETE.json").is_file())
        preserved = (self.state / "cache/prott5/P1.npy").read_bytes()

        pairs = self.root / "retry.tsv"
        self.write_pairs(
            pairs,
            [("P3", "sequence"), ("P2", "text"), ("P2", "structure")],
        )
        second = self.root / "generated-second"
        self.save(second, "prott5", "P3", 4, 3)
        self.save(second, "exp_text_embeddings_temporal", "P2", 3, 2)
        self.save(second, "IF1", "P2", 2, 2)
        result = self.run_state(
            "merge",
            "--state-root",
            str(self.state),
            "--generated-cache-root",
            str(second),
            "--attempt-id",
            "attempt-2",
            "--requested-pairs",
            str(pairs),
        )
        summary = json.loads(result.stdout)
        self.assertTrue(summary["embedding_gate_passed"])
        self.assertTrue((self.state / "EMBEDDING_GATE_PASSED.json").is_file())
        self.assertFalse((self.state / "GENERATION_INCOMPLETE.json").exists())
        self.assertEqual((self.state / "cache/prott5/P1.npy").read_bytes(), preserved)

        with (self.state / "needs_retry.tsv").open(encoding="utf-8") as handle:
            retry_rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertIn(("P2", "ppi"), {(row["protein_id"], row["modality"]) for row in retry_rows})
        self.assertNotIn(("P2", "text"), {(row["protein_id"], row["modality"]) for row in retry_rows})

    def test_contract_drift_is_rejected(self) -> None:
        self.initialize()
        contract = json.loads((self.state / "contract.json").read_text(encoding="utf-8"))
        self.assertNotIn("baseline", contract)
        result = self.run_state(
            "initialize",
            "--state-root",
            str(self.state),
            "--benchmark-id",
            "different",
            "--benchmark-dir",
            str(self.benchmark),
            "--data-dir",
            str(self.data),
            "--policy",
            str(self.policy),
            "--pfp-commit",
            "1" * 40,
            "--framework-commit",
            "2" * 40,
            "--environment-report",
            str(self.environment),
            "--source-file",
            f"fixture={self.source}",
            "--runtime-value",
            "text_cutoff_date=2016-02-17",
            check=False,
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("contract mismatch", result.stderr)

    def test_framework_commit_drift_is_permissive_by_default_without_weakening_contract(self) -> None:
        self.initialize()
        contract_path = self.state / "contract.json"
        original_contract = contract_path.read_bytes()
        common_arguments = (
            "initialize",
            "--state-root",
            str(self.state),
            "--benchmark-id",
            "fixture",
            "--benchmark-dir",
            str(self.benchmark),
            "--data-dir",
            str(self.data),
            "--policy",
            str(self.policy),
            "--pfp-commit",
            "1" * 40,
            "--framework-commit",
            "3" * 40,
            "--environment-report",
            str(self.environment),
            "--source-file",
            f"fixture={self.source}",
            "--runtime-value",
            "text_cutoff_date=2016-02-17",
        )

        accepted = self.run_state(*common_arguments)
        self.assertIn("every scientific contract field matched exactly", accepted.stderr)
        self.assertEqual(contract_path.read_bytes(), original_contract)
        self.assertEqual(
            json.loads(accepted.stdout)["contract_sha256"],
            json.loads(original_contract)["contract_sha256"],
        )

        rejected_arguments = list(common_arguments)
        benchmark_index = rejected_arguments.index("fixture")
        rejected_arguments[benchmark_index] = "different"
        rejected = self.run_state(*rejected_arguments, check=False)
        self.assertNotEqual(rejected.returncode, 0)
        self.assertIn("contract mismatch", rejected.stderr)
        self.assertIn("differing_fields=$.benchmark_id:value", rejected.stderr)

        strict = self.run_state(
            *common_arguments,
            "--strict-framework-commit",
            check=False,
        )
        self.assertNotEqual(strict.returncode, 0)
        self.assertIn("contract mismatch", strict.stderr)

        self.source.write_text("changed-scientific-source\n", encoding="utf-8")
        source_drift = self.run_state(*common_arguments, check=False)
        self.assertNotEqual(source_drift.returncode, 0)
        self.assertIn("contract mismatch", source_drift.stderr)

    def test_invalid_array_stays_needs_retry(self) -> None:
        self.initialize()
        generated = self.root / "generated"
        self.save(generated, "prott5", "P1", 3, 1)
        pairs = self.root / "pairs.tsv"
        self.write_pairs(pairs, [("P1", "sequence")])
        self.run_state(
            "merge",
            "--state-root",
            str(self.state),
            "--generated-cache-root",
            str(generated),
            "--attempt-id",
            "invalid",
            "--requested-pairs",
            str(pairs),
        )
        self.assertFalse((self.state / "cache/prott5/P1.npy").exists())
        ledger = (self.state / "failure_ledger.tsv").read_text(encoding="utf-8")
        self.assertIn("invalid_generated_array", ledger)

    def test_archive_baseline_and_delta_share_one_logical_state(self) -> None:
        reuse_table = self.root / "reuse.tsv"
        regenerate_table = self.root / "regenerate.tsv"
        with reuse_table.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["protein_id", "sequence", "sequence_sha256"])
            writer.writerow(["P1", "ACDE", hashlib.sha256(b"ACDE").hexdigest()])
        with regenerate_table.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["protein_id", "sequence", "sequence_sha256"])
            writer.writerow(["P2", "FGHI", hashlib.sha256(b"FGHI").hexdigest()])
            writer.writerow(["P3", "KLMN", hashlib.sha256(b"KLMN").hexdigest()])

        package = self.root / "package/data/embedding_cache"
        specifications = {
            "prott5": ("prott5", "sequence", 4),
            "text": ("exp_text_embeddings_temporal", "text", 3),
            "structure": ("IF1", "structure", 2),
            "ppi": ("ppi", "ppi", 2),
        }
        for directory, _, dimension in specifications.values():
            self.save(package, directory, "P1", dimension, 1)
        self.save(package, "prott5", "P2", 4, 2)

        report = self.root / "embedding_assembly.tsv.gz"
        with gzip.open(report, "wt", encoding="utf-8", newline="") as handle:
            fields = ["protein_id", "modality", "status", "dimension"]
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            for protein_id in ("P1", "P2", "P3"):
                for report_modality, (_, _, dimension) in specifications.items():
                    available = protein_id == "P1" or (
                        protein_id == "P2" and report_modality == "prott5"
                    )
                    writer.writerow(
                        {
                            "protein_id": protein_id,
                            "modality": report_modality,
                            "status": "available" if available else "missing",
                            "dimension": dimension,
                        }
                    )
        archive = self.root / "baseline.tar.gz"
        with tarfile.open(archive, "w:gz") as handle:
            handle.add(self.root / "package/data", arcname="data")

        result = self.run_state(
            "initialize",
            "--state-root",
            str(self.state),
            "--benchmark-id",
            "archive-fixture",
            "--benchmark-dir",
            str(self.benchmark),
            "--target-table",
            str(reuse_table),
            "--target-table",
            str(regenerate_table),
            "--policy",
            str(self.policy),
            "--pfp-commit",
            "1" * 40,
            "--framework-commit",
            "2" * 40,
            "--baseline-archive",
            str(archive),
            "--baseline-assembly-report",
            str(report),
        )
        summary = json.loads(result.stdout)
        self.assertEqual(summary["coverage"]["sequence"]["accepted"], 2)
        self.assertEqual(summary["coverage"]["text"]["accepted"], 1)
        self.assertEqual(
            len(list((self.state / "cache/prott5").glob("*.npy"))), 0
        )
        with (self.state / "pair_status.tsv").open(
            encoding="utf-8", newline=""
        ) as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        accepted = [row for row in rows if row["state"] == "accepted"]
        self.assertTrue(accepted)
        self.assertTrue(all(len(row["embedding_sha256"]) == 64 for row in accepted))

        controls = self.root / "archive-controls.tsv"
        self.write_pairs(controls, [("P1", "sequence")])
        materialized = self.root / "materialized"
        self.run_state(
            "materialize",
            "--state-root",
            str(self.state),
            "--pairs",
            str(controls),
            "--output-cache-root",
            str(materialized),
        )
        np.testing.assert_array_equal(
            np.load(materialized / "prott5/P1.npy"),
            np.full(4, 1, dtype=np.float32),
        )

        retry = self.root / "archive-retry.tsv"
        self.write_pairs(retry, [("P3", "sequence")])
        generated = self.root / "archive-generated"
        self.save(generated, "prott5", "P3", 4, 3)
        self.run_state(
            "merge",
            "--state-root",
            str(self.state),
            "--generated-cache-root",
            str(generated),
            "--attempt-id",
            "archive-attempt",
            "--requested-pairs",
            str(retry),
        )
        self.assertTrue((self.state / "cache/prott5/P3.npy").is_file())

        contract_before = (self.state / "contract.json").read_bytes()
        coverage_before = json.loads(
            (self.state / "coverage.json").read_text(encoding="utf-8")
        )["coverage"]
        baseline_path = self.state / "baseline_accepted.tsv"
        with baseline_path.open(encoding="utf-8", newline="") as handle:
            baseline_rows = list(csv.DictReader(handle, delimiter="\t"))
        legacy_lines = ["protein_id\tmodality\tarchive_member"]
        legacy_lines.extend(
            "\t".join(
                (row["protein_id"], row["modality"], row["archive_member"])
            )
            for row in baseline_rows
        )
        legacy_index = "\n".join(legacy_lines) + "\n"
        baseline_path.write_text(legacy_index, encoding="utf-8")

        tampered_mapping = legacy_index.replace(
            "P1\tsequence\t", "P2\tsequence\t", 1
        )
        baseline_path.write_text(tampered_mapping, encoding="utf-8")
        rejected_mapping = self.run_state(
            "upgrade-evidence-hashes",
            "--state-root",
            str(self.state),
            check=False,
        )
        self.assertNotEqual(rejected_mapping.returncode, 0)
        self.assertIn("mapping differs", rejected_mapping.stderr)
        self.assertEqual(baseline_path.read_text(encoding="utf-8"), tampered_mapping)
        baseline_path.write_text(legacy_index, encoding="utf-8")

        upgrade_report = self.root / "evidence-upgrade.json"
        result = self.run_state(
            "upgrade-evidence-hashes",
            "--state-root",
            str(self.state),
            "--report",
            str(upgrade_report),
        )
        upgrade = json.loads(result.stdout)
        self.assertTrue(upgrade["accepted_membership_unchanged"])
        self.assertTrue(upgrade["baseline_index_changed"])
        self.assertEqual(upgrade["baseline_pairs_hashed"], 5)
        self.assertEqual(
            upgrade["accepted_counts_before"],
            {key: value["accepted"] for key, value in coverage_before.items()},
        )
        self.assertEqual(
            upgrade["accepted_counts_after"], upgrade["accepted_counts_before"]
        )
        self.assertEqual(
            (self.state / "contract.json").read_bytes(), contract_before
        )
        self.assertTrue((self.state / "cache/prott5/P3.npy").is_file())
        self.assertTrue(
            (self.state / "EVIDENCE_HASHES_COMPLETE.json").is_file()
        )
        with baseline_path.open(encoding="utf-8", newline="") as handle:
            upgraded_rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertTrue(
            all(len(row["embedding_sha256"]) == 64 for row in upgraded_rows)
        )
        with (self.state / "pair_status.tsv").open(
            encoding="utf-8", newline=""
        ) as handle:
            status_rows = list(csv.DictReader(handle, delimiter="\t"))
        accepted_rows = [row for row in status_rows if row["state"] == "accepted"]
        self.assertTrue(
            all(len(row["embedding_sha256"]) == 64 for row in accepted_rows)
        )

        repeated = json.loads(
            self.run_state(
                "upgrade-evidence-hashes", "--state-root", str(self.state)
            ).stdout
        )
        self.assertFalse(repeated["baseline_index_changed"])
        self.assertEqual(
            repeated["accepted_counts_after"], upgrade["accepted_counts_after"]
        )

        self.run_state("summary", "--state-root", str(self.state))
        self.assertFalse(
            (self.state / "EVIDENCE_HASHES_COMPLETE.json").exists()
        )
        self.run_state("upgrade-evidence-hashes", "--state-root", str(self.state))
        self.assertTrue(
            (self.state / "EVIDENCE_HASHES_COMPLETE.json").is_file()
        )

        baseline_before_tamper = baseline_path.read_bytes()
        pair_status_before_tamper = (self.state / "pair_status.tsv").read_bytes()
        marker_before_tamper = (
            self.state / "EVIDENCE_HASHES_COMPLETE.json"
        ).read_bytes()
        archive_before_tamper = archive.read_bytes()
        with archive.open("ab") as handle:
            handle.write(b"tamper")
        rejected = self.run_state(
            "upgrade-evidence-hashes",
            "--state-root",
            str(self.state),
            check=False,
        )
        self.assertNotEqual(rejected.returncode, 0)
        self.assertIn("archive size changed", rejected.stderr)
        self.assertEqual(baseline_path.read_bytes(), baseline_before_tamper)
        self.assertEqual(
            (self.state / "pair_status.tsv").read_bytes(),
            pair_status_before_tamper,
        )
        self.assertEqual(
            (self.state / "EVIDENCE_HASHES_COMPLETE.json").read_bytes(),
            marker_before_tamper,
        )
        archive.write_bytes(archive_before_tamper)

        hydrated = self.root / "hydrated"
        self.run_state(
            "hydrate",
            "--state-root",
            str(self.state),
            "--output-cache-root",
            str(hydrated),
            "--preserve-evidence",
        )
        self.assertTrue((hydrated / "prott5/P1.npy").is_file())
        self.assertTrue((hydrated / "prott5/P2.npy").is_file())
        self.assertTrue((hydrated / "prott5/P3.npy").is_file())
        self.assertTrue(
            (self.state / "EVIDENCE_HASHES_COMPLETE.json").is_file()
        )

    def test_retry_workspace_and_equivalence_use_only_selected_ids(self) -> None:
        self.initialize()
        for aspect in ("BPO", "CCO", "MFO"):
            for split in ("train", "valid", "test"):
                names = np.asarray(["P1", "P2", "P3"], dtype=object)
                np.save(self.data / f"{aspect}_{split}_names.npy", names)
                (self.data / f"{aspect}_{split}_sequences.json").write_text(
                    json.dumps({"P1": "ACDE", "P2": "FGHI", "P3": "KLMN"}),
                    encoding="utf-8",
                )
        requested = self.root / "requested.tsv"
        controls = self.root / "controls.tsv"
        self.write_pairs(requested, [("P2", "sequence")])
        self.write_pairs(controls, [("P1", "sequence")])
        subprocess.run(
            [
                sys.executable,
                str(RETRY_WORKSPACE),
                "--data-dir",
                str(self.data),
                "--requested-pairs",
                str(requested),
                "--control-pairs",
                str(controls),
                "--modality",
                "sequence",
                "--report",
                str(self.root / "workspace.json"),
            ],
            check=True,
        )
        names = np.load(self.data / "BPO_train_names.npy", allow_pickle=True)
        self.assertEqual(set(names), {"P1", "P2"})

        reference = self.state / "cache/prott5"
        generated = self.root / "equivalence-generated/prott5"
        reference.mkdir(parents=True, exist_ok=True)
        generated.mkdir(parents=True)
        array = np.arange(4, dtype=np.float32)
        np.save(reference / "P1.npy", array)
        np.save(generated / "P1.npy", array + 1e-7)
        report = self.root / "equivalence.json"
        result = subprocess.run(
            [
                sys.executable,
                str(EQUIVALENCE),
                "--state-root",
                str(self.state),
                "--generated-cache-root",
                str(self.root / "equivalence-generated"),
                "--control-pairs",
                str(controls),
                "--modality",
                "sequence",
                "--minimum-compared",
                "1",
                "--report",
                str(report),
            ],
            check=False,
        )
        self.assertEqual(result.returncode, 0)
        self.assertEqual(json.loads(report.read_text())["failed"], 0)

    def test_diagnostic_controls_are_balanced_across_global_splits(self) -> None:
        self.initialize()
        generated = self.root / "balanced-generated"
        for protein_id in ("P1", "P2", "P3"):
            self.save(generated, "exp_text_embeddings_temporal", protein_id, 3, 1)
        self.run_state(
            "merge",
            "--state-root",
            str(self.state),
            "--generated-cache-root",
            str(generated),
            "--attempt-id",
            "balanced-controls",
        )

        plan = self.root / "balanced-plan"
        plan.mkdir()
        fields = ["protein_id", "target_memberships"]
        tables = {
            "reuse_proteins.tsv": [
                {"protein_id": "P1", "target_memberships": '["bp-training.csv"]'},
                {"protein_id": "P2", "target_memberships": '["cc-validation.csv"]'},
            ],
            "regenerate_proteins.tsv": [
                {"protein_id": "P3", "target_memberships": '["mf-test.csv"]'},
            ],
        }
        for name, rows in tables.items():
            with (plan / name).open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
                writer.writeheader()
                writer.writerows(rows)

        controls = self.root / "balanced-controls.tsv"
        result = self.run_state(
            "controls",
            "--state-root",
            str(self.state),
            "--modality",
            "text",
            "--count",
            "3",
            "--plan-dir",
            str(plan),
            "--balance-global-splits",
            "--output",
            str(controls),
        )
        report = json.loads(result.stdout)
        self.assertEqual(
            report["global_split_counts"],
            {"training": 1, "validation": 1, "test": 1},
        )
        with controls.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual([row["protein_id"] for row in rows], ["P1", "P2", "P3"])

    def test_regenerated_modalities_use_compatibility_tolerance(self) -> None:
        self.initialize()
        controls = self.root / "tolerance-controls.tsv"
        generated_root = self.root / "tolerance-generated"

        def verify(
            modality: str, directory: str, dimension: int
        ) -> tuple[subprocess.CompletedProcess, dict]:
            self.write_pairs(controls, [("P1", modality)])
            reference = self.state / "cache" / directory
            generated = generated_root / directory
            reference.mkdir(parents=True, exist_ok=True)
            generated.mkdir(parents=True, exist_ok=True)
            np.save(reference / "P1.npy", np.zeros(dimension, dtype=np.float32))
            np.save(
                generated / "P1.npy",
                np.asarray([7.1e-5] + [0.0] * (dimension - 1), dtype=np.float32),
            )
            report_path = self.root / f"{modality}-tolerance.json"
            result = subprocess.run(
                [
                    sys.executable,
                    str(EQUIVALENCE),
                    "--state-root",
                    str(self.state),
                    "--generated-cache-root",
                    str(generated_root),
                    "--control-pairs",
                    str(controls),
                    "--modality",
                    modality,
                    "--minimum-compared",
                    "1",
                    "--report",
                    str(report_path),
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            return result, json.loads(report_path.read_text(encoding="utf-8"))

        structure_result, structure_report = verify("structure", "IF1", 2)
        self.assertEqual(structure_result.returncode, 0, structure_result.stderr)
        self.assertEqual(structure_report["atol"], 1e-4)
        self.assertEqual(structure_report["failed"], 0)

        text_result, text_report = verify("text", "exp_text_embeddings_temporal", 3)
        self.assertEqual(text_result.returncode, 0, text_result.stderr)
        self.assertEqual(text_report["atol"], 1e-4)
        self.assertEqual(text_report["failed"], 0)

        ppi_result, ppi_report = verify("ppi", "ppi", 2)
        self.assertEqual(ppi_result.returncode, 0, ppi_result.stderr)
        self.assertEqual(ppi_report["atol"], 1e-4)
        self.assertEqual(ppi_report["failed"], 0)

        sequence_result, sequence_report = verify("sequence", "prott5", 4)
        self.assertNotEqual(sequence_result.returncode, 0)
        self.assertEqual(sequence_report["atol"], 1e-6)
        self.assertEqual(sequence_report["failed"], 1)

    def test_alphafold_prefetch_reuses_valid_persistent_pdb_without_api(self) -> None:
        pfp = self.root / "pfp"
        scripts = pfp / "scripts"
        scripts.mkdir(parents=True)
        (scripts / "check_alphafold_coverage.py").write_text(
            """
def get_all_cafa_proteins(data_dir):
    return {"P1"}
def build_cafa_to_accession_mapping(*args, **kwargs):
    raise AssertionError("API mapping should not run for a cached PDB")
def check_alphafold_coverage(*args, **kwargs):
    raise AssertionError("API check should not run for a cached PDB")
""",
            encoding="utf-8",
        )
        cafa = self.root / "cafa"
        cafa.mkdir()
        cache = self.root / "pdb-cache"
        cache.mkdir()
        (cache / "P1.pdb").write_text(
            "HEADER cached fixture\n"
            "ATOM      1  N   ALA A   1      11.000  12.000  13.000  1.00 90.00           N\n"
            "ATOM      2  CA  ALA A   1      12.000  13.000  14.000  1.00 90.00           C\n",
            encoding="ascii",
        )
        pdb_path = cache / "P1.pdb"
        import hashlib

        pdb_sha = hashlib.sha256(pdb_path.read_bytes()).hexdigest()
        (cache / "alphafold_source_manifest.tsv").write_text(
            "protein_id\tsha256\tsize_bytes\tpdb_url\talphafold_version\tresolved_accession\tacquired_at\n"
            f"P1\t{pdb_sha}\t{pdb_path.stat().st_size}\thttps://example/P1.pdb\t4\tP1\tfixture\n",
            encoding="utf-8",
        )
        workspace = self.root / "workspace-pdb"
        report = self.root / "prefetch.json"
        coverage = self.root / "coverage.txt"
        subprocess.run(
            [
                sys.executable,
                str(PREFETCH),
                "--pfp-root",
                str(pfp),
                "--cafa-assessment-dir",
                str(cafa),
                "--data-dir",
                str(self.data),
                "--persistent-cache-dir",
                str(cache),
                "--workspace-pdb-dir",
                str(workspace),
                "--coverage-report",
                str(coverage),
                "--report",
                str(report),
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        payload = json.loads(report.read_text(encoding="utf-8"))
        self.assertEqual(payload["cached_before"], 1)
        self.assertEqual(payload["missing_before"], 0)
        self.assertEqual(payload["available_for_if1"], 1)
        self.assertTrue((workspace / "P1.pdb").is_file())


if __name__ == "__main__":
    unittest.main()
