from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parents[1]
ACTION_COLUMNS = (
    "protein_id",
    "sequence",
    "sequence_sha256",
    "action",
    "reason",
    "matching_embedded_benchmarks",
    "embedded_benchmark_memberships",
    "target_memberships",
    "regenerate_modalities",
)
CSV_NAMES = tuple(
    f"{aspect}-{split}.csv"
    for aspect in ("bp", "cc", "mf")
    for split in ("training", "validation", "test")
)


def load_script(name: str):
    path = SCRIPT_DIR / name
    specification = importlib.util.spec_from_file_location(path.stem, path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"Cannot import {path}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def digest_text(value: str) -> str:
    return hashlib.sha256(value.encode("ascii")).hexdigest()


def action_row(protein_id: str, sequence: str, action: str, membership: str) -> dict[str, str]:
    return {
        "protein_id": protein_id,
        "sequence": sequence,
        "sequence_sha256": digest_text(sequence),
        "action": action,
        "reason": "fixture",
        "matching_embedded_benchmarks": "[]",
        "embedded_benchmark_memberships": "[]",
        "target_memberships": json.dumps([membership]),
        "regenerate_modalities": json.dumps(
            ["prott5", "text", "structure", "ppi"] if action == "regenerate" else []
        ),
    }


def write_action_table(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=ACTION_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_plan(
    plan_dir: Path,
    benchmark_dir: Path,
    reuse_rows: list[dict[str, str]],
    regenerate_rows: list[dict[str, str]],
) -> None:
    plan_dir.mkdir()
    benchmark_dir.mkdir()
    input_csvs = []
    for name in CSV_NAMES:
        contents = "proteins,sequences,GO:0000001\n"
        path = benchmark_dir / name
        path.write_text(contents, encoding="utf-8")
        input_csvs.append(
            {"relative_path": name, "sha256": hashlib.sha256(contents.encode()).hexdigest()}
        )
    write_action_table(plan_dir / "reuse_proteins.tsv", reuse_rows)
    write_action_table(plan_dir / "regenerate_proteins.tsv", regenerate_rows)
    (plan_dir / "summary.json").write_text(
        json.dumps(
            {
                "counts": {
                    "target_proteins": len(reuse_rows) + len(regenerate_rows),
                    "reuse": len(reuse_rows),
                    "regenerate": len(regenerate_rows),
                }
            }
        ),
        encoding="utf-8",
    )
    (plan_dir / "run_manifest.json").write_text(
        json.dumps({"benchmarks": {"target": {"input_csvs": input_csvs}}}),
        encoding="utf-8",
    )


class WorkspaceTests(unittest.TestCase):
    def test_workspace_is_bound_to_target_csvs_and_preserves_global_splits(self) -> None:
        module = load_script("prepare_regeneration_workspace.py")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            plan = root / "plan"
            benchmark = root / "benchmark"
            data = root / "data"
            rows = [
                action_row("TRAIN", "AAAA", "regenerate", "bp-training.csv"),
                action_row("VALID", "BBBB", "regenerate", "cc-validation.csv"),
                action_row("TEST", "CCCC", "regenerate", "mf-test.csv"),
            ]
            write_plan(plan, benchmark, [], rows)

            checksums = module.validate_target_benchmark(plan, benchmark)
            selected = module.select_rows(module.load_regeneration_rows(plan), 1)
            manifest = module.write_workspace(selected, data)

            self.assertEqual(set(checksums), set(CSV_NAMES))
            self.assertEqual(manifest["protein_count"], 3)
            self.assertEqual(
                manifest["global_split_counts"],
                {"test": 1, "training": 1, "validation": 1},
            )
            self.assertEqual(
                (data / "proteins.fasta").read_text(encoding="ascii"),
                ">TEST\nCCCC\n>TRAIN\nAAAA\n>VALID\nBBBB\n",
            )
            np.testing.assert_array_equal(
                np.load(data / "BPO_train_names.npy", allow_pickle=True), ["TRAIN"]
            )

            (benchmark / "bp-test.csv").write_text("changed\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "does not match the reuse plan"):
                module.validate_target_benchmark(plan, benchmark)

            write_action_table(
                plan / "regenerate_proteins.tsv",
                [action_row("BAD/ID", "AAAA", "regenerate", "bp-training.csv")],
            )
            with self.assertRaisesRegex(ValueError, "Unsafe protein ID"):
                module.load_regeneration_rows(plan)

    def test_retry_workspace_selects_requested_and_control_from_both_actions(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            plan = root / "plan"
            benchmark = root / "benchmark"
            data = root / "data"
            write_plan(
                plan,
                benchmark,
                [action_row("CONTROL", "AAAA", "reuse", "bp-training.csv")],
                [action_row("REQUEST", "BBBB", "regenerate", "mf-test.csv")],
            )
            requested = root / "requested.tsv"
            controls = root / "controls.tsv"
            for path, protein_id in ((requested, "REQUEST"), (controls, "CONTROL")):
                with path.open("w", encoding="utf-8", newline="") as handle:
                    writer = csv.writer(handle, delimiter="\t")
                    writer.writerow(["protein_id", "modality"])
                    writer.writerow([protein_id, "structure"])
            report = root / "report.json"
            subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT_DIR / "prepare_contemporary_retry_workspace.py"),
                    "--plan-dir",
                    str(plan),
                    "--target-benchmark-dir",
                    str(benchmark),
                    "--data-dir",
                    str(data),
                    "--requested-pairs",
                    str(requested),
                    "--control-pairs",
                    str(controls),
                    "--modality",
                    "structure",
                    "--report",
                    str(report),
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            payload = json.loads(report.read_text(encoding="utf-8"))
            self.assertEqual(payload["requested_count"], 1)
            self.assertEqual(payload["control_count"], 1)
            self.assertEqual(payload["protein_count"], 2)
            self.assertEqual(
                set(json.loads((data / "BPO_train_sequences.json").read_text())),
                {"CONTROL"},
            )
            self.assertEqual(
                set(json.loads((data / "MFO_test_sequences.json").read_text())),
                {"REQUEST"},
            )


class TextReductionTests(unittest.TestCase):
    def test_reducer_keeps_exact_cls_and_skips_producer_temporary_files(self) -> None:
        module = load_script("reduce_text_embeddings_to_cls.py")
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            raw = np.arange(1 * 4 * 768, dtype=np.float64).reshape(1, 4, 768)
            np.save(directory / "P1.npy", raw)
            np.save(directory / "P2.tmp.npy", raw)
            processed: set[str] = set()

            counts = module.sweep(directory, processed)
            self.assertEqual(counts, {"reduced": 1, "already_cls": 0})
            self.assertEqual(processed, {"P1.npy"})
            reduced = np.load(directory / "P1.npy", allow_pickle=False)
            self.assertEqual(reduced.shape, (768,))
            self.assertEqual(reduced.dtype, np.float32)
            np.testing.assert_array_equal(reduced, raw[0, 0, :].astype(np.float32))
            self.assertEqual(
                module.sweep(directory, processed), {"reduced": 0, "already_cls": 0}
            )


class PpiCompatibilityTests(unittest.TestCase):
    def test_compatibility_copy_does_not_modify_source(self) -> None:
        module = load_script("build_pfp_ppi_compat_copy.py")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source.py"
            output = root / "compat.py"
            report = root / "report.json"
            original = 'print(f"{mapped_cafa3/len(cafa3_ids)*100:.1f}")\n'
            source.write_text(original, encoding="utf-8")

            with mock.patch.object(
                sys,
                "argv",
                ["compat", "--source", str(source), "--output", str(output), "--report", str(report)],
            ):
                self.assertEqual(module.main(), 0)

            self.assertEqual(source.read_text(encoding="utf-8"), original)
            self.assertIn("max(1, len(cafa3_ids))", output.read_text(encoding="utf-8"))
            self.assertFalse(json.loads(report.read_text())["upstream_source_modified"])


class If1CompatibilityTests(unittest.TestCase):
    def test_compatibility_copy_is_cuda_safe_and_keeps_source_unchanged(self) -> None:
        module = load_script("build_pfp_if1_compat_copy.py")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source.py"
            output = root / "compat.py"
            report = root / "report.json"
            original = (
                "def extract():\n"
                + module.IMPORT_OLD
                + module.ENCODER_OLD
                + module.SUMMARY_OLD
            )
            source.write_text(original, encoding="utf-8")

            with mock.patch.object(
                sys,
                "argv",
                ["compat", "--source", str(source), "--output", str(output), "--report", str(report)],
            ):
                self.assertEqual(module.main(), 0)

            compat = output.read_text(encoding="utf-8")
            self.assertEqual(source.read_text(encoding="utf-8"), original)
            self.assertIn("batch_converter(batch, device=device)", compat)
            self.assertIn('encoder_out["encoder_out"][0][1:-1, 0]', compat)
            self.assertIn("failed for all", compat)
            payload = json.loads(report.read_text(encoding="utf-8"))
            self.assertFalse(payload["upstream_source_modified"])
            self.assertFalse(payload["scientific_output_change"])

    def test_compatibility_copy_rejects_unvalidated_source_drift(self) -> None:
        module = load_script("build_pfp_if1_compat_copy.py")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source.py"
            source.write_text(module.IMPORT_OLD, encoding="utf-8")
            with mock.patch.object(
                sys,
                "argv",
                [
                    "compat",
                    "--source", str(source),
                    "--output", str(root / "compat.py"),
                    "--report", str(root / "report.json"),
                ],
            ):
                with self.assertRaisesRegex(SystemExit, "IF1 encoder-output"):
                    module.main()


class AssemblyTests(unittest.TestCase):
    def write_array(self, root: Path, directory: str, protein_id: str, dimension: int, value: float) -> None:
        path = root / directory / f"{protein_id}.npy"
        path.parent.mkdir(parents=True, exist_ok=True)
        np.save(path, np.full(dimension, value, dtype=np.float32))

    def test_assembly_uses_action_specific_sources_and_logs_missing_modalities(self) -> None:
        module = load_script("assemble_contemporary_embedding_cache.py")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            plan = root / "plan"
            benchmark = root / "benchmark"
            published = root / "published"
            generated = root / "generated"
            output = root / "output"
            reports = root / "reports"
            reuse = [action_row("OLD", "AAAA", "reuse", "bp-training.csv")]
            regenerate = [
                action_row("NEW1", "BBBB", "regenerate", "cc-validation.csv"),
                action_row("NEW2", "CCCC", "regenerate", "mf-test.csv"),
            ]
            write_plan(plan, benchmark, reuse, regenerate)

            for specification in module.MODALITIES.values():
                self.write_array(
                    published, specification["directory"], "OLD", specification["dimension"], 1.0
                )
                self.write_array(
                    generated, specification["directory"], "NEW1", specification["dimension"], 2.0
                )
            self.write_array(published, "prott5", "NEW1", 1024, 99.0)
            self.write_array(generated, "prott5", "NEW2", 1024, 3.0)

            arguments = [
                "assemble",
                "--plan-dir", str(plan),
                "--published-cache", str(published),
                "--generated-cache", str(generated),
                "--output-cache", str(output),
                "--report-dir", str(reports),
            ]
            with mock.patch.object(sys, "argv", arguments):
                self.assertEqual(module.main(), 0)

            self.assertEqual(float(np.load(output / "prott5/OLD.npy")[0]), 1.0)
            self.assertEqual(float(np.load(output / "prott5/NEW1.npy")[0]), 2.0)
            self.assertEqual(float(np.load(output / "prott5/NEW2.npy")[0]), 3.0)
            self.assertFalse((output / "exp_text_embeddings_temporal/NEW2.npy").exists())
            self.assertEqual(
                (reports / "regenerate_missing_text.txt").read_text(encoding="utf-8"), "NEW2\n"
            )
            summary = json.loads((reports / "assembly_summary.json").read_text())
            self.assertEqual(summary["modalities"]["prott5"]["combined"]["available"], 3)

    def test_assembly_rejects_generated_proteins_outside_plan(self) -> None:
        module = load_script("assemble_contemporary_embedding_cache.py")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            generated = root / "generated"
            self.write_array(generated, "prott5", "EXTRA", 1024, 1.0)
            with self.assertRaisesRegex(ValueError, "outside the regenerate partition"):
                module.validate_generated_scope(generated, {"EXPECTED"})


class WorkflowScriptTests(unittest.TestCase):
    def test_completion_marker_does_not_invoke_git_inside_python_runtime(self) -> None:
        workflow = (SCRIPT_DIR / "run_contemporary_embedding_generation.sh").read_text(
            encoding="utf-8"
        )
        wrapper = (
            SCRIPT_DIR.parents[1]
            / "hpc_jobs/active/hpc_contemporary_embedding_generation.sh"
        ).read_text(encoding="utf-8")
        self.assertNotIn("subprocess.check_output", workflow)
        self.assertNotIn("import subprocess", workflow)
        self.assertIn('os.environ.get("PFP_COMMIT", "unknown")', workflow)
        self.assertIn('PFP_COMMIT="$PFP_COMMIT"', wrapper)
        self.assertIn('FRAMEWORK_COMMIT="$FRAMEWORK_COMMIT"', wrapper)


if __name__ == "__main__":
    unittest.main()
