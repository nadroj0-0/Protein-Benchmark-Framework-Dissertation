from __future__ import annotations

import ast
import contextlib
import hashlib
import io
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
sys.path.insert(0, str(SRC))

from pfp_benchmark_reuse.benchmark import BenchmarkError, parse_benchmark  # noqa: E402
from pfp_benchmark_reuse.cli import main  # noqa: E402
from pfp_benchmark_reuse.planner import PlanningError, build_plan  # noqa: E402
from pfp_benchmark_reuse.reports import (  # noqa: E402
    ReportError,
    _validate_output_manifest,
    write_reports,
)

from helpers import read_tsv, rows_in, write_benchmark  # noqa: E402


class ReportAndCliTests(unittest.TestCase):
    def test_self_comparison_outputs_exact_empty_regenerate_artifacts(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            directory = write_benchmark(
                root / "benchmark",
                rows_in("bp-training.csv", ("P2", "BBBB"), ("P1", "AAAA")),
            )
            plan = build_plan(
                (parse_benchmark("source", directory),),
                parse_benchmark("target", directory),
            )
            output = write_reports(plan, root / "output")

            self.assertEqual([row["protein_id"] for row in read_tsv(output / "reuse_proteins.tsv")], ["P1", "P2"])
            self.assertEqual(read_tsv(output / "regenerate_proteins.tsv"), [])
            self.assertEqual((output / "regenerate_proteins.txt").read_text(), "")
            self.assertEqual((output / "regenerate_proteins.fasta").read_text(), "")

    def test_outputs_partition_txt_fasta_json_and_manifests_exactly(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = parse_benchmark(
                "source",
                write_benchmark(
                    root / "source",
                    rows_in("bp-training.csv", ("P1", "AAAA"), ("OLD", "OOOO")),
                ),
            )
            target = parse_benchmark(
                "target",
                write_benchmark(
                    root / "target",
                    rows_in(
                        "mf-test.csv",
                        ("P4", "DDDD"),
                        ("P1", "AAAA"),
                        ("P2", "BBBB"),
                    ),
                ),
            )
            output = write_reports(build_plan((source,), target), root / "output")

            expected_payloads = (
                "known_embedded_proteins.tsv",
                "regenerate_proteins.fasta",
                "regenerate_proteins.tsv",
                "regenerate_proteins.txt",
                "reuse_proteins.tsv",
                "reuse_proteins.txt",
                "run_manifest.json",
                "summary.json",
                "summary.md",
            )
            expected_files = set(expected_payloads) | {"output_manifest.json", "RUN_COMPLETE.json"}
            self.assertEqual({path.name for path in output.iterdir()}, expected_files)
            reuse = read_tsv(output / "reuse_proteins.tsv")
            regenerate = read_tsv(output / "regenerate_proteins.tsv")
            self.assertEqual(
                tuple(reuse[0]),
                (
                    "protein_id",
                    "sequence",
                    "sequence_sha256",
                    "action",
                    "reason",
                    "matching_embedded_benchmarks",
                    "embedded_benchmark_memberships",
                    "target_memberships",
                    "regenerate_modalities",
                ),
            )
            self.assertEqual([row["protein_id"] for row in reuse], ["P1"])
            self.assertEqual([row["protein_id"] for row in regenerate], ["P2", "P4"])
            self.assertEqual(json.loads(regenerate[0]["regenerate_modalities"]), ["ppi", "prott5", "structure", "text"])
            self.assertEqual(json.loads(reuse[0]["matching_embedded_benchmarks"]), ["source"])
            self.assertEqual(json.loads(regenerate[0]["target_memberships"]), ["mf-test.csv"])
            self.assertEqual(
                regenerate[0]["sequence_sha256"], hashlib.sha256(b"BBBB").hexdigest()
            )
            self.assertEqual((output / "reuse_proteins.txt").read_text(), "P1\n")
            self.assertEqual((output / "regenerate_proteins.txt").read_text(), "P2\nP4\n")
            self.assertEqual(
                (output / "regenerate_proteins.fasta").read_text(),
                ">P2\nBBBB\n>P4\nDDDD\n",
            )

            known_ids = [row["protein_id"] for row in read_tsv(output / "known_embedded_proteins.tsv")]
            self.assertEqual(known_ids, ["OLD", "P1"])
            for name in ("summary.json", "run_manifest.json", "output_manifest.json", "RUN_COMPLETE.json"):
                self.assertIsInstance(json.loads((output / name).read_text()), dict)

            run_manifest = json.loads((output / "run_manifest.json").read_text())
            self.assertEqual(len(run_manifest["benchmarks"]["embedded"][0]["input_csvs"]), 9)
            self.assertEqual(len(run_manifest["benchmarks"]["target"]["input_csvs"]), 9)
            self.assertEqual(run_manifest["counts"]["regenerate"], 2)
            self.assertTrue(
                all(
                    len(item["sha256"]) == 64
                    for item in run_manifest["benchmarks"]["target"]["input_csvs"]
                )
            )

            manifest = json.loads((output / "output_manifest.json").read_text())
            self.assertEqual([item["path"] for item in manifest["files"]], list(expected_payloads))
            self.assertNotIn("output_manifest.json", {item["path"] for item in manifest["files"]})
            self.assertNotIn("RUN_COMPLETE.json", {item["path"] for item in manifest["files"]})
            for item in manifest["files"]:
                data = (output / item["path"]).read_bytes()
                self.assertEqual(item["size_bytes"], len(data))
                self.assertEqual(item["sha256"], hashlib.sha256(data).hexdigest())

            complete = json.loads((output / "RUN_COMPLETE.json").read_text())
            manifest_bytes = (output / "output_manifest.json").read_bytes()
            self.assertTrue(complete["complete"])
            self.assertEqual(complete["output_manifest"]["size_bytes"], len(manifest_bytes))
            self.assertEqual(
                complete["output_manifest"]["sha256"],
                hashlib.sha256(manifest_bytes).hexdigest(),
            )

    def test_reference_argument_order_produces_byte_identical_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            alpha = parse_benchmark(
                "alpha",
                write_benchmark(root / "alpha", rows_in("bp-training.csv", ("P1", "AAAA"))),
            )
            zeta = parse_benchmark(
                "zeta",
                write_benchmark(root / "zeta", rows_in("bp-training.csv", ("P2", "BBBB"))),
            )
            target = parse_benchmark(
                "target",
                write_benchmark(
                    root / "target",
                    rows_in("bp-training.csv", ("P1", "AAAA"), ("P2", "BBBB")),
                ),
            )
            output_path = root / "deterministic"
            first = write_reports(build_plan((zeta, alpha), target), output_path)
            first_bytes = {path.name: path.read_bytes() for path in first.iterdir()}
            shutil.rmtree(first)
            second = write_reports(build_plan((alpha, zeta), target), output_path)
            second_bytes = {path.name: path.read_bytes() for path in second.iterdir()}
            self.assertEqual(first_bytes, second_bytes)

            manifest = json.loads((second / "run_manifest.json").read_text())
            self.assertEqual(
                [item["name"] for item in manifest["benchmarks"]["embedded"]],
                ["alpha", "zeta"],
            )

    def test_interrupt_during_report_generation_leaves_no_destination_or_stage(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = parse_benchmark("source", write_benchmark(root / "source"))
            target = parse_benchmark("target", write_benchmark(root / "target"))
            output = root / "interrupted"
            with mock.patch(
                "pfp_benchmark_reuse.reports._write_summary_markdown",
                side_effect=KeyboardInterrupt,
            ):
                with self.assertRaises(KeyboardInterrupt):
                    write_reports(build_plan((source,), target), output)
            self.assertFalse(os.path.lexists(output))
            self.assertEqual(list(root.glob(".interrupted.staging-*")), [])

    def test_input_change_before_publication_cleans_stage_and_destination(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source_dir = write_benchmark(root / "source")
            target_dir = write_benchmark(root / "target")
            source = parse_benchmark("source", source_dir)
            target = parse_benchmark("target", target_dir)
            plan = build_plan((source,), target)
            with (target_dir / "bp-training.csv").open("a", encoding="utf-8") as handle:
                handle.write("P2,BBBB,1\n")
            output = root / "changed-input"
            with self.assertRaisesRegex(BenchmarkError, "changed after planning"):
                write_reports(plan, output)
            self.assertFalse(output.exists())
            self.assertEqual(list(root.glob(".changed-input.staging-*")), [])

    def test_existing_directory_file_and_broken_symlink_are_refused(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = parse_benchmark("source", write_benchmark(root / "source"))
            target = parse_benchmark("target", write_benchmark(root / "target"))
            plan = build_plan((source,), target)

            existing_directory = root / "existing-dir"
            existing_directory.mkdir()
            existing_file = root / "existing-file"
            existing_file.write_text("keep")
            for destination in (existing_directory, existing_file):
                with self.subTest(destination=destination.name):
                    with self.assertRaisesRegex(ReportError, "Refusing to overwrite"):
                        write_reports(plan, destination)

            broken = root / "broken-link"
            try:
                broken.symlink_to(root / "missing-target", target_is_directory=True)
            except OSError as exc:
                self.skipTest("Cannot create a symlink: %s" % exc)
            with self.assertRaisesRegex(ReportError, "Refusing to overwrite"):
                write_reports(plan, broken)

    def test_manifest_validator_detects_tampering(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = parse_benchmark("source", write_benchmark(root / "source"))
            target = parse_benchmark("target", write_benchmark(root / "target"))
            output = write_reports(build_plan((source,), target), root / "output")
            (output / "summary.md").write_text("tampered\n", encoding="utf-8")
            with self.assertRaisesRegex(ReportError, "Manifest identity does not match"):
                _validate_output_manifest(output)

    def test_forged_plan_record_cannot_publish_target_mismatched_fasta(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = parse_benchmark("source", write_benchmark(root / "source"))
            target = parse_benchmark(
                "target",
                write_benchmark(root / "target", rows_in("bp-training.csv", ("P2", "BBBB"))),
            )
            plan = build_plan((source,), target)
            forged = replace(
                plan.records[0], sequence="ZZZZ", sequence_sha256=hashlib.sha256(b"ZZZZ").hexdigest()
            )
            output = root / "must-not-publish"
            with self.assertRaisesRegex(PlanningError, "do not exactly match target"):
                write_reports(replace(plan, records=(forged,)), output)
            self.assertFalse(output.exists())

    def test_module_cli_supports_multiple_references_and_paths_with_spaces_and_equals(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            first = write_benchmark(
                root / "embedded benchmark=one",
                rows_in("bp-training.csv", ("P1", "AAAA")),
            )
            second = write_benchmark(
                root / "embedded benchmark two",
                rows_in("bp-training.csv", ("P2", "BBBB")),
            )
            target = write_benchmark(
                root / "target benchmark=three",
                rows_in("mf-test.csv", ("P1", "AAAA"), ("P2", "BBBB")),
            )
            output = root / "output plan"
            environment = dict(os.environ)
            environment["PYTHONPATH"] = str(SRC)
            environment["PYTHONDONTWRITEBYTECODE"] = "1"
            result = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "pfp_benchmark_reuse",
                    "plan",
                    "--embedded-benchmark",
                    "zeta=%s" % second,
                    "--embedded-benchmark",
                    "alpha=%s" % first,
                    "--target-benchmark",
                    "target=%s" % target,
                    "--output-dir",
                    str(output),
                ],
                cwd=ROOT,
                env=environment,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=False,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual([row["protein_id"] for row in read_tsv(output / "reuse_proteins.tsv")], ["P1", "P2"])
            self.assertIn("reuse=2, regenerate=0", result.stdout)

    def test_cli_rejects_repeated_target_duplicate_names_and_unsafe_names(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = write_benchmark(root / "source")
            target = write_benchmark(root / "target")
            cases = (
                [
                    "plan",
                    "--embedded-benchmark",
                    "source=%s" % source,
                    "--target-benchmark",
                    "target=%s" % target,
                    "--target-benchmark",
                    "other=%s" % target,
                    "--output-dir",
                    str(root / "out-one"),
                ],
                [
                    "plan",
                    "--embedded-benchmark",
                    "same=%s" % source,
                    "--target-benchmark",
                    "same=%s" % target,
                    "--output-dir",
                    str(root / "out-two"),
                ],
                [
                    "plan",
                    "--embedded-benchmark",
                    "bad name=%s" % source,
                    "--target-benchmark",
                    "target=%s" % target,
                    "--output-dir",
                    str(root / "out-three"),
                ],
                [
                    "plan",
                    "--embedded-benchmark",
                    "same=%s" % source,
                    "--embedded-benchmark",
                    "same=%s" % source,
                    "--target-benchmark",
                    "target=%s" % target,
                    "--output-dir",
                    str(root / "out-four"),
                ],
                [
                    "plan",
                    "--embedded-benchmark",
                    "=%s" % source,
                    "--target-benchmark",
                    "target=%s" % target,
                    "--output-dir",
                    str(root / "out-five"),
                ],
            )
            for arguments in cases:
                with self.subTest(arguments=arguments):
                    stderr = io.StringIO()
                    with contextlib.redirect_stderr(stderr):
                        self.assertEqual(main(arguments), 2)
                    self.assertIn("error:", stderr.getvalue())

    def test_cli_required_benchmark_arguments_fail_nonzero(self) -> None:
        cases = (
            ["plan", "--target-benchmark", "target=/tmp/target", "--output-dir", "/tmp/out"],
            ["plan", "--embedded-benchmark", "source=/tmp/source", "--output-dir", "/tmp/out"],
        )
        for arguments in cases:
            with self.subTest(arguments=arguments):
                with contextlib.redirect_stderr(io.StringIO()):
                    with self.assertRaises(SystemExit) as raised:
                        main(arguments)
                self.assertNotEqual(raised.exception.code, 0)

    def test_decoy_array_is_never_opened_and_source_has_no_inventory_or_numpy_import(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source_dir = write_benchmark(root / "source")
            target_dir = write_benchmark(root / "target")
            (source_dir / "must-not-open.npy").write_bytes(b"not an embedding")
            original_open = Path.open
            original_builtin_open = open

            def guarded_open(path: Path, *args: object, **kwargs: object):
                if path.suffix == ".npy":
                    raise AssertionError("Planner attempted to open an array")
                return original_open(path, *args, **kwargs)

            def guarded_builtin_open(file: object, *args: object, **kwargs: object):
                try:
                    suffix = Path(file).suffix  # type: ignore[arg-type]
                except TypeError:
                    suffix = ""
                if suffix == ".npy":
                    raise AssertionError("Planner attempted to open an array")
                return original_builtin_open(file, *args, **kwargs)

            with mock.patch("pathlib.Path.open", new=guarded_open), mock.patch(
                "builtins.open", new=guarded_builtin_open
            ):
                source = parse_benchmark("source", source_dir)
                target = parse_benchmark("target", target_dir)
                write_reports(build_plan((source,), target), root / "output")

        forbidden = {"embedding_inventory", "pfp_embedding_inventory", "numpy"}
        for path in sorted(SRC.rglob("*.py")):
            source_text = path.read_text(encoding="utf-8")
            tree = ast.parse(source_text, filename=str(path))
            imported = set()
            for node in ast.walk(tree):
                if isinstance(node, ast.Import):
                    imported.update(alias.name.split(".")[0] for alias in node.names)
                elif isinstance(node, ast.ImportFrom) and node.module:
                    imported.add(node.module.split(".")[0])
            self.assertFalse(imported & forbidden, "%s imports %s" % (path, imported & forbidden))
            for traversal_marker in (".npy", ".rglob(", ".glob(", "os.walk("):
                self.assertNotIn(traversal_marker, source_text)


if __name__ == "__main__":
    unittest.main()
