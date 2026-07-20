from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix, save_npz


DIAGNOSTICS = Path(__file__).parents[1]
AUDIT = DIAGNOSTICS / "audit_pfp_label_space.py"
COMPARE = DIAGNOSTICS / "compare_pfp_label_audits.py"
ASPECTS = {
    "BPO": ("bp", "GO:0008150", "GO:0009987", "biological_process"),
    "CCO": ("cc", "GO:0005575", "GO:0005622", "cellular_component"),
    "MFO": ("mf", "GO:0003674", "GO:0003824", "molecular_function"),
}
SPLITS = {"training": "train", "validation": "valid", "test": "test"}


class LabelSpaceDiagnosticsTests(unittest.TestCase):
    def make_fixture(self, root: Path, benchmark_id: str, root_only_test: bool) -> tuple[Path, Path, Path, Path]:
        benchmark = root / f"benchmark-{benchmark_id}"
        prepared = root / f"prepared-{benchmark_id}"
        benchmark.mkdir()
        prepared.mkdir()
        obo = root / f"{benchmark_id}.obo"
        stanzas = ["format-version: 1.2", ""]
        for _, (_, ontology_root, child, namespace) in ASPECTS.items():
            stanzas.extend(
                [
                    "[Term]",
                    f"id: {ontology_root}",
                    f"namespace: {namespace}",
                    "",
                    "[Term]",
                    f"id: {child}",
                    f"namespace: {namespace}",
                    f"is_a: {ontology_root} ! root",
                    "",
                ]
            )
        obo.write_text("\n".join(stanzas), encoding="utf-8")

        rows_by_split = {
            "training": [("TRAIN1", "AAAA", [1, 1]), ("TRAIN2", "BBBB", [1, 0])],
            "validation": [("VALID1", "CCCC", [1, 1])],
            "test": [
                ("TEST1", "DDDD", [1, 0 if root_only_test else 1]),
                ("TEST2", "EEEE", [1, 1]),
            ],
        }
        for aspect, (prefix, ontology_root, child, _) in ASPECTS.items():
            terms = [ontology_root, child]
            (prepared / f"{aspect}_go_terms.json").write_text(
                json.dumps(terms), encoding="utf-8"
            )
            for csv_split, pfp_split in SPLITS.items():
                rows = rows_by_split[csv_split]
                first_header = "protein" if aspect == "MFO" and csv_split == "test" else "proteins"
                with (benchmark / f"{prefix}-{csv_split}.csv").open(
                    "w", newline="", encoding="utf-8"
                ) as handle:
                    writer = csv.writer(handle, lineterminator="\n")
                    writer.writerow([first_header, "sequences", *terms])
                    for protein_id, sequence, labels in rows:
                        writer.writerow([protein_id, sequence, *labels])
                names = np.asarray([row[0] for row in rows], dtype=object)
                labels = csr_matrix(np.asarray([row[2] for row in rows], dtype=np.float32))
                sequences = {row[0]: row[1] for row in rows}
                np.save(prepared / f"{aspect}_{pfp_split}_names.npy", names)
                save_npz(prepared / f"{aspect}_{pfp_split}_labels.npz", labels)
                (prepared / f"{aspect}_{pfp_split}_sequences.json").write_text(
                    json.dumps(sequences), encoding="utf-8"
                )
        config = root / f"{benchmark_id}.json"
        config.write_text(
            json.dumps(
                {
                    "schema_version": 1,
                    "name": benchmark_id,
                    "benchmark_contract": {
                        "allow_legacy_singular_protein_header": True,
                        "allow_all_zero_rows": False,
                    },
                }
            ),
            encoding="utf-8",
        )
        return benchmark, prepared, obo, config

    def run_audit(
        self,
        benchmark_id: str,
        benchmark: Path,
        prepared: Path,
        obo: Path,
        config: Path,
        output: Path,
        ia_dir: Path | None = None,
    ) -> subprocess.CompletedProcess[str]:
        command = [
            sys.executable,
            str(AUDIT),
            "--benchmark-id",
            benchmark_id,
            "--benchmark-dir",
            str(benchmark),
            "--obo-file",
            str(obo),
            "--config",
            str(config),
            "--prepared-data",
            f"fixture={prepared}",
            "--metadata",
            "identity=30",
            "--output-dir",
            str(output),
        ]
        if ia_dir is not None:
            command.extend(("--ia-file-dir", str(ia_dir)))
        return subprocess.run(
            command,
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )

    def test_audit_is_generic_and_verifies_prepared_data(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            benchmark, prepared, obo, config = self.make_fixture(root, "fixture", True)
            output = root / "audit"
            result = self.run_audit(
                "fixture", benchmark, prepared, obo, config, output
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            report = json.loads((output / "label_space_audit.json").read_text())
            self.assertEqual(report["status"], "passed")
            self.assertEqual(report["metadata"]["identity"], "30")
            self.assertTrue(report["prepared_data_verification"][0]["passed"])
            self.assertEqual(
                report["header_compatibility_aliases"],
                [{"file": "mf-test.csv", "normalized": "proteins", "source": "protein"}],
            )
            for aspect in ASPECTS:
                test = report["files"][f"{aspect}:test"]
                self.assertEqual(test["root_only_rows"], 1)
                self.assertEqual(test["root_only_fraction"], 0.5)
                self.assertAlmostEqual(
                    test["root_only_diagnostic_baseline"]["macro_f"],
                    6.0 / 7.0,
                )
            self.assertTrue((output / "RUN_COMPLETE.json").is_file())
            self.assertTrue((output / "root_only_targets.tsv").is_file())

    def test_prepared_label_mismatch_fails_without_publication(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            benchmark, prepared, obo, config = self.make_fixture(root, "bad", True)
            labels = csr_matrix(np.asarray([[1, 1], [1, 1]], dtype=np.float32))
            save_npz(prepared / "BPO_test_labels.npz", labels)
            output = root / "audit"
            result = self.run_audit("bad", benchmark, prepared, obo, config, output)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("Prepared-data ordered_labels_sha256 differs", result.stderr)
            self.assertFalse(output.exists())

    def test_ia_inputs_are_bound_into_the_benchmark_fingerprint(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            benchmark, prepared, obo, config = self.make_fixture(root, "fixture", True)
            ia_dir = root / "ia"
            ia_dir.mkdir()
            for aspect, (_, ontology_root, child, _) in ASPECTS.items():
                (ia_dir / f"{aspect}_ia.txt").write_text(
                    f"{ontology_root}\t0\n{child}\t1\n", encoding="utf-8"
                )
            first_output = root / "audit-first"
            first = self.run_audit(
                "fixture", benchmark, prepared, obo, config, first_output, ia_dir
            )
            self.assertEqual(first.returncode, 0, first.stderr)
            first_report = json.loads(
                (first_output / "label_space_audit.json").read_text()
            )
            (ia_dir / "BPO_ia.txt").write_text(
                "GO:0008150\t0\nGO:0009987\t2\n", encoding="utf-8"
            )
            second_output = root / "audit-second"
            second = self.run_audit(
                "fixture", benchmark, prepared, obo, config, second_output, ia_dir
            )
            self.assertEqual(second.returncode, 0, second.stderr)
            second_report = json.loads(
                (second_output / "label_space_audit.json").read_text()
            )
            self.assertNotEqual(
                first_report["benchmark_fingerprint"],
                second_report["benchmark_fingerprint"],
            )
            self.assertNotEqual(
                first_report["inputs"]["ia_files"]["BPO"]["sha256"],
                second_report["inputs"]["ia_files"]["BPO"]["sha256"],
            )

    def test_child_positive_root_negative_row_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            benchmark, prepared, obo, config = self.make_fixture(root, "bad", True)
            path = benchmark / "bp-test.csv"
            with path.open(newline="", encoding="utf-8") as handle:
                rows = list(csv.reader(handle))
            rows[2][2:] = ["0", "1"]
            with path.open("w", newline="", encoding="utf-8") as handle:
                csv.writer(handle, lineterminator="\n").writerows(rows)
            output = root / "audit"
            result = self.run_audit("bad", benchmark, prepared, obo, config, output)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("labels are not ancestor-closed", result.stderr)
            self.assertFalse(output.exists())

    def test_disconnected_zero_support_header_term_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            benchmark, prepared, obo, config = self.make_fixture(root, "bad", True)
            with obo.open("a", encoding="utf-8") as handle:
                handle.write(
                    "\n[Term]\nid: GO:0099999\nnamespace: biological_process\n"
                )
            path = benchmark / "bp-test.csv"
            with path.open(newline="", encoding="utf-8") as handle:
                rows = list(csv.reader(handle))
            rows[0].append("GO:0099999")
            for row in rows[1:]:
                row.append("0")
            with path.open("w", newline="", encoding="utf-8") as handle:
                csv.writer(handle, lineterminator="\n").writerows(rows)
            output = root / "audit"
            result = self.run_audit("bad", benchmark, prepared, obo, config, output)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("term disconnected", result.stderr)
            self.assertFalse(output.exists())

    def test_comparator_accepts_distinct_benchmark_reports(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            reports = []
            for benchmark_id, root_only in (("first", True), ("second", False)):
                benchmark, prepared, obo, config = self.make_fixture(
                    root, benchmark_id, root_only
                )
                output = root / f"audit-{benchmark_id}"
                result = self.run_audit(
                    benchmark_id, benchmark, prepared, obo, config, output
                )
                self.assertEqual(result.returncode, 0, result.stderr)
                reports.append(output / "label_space_audit.json")
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
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            value = json.loads(
                (comparison / "label_space_comparison.json").read_text()
            )
            self.assertEqual(len(value["rows"]), 6)
            first_bpo = next(
                row
                for row in value["rows"]
                if row["benchmark_id"] == "first" and row["aspect"] == "BPO"
            )
            second_bpo = next(
                row
                for row in value["rows"]
                if row["benchmark_id"] == "second" and row["aspect"] == "BPO"
            )
            self.assertEqual(first_bpo["root_only_fraction"], 0.5)
            self.assertEqual(second_bpo["root_only_fraction"], 0.0)


if __name__ == "__main__":
    unittest.main()
