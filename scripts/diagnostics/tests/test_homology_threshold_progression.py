from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest


SCRIPT = Path(__file__).resolve().parents[1] / "audit_homology_threshold_progression.py"
SPEC = importlib.util.spec_from_file_location("homology_threshold_progression", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def _write_benchmark(
    root: Path,
    assignments: dict[str, tuple[str, str]],
    *,
    sequence_override: dict[str, str] | None = None,
) -> None:
    root.mkdir()
    sequence_override = sequence_override or {}
    with (root / "protein_cluster_assignments.tsv").open("w", encoding="utf-8") as handle:
        handle.write("uniprot_accession\tmmseqs_cluster_id\tsplit\n")
        for protein, (cluster, split_name) in assignments.items():
            handle.write(f"{protein}\t{cluster}\t{split_name}\n")
    with (root / "cluster_split_assignments.tsv").open("w", encoding="utf-8") as handle:
        handle.write("mmseqs_cluster_id\tsplit\n")
        seen = set()
        for cluster, split_name in assignments.values():
            if cluster and cluster not in seen:
                handle.write(f"{cluster}\t{split_name}\n")
                seen.add(cluster)
    for aspect in MODULE.ASPECTS:
        for split_name in MODULE.SPLITS:
            with (root / f"{aspect}-{split_name}.csv").open("w", encoding="ascii") as handle:
                handle.write("proteins,sequences,GO:0000001\n")
                for protein, (_, assigned_split) in assignments.items():
                    if assigned_split == split_name:
                        sequence = sequence_override.get(protein, f"M{protein}")
                        handle.write(f"{protein},{sequence},1\n")


class HomologyThresholdProgressionTests(unittest.TestCase):
    def _args(self, root: Path, labels: tuple[str, str] = ("30", "25")):
        return MODULE.argparse.Namespace(
            benchmark=[
                f"{labels[0]}={root / labels[0]}",
                f"{labels[1]}={root / labels[1]}",
            ],
            output_dir=root / "out",
        )

    def test_identical_benchmarks_are_reported_exactly(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            assignments = {
                "P1": ("C1", "training"),
                "P2": ("C1", "training"),
                "P3": ("C2", "test"),
                "P4": ("", ""),
            }
            _write_benchmark(root / "30", assignments)
            _write_benchmark(root / "25", assignments)
            payload = MODULE.run(self._args(root))
            comparison = payload["adjacent_comparisons"][0]
            self.assertTrue(comparison["all_required_files_byte_identical"])
            self.assertTrue(comparison["retained_partition"]["partitions_exactly_identical"])
            self.assertEqual(comparison["global_split"]["changed_state_proteins"], 0)
            self.assertEqual(comparison["aspects"]["bp"]["changed_state_proteins"], 0)
            self.assertTrue((root / "out" / "RUN_COMPLETE.json").is_file())
            self.assertIn("30 to 25", (root / "out" / "summary.md").read_text())

    def test_detects_partition_and_split_changes(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            _write_benchmark(root / "30", {
                "P1": ("C1", "training"),
                "P2": ("C1", "training"),
                "P3": ("C2", "test"),
            })
            _write_benchmark(root / "25", {
                "P1": ("D1", "training"),
                "P2": ("D2", "test"),
                "P3": ("D2", "test"),
            })
            comparison = MODULE.run(self._args(root))["adjacent_comparisons"][0]
            self.assertFalse(comparison["retained_partition"]["partitions_exactly_identical"])
            self.assertEqual(comparison["global_split"]["changed_state_proteins"], 1)
            self.assertEqual(comparison["aspects"]["mf"]["transition_matrix"]["training->test"], 1)
            self.assertEqual(comparison["aspects"]["cc"]["common_complete_row_disagreements"], 0)

    def test_rejects_sequence_disagreement(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            assignments = {"P1": ("C1", "training")}
            _write_benchmark(root / "30", assignments)
            _write_benchmark(root / "25", assignments, sequence_override={"P1": "MOTHER"})
            with self.assertRaisesRegex(ValueError, "Sequence content differs"):
                MODULE.run(self._args(root))

    def test_rejects_csv_split_disagreement(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            assignments = {"P1": ("C1", "training")}
            _write_benchmark(root / "30", assignments)
            _write_benchmark(root / "25", assignments)
            csv_path = root / "25" / "bp-training.csv"
            row = csv_path.read_text().splitlines()[1]
            with (root / "25" / "bp-test.csv").open("a", encoding="ascii") as handle:
                handle.write(row + "\n")
            with self.assertRaisesRegex(ValueError, "multiple splits"):
                MODULE.run(self._args(root))


if __name__ == "__main__":
    unittest.main()
