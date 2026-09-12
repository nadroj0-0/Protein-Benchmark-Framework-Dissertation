from __future__ import annotations

from pathlib import Path
import sys
import tempfile
import unittest


DIAGNOSTICS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(DIAGNOSTICS))

import audit_consecutive_homology_partitions as module  # noqa: E402


class ConsecutiveHomologyPartitionTests(unittest.TestCase):
    def test_typed_pair_counts(self) -> None:
        result = module.typed_pair_counts([(1, 1), (2, 1)])
        self.assertEqual(result, {
            "linked_linked": 1,
            "linked_background": 3,
            "background_background": 0,
        })

    def test_lineage_categories(self) -> None:
        self.assertEqual(module.lineage_category(1, 1, 0, 0), "unchanged")
        self.assertEqual(module.lineage_category(1, 1, 0, 3), "expanded")
        self.assertEqual(module.lineage_category(2, 1, 0, 0), "merged")
        self.assertEqual(module.lineage_category(1, 2, 0, 0), "split")
        self.assertEqual(module.lineage_category(2, 2, 0, 0), "complex")

    def test_member_reader_rejects_duplicates(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            path = Path(raw) / "members.tsv"
            path.write_text("M1\tC1\ttraining\t1\nM1\tC2\ttest\t0\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                list(module.member_rows(path))

    def test_partition_metrics_detect_identical_partition(self) -> None:
        result = module.partition_metrics([2, 1], [2, 1], [2, 1], 3)
        self.assertAlmostEqual(result["adjusted_rand_index"], 1.0)
        self.assertAlmostEqual(result["pair_jaccard"], 1.0)

    def test_pair_comparison_separates_labelled_and_background_merges(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            left_path = root / "left.tsv"
            right_path = root / "right.tsv"
            left_path.write_text(
                "M1\tA\ttraining\t1\n"
                "M2\tA\ttraining\t0\n"
                "M3\tB\ttest\t1\n"
                "M4\tB\ttest\t0\n",
                encoding="utf-8",
            )
            right_path.write_text(
                "M1\tX\ttest\t1\n"
                "M2\tX\ttest\t0\n"
                "M3\tX\ttest\t1\n"
                "M4\tX\ttest\t0\n",
                encoding="utf-8",
            )
            result = module.compare_pair(
                "30",
                "25",
                left_path,
                right_path,
                {
                    "A": {"split": "training", "members": 2, "proteins": 1},
                    "B": {"split": "test", "members": 2, "proteins": 1},
                },
                {"X": {"split": "test", "members": 4, "proteins": 2}},
                {
                    "P1": ("A", "training", "M1"),
                    "P2": ("B", "test", "M3"),
                },
                {
                    "P1": ("X", "test", "M1"),
                    "P2": ("X", "test", "M3"),
                },
            )
            self.assertEqual(result["summary"]["supervised_changed_split"], 1)
            changes = {
                (row["population"], row["pair_type"]): row
                for row in result["pair_changes"]
            }
            self.assertEqual(
                changes[("canonical_supervised_proteins", "labelled_labelled")]["gained_same_cluster_pairs"],
                1,
            )
            self.assertEqual(
                changes[("retained_uniref50_members", "linked_background")]["gained_same_cluster_pairs"],
                2,
            )
            self.assertEqual(result["category_summary"][0]["category"], "merged")
            self.assertEqual(result["destinations"][0]["background_from_other_predecessors"], 1)


if __name__ == "__main__":
    unittest.main()
