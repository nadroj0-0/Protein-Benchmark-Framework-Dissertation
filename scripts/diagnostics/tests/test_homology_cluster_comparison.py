from __future__ import annotations

import gzip
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest


SCRIPT = Path(__file__).resolve().parents[1] / "compare_homology_cluster_assignments.py"
SPEC = importlib.util.spec_from_file_location("compare_homology_clusters", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def _write(path: Path, rows: list[tuple[str, str]]) -> None:
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        for cluster, member in rows:
            handle.write(f"{cluster}\t{member}\n")


class HomologyClusterComparisonTests(unittest.TestCase):
    def _args(self, root: Path, left: Path, right: Path):
        return MODULE.argparse.Namespace(
            left=left,
            right=right,
            left_label="left",
            right_label="right",
            output_dir=root / "out",
            scratch_dir=root / "scratch",
            sort_binary="sort",
            sort_parallel=1,
            sort_memory="64M",
        )

    def test_detects_label_independent_exact_partitions(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            left = root / "left.tsv.gz"
            right = root / "right.tsv.gz"
            _write(left, [("a", "a"), ("a", "b"), ("c", "c")])
            _write(right, [("b", "a"), ("b", "b"), ("c", "c")])
            payload = MODULE.compare(self._args(root, left, right))
            self.assertTrue(payload["partitions_exactly_identical"])
            self.assertEqual(payload["exact_partition_blocks"]["exact_cluster_blocks"], 2)
            self.assertEqual(payload["exact_partition_blocks"]["members_in_exact_cluster_blocks"], 3)
            self.assertEqual(payload["member_universe"]["members_with_same_raw_representative"], 1)

    def test_reports_splits_and_merges(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            left = root / "left.tsv.gz"
            right = root / "right.tsv.gz"
            _write(left, [("a", "a"), ("a", "b"), ("c", "c"), ("d", "d"), ("d", "e")])
            _write(right, [("b", "a"), ("b", "b"), ("c", "c"), ("c", "d"), ("e", "e")])
            payload = MODULE.compare(self._args(root, left, right))
            self.assertFalse(payload["partitions_exactly_identical"])
            self.assertEqual(payload["left_partition"]["divergent_clusters"], 1)
            self.assertEqual(payload["right_partition"]["divergent_clusters"], 1)
            self.assertEqual(payload["exact_partition_blocks"]["exact_cluster_blocks"], 1)
            self.assertEqual(payload["exact_partition_blocks"]["members_in_exact_cluster_blocks"], 2)
            parsed = json.loads((root / "out" / "comparison.json").read_text())
            self.assertIn("adjusted_rand_index", parsed["pairwise_agreement"])

    def test_rejects_duplicate_members(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            left = root / "left.tsv.gz"
            right = root / "right.tsv.gz"
            _write(left, [("a", "a"), ("b", "a")])
            _write(right, [("a", "a")])
            with self.assertRaisesRegex(ValueError, "duplicate member"):
                MODULE.compare(self._args(root, left, right))

    def test_all_singletons_have_defined_report(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            left = root / "left.tsv.gz"
            right = root / "right.tsv.gz"
            _write(left, [("a", "a"), ("b", "b")])
            _write(right, [("a", "a"), ("b", "b")])
            payload = MODULE.compare(self._args(root, left, right))
            self.assertTrue(payload["partitions_exactly_identical"])
            self.assertIsNone(payload["pairwise_agreement"]["pair_jaccard"])
            self.assertIn("Pair Jaccard: n/a", (root / "out" / "summary.md").read_text())


if __name__ == "__main__":
    unittest.main()
