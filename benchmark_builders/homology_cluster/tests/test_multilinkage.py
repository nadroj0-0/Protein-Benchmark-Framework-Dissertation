from __future__ import annotations

import gzip
from pathlib import Path
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from homology_cluster_benchmark.multilinkage import build_multilinkage_index


class MultilinkageIndexTests(unittest.TestCase):
    def test_unordered_rows_duplicates_and_uniparc_members_are_preserved(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "cliques.csv.gz"
            with gzip.open(source, "wt", encoding="utf-8", newline="") as handle:
                handle.write(
                    "2,UniRef50_B,P2\n"
                    "1,UniRef50_A,UPI0001\n"
                    "2,UniRef50_B,UPI0002\n"
                    "1,UniRef50_A,P1\n"
                    "1,UniRef50_A,P1\n"
                )
            index, stats = build_multilinkage_index(source, root / "cliques.sqlite")
            self.assertEqual(stats.raw_rows, 5)
            self.assertEqual(stats.duplicate_rows, 1)
            self.assertEqual(stats.members, 4)
            self.assertEqual(stats.clusters, 2)
            self.assertEqual(stats.uniparc_members, 2)
            self.assertEqual(index.cluster_for("P1"), "1")
            self.assertEqual(
                list(index.iter_assignments()),
                [("1", "P1"), ("1", "UPI0001"), ("2", "P2"), ("2", "UPI0002")],
            )

    def test_member_cannot_appear_in_two_cliques(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "cliques.csv"
            source.write_text("1,R1,P1\n2,R2,P1\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "members_in_multiple_clusters=1"):
                build_multilinkage_index(source, root / "cliques.sqlite")

    def test_clique_cannot_have_two_representatives(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "cliques.csv"
            source.write_text("1,R1,P1\n1,R2,P2\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "clusters_with_multiple_representatives=1"):
                build_multilinkage_index(source, root / "cliques.sqlite")

    def test_malformed_row_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "cliques.csv"
            source.write_text("1,R1\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Malformed multi-linkage row"):
                build_multilinkage_index(source, root / "cliques.sqlite")

            source.write_text("cluster,R1,P1\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "cluster_id must be numeric"):
                build_multilinkage_index(source, root / "cliques.sqlite")


if __name__ == "__main__":
    unittest.main()
