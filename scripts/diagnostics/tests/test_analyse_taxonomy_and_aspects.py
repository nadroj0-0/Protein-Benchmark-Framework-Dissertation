from __future__ import annotations

import csv
from pathlib import Path
import sys
import tempfile
import unittest


DIAGNOSTICS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(DIAGNOSTICS))

import analyse_taxonomy_and_aspects as module  # noqa: E402


class TaxonomyAndAspectTests(unittest.TestCase):
    def make_taxonomy(self, root: Path) -> None:
        (root / "nodes.dmp").write_text(
            "1\t|\t1\t|\tno rank\t|\n"
            "131567\t|\t1\t|\tno rank\t|\n"
            "2\t|\t131567\t|\tsuperkingdom\t|\n"
            "2759\t|\t131567\t|\tsuperkingdom\t|\n"
            "9606\t|\t2759\t|\tspecies\t|\n"
            "83333\t|\t2\t|\tstrain\t|\n",
            encoding="utf-8",
        )
        (root / "names.dmp").write_text(
            "1\t|\troot\t|\t\t|\tscientific name\t|\n"
            "2\t|\tBacteria\t|\t\t|\tscientific name\t|\n"
            "2759\t|\tEukaryota\t|\t\t|\tscientific name\t|\n"
            "9606\t|\tHomo sapiens\t|\t\t|\tscientific name\t|\n"
            "83333\t|\tEscherichia coli K-12\t|\t\t|\tscientific name\t|\n",
            encoding="utf-8",
        )
        (root / "merged.dmp").write_text("999\t|\t9606\t|\n", encoding="utf-8")
        (root / "delnodes.dmp").write_text("888\t|\n", encoding="utf-8")

    def test_resolves_domains_and_lifecycle(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            self.make_taxonomy(root)
            taxonomy = module.Taxonomy(root)
            self.assertEqual(taxonomy.resolve("9606")["domain"], "Eukaryota")
            self.assertEqual(taxonomy.resolve("83333")["domain"], "Bacteria")
            merged = taxonomy.resolve("999")
            self.assertEqual(merged["action"], "merged")
            self.assertEqual(merged["resolved_taxid"], "9606")
            self.assertEqual(taxonomy.resolve("888")["domain"], "Unresolved")

    def test_aspect_pattern_is_ordered(self) -> None:
        memberships = {"BPO": {"P1"}, "CCO": set(), "MFO": {"P1"}}
        self.assertEqual(module.pattern_for("P1", memberships), "BPO+MFO")

    def test_temporal_alias_conflict_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            path = Path(raw) / "flow.tsv"
            with path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle, fieldnames=["t0_id", "t1_id", "taxon_id"], delimiter="\t"
                )
                writer.writeheader()
                writer.writerow({"t0_id": "P1", "t1_id": "P2", "taxon_id": "9606"})
                writer.writerow({"t0_id": "P1", "t1_id": "P3", "taxon_id": "83333"})
            with self.assertRaisesRegex(ValueError, "Conflicting taxa"):
                module.temporal_taxa(path)


if __name__ == "__main__":
    unittest.main()
