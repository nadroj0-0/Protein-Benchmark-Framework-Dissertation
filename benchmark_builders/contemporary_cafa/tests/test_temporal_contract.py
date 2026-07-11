from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from cafa_benchmark_builder.ontology import Ontology
from cafa_benchmark_builder.goa import load_normalized_annotation_map
from cafa_benchmark_builder.models import IdentityMatch, ProteinCatalog, ProteinRecord
from cafa_benchmark_builder.parsers import load_protein_catalog
from cafa_benchmark_builder.snapshot import (
    _build_identity_crosswalk,
    _build_test_dataframe,
    _drop_protein_binding_only,
)


def dat_record(entry: str, accessions: list[str], sequence: str) -> str:
    return (
        f"ID   {entry} Reviewed;\n"
        f"AC   {'; '.join(accessions)};\n"
        "OX   NCBI_TaxID=9606;\n"
        f"SQ   SEQUENCE   {len(sequence)} AA;\n"
        f"     {sequence}\n"
        "//\n"
    )


class TemporalContractTest(unittest.TestCase):
    @staticmethod
    def _catalog(protein_id: str = "P00001") -> ProteinCatalog:
        record = ProteinRecord(
            protein_id=protein_id,
            sequence="MAAAA",
            taxon_id="9606",
            reviewed=True,
            entry_name="P_HUMAN",
            accessions=(protein_id,),
        )
        return ProteinCatalog(
            records={protein_id: record},
            alias_to_primary={protein_id: protein_id},
        )

    def test_secondary_accession_maps_cross_release_identity(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            t0 = root / "t0.dat"
            t1 = root / "t1.dat"
            t0.write_text(dat_record("OLD_HUMAN", ["POLD", "SOLD"], "MAAAA"))
            t1.write_text(dat_record("NEW_HUMAN", ["PNEW", "POLD"], "MAAAA"))

            old = load_protein_catalog((t0,))
            new = load_protein_catalog((t1,))
            matches, reverse = _build_identity_crosswalk(old, new, "exclude")

            self.assertEqual(reverse, {"PNEW": "POLD"})
            self.assertEqual(matches[0].status, "matched")
            self.assertEqual(matches[0].t1_id, "PNEW")

    def test_changed_sequence_is_excluded_but_can_explicitly_use_t0(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            t0 = root / "t0.dat"
            t1 = root / "t1.dat"
            t0.write_text(dat_record("P_HUMAN", ["P00001"], "MAAAA"))
            t1.write_text(dat_record("P_HUMAN", ["P00001"], "MBBBB"))
            old = load_protein_catalog((t0,))
            new = load_protein_catalog((t1,))

            excluded, reverse = _build_identity_crosswalk(old, new, "exclude")
            self.assertEqual(reverse, {})
            self.assertEqual(excluded[0].reason, "sequence_changed")

            retained, reverse = _build_identity_crosswalk(old, new, "use-t0")
            self.assertEqual(reverse, {"P00001": "P00001"})
            self.assertTrue(retained[0].sequence_changed)

    def test_many_to_one_accession_merge_is_not_guessed(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            t0 = root / "t0.dat"
            t1 = root / "t1.dat"
            t0.write_text(
                dat_record("A_HUMAN", ["P00001"], "MAAAA")
                + dat_record("B_HUMAN", ["P00002"], "MBBBB")
            )
            t1.write_text(dat_record("AB_HUMAN", ["PNEW", "P00001", "P00002"], "MCCCC"))
            old = load_protein_catalog((t0,))
            new = load_protein_catalog((t1,))

            matches, reverse = _build_identity_crosswalk(old, new, "use-t0")
            self.assertEqual(reverse, {})
            self.assertEqual(
                {match.reason for match in matches},
                {"many_t0_ids_map_to_one_t1"},
            )

    def test_go_alt_ids_and_single_replacements_are_canonicalised(self):
        with tempfile.TemporaryDirectory() as tmp:
            obo = Path(tmp) / "go.obo"
            obo.write_text(
                "format-version: 1.2\n\n"
                "[Term]\n"
                "id: GO:0000001\n"
                "alt_id: GO:0000009\n"
                "name: live\n"
                "namespace: biological_process\n\n"
                "[Term]\n"
                "id: GO:0000008\n"
                "name: old\n"
                "namespace: biological_process\n"
                "is_obsolete: true\n"
                "replaced_by: GO:0000001\n"
            )
            go = Ontology(obo)
            self.assertEqual(go.resolve_term("GO:0000009"), "GO:0000001")
            self.assertEqual(go.resolve_term("GO:0000008"), "GO:0000001")

    def test_protein_binding_only_mf_label_is_removed_without_dropping_bp(self):
        with tempfile.TemporaryDirectory() as tmp:
            obo = Path(tmp) / "go.obo"
            obo.write_text(
                "format-version: 1.2\n\n"
                "[Term]\n"
                "id: GO:0003674\n"
                "name: molecular_function\n"
                "namespace: molecular_function\n\n"
                "[Term]\n"
                "id: GO:0005515\n"
                "name: protein binding\n"
                "namespace: molecular_function\n"
                "is_a: GO:0003674 ! molecular_function\n\n"
                "[Term]\n"
                "id: GO:0008150\n"
                "name: biological_process\n"
                "namespace: biological_process\n"
            )
            go = Ontology(obo)
            retained = _drop_protein_binding_only(
                {"GO:0005515", "GO:0008150"}, go, "drop-mf-protein-binding-only"
            )
            self.assertEqual(retained, {"GO:0008150"})

    def test_cafa_ontology_policy_keeps_limited_knowledge_target(self):
        go = Ontology(ROOT / "tests" / "fixtures" / "go-mini.obo")
        catalog = self._catalog()
        matches = [IdentityMatch("P00001", "P00001", "matched", "matched", False)]
        t0_annots = {"P00001": {"GO:0005488"}}
        t1_annots = {"P00001": {"GO:0009987", "GO:0005488"}}

        cafa, _ = _build_test_dataframe(
            go, catalog, catalog, matches, {"P00001": "P00001"},
            t0_annots, t1_annots, "keep", "ontology-no-knowledge",
        )
        supervisor, _ = _build_test_dataframe(
            go, catalog, catalog, matches, {"P00001": "P00001"},
            t0_annots, t1_annots, "keep", "global-no-knowledge",
        )

        self.assertEqual(cafa["proteins"].tolist(), ["P00001"])
        self.assertIn("GO:0009987", set(cafa.loc[0, "annotations"]))
        self.assertNotIn("GO:0005488", set(cafa.loc[0, "annotations"]))
        self.assertTrue(supervisor.empty)

    def test_t1_annotation_dates_are_bounded_at_both_ends(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            gaf = root / "annotations.gaf"
            gaf.write_text(
                "!gaf-version: 2.2\n"
                "UniProtKB\tP00001\tP\t\tGO:0009987\tPMID:1\tIDA\t\tP\tProtein\t\tprotein\t"
                "taxon:9606\t20170101\tUniProt\t\t\n"
                "UniProtKB\tP00001\tP\t\tGO:0009987\tPMID:2\tIDA\t\tP\tProtein\t\tprotein\t"
                "taxon:9606\t20170601\tUniProt\t\t\n"
                "UniProtKB\tP00001\tP\t\tGO:0009987\tPMID:3\tIDA\t\tP\tProtein\t\tprotein\t"
                "taxon:9606\t20171201\tUniProt\t\t\n"
            )
            go = Ontology(ROOT / "tests" / "fixtures" / "go-mini.obo")
            result = load_normalized_annotation_map(
                gaf,
                alias_to_primary={"P00001": "P00001"},
                source_ontology=go,
                benchmark_ontology=go,
                exclude_on_or_before="20170213",
                include_on_or_before="20171115",
            )

            self.assertEqual(result.annotations, {"P00001": {"GO:0009987"}})
            self.assertEqual(result.counters["skipped_backfill"], 1)
            self.assertEqual(result.counters["skipped_after_cutoff"], 1)


if __name__ == "__main__":
    unittest.main()
