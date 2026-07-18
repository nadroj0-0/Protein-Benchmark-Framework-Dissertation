from __future__ import annotations

import gzip
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from homology_cluster_benchmark.goa import load_goa
from homology_cluster_benchmark.config import SUPERVISOR_EVIDENCE_CODES
from homology_cluster_benchmark.idmapping import load_uniref90_mappings
from homology_cluster_benchmark.mapping import (
    canonicalize_goa_accessions,
    load_requested_proteins,
)
from homology_cluster_benchmark.ontology import Ontology
from homology_cluster_benchmark.uniref import UniRefIndex, iter_fasta

from tests.helpers import FIXTURES


class ParserTests(unittest.TestCase):
    def test_ontology_alt_obsolete_propagation_and_all_root_retention(self):
        ontology = Ontology(FIXTURES / "go-mini.obo", include_relationships=True)
        self.assertEqual(ontology.data_version, "releases/2026-06-15")
        self.assertEqual(ontology.resolve("GO:9990002"), "GO:0005488")
        self.assertEqual(ontology.resolve("GO:9990001"), "GO:0009987")
        self.assertIsNone(ontology.resolve("GO:9990003"))
        self.assertEqual(
            ontology.ancestors("GO:0006355"),
            {"GO:0006355", "GO:0009987", "GO:0008150"},
        )
        first = ontology.ancestors("GO:0006355")
        first.clear()
        self.assertEqual(
            ontology.ancestors("GO:0006355"),
            {"GO:0006355", "GO:0009987", "GO:0008150"},
        )
        self.assertEqual(
            ontology.ancestors("GO:0005886"), {"GO:0005886", "GO:0005575"}
        )
        self.assertEqual(
            ontology.ancestors("GO:0005488"), {"GO:0005488", "GO:0003674"}
        )
        self.assertIn("part_of", ontology.relationship_types)

    def test_goa_exact_policy_and_diagnostics(self):
        ontology = Ontology(FIXTURES / "go-mini.obo")
        result = load_goa(FIXTURES / "goa.gaf", ontology)
        self.assertEqual(result.counters["kept_rows"], 17)
        self.assertEqual(result.counters["qualifying_proteins"], 17)
        self.assertNotIn("PNOT", result.annotations)
        self.assertNotIn("PLOW", result.annotations)
        self.assertEqual(result.annotations["P2MF"], {"GO:0005488"})
        self.assertEqual(result.annotations["P5BP"], {"GO:0009987"})
        reasons = {item.rejection_reason for item in result.excluded}
        self.assertTrue({
            "not_qualifier", "evidence_code", "isoform_specific", "unknown",
            "namespace_mismatch", "object_type", "database",
        }.issubset(reasons))
        self.assertEqual(set(result.evidence_counts), {
            "EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "HTP", "HDA", "HMP", "HGI",
            "HEP", "TAS", "NAS", "IGC", "RCA", "ND", "IC",
        })
        first = result.records[0]
        self.assertEqual(first.database, "UniProtKB")
        self.assertEqual(first.reference, "PMID:1")
        self.assertEqual(first.assigned_date, "20260617")
        self.assertEqual(first.assigned_by, "UniProt")

    def test_every_supervisor_evidence_code_is_accepted_literally(self):
        ontology = Ontology(FIXTURES / "go-mini.obo")
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for code in sorted(SUPERVISOR_EVIDENCE_CODES):
                with self.subTest(code=code):
                    path = root / f"{code}.gaf"
                    path.write_text(
                        "!gaf-version: 2.2\n"
                        f"UniProtKB\tP{code}\tP{code}\tinvolved_in\tGO:0009987\t"
                        f"PMID:1\t{code}\t\tP\t\t\tprotein\ttaxon:9606\t20260617\tUniProt\t\t\n"
                    )
                    result = load_goa(path, ontology)
                    self.assertEqual(result.counters["kept_rows"], 1)
                    self.assertEqual(set(result.evidence_counts), {code})

    def test_compressed_goa_and_uniref_are_streamed(self):
        ontology = Ontology(FIXTURES / "go-mini.obo")
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            goa_gz = root / "goa.gaf.gz"
            fasta_gz = root / "uniref90.fasta.gz"
            with gzip.open(goa_gz, "wt") as handle:
                handle.write((FIXTURES / "goa.gaf").read_text())
            with gzip.open(fasta_gz, "wt") as handle:
                handle.write((FIXTURES / "uniref90.fasta").read_text())
            self.assertEqual(load_goa(goa_gz, ontology).counters["kept_rows"], 17)
            self.assertEqual(len(list(iter_fasta(fasta_gz))), 7)

    def test_malformed_gaf_is_distinguished_and_strict_mode_fails(self):
        ontology = Ontology(FIXTURES / "go-mini.obo")
        with tempfile.TemporaryDirectory() as tmp:
            malformed = Path(tmp) / "bad.gaf"
            malformed.write_text("!gaf-version: 2.2\nUniProtKB\tP1\n")
            result = load_goa(malformed, ontology, strict_malformed=False)
            self.assertEqual(result.counters["malformed"], 1)
            self.assertEqual(result.excluded[0].rejection_reason, "malformed")
            with self.assertRaisesRegex(ValueError, "malformed"):
                load_goa(malformed, ontology, strict_malformed=True)

            fifteen = Path(tmp) / "fifteen.gaf"
            fifteen.write_text(
                "!gaf-version: 2.2\n"
                "UniProtKB\tP1\tP1\tinvolved_in\tGO:0009987\tPMID:1\tEXP\t\tP\t\t\t"
                "protein\ttaxon:9606\t20260617\tUniProt\n"
            )
            with self.assertRaisesRegex(ValueError, "malformed"):
                load_goa(fifteen, ontology)

    def test_isoform_accession_is_explicitly_excluded_without_suffix_stripping(self):
        ontology = Ontology(FIXTURES / "go-mini.obo")
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "isoform.gaf"
            path.write_text(
                "!gaf-version: 2.2\n"
                "UniProtKB\tP12345-2\tISO\tinvolved_in\tGO:0009987\tPMID:1\tEXP\t\tP\t\t\t"
                "protein\ttaxon:9606\t20260617\tUniProt\t\t\n"
            )
            result = load_goa(path, ontology)
        self.assertNotIn("P12345", result.annotations)
        self.assertEqual(result.excluded[0].rejection_reason, "isoform_accession")

    def test_mapping_chain_reports_ambiguous_and_unmapped(self):
        ontology = Ontology(FIXTURES / "go-mini.obo")
        goa = load_goa(FIXTURES / "goa.gaf", ontology)
        catalog = load_requested_proteins(FIXTURES / "uniprot.fasta", set(goa.annotations))
        canonicalize_goa_accessions(goa, catalog)
        with tempfile.TemporaryDirectory() as tmp:
            index = UniRefIndex.build(FIXTURES / "uniref90.fasta", Path(tmp) / "uniref.sqlite")
            decisions = load_uniref90_mappings(
                FIXTURES / "idmapping_selected.tab", set(record.raw_accession for record in goa.records),
                catalog, index,
            )
        statuses = {item.raw_accession: item.status for item in decisions}
        by_accession = {item.raw_accession: item for item in decisions}
        self.assertEqual(statuses["PAMB"], "ambiguous")
        self.assertIsNone(by_accession["PAMB"].exists_in_fasta)
        self.assertIn("present_in_fasta=UniRef90_U2;UniRef90_U3", by_accession["PAMB"].detail)
        self.assertEqual(statuses["PUN"], "unmapped-absent")
        self.assertIsNone(by_accession["PUN"].exists_in_fasta)
        self.assertEqual(statuses["P1BP"], "mapped")
        self.assertEqual(
            {item.uniref90_id for item in decisions if item.raw_accession.startswith("P1")},
            {"UniRef90_U1"},
        )

    def test_dat_secondary_accession_is_explicitly_canonicalized(self):
        with tempfile.TemporaryDirectory() as tmp:
            dat = Path(tmp) / "mini.dat"
            dat.write_text(
                "ID   PRIMARY_HUMAN Reviewed; 4 AA.\n"
                "AC   PPRIMARY; SSECOND;\n"
                "OX   NCBI_TaxID=9606;\n"
                "SQ   SEQUENCE   4 AA;\n"
                "     MAAA\n"
                "//\n"
            )
            catalog = load_requested_proteins(dat, {"SSECOND"})
        self.assertEqual(catalog.alias_to_primary["SSECOND"], "PPRIMARY")
        self.assertEqual(catalog.records["PPRIMARY"].sequence, "MAAA")

    def test_truncated_dat_record_fails_loudly(self):
        with tempfile.TemporaryDirectory() as tmp:
            dat = Path(tmp) / "truncated.dat"
            dat.write_text(
                "ID   PRIMARY_HUMAN Reviewed; 4 AA.\n"
                "AC   PPRIMARY;\n"
                "SQ   SEQUENCE   4 AA;\n"
                "     MAAA\n"
            )
            with self.assertRaisesRegex(ValueError, "Unterminated"):
                load_requested_proteins(dat, {"PPRIMARY"})

    def test_duplicate_requested_secondary_with_conflicting_sequences_is_ambiguous(self):
        with tempfile.TemporaryDirectory() as tmp:
            dat = Path(tmp) / "conflict.dat"
            dat.write_text(
                "ID   FIRST Reviewed; 4 AA.\nAC   P11111; SSECOND;\n"
                "SQ   SEQUENCE   4 AA;\n     MAAA\n//\n"
                "ID   SECOND Reviewed; 4 AA.\nAC   P22222; SSECOND;\n"
                "SQ   SEQUENCE   4 AA;\n     MCCC\n//\n"
            )
            catalog = load_requested_proteins(dat, {"SSECOND"})
            self.assertEqual(set(catalog.records), {"P11111", "P22222"})
            self.assertIn("SSECOND", catalog.ambiguous_aliases)
            self.assertNotIn("SSECOND", catalog.alias_to_primary)
            self.assertEqual(
                catalog.collision_counts["ambiguous-secondary-conflicting"], 1
            )


if __name__ == "__main__":
    unittest.main()
