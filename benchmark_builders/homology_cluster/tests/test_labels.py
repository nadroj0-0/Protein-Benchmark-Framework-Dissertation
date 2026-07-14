from __future__ import annotations

from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from homology_cluster_benchmark.labels import build_labels
from homology_cluster_benchmark.models import (
    AnnotationRecord,
    GoaLoadResult,
    MappingDecision,
    ProteinCatalog,
    ProteinRecord,
    SplitAssignment,
)
from homology_cluster_benchmark.ontology import Ontology

from tests.helpers import FIXTURES


def annotation(raw: str, protein: str, go_id: str) -> AnnotationRecord:
    return AnnotationRecord(
        database="UniProtKB", raw_accession=raw, protein_id=protein, symbol=raw,
        raw_go_id=go_id, go_id=go_id, namespace="biological_process", aspect="P",
        evidence="EXP", qualifier="involved_in", reference="PMID:1", with_from="",
        taxon_id="9606", assigned_date="20260617", assigned_by="UniProt",
        annotation_extension="", gene_product_form="", line_number=1,
        term_action="primary",
    )


class LabelTests(unittest.TestCase):
    def test_term_universe_uses_complete_development_and_retains_roots(self):
        ontology = Ontology(FIXTURES / "go-mini.obo")
        goa = GoaLoadResult(
            records=[
                annotation("PTRAIN", "PTRAIN", "GO:0009987"),
                annotation("PVALID", "PVALID", "GO:0006355"),
                annotation("PTEST", "PTEST", "GO:0006355"),
            ], excluded=[],
            annotations={
                "PTRAIN": {"GO:0009987"},
                "PVALID": {"GO:0006355"},
                "PTEST": {"GO:0006355"},
            },
        )
        catalog = ProteinCatalog(records={
            "PTRAIN": ProteinRecord("PTRAIN", "MAAA"),
            "PVALID": ProteinRecord("PVALID", "MBBB"),
            "PTEST": ProteinRecord("PTEST", "MCCC"),
        })
        mappings = [
            MappingDecision(protein, protein, "exact", f"UniRef90_{cluster}", "mapped", "", True, True, cluster)
            for protein, cluster in (("PTRAIN", "C1"), ("PVALID", "C2"), ("PTEST", "C3"))
        ]
        assignments = {
            "C1": SplitAssignment("C1", "training", 1, 1, "test"),
            "C2": SplitAssignment("C2", "validation", 1, 1, "test"),
            "C3": SplitAssignment("C3", "test", 1, 1, "test"),
        }
        labels = build_labels(ontology, goa, catalog, mappings, assignments, min_count=1)
        self.assertEqual(
            labels.term_universe, ("GO:0006355", "GO:0008150", "GO:0009987")
        )
        self.assertEqual(
            labels.frames["validation"].iloc[0].annotations,
            ("GO:0006355", "GO:0008150", "GO:0009987"),
        )
        self.assertEqual(
            labels.restricted_annotations["PVALID"],
            ("GO:0006355", "GO:0008150", "GO:0009987"),
        )
        self.assertIn("GO:0008150", labels.term_universe)

    def test_ambiguous_raw_alias_cannot_leak_its_annotation_to_mapped_primary(self):
        ontology = Ontology(FIXTURES / "go-mini.obo")
        goa = GoaLoadResult(
            records=[
                annotation("PPRIMARY", "PPRIMARY", "GO:0009987"),
                annotation("SSECOND", "PPRIMARY", "GO:0006355"),
            ],
            excluded=[],
            annotations={"PPRIMARY": {"GO:0009987", "GO:0006355"}},
            candidate_accessions={"PPRIMARY", "SSECOND"},
        )
        catalog = ProteinCatalog(records={
            "PPRIMARY": ProteinRecord("PPRIMARY", "MAAA")
        })
        mappings = [
            MappingDecision(
                "PPRIMARY", "PPRIMARY", "exact", "UniRef90_U1", "mapped", "",
                True, True, "C1", "training",
            ),
            MappingDecision(
                "SSECOND", "PPRIMARY", "secondary-to-primary", status="ambiguous",
                detail="multiple UniRef90 mappings",
            ),
        ]
        assignments = {
            "C1": SplitAssignment("C1", "training", 1, 1, "fixture")
        }
        labels = build_labels(ontology, goa, catalog, mappings, assignments, min_count=1)
        self.assertEqual(
            labels.unrestricted_annotations["PPRIMARY"], ("GO:0008150", "GO:0009987")
        )
        self.assertEqual(labels.annotation_exclusion_counts["mapping_status:ambiguous"], 1)
        self.assertEqual(sum(labels.row_attrition_counts.values()), 2)
        self.assertEqual(sum(labels.protein_attrition_counts.values()), 2)
        self.assertEqual(labels.protein_attrition_counts["mapping_status:ambiguous"], 1)

    def test_same_cluster_proteins_keep_only_their_own_disjoint_labels(self):
        ontology = Ontology(FIXTURES / "go-mini.obo")
        goa = GoaLoadResult(
            records=[
                annotation("PBP", "PBP", "GO:0006355"),
                annotation("PCC", "PCC", "GO:0005886"),
            ],
            excluded=[],
            annotations={"PBP": {"GO:0006355"}, "PCC": {"GO:0005886"}},
            candidate_accessions={"PBP", "PCC"},
        )
        catalog = ProteinCatalog(records={
            "PBP": ProteinRecord("PBP", "MAAA"),
            "PCC": ProteinRecord("PCC", "MBBB"),
        })
        mappings = [
            MappingDecision(
                protein, protein, "exact", "UniRef90_SHARED", "mapped", "", True,
                True, "C_SHARED", "training",
            )
            for protein in ("PBP", "PCC")
        ]
        assignments = {
            "C_SHARED": SplitAssignment("C_SHARED", "training", 2, 2, "fixture")
        }
        labels = build_labels(ontology, goa, catalog, mappings, assignments, min_count=1)
        self.assertIn("GO:0006355", labels.unrestricted_annotations["PBP"])
        self.assertNotIn("GO:0005886", labels.unrestricted_annotations["PBP"])
        self.assertIn("GO:0005886", labels.unrestricted_annotations["PCC"])
        self.assertNotIn("GO:0006355", labels.unrestricted_annotations["PCC"])

    def test_fully_rejected_and_missing_sequence_proteins_have_distinct_outcomes(self):
        ontology = Ontology(FIXTURES / "go-mini.obo")
        goa = GoaLoadResult(
            records=[annotation("PMISSING", "PMISSING", "GO:0009987")],
            excluded=[],
            annotations={"PMISSING": {"GO:0009987"}},
            candidate_accessions={"PMISSING", "PREJECTED"},
        )
        mappings = [MappingDecision(
            "PMISSING", "PMISSING", "exact", "UniRef90_M", "mapped", "", True,
            False, "C1", "training",
        )]
        labels = build_labels(
            ontology,
            goa,
            ProteinCatalog(records={}),
            mappings,
            {"C1": SplitAssignment("C1", "training", 1, 1, "fixture")},
            min_count=1,
        )
        self.assertEqual(labels.protein_attrition_counts["missing_canonical_sequence"], 1)
        self.assertEqual(
            labels.protein_attrition_counts["ontology_resolution_or_namespace_rejection"], 1
        )
        self.assertEqual(sum(labels.protein_attrition_counts.values()), 2)

    def test_min_count_uses_49_plus_1_development_but_never_test(self):
        ontology = Ontology(FIXTURES / "go-mini.obo")
        records = []
        catalog_records = {}
        mappings = []
        assignments = {}
        annotations = {}
        for index in range(100):
            protein = f"P{index:03d}"
            cluster = f"C{index:03d}"
            if index < 49:
                split = "training"
            elif index == 49:
                split = "validation"
            else:
                split = "test"
            terms = ["GO:0006355"] if index < 50 else ["GO:0005488"]
            if index < 49:
                terms.append("GO:0005886")
            annotations[protein] = set(terms)
            for term in terms:
                namespace = ontology.namespace(term)
                aspect = {
                    "biological_process": "P",
                    "cellular_component": "C",
                    "molecular_function": "F",
                }[namespace]
                records.append(AnnotationRecord(
                    database="UniProtKB", raw_accession=protein, protein_id=protein,
                    symbol=protein, raw_go_id=term, go_id=term, namespace=namespace,
                    aspect=aspect, evidence="EXP", qualifier="", reference="PMID:1",
                    with_from="", taxon_id="9606", assigned_date="20260617",
                    assigned_by="UniProt", annotation_extension="", gene_product_form="",
                    line_number=index + 1, term_action="primary",
                ))
            catalog_records[protein] = ProteinRecord(protein, "MAAA")
            mappings.append(MappingDecision(
                raw_accession=protein, protein_id=protein, accession_action="exact",
                uniref90_id=f"UniRef90_{cluster}", status="mapped", exists_in_fasta=True,
                canonical_sequence_available=True, mmseqs_cluster_id=cluster, split=split,
            ))
            assignments[cluster] = SplitAssignment(cluster, split, 1, 1, "fixture")
        labels = build_labels(
            ontology,
            GoaLoadResult(
                records=records, excluded=[], annotations=annotations,
                candidate_accessions=set(annotations),
            ),
            ProteinCatalog(records=catalog_records),
            mappings,
            assignments,
            min_count=50,
        )
        self.assertIn("GO:0006355", labels.term_universe)
        self.assertIn("GO:0008150", labels.term_universe)
        self.assertNotIn("GO:0005886", labels.term_universe)
        self.assertNotIn("GO:0005575", labels.term_universe)
        self.assertNotIn("GO:0005488", labels.term_universe)
        self.assertNotIn("GO:0003674", labels.term_universe)
        self.assertEqual(list(labels.term_universe), sorted(labels.term_universe))


if __name__ == "__main__":
    unittest.main()
