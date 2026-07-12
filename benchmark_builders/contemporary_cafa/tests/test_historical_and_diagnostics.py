from __future__ import annotations

# ruff: noqa: E402 - tests add the local src tree before importing the package.

import json
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from cafa_benchmark_builder.builder import build_benchmark
from cafa_benchmark_builder.config import BuildConfig
from cafa_benchmark_builder.goa import load_normalized_annotation_map
from cafa_benchmark_builder.official_targets import load_official_target_catalog
from cafa_benchmark_builder.ontology import Ontology
from cafa_benchmark_builder.parsers import load_protein_catalog


FIXTURES = ROOT / "tests" / "fixtures"


def dat_record(entry: str, accession: str, sequence: str) -> str:
    return (
        f"ID   {entry} Reviewed;\n"
        f"AC   {accession};\n"
        "OX   NCBI_TaxID=9606;\n"
        f"SQ   SEQUENCE   {len(sequence)} AA;\n"
        f"     {sequence}\n"
        "//\n"
    )


def gaf_row(accession: str, go_id: str, date: str = "20260101") -> str:
    return (
        f"UniProtKB\t{accession}\tP\t\t{go_id}\tPMID:1\tIDA\t\tP\tProtein\t\t"
        f"protein\ttaxon:9606\t{date}\tUniProt\t\t\n"
    )


class HistoricalTargetModeTest(unittest.TestCase):
    def test_historical_runner_defaults_to_pre_freeze_release_and_official_targets(self):
        script = (
            ROOT.parents[1] / "scripts" / "validation"
            / "run_cafa3_historical_validation.sh"
        ).read_text()
        self.assertIn(
            'HISTORICAL_TRAINING_SNAPSHOT="${HISTORICAL_TRAINING_SNAPSHOT:-september-2016}"',
            script,
        )
        self.assertIn(
            'TARGET_UNIVERSE_POLICY="${TARGET_UNIVERSE_POLICY:-official-cafa3-targets}"',
            script,
        )
        self.assertIn("release-2016_08", script)
        self.assertIn("bad05b535790955ddd5e0d1833915f9f", script)

    def test_mapping_failure_preserves_official_target(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            reference = root / "reference.dat"
            reference.write_text(dat_record("KNOWN_HUMAN", "P00001", "MAAAA"))
            fasta = root / "targets.fasta"
            fasta.write_text(">T960600000001 UNKNOWN_HUMAN\nMZZZZ\n")
            mappings = root / "mappings"
            mappings.mkdir()
            (mappings / "sp_species.9606.map").write_text(
                "T960600000001\tUNKNOWN_HUMAN\n"
            )

            result = load_official_target_catalog(
                (fasta,), mappings, load_protein_catalog((reference,)),
                frozenset({"9606"}), "t0",
            )

            self.assertIn("T960600000001", result.catalog.records)
            self.assertEqual(result.catalog.records["T960600000001"].sequence, "MZZZZ")
            self.assertEqual(result.rows[0]["status"], "unmapped")
            self.assertEqual(result.rows[0]["reason"], "no_uniprot_mapping")
            self.assertEqual(result.rows[0]["present_in_snapshot"], 0)

    def test_official_target_id_and_training_snapshot_reach_outputs_and_manifest(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            uniprot = root / "uniprot.dat"
            uniprot.write_text(
                dat_record("TRAIN_HUMAN", "PTRAIN", "MAAAAA")
                + dat_record("TARGET_HUMAN", "PTEST", "MBBBBB")
            )
            targets = root / "targets.fasta"
            targets.write_text(">T960600000001 TARGET_HUMAN\nMBBBBB\n")
            mappings = root / "mappings"
            mappings.mkdir()
            (mappings / "sp_species.9606.map").write_text(
                "T960600000001\tTARGET_HUMAN\n"
            )
            released_labels = root / "uniprot_sprot_exp.txt"
            released_labels.write_text(
                "PTRAIN\tGO:0009987\tP\n"
                "PTRAIN\tGO:0005488\tF\n"
                "PTRAIN\tGO:0005886\tC\n"
            )
            t0_gaf = root / "t0.gaf"
            t0_gaf.write_text("!gaf-version: 2.2\n")
            t1_gaf = root / "t1.gaf"
            t1_gaf.write_text("!gaf-version: 2.2\n" + gaf_row("PTEST", "GO:0009987"))
            output = root / "output"
            reports = root / "reports"

            written = build_benchmark(BuildConfig(
                uniprot_t0=(uniprot,),
                uniprot_t1=(uniprot,),
                target_uniprot_t0=(uniprot,),
                target_uniprot_t1=(uniprot,),
                goa_t0=t0_gaf,
                goa_t1=t1_gaf,
                go_obo=FIXTURES / "go-mini.obo",
                output_dir=output,
                report_dir=reports,
                target_universe_policy="official-cafa3-targets",
                official_target_fastas=(targets,),
                official_target_mapping_dir=mappings,
                training_annotations_file=released_labels,
                training_snapshot_id="UniProtKB-2016_08",
                training_snapshot_date="07-Sep-2016",
                profile_name="cafa3-reconstructed",
                target_taxa=frozenset({"9606"}),
                min_count=1,
                split=0.5,
                t0_cutoff="20250131",
                t1_cutoff="20260131",
                strict_qc=False,
                write_checksums=False,
            ))

            test = pd.read_csv(written["bp-test"])
            self.assertEqual(test["proteins"].tolist(), ["T960600000001"])
            mapping_rows = (reports / "official_target_mapping.tsv").read_text()
            self.assertIn("source-and-exact-sequence", mapping_rows)
            manifest = json.loads((reports / "build_manifest.json").read_text())
            self.assertEqual(manifest["training_snapshot_id"], "UniProtKB-2016_08")
            self.assertEqual(manifest["target_universe_policy"], "official-cafa3-targets")


class OntologyDiagnosticTest(unittest.TestCase):
    def test_frozen_graph_resolves_nearest_source_snapshot_mismatch(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            frozen_path = root / "frozen.obo"
            frozen_path.write_text(
                "format-version: 1.2\ndata-version: releases/2025-02-06\n\n"
                "[Term]\nid: GO:1234567\nname: old live term\n"
                "namespace: biological_process\n"
            )
            source_path = root / "source.obo"
            source_path.write_text(
                "format-version: 1.2\ndata-version: releases/2025-03-16\n\n"
                "[Term]\nid: GO:1234567\nname: obsolete term\n"
                "namespace: biological_process\nis_obsolete: true\n"
                "consider: GO:7654321\nconsider: GO:1111111\n"
            )
            gaf = root / "annotations.gaf"
            gaf.write_text("!gaf-version: 2.2\n" + gaf_row("P00001", "GO:1234567"))
            frozen = Ontology(frozen_path)
            source = Ontology(source_path)

            handled = load_normalized_annotation_map(
                gaf,
                alias_to_primary={"P00001": "P00001"},
                source_ontology=source,
                benchmark_ontology=frozen,
                snapshot="t0",
                allow_frozen_source_fallback=True,
            )
            rejected = load_normalized_annotation_map(
                gaf,
                alias_to_primary={"P00001": "P00001"},
                source_ontology=source,
                benchmark_ontology=frozen,
                snapshot="t0",
                allow_frozen_source_fallback=False,
            )

            self.assertEqual(handled.annotations, {"P00001": {"GO:1234567"}})
            self.assertEqual(handled.source_diagnostics[0]["final_action"], "use_frozen_term")
            self.assertEqual(handled.source_diagnostics[0]["consider"], "GO:7654321|GO:1111111")
            self.assertEqual(rejected.annotations, {})
            self.assertEqual(rejected.unmapped_terms, {"GO:1234567": 1})
            self.assertEqual(rejected.source_diagnostics[0]["final_action"], "fail_strict_qc")

    def test_valid_t1_only_term_is_reported_outside_frozen_graph(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source_path = root / "t1.obo"
            source_path.write_text(
                "format-version: 1.2\ndata-version: releases/2026-06-15\n\n"
                "[Term]\nid: GO:1234567\nname: new term\nnamespace: biological_process\n"
            )
            gaf = root / "annotations.gaf"
            gaf.write_text("!gaf-version: 2.2\n" + gaf_row("P00001", "GO:1234567"))
            result = load_normalized_annotation_map(
                gaf,
                alias_to_primary={"P00001": "P00001"},
                source_ontology=Ontology(source_path),
                benchmark_ontology=Ontology(FIXTURES / "go-mini.obo"),
                snapshot="t1",
            )

            self.assertEqual(result.annotations, {})
            self.assertFalse(result.unmapped_terms)
            self.assertEqual(result.out_of_benchmark_terms, {"GO:1234567": 1})
            self.assertEqual(
                result.outside_frozen_diagnostics[0]["final_action"],
                "exclude_from_frozen_label_space",
            )

    def test_strict_failure_persists_unresolved_row_and_preflight_report(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            t0_gaf = root / "t0.gaf"
            t0_gaf.write_text(
                "!gaf-version: 2.2\n" + gaf_row("P00001", "GO:9999999", "20250101")
            )
            t1_gaf = root / "t1.gaf"
            t1_gaf.write_text("!gaf-version: 2.2\n")
            reports = root / "reports"
            config = BuildConfig(
                uniprot_t0=(FIXTURES / "uniprot-t0.fasta",),
                uniprot_t1=(FIXTURES / "uniprot-t1.fasta",),
                goa_t0=t0_gaf,
                goa_t1=t1_gaf,
                go_obo=FIXTURES / "go-mini.obo",
                output_dir=root / "output",
                report_dir=reports,
                training_taxa=frozenset({"9606"}),
                target_taxa=frozenset({"9606"}),
                t0_cutoff="20250131",
                t1_cutoff="20260131",
                strict_qc=True,
                write_checksums=False,
            )

            with self.assertRaisesRegex(ValueError, "cannot be resolved"):
                build_benchmark(config)

            unresolved = reports / "unresolved_source_go_annotations.tsv"
            self.assertTrue(unresolved.is_file())
            self.assertIn("GO:9999999", unresolved.read_text())
            report = (reports / "benchmark_build_report.md").read_text()
            self.assertIn("GO:9999999", report)
            self.assertIn("strict-QC failures: 1", report)


if __name__ == "__main__":
    unittest.main()
