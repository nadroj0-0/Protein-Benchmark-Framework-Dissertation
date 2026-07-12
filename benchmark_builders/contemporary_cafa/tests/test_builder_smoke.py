from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path
import json

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from cafa_benchmark_builder.builder import (
    build_benchmark,
    export_from_deepgoplus_pickles,
    generate_deepgoplus_pickles_from_cafa_files,
)
from cafa_benchmark_builder.config import BENCHMARK_PROFILES, BuildConfig, EVIDENCE_POLICIES
from cafa_benchmark_builder.goa import load_annotation_map
from cafa_benchmark_builder.parsers import iter_uniprot


FIXTURES = Path(__file__).resolve().parent / "fixtures"


class BenchmarkBuilderSmokeTest(unittest.TestCase):
    def test_parsers_keep_cafa3_final_evidence_and_remove_not(self):
        proteins = {rec.protein_id: rec for rec in iter_uniprot(FIXTURES / "uniprot-t0.fasta")}
        self.assertEqual(proteins["P00001"].taxon_id, "9606")
        self.assertTrue(proteins["P00001"].reviewed)

        annots = load_annotation_map(FIXTURES / "goa-t0.gaf", target_taxa=frozenset({"9606"}))
        self.assertIn("GO:0005488", annots["P00002"])  # TAS kept by final CAFA3 policy
        self.assertIn("GO:0005886", annots["P00003"])  # IC kept by final CAFA3 policy
        self.assertNotIn("P00004", annots)             # NOT removed

        filtered = load_annotation_map(
            FIXTURES / "goa-t0.gaf",
            target_taxa=frozenset({"9606"}),
            allowed_proteins={"P00002"},
        )
        self.assertEqual(set(filtered), {"P00002"})

    def test_builds_nine_pfp_compatible_csvs(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            written = build_benchmark(BuildConfig(
                uniprot_t0=(FIXTURES / "uniprot-t0.fasta",),
                uniprot_t1=(FIXTURES / "uniprot-t1.fasta",),
                goa_t0=FIXTURES / "goa-t0.gaf",
                goa_t1=FIXTURES / "goa-t1.gaf",
                go_obo=FIXTURES / "go-mini.obo",
                output_dir=out,
                training_taxa=frozenset({"9606"}),
                target_taxa=frozenset({"9606"}),
                min_count=1,
                t0_cutoff="20250131",
                write_checksums=False,
                strict_qc=False,
            ))

            expected = {
                "bp-training.csv", "bp-validation.csv", "bp-test.csv",
                "cc-training.csv", "cc-validation.csv", "cc-test.csv",
                "mf-training.csv", "mf-validation.csv", "mf-test.csv",
            }
            self.assertTrue(expected.issubset({p.name for p in written.values()}))
            for name in expected:
                df = pd.read_csv(out / name)
                self.assertIn("proteins", df.columns)
                self.assertIn("sequences", df.columns)

            bp_test = pd.read_csv(out / "bp-test.csv")
            self.assertEqual(set(bp_test["proteins"]), {"P00004"})
            self.assertEqual(bp_test.loc[0, "sequences"], "MEEEEE")  # always the t0 sequence
            self.assertIn("GO:0009987", bp_test.columns)
            self.assertIn("GO:0008150", bp_test.columns)

            flow = pd.read_csv(out / "reports" / "protein_flow.tsv", sep="\t")
            future_only = flow.loc[flow["t1_id"] == "P00005"].iloc[0]
            self.assertEqual(future_only["reason"], "not_present_at_t0")
            self.assertNotIn("P00005", set(bp_test["proteins"]))
            stats = json.loads((out / "reports" / "benchmark_statistics.json").read_text())
            self.assertEqual(stats["t1_goa"]["skipped_backfill"], 1)

            all_train_valid = []
            for name in ["bp-training.csv", "bp-validation.csv", "cc-training.csv",
                         "cc-validation.csv", "mf-training.csv", "mf-validation.csv"]:
                all_train_valid.extend(pd.read_csv(out / name)["proteins"].tolist())
            self.assertIn("P00002", all_train_valid)  # TAS survived
            self.assertIn("P00003", all_train_valid)  # IC survived
            self.assertNotIn("P00004", all_train_valid)

    def test_exports_from_deepgoplus_pickles(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "deepgoplus"
            out = root / "out"
            source.mkdir()

            pd.DataFrame({
                "proteins": ["P00001", "P00002"],
                "sequences": ["MAAA", "MBBB"],
                "annotations": [{"GO:0009987", "GO:0008150"}, {"GO:0005488", "GO:0003674"}],
            }).to_pickle(source / "train_data.pkl")
            pd.DataFrame({
                "proteins": ["P00001"],
                "sequences": ["MAAA"],
                "annotations": [{"GO:0009987", "GO:0008150"}],
            }).to_pickle(source / "train_data_train.pkl")
            pd.DataFrame({
                "proteins": ["P00002"],
                "sequences": ["MBBB"],
                "annotations": [{"GO:0005488", "GO:0003674"}],
                "preds": [None],
            }).to_pickle(source / "train_data_valid.pkl")
            pd.DataFrame({
                "proteins": ["T96060000001"],
                "sequences": ["MCCC"],
                "annotations": [{"GO:0005886", "GO:0005575"}],
            }).to_pickle(source / "test_data.pkl")
            pd.DataFrame({
                "terms": ["GO:0009987", "GO:0008150", "GO:0005488", "GO:0003674", "GO:0005886", "GO:0005575"],
            }).to_pickle(source / "terms.pkl")

            written = export_from_deepgoplus_pickles(
                deepgoplus_dir=source,
                go_obo=FIXTURES / "go-mini.obo",
                output_dir=out,
            )

            expected = {
                "bp-training.csv", "bp-validation.csv", "bp-test.csv",
                "cc-training.csv", "cc-validation.csv", "cc-test.csv",
                "mf-training.csv", "mf-validation.csv", "mf-test.csv",
                "train_data.pkl", "train_data_train.pkl", "train_data_valid.pkl",
                "test_data.pkl", "terms.pkl",
            }
            self.assertTrue(expected.issubset({p.name for p in written.values()}))
            bp_train = pd.read_csv(out / "bp-training.csv")
            self.assertEqual(bp_train["proteins"].tolist(), ["P00001"])
            cc_test = pd.read_csv(out / "cc-test.csv")
            self.assertEqual(cc_test["proteins"].tolist(), ["T96060000001"])

    def test_generates_deepgoplus_pickles_from_cafa_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "out"
            train_fasta = root / "train.fasta"
            train_annots = root / "train.tsv"
            test_fasta = root / "test.fasta"
            test_annots = root / "test.tsv"

            train_fasta.write_text(">P00001\nMAAA\n>P00002 extra text\nMBBB\n>P99999\nMXXX\n")
            train_annots.write_text("P00001\tGO:0009987\tP\nP00002\tGO:0005488\tF\n")
            test_fasta.write_text(">T96060000001 1433B_MOUSE\nMCCC\n>T00000000002\nMDDD\n")
            test_annots.write_text("T96060000001\tGO:0005886\n")

            written = generate_deepgoplus_pickles_from_cafa_files(
                go_obo=FIXTURES / "go-mini.obo",
                train_sequences_file=train_fasta,
                train_annotations_file=train_annots,
                test_sequences_file=test_fasta,
                test_annotations_file=test_annots,
                output_dir=out,
                min_count=1,
            )

            self.assertEqual({p.name for p in written.values()}, {"train_data.pkl", "test_data.pkl", "terms.pkl"})
            train_df = pd.read_pickle(out / "train_data.pkl")
            test_df = pd.read_pickle(out / "test_data.pkl")
            terms_df = pd.read_pickle(out / "terms.pkl")

            self.assertEqual(train_df["proteins"].tolist(), ["P00001", "P00002"])
            self.assertEqual(test_df["proteins"].tolist(), ["T96060000001"])
            self.assertEqual(set(train_df.loc[0, "annotations"]), {"GO:0009987", "GO:0008150"})
            self.assertEqual(set(test_df.loc[0, "annotations"]), {"GO:0005886", "GO:0005575"})
            self.assertEqual(
                set(terms_df["terms"]),
                {"GO:0009987", "GO:0008150", "GO:0005488", "GO:0003674"},
            )

    def test_snapshot_build_can_use_released_cafa3_test_groundtruth(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "out"
            mapping_dir = root / "mappings"
            mapping_dir.mkdir()
            training_annotations = root / "training.tsv"
            test_fasta = root / "targets.fasta"
            test_annotations = root / "leafonly_all.txt"
            training_annotations.write_text(
                "P00001\tGO:0009987\tP\n"
                "P00002\tGO:0005488\tF\n"
                "P00003\tGO:0005886\tC\n"
            )
            test_fasta.write_text(">T96060000001 TARGET_HUMAN\nMCCC\n")
            test_annotations.write_text("T96060000001\tGO:0005886\n")

            written = build_benchmark(BuildConfig(
                uniprot_t0=(FIXTURES / "uniprot-t0.fasta",),
                uniprot_t1=(),
                goa_t0=None,
                goa_t1=None,
                go_obo=FIXTURES / "go-mini.obo",
                output_dir=out,
                target_universe_policy="official-cafa3-targets",
                official_target_fastas=(test_fasta,),
                official_target_mapping_dir=mapping_dir,
                training_annotations_file=training_annotations,
                test_annotations_file=test_annotations,
                profile_name="cafa3-reconstructed",
                target_taxa=frozenset({"9606"}),
                min_count=1,
                write_checksums=False,
                strict_qc=False,
            ))

            test_df = pd.read_pickle(out / "test_data.pkl")
            self.assertEqual(test_df["proteins"].tolist(), ["T96060000001"])
            self.assertEqual(test_df["sequences"].tolist(), ["MCCC"])
            self.assertEqual(set(test_df.loc[0, "annotations"]), {"GO:0005886", "GO:0005575"})
            stats = json.loads((out / "reports" / "benchmark_statistics.json").read_text())
            self.assertEqual(stats["test_annotation_source"], "released_official_groundtruth")
            self.assertIn("cc-test", written)

    def test_named_evidence_policies_are_available(self):
        self.assertIn("TAS", EVIDENCE_POLICIES["cafa3-final"])
        self.assertIn("IC", EVIDENCE_POLICIES["cafa3-final"])
        self.assertIn("NAS", EVIDENCE_POLICIES["supervisor"])
        self.assertIn("ND", EVIDENCE_POLICIES["supervisor"])
        self.assertNotIn("TAS", EVIDENCE_POLICIES["cafa3-public-python"])

    def test_profiles_separate_training_and_target_policy(self):
        cafa = BENCHMARK_PROFILES["contemporary-cafa3-style"]
        supervisor = BENCHMARK_PROFILES["supervisor"]
        self.assertEqual(cafa.training_taxon_policy, "all")
        self.assertEqual(cafa.target_taxon_policy, "cafa3-targets")
        self.assertEqual(cafa.test_eligibility_policy, "ontology-no-knowledge")
        self.assertTrue(cafa.training_reviewed_only)
        self.assertFalse(cafa.target_reviewed_only)
        self.assertEqual(supervisor.training_taxon_policy, "cafa3-targets")
        self.assertEqual(supervisor.evidence_policy, "supervisor")
        self.assertEqual(supervisor.test_eligibility_policy, "global-no-knowledge")

    def test_snapshot_outputs_are_deterministic(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            outputs = []
            for name in ("first", "second"):
                out = root / name
                build_benchmark(BuildConfig(
                    uniprot_t0=(FIXTURES / "uniprot-t0.fasta",),
                    uniprot_t1=(FIXTURES / "uniprot-t1.fasta",),
                    goa_t0=FIXTURES / "goa-t0.gaf",
                    goa_t1=FIXTURES / "goa-t1.gaf",
                    go_obo=FIXTURES / "go-mini.obo",
                    output_dir=out,
                    training_taxa=frozenset({"9606"}),
                    target_taxa=frozenset({"9606"}),
                    min_count=1,
                    t0_cutoff="20250131",
                    write_checksums=False,
                    strict_qc=False,
                ))
                outputs.append(out)

            core_files = [
                "train_data.pkl", "train_data_train.pkl", "train_data_valid.pkl",
                "test_data.pkl", "terms.pkl",
                "bp-training.csv", "bp-validation.csv", "bp-test.csv",
                "cc-training.csv", "cc-validation.csv", "cc-test.csv",
                "mf-training.csv", "mf-validation.csv", "mf-test.csv",
            ]
            for filename in core_files:
                self.assertEqual(
                    (outputs[0] / filename).read_bytes(),
                    (outputs[1] / filename).read_bytes(),
                    filename,
                )


if __name__ == "__main__":
    unittest.main()
