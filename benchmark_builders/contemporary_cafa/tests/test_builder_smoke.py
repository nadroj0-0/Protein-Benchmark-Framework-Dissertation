from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from cafa_benchmark_builder.builder import build_benchmark
from cafa_benchmark_builder.config import BuildConfig
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
                target_taxa=frozenset({"9606"}),
                min_count=1,
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
            self.assertEqual(set(bp_test["proteins"]), {"P00005"})
            self.assertIn("GO:0009987", bp_test.columns)
            self.assertIn("GO:0008150", bp_test.columns)

            all_train_valid = []
            for name in ["bp-training.csv", "bp-validation.csv", "cc-training.csv",
                         "cc-validation.csv", "mf-training.csv", "mf-validation.csv"]:
                all_train_valid.extend(pd.read_csv(out / name)["proteins"].tolist())
            self.assertIn("P00002", all_train_valid)  # TAS survived
            self.assertIn("P00003", all_train_valid)  # IC survived
            self.assertNotIn("P00004", all_train_valid)


if __name__ == "__main__":
    unittest.main()
