from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT.parents[1] / "embedding_inventory" / "src"))

from cafa_benchmark_builder.config import SUPERVISOR_EXP_CODES  # noqa: E402
from cafa_benchmark_builder.ontology import Ontology  # noqa: E402
from cafa_benchmark_builder.random_baseline import (  # noqa: E402
    _relabel_at_snapshot,
    build_random_baseline,
)
from pfp_embedding_inventory.benchmark import parse_benchmark  # noqa: E402
from pfp_embedding_inventory.models import BenchmarkContract  # noqa: E402


FIXTURES = Path(__file__).resolve().parent / "fixtures"
ANNOTATIONS = {
    "GO:0009987", "GO:0008150",
    "GO:0005488", "GO:0003674",
    "GO:0005886", "GO:0005575",
}


class RandomBaselineTest(unittest.TestCase):
    def _source(self, root: Path) -> Path:
        source = root / "source"
        source.mkdir()
        proteins = [f"P{index:05d}" for index in range(1, 10)]
        sequences = {
            protein: "M" + chr(64 + index) * 5
            for index, protein in enumerate(proteins, start=1)
        }
        frames = {
            "train_data_train.pkl": proteins[:5],
            "train_data_valid.pkl": proteins[5:7],
            "test_data.pkl": proteins[7:],
        }
        terms_by_aspect = {
            "bp": ["GO:0009987", "GO:0008150"],
            "cc": ["GO:0005886", "GO:0005575"],
            "mf": ["GO:0005488", "GO:0003674"],
        }
        for filename, ids in frames.items():
            pd.DataFrame({
                "proteins": ids,
                "sequences": [sequences[protein] for protein in ids],
                "annotations": [ANNOTATIONS for _ in ids],
            }).to_pickle(source / filename)
        pd.DataFrame({"terms": sorted(ANNOTATIONS)}).to_pickle(source / "terms.pkl")
        for filename, ids in frames.items():
            split = {
                "train_data_train.pkl": "training",
                "train_data_valid.pkl": "validation",
                "test_data.pkl": "test",
            }[filename]
            for prefix, terms in terms_by_aspect.items():
                data = {
                    "proteins": ids,
                    "sequences": [sequences[protein] for protein in ids],
                    **{term: [1] * len(ids) for term in terms},
                }
                pd.DataFrame(data).to_csv(source / f"{prefix}-{split}.csv", index=False)
        return source

    def test_prepared_mode_preserves_population_and_writes_evidence(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = self._source(root)
            output = root / "output"
            build_random_baseline(
                mode="prepared",
                benchmark_id="fixture-random-h",
                source_dir=source,
                go_obo=FIXTURES / "go-mini.obo",
                output_dir=output,
            )

            manifest = json.loads((output / "build_manifest.json").read_text())
            validation = json.loads((output / "validation_report.json").read_text())
            self.assertEqual(manifest["source_split_counts"], {"training": 5, "validation": 2, "test": 2})
            self.assertEqual(manifest["source_split_counts"], manifest["output_split_counts"])
            self.assertTrue(validation["valid"])
            self.assertTrue((output / "RUN_COMPLETE.json").is_file())
            self.assertEqual(len(list(output.glob("*-*.csv"))), 9)
            self.assertEqual(len(pd.read_pickle(output / "train_data.pkl")), 7)
            self.assertEqual(len(pd.read_pickle(output / "test_data.pkl")), 2)
            parsed = parse_benchmark(
                output,
                BenchmarkContract(
                    id_overlap="global-disjoint",
                    sequence_overlap="global-disjoint",
                    protein_id_pattern=r"^[^\s/\\]+$",
                    sequence_pattern=r"^[A-Za-z*.-]+$",
                ),
            )
            self.assertEqual(len(parsed.proteins), 9)

    def test_t1_relabelling_uses_current_snapshot_labels(self):
        source = pd.DataFrame({
            "proteins": ["P00001", "P00004", "P00005"],
            "sequences": ["OLD1", "OLD4", "OLD5"],
            "annotations": [{"GO:0008150"}] * 3,
        })
        benchmark_go = Ontology(FIXTURES / "go-mini.obo", with_rels=True)
        relabelled, stats = _relabel_at_snapshot(
            source,
            (FIXTURES / "uniprot-t1.fasta",),
            FIXTURES / "goa-t1.gaf",
            benchmark_go,
            benchmark_go,
            SUPERVISOR_EXP_CODES,
        )

        by_id = relabelled.set_index("proteins")
        self.assertEqual(by_id.loc["P00004", "sequences"], "MEEEEE")
        self.assertIn("GO:0005488", set(by_id.loc["P00004", "annotations"]))
        self.assertIn("GO:0005886", set(by_id.loc["P00004", "annotations"]))
        self.assertEqual(stats["qualifying_t1_proteins"], 3)
        self.assertNotIn("TAS", stats["goa_evidence_counts"])
        self.assertNotIn("IC", stats["goa_evidence_counts"])

    def test_t1_relabelling_resolves_split_alias_by_exact_sequence(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            uniprot = root / "uniprot.dat"
            uniprot.write_text(
                "ID   FIRST_HUMAN Reviewed; 6 AA.\n"
                "AC   PNEW01; POLD01;\n"
                "OX   NCBI_TaxID=9606;\n"
                "SQ   SEQUENCE 6 AA;\n"
                "     MAAAAA\n"
                "//\n"
                "ID   SECOND_HUMAN Reviewed; 7 AA.\n"
                "AC   PNEW02; POLD01;\n"
                "OX   NCBI_TaxID=9606;\n"
                "SQ   SEQUENCE 7 AA;\n"
                "     MMAAAAA\n"
                "//\n",
                encoding="utf-8",
            )
            goa = root / "goa.gaf"
            goa.write_text(
                "!gaf-version: 2.2\n"
                "UniProtKB\tPNEW01\tP1\t\tGO:0009987\tPMID:1\tIDA\t\tP\t"
                "Protein one\t\tprotein\ttaxon:9606\t20260101\tUniProt\t\t\n",
                encoding="utf-8",
            )
            source = pd.DataFrame({
                "proteins": ["POLD01"],
                "sequences": ["MAAAAA"],
                "annotations": [{"GO:0008150"}],
            })
            benchmark_go = Ontology(FIXTURES / "go-mini.obo", with_rels=True)

            relabelled, stats = _relabel_at_snapshot(
                source,
                (uniprot,),
                goa,
                benchmark_go,
                benchmark_go,
                SUPERVISOR_EXP_CODES,
            )

            self.assertEqual(relabelled.loc[0, "proteins"], "POLD01")
            self.assertEqual(relabelled.loc[0, "sequences"], "MAAAAA")
            self.assertEqual(stats["sequence_disambiguated_source_proteins"], 1)


if __name__ == "__main__":
    unittest.main()
