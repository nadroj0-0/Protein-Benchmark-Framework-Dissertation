from __future__ import annotations

# ruff: noqa: E402 - tests add the local src tree before importing the package.

import json
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from cafa_benchmark_builder.organizer_replay import (
    OrganizerReplayConfig,
    run_organizer_replay,
)


def dat_record(entry: str, accession: str, sequence: str) -> str:
    return (
        f"ID   {entry} Reviewed;\n"
        f"AC   {accession};\n"
        "OX   NCBI_TaxID=9606;\n"
        f"SQ   SEQUENCE   {len(sequence)} AA;\n"
        f"     {sequence}\n"
        "//\n"
    )


def gaf_row(
    accession: str,
    go_id: str,
    aspect: str,
    date: str,
    qualifier: str = "",
) -> str:
    return (
        f"UniProtKB\t{accession}\tP\t{qualifier}\t{go_id}\tPMID:1\tIDA\t\t"
        f"{aspect}\tProtein\t\tprotein\ttaxon:9606\t{date}\tUniProt\t\t\n"
    )


OBO = """format-version: 1.2

[Term]
id: GO:0008150
name: biological_process
namespace: biological_process

[Term]
id: GO:0009987
name: cellular process
namespace: biological_process
is_a: GO:0008150 ! biological_process

[Term]
id: GO:1234567
name: child process
namespace: biological_process
is_a: GO:0009987 ! cellular process

[Term]
id: GO:0003674
name: molecular_function
namespace: molecular_function

[Term]
id: GO:0005488
name: binding
namespace: molecular_function
is_a: GO:0003674 ! molecular_function

[Term]
id: GO:0005515
name: protein binding
namespace: molecular_function
is_a: GO:0005488 ! binding

[Term]
id: GO:0003824
name: catalytic activity
namespace: molecular_function
is_a: GO:0003674 ! molecular_function

[Term]
id: GO:0005575
name: cellular_component
namespace: cellular_component

[Term]
id: GO:0005886
name: plasma membrane
namespace: cellular_component
is_a: GO:0005575 ! cellular_component
"""


class OrganizerReplayTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.obo = self.root / "go.obo"
        self.obo.write_text(OBO)

        records = []
        fasta = []
        mappings = []
        for index in range(1, 9):
            accession = f"P{index}"
            entry = f"ENTRY{index}_HUMAN"
            sequence = "M" + chr(64 + index) * 5
            target = f"T96060000000{index}"
            records.append(dat_record(entry, accession, sequence))
            fasta.append(f">{target} {entry}\n{sequence}\n")
            mappings.append(f"{target}\t{entry}\n")

        self.uniprot_t0 = self.root / "t0.dat"
        self.uniprot_t1 = self.root / "t1.dat"
        self.uniprot_t0.write_text("".join(records))
        self.uniprot_t1.write_text("".join(records))
        self.targets = self.root / "targets.fasta"
        self.targets.write_text("".join(fasta))
        self.mapping_dir = self.root / "mappings"
        self.mapping_dir.mkdir()
        (self.mapping_dir / "sp_species.9606.map").write_text("".join(mappings))

        self.goa_t0 = self.root / "t0.gaf"
        self.goa_t0.write_text(
            "!gaf-version: 2.2\n"
            + gaf_row("P2", "GO:0003824", "F", "20170101")
            + gaf_row("P3", "GO:0008150", "P", "20170101")
            + gaf_row("P4", "GO:0005515", "F", "20170101")
            + gaf_row("P8", "GO:0005886", "C", "20170101")
        )
        self.goa_t1 = self.root / "t1.gaf"
        self.goa_t1.write_text(
            "!gaf-version: 2.2\n"
            + gaf_row("P1", "GO:1234567", "P", "20170301")
            + gaf_row("P2", "GO:1234567", "P", "20170301")
            + gaf_row("P3", "GO:1234567", "P", "20170301")
            + gaf_row("P4", "GO:1234567", "P", "20170301")
            + gaf_row("P4", "GO:0003824", "F", "20170301")
            + gaf_row("P5", "GO:0005515", "F", "20170301")
            + gaf_row("P6", "GO:1234567", "P", "20170101")
            + gaf_row("P7", "GO:0005886", "C", "20170301", "NOT")
            + gaf_row("P8", "GO:0003824", "F", "20170301")
        )

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def _official_reference(self) -> Path:
        root = self.root / "official" / "benchmark20171115"
        (root / "lists").mkdir(parents=True)
        (root / "groundtruth").mkdir()
        lists = {
            "bpo_all_type1.txt": ["T960600000001", "T960600000004"],
            "bpo_all_type2.txt": ["T960600000002"],
            "bpo_all_typex.txt": [
                "T960600000001", "T960600000002", "T960600000004",
            ],
            "cco_all_type1.txt": [],
            "cco_all_type2.txt": [],
            "cco_all_typex.txt": [],
            "mfo_all_type1.txt": ["T960600000004"],
            "mfo_all_type2.txt": ["T960600000008"],
            "mfo_all_typex.txt": ["T960600000004", "T960600000008"],
        }
        for name, values in lists.items():
            (root / "lists" / name).write_text("".join(f"{value}\n" for value in values))
        (root / "lists" / "xxo_all_typex.txt").write_text(
            "T960600000001\n"
            "T960600000002\n"
            "T960600000004\n"
            "T960600000004\n"
            "T960600000008\n"
        )
        (root / "groundtruth" / "leafonly_BPO.txt").write_text(
            "T960600000001\tGO:1234567\n"
            "T960600000002\tGO:1234567\n"
            "T960600000004\tGO:1234567\n"
        )
        (root / "groundtruth" / "leafonly_CCO.txt").write_text("")
        (root / "groundtruth" / "leafonly_MFO.txt").write_text(
            "T960600000004\tGO:0003824\n"
            "T960600000008\tGO:0003824\n"
        )
        return root.parent

    def test_cafa2_aligned_membership_leaf_truth_and_comparison(self) -> None:
        output = self.root / "replay"
        written = run_organizer_replay(OrganizerReplayConfig(
            uniprot_t0=(self.uniprot_t0,),
            uniprot_t1=(self.uniprot_t1,),
            goa_t0=self.goa_t0,
            goa_t1=self.goa_t1,
            go_obo=self.obo,
            official_target_fastas=(self.targets,),
            official_target_mapping_dir=self.mapping_dir,
            official_benchmark_dir=self._official_reference(),
            output_dir=output,
            write_input_checksums=False,
        ))

        self.assertEqual(
            (output / "benchmark20171115/lists/bpo_all_type1.txt").read_text().splitlines(),
            ["T960600000001", "T960600000004"],
        )
        self.assertEqual(
            (output / "benchmark20171115/lists/bpo_all_type2.txt").read_text().splitlines(),
            ["T960600000002"],
        )
        self.assertEqual(
            (output / "benchmark20171115/lists/mfo_all_type1.txt").read_text().splitlines(),
            ["T960600000004"],
        )
        self.assertNotIn(
            "T960600000005",
            (output / "benchmark20171115/lists/mfo_all_typex.txt").read_text(),
        )
        self.assertEqual(
            (output / "benchmark20171115/groundtruth/leafonly_BPO.txt")
            .read_text().splitlines(),
            [
                "T960600000001\tGO:1234567",
                "T960600000002\tGO:1234567",
                "T960600000004\tGO:1234567",
            ],
        )
        self.assertEqual(
            (output / "benchmark20171115/lists/xxo_all_typex.txt").read_text().count(
                "T960600000004\n"
            ),
            2,
        )
        self.assertIn(
            "T960600000004\tmfo\tGO:0005515",
            (output / "reports/stage_01_t0_filtered_direct_annotations.tsv").read_text(),
        )
        self.assertNotIn(
            "T960600000004\tmfo\tGO:0005515",
            (output / "reports/stage_03_t0_after_protein_binding_filter.tsv").read_text(),
        )

        comparisons = json.loads(
            (output / "reports/official_comparison.json").read_text()
        )
        self.assertTrue(all(value["jaccard"] == 1.0 for value in comparisons.values()))
        manifest = json.loads((output / "reports/replay_manifest.json").read_text())
        self.assertEqual(manifest["membership_counts"]["bpo"]["typex"], 3)
        self.assertEqual(
            manifest["policy"]["protein_binding_policy"],
            "CAFA2 pre-classification removal at t0 and t1",
        )
        self.assertTrue(written["completion"].is_file())

    def test_validation_happens_before_output_creation(self) -> None:
        output = self.root / "bad-output"
        with self.assertRaises(FileNotFoundError):
            run_organizer_replay(OrganizerReplayConfig(
                uniprot_t0=(self.root / "missing.dat",),
                uniprot_t1=(self.uniprot_t1,),
                goa_t0=self.goa_t0,
                goa_t1=self.goa_t1,
                go_obo=self.obo,
                official_target_fastas=(self.targets,),
                official_target_mapping_dir=self.mapping_dir,
                output_dir=output,
            ))
        self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
