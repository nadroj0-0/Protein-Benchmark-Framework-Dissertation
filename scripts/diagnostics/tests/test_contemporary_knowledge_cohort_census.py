from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


FRAMEWORK = Path(__file__).parents[3]
SCRIPT = FRAMEWORK / "scripts" / "diagnostics" / "build_contemporary_knowledge_cohort_census.py"


def write_dat(path: Path) -> None:
    records = (
        ("P1_HUMAN", "P1", "AAAA"),
        ("P2_HUMAN", "P2", "CCCC"),
        ("P3_HUMAN", "P3", "GGGG"),
    )
    with path.open("w", encoding="utf-8") as handle:
        for entry, accession, sequence in records:
            handle.write(
                f"ID   {entry} Reviewed;\n"
                f"AC   {accession};\n"
                "OX   NCBI_TaxID=9606;\n"
                f"SQ   SEQUENCE   {len(sequence)} AA;\n"
                f"     {sequence}\n"
                "//\n"
            )


def write_obo(path: Path) -> None:
    terms = (
        ("GO:0008150", "biological_process", None, None),
        ("GO:1000001", "biological_process", "GO:0008150", None),
        ("GO:1000002", "biological_process", None, "GO:1000001"),
        ("GO:0005575", "cellular_component", None, None),
        ("GO:2000001", "cellular_component", "GO:0005575", None),
        ("GO:0003674", "molecular_function", None, None),
        ("GO:3000001", "molecular_function", "GO:0003674", None),
    )
    with path.open("w", encoding="utf-8") as handle:
        handle.write("format-version: 1.2\ndata-version: fixture\n")
        for term, namespace, parent, part_of in terms:
            handle.write(f"\n[Term]\nid: {term}\nname: {term}\nnamespace: {namespace}\n")
            if parent:
                handle.write(f"is_a: {parent} ! parent\n")
            if part_of:
                handle.write(f"relationship: part_of {part_of} ! parent\n")


def write_gaf(path: Path, rows: list[tuple[str, str, str]]) -> None:
    aspect_code = {"BPO": "P", "CCO": "C", "MFO": "F"}
    with path.open("w", encoding="utf-8") as handle:
        handle.write("!gaf-version: 2.2\n")
        for protein_id, aspect, term in rows:
            values = [
                "UniProtKB",
                protein_id,
                protein_id,
                "",
                term,
                "PMID:1",
                "EXP",
                "",
                aspect_code[aspect],
                "",
                "",
                "protein",
                "taxon:9606",
                "20250101",
                "TEST",
                "",
                "",
            ]
            handle.write("\t".join(values) + "\n")


def write_csv(
    path: Path,
    terms: list[str],
    rows: list[tuple[str, str, set[str]]],
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["proteins", "sequences", *terms])
        for protein_id, sequence, positives in rows:
            writer.writerow(
                [protein_id, sequence, *(int(term in positives) for term in terms)]
            )


class ContemporaryKnowledgeCohortCensusTests(unittest.TestCase):
    def test_end_to_end_census_reconstructs_nk_lk_and_pk(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            benchmark = root / "benchmark"
            benchmark.mkdir()
            t0_sprot = root / "t0_sprot.dat"
            t1_sprot = root / "t1_sprot.dat"
            t0_trembl = root / "t0_trembl.dat"
            t1_trembl = root / "t1_trembl.dat"
            write_dat(t0_sprot)
            write_dat(t1_sprot)
            t0_trembl.write_text("# fixture has no TrEMBL records\n", encoding="utf-8")
            t1_trembl.write_text("# fixture has no TrEMBL records\n", encoding="utf-8")
            obo = root / "go.obo"
            write_obo(obo)
            t0_gaf = root / "t0.gaf"
            t1_gaf = root / "t1.gaf"
            write_gaf(
                t0_gaf,
                [
                    ("P2", "CCO", "GO:2000001"),
                    ("P3", "BPO", "GO:1000001"),
                ],
            )
            write_gaf(
                t1_gaf,
                [
                    ("P1", "BPO", "GO:1000002"),
                    ("P2", "CCO", "GO:2000001"),
                    ("P2", "BPO", "GO:1000002"),
                    ("P3", "BPO", "GO:1000001"),
                    ("P3", "BPO", "GO:1000002"),
                ],
            )
            manifest = {
                "profile": "supervisor",
                "test_eligibility_policy": "global-no-knowledge",
                "target_universe_policy": "reconstructed-all-qualifying",
                "training_reviewed_only": False,
                "target_reviewed_only": False,
                "exclude_t1_backfill": False,
                "t1_endpoint_policy": "snapshot-membership",
                "require_t0_presence": True,
                "sequence_change_policy": "exclude",
                "protein_binding_policy": "drop-mf-protein-binding-only",
                "include_relationships": True,
                "target_taxa": ["9606"],
                "training_taxa": ["9606"],
                "evidence_codes": ["EXP"],
                "allow_frozen_source_fallback": True,
            }
            (benchmark / "build_manifest.json").write_text(
                json.dumps(manifest), encoding="utf-8"
            )
            terms = {
                "bp": ["GO:0008150", "GO:1000001", "GO:1000002"],
                "cc": ["GO:0005575", "GO:2000001"],
                "mf": ["GO:0003674", "GO:3000001"],
            }
            for prefix in terms:
                write_csv(benchmark / f"{prefix}-training.csv", terms[prefix], [])
                write_csv(benchmark / f"{prefix}-validation.csv", terms[prefix], [])
                write_csv(benchmark / f"{prefix}-test.csv", terms[prefix], [])
            write_csv(
                benchmark / "cc-training.csv",
                terms["cc"],
                [("P2", "CCCC", {"GO:0005575", "GO:2000001"})],
            )
            write_csv(
                benchmark / "bp-validation.csv",
                terms["bp"],
                [("P3", "GGGG", {"GO:0008150", "GO:1000001"})],
            )
            write_csv(
                benchmark / "bp-test.csv",
                terms["bp"],
                [
                    (
                        "P1",
                        "AAAA",
                        {"GO:0008150", "GO:1000001", "GO:1000002"},
                    )
                ],
            )
            output = root / "output"
            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--accepted-benchmark-dir",
                    str(benchmark),
                    "--t0-sprot",
                    str(t0_sprot),
                    "--t0-trembl",
                    str(t0_trembl),
                    "--t1-sprot",
                    str(t1_sprot),
                    "--t1-trembl",
                    str(t1_trembl),
                    "--goa-t0",
                    str(t0_gaf),
                    "--goa-t1",
                    str(t1_gaf),
                    "--benchmark-obo",
                    str(obo),
                    "--t0-source-obo",
                    str(obo),
                    "--t1-source-obo",
                    str(obo),
                    "--output-dir",
                    str(output),
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            report = json.loads(
                (output / "temporal_annotation_ledger.json").read_text()
            )
            gains = report["gainer_cohort_counts"]["BPO"]
            self.assertEqual(gains["no_qualifying"], 1)
            self.assertEqual(gains["cross_ontology_known"], 1)
            self.assertEqual(gains["same_aspect_partial"], 1)
            self.assertTrue(
                report["preparation"]["accepted_test_alignment"]["BPO"][
                    "t1_truth_exact_match_verified"
                ]
            )
            self.assertTrue((output / "cohort_census.tsv").is_file())
            self.assertTrue((output / "RUN_COMPLETE.json").is_file())


if __name__ == "__main__":
    unittest.main()
