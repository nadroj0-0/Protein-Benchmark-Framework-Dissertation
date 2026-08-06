from __future__ import annotations

import csv
import gzip
import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest


SCRIPT = Path(__file__).resolve().parents[1] / "audit_homology_evidence_policy.py"
SPEC = importlib.util.spec_from_file_location("homology_evidence_policy", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def _write_benchmark(root: Path) -> None:
    root.mkdir()
    for aspect in MODULE.ASPECTS:
        for split_name in MODULE.SPLITS:
            (root / f"{aspect}-{split_name}.csv").write_text(
                "proteins,sequences,GO:0000001\nP1,MA,1\n", encoding="ascii"
            )
    with gzip.open(root / "qualifying_annotations.tsv.gz", "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["uniprot_accession", "aspect", "evidence_code", "split"],
            delimiter="\t",
        )
        writer.writeheader()
        for code in sorted(MODULE.EXPECTED_CODES):
            writer.writerow({
                "uniprot_accession": "P1", "aspect": "P",
                "evidence_code": code, "split": "training",
            })


class HomologyEvidencePolicyTests(unittest.TestCase):
    def test_separates_experimental_and_non_experimental_sources(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            _write_benchmark(root / "benchmark")
            args = MODULE.argparse.Namespace(
                benchmark=[f"fixture={root / 'benchmark'}"], output_dir=root / "out"
            )
            payload = MODULE.run(args)
            audit = payload["benchmarks"][0]
            categories = {
                row["category"]: row for row in audit["category_rows"]
                if row["aspect"] == "bp" and row["split"] == "training"
            }
            self.assertEqual(categories["experimental"]["proteins"], 1)
            self.assertEqual(categories["no_biological_data"]["proteins"], 1)
            self.assertTrue((root / "out" / "RUN_COMPLETE.json").is_file())


if __name__ == "__main__":
    unittest.main()
