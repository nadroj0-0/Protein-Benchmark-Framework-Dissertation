from __future__ import annotations

import gzip
import io
import subprocess
import sys
import tarfile
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
FILTER = ROOT / "scripts" / "benchmark_generation" / "filter_uniprot_dat.py"
EXTRACT = ROOT / "scripts" / "benchmark_generation" / "extract_tar_member.py"


def dat_record(accession: str, taxon: str) -> str:
    return (
        f"ID   {accession}_ENTRY Unreviewed;\n"
        f"AC   {accession};\n"
        f"OX   NCBI_TaxID={taxon};\n"
        "SQ   SEQUENCE   4 AA;\n"
        "     MAAA\n"
        "//\n"
    )


class AcquisitionHelperTest(unittest.TestCase):
    def test_uniprot_filter_keeps_only_requested_taxa(self):
        with tempfile.TemporaryDirectory() as tmp:
            taxa = Path(tmp) / "taxa.txt"
            taxa.write_text("9606\n")
            source = dat_record("PKEEP", "9606") + dat_record("PDROP", "10090")
            result = subprocess.run(
                [sys.executable, str(FILTER), "--taxa-file", str(taxa)],
                input=source,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=True,
            )
            self.assertIn("PKEEP", result.stdout)
            self.assertNotIn("PDROP", result.stdout)
            self.assertIn("processed=2 kept=1", result.stderr)

    def test_tar_member_extractor_streams_nested_gzip_member(self):
        with tempfile.TemporaryDirectory() as tmp:
            archive = Path(tmp) / "knowledgebase.tar.gz"
            payload = gzip.compress(dat_record("PKEEP", "9606").encode())
            with tarfile.open(archive, "w:gz") as tar:
                info = tarfile.TarInfo("release/uniprot_trembl.dat.gz")
                info.size = len(payload)
                tar.addfile(info, io.BytesIO(payload))
                other = b"ignored"
                info = tarfile.TarInfo("docs/readme.txt")
                info.size = len(other)
                tar.addfile(info, io.BytesIO(other))

            result = subprocess.run(
                [
                    sys.executable,
                    str(EXTRACT),
                    "--archive",
                    str(archive),
                    "--suffix",
                    "uniprot_trembl.dat.gz",
                ],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=True,
            )
            self.assertEqual(gzip.decompress(result.stdout).decode(), dat_record("PKEEP", "9606"))
            self.assertIn(b"release/uniprot_trembl.dat.gz", result.stderr)


if __name__ == "__main__":
    unittest.main()
