from __future__ import annotations

from contextlib import redirect_stdout, redirect_stderr
import io
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from homology_cluster_benchmark.cli import _identities, main
from homology_cluster_benchmark.config import SUPPORTED_IDENTITIES

class CLITests(unittest.TestCase):
    def test_all_expands_to_exact_six_locked_thresholds(self):
        self.assertEqual(_identities(["all"]), SUPPORTED_IDENTITIES)

    def test_dry_run_prints_exact_commands_without_resolving_inputs(self):
        with tempfile.TemporaryDirectory() as tmp, redirect_stdout(io.StringIO()) as output:
            status = main([
                "build", "--identity", "5", "--output-dir", str(Path(tmp) / "out"),
                "--uniref90-fasta-url", "https://example.invalid/uniref90.fasta.gz",
                "--idmapping-url", "https://example.invalid/idmapping.gz",
                "--uniprot-sequences-url", "https://example.invalid/uniprot.fasta.gz",
                "--goa-url", "https://example.invalid/goa.gaf.gz",
                "--go-obo-url", "https://example.invalid/go.obo",
                "--dry-run",
            ])
            self.assertEqual(status, 0)
            text = output.getvalue()
            self.assertIn("--min-seq-id 0.05", text)
            self.assertIn("-c 0.8", text)
            self.assertIn("--cov-mode 0", text)
            self.assertIn("seed_0", text)
            self.assertIn("min_count_50", text)

    def test_cli_rejects_unlocked_identity(self):
        with tempfile.TemporaryDirectory() as tmp, redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit) as caught:
                main(["build", "--identity", "40", "--output-dir", str(Path(tmp) / "out"), "--dry-run"])
            self.assertEqual(caught.exception.code, 2)

    def test_cli_rejects_all_members_before_input_resolution(self):
        with tempfile.TemporaryDirectory() as tmp, redirect_stderr(io.StringIO()) as errors:
            with self.assertRaises(SystemExit) as caught:
                main([
                    "build", "--identity", "30", "--output-dir", str(Path(tmp) / "out"),
                    "--training-population", "all-cluster-members", "--dry-run",
                ])
            self.assertEqual(caught.exception.code, 2)
            self.assertIn("intentionally unsupported", errors.getvalue())


if __name__ == "__main__":
    unittest.main()
