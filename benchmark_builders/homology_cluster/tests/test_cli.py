from __future__ import annotations

from contextlib import redirect_stdout, redirect_stderr
import io
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

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
                "--uniprot-source-scope", "sprot-only",
                "--uniprot-sprot-sequences-url", "https://example.invalid/uniprot.dat.gz",
                "--goa-url", "https://example.invalid/goa.gaf.gz",
                "--go-obo-url", "https://example.invalid/go.obo",
                "--fixture-mode",
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

    def test_authorize_array_validates_the_complete_pilot_publication(self):
        with tempfile.TemporaryDirectory() as tmp, redirect_stdout(io.StringIO()):
            root = Path(tmp)
            pilot = root / "pilot"
            pilot.mkdir()
            marker = pilot / "RUN_COMPLETE.json"
            attrition = pilot / "attrition_report.json"
            manifest = root / "frozen.json"
            for path in (marker, attrition, manifest):
                path.write_text("{}\n")
            arguments = [
                "authorize-array",
                "--attrition-policy", str(root / "policy.json"),
                "--pilot-approval", str(root / "approval.json"),
                "--pilot-completion-marker", str(marker),
                "--pilot-attrition-report", str(attrition),
                "--pilot-run-dir", str(pilot),
                "--pilot-task-context", str(root / "context.json"),
                "--pilot-measurement-evidence", str(root / "measurements.json"),
                "--frozen-input-manifest", str(manifest),
                "--framework-revision", "a" * 40,
                "--uniprot-source-scope", "sprot-only",
                "--split-policy", "sequence-balanced",
                "--training-population", "annotated-only",
                "--expected-mmseqs-version", "15-6f452",
            ]
            with (
                mock.patch(
                    "homology_cluster_benchmark.cli.validate_publication"
                ) as validate,
                mock.patch(
                    "homology_cluster_benchmark.cli.load_attrition_policy",
                    return_value=({}, "e" * 64),
                ),
                mock.patch(
                    "homology_cluster_benchmark.cli.validate_pilot_approval"
                ) as approve,
            ):
                self.assertEqual(main(arguments), 0)
            validate.assert_called_once_with(pilot.resolve())
            self.assertEqual(
                approve.call_args.kwargs["reviewed_attrition_policy_sha256"],
                "e" * 64,
            )

            wrong_marker = root / "detached-marker.json"
            wrong_marker.write_text("{}\n")
            arguments[arguments.index(str(marker))] = str(wrong_marker)
            with redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as caught:
                main(arguments)
            self.assertEqual(caught.exception.code, 2)


if __name__ == "__main__":
    unittest.main()
