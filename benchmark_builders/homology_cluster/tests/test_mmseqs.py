from __future__ import annotations

from pathlib import Path
import contextlib
import io
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from homology_cluster_benchmark.config import SUPPORTED_IDENTITIES, parse_identity
from homology_cluster_benchmark.mmseqs import (
    ClusterIndex,
    CommandSpec,
    build_mmseqs_commands,
    execute_commands,
    resolve_mmseqs_runtime,
    validate_exact_mmseqs_version,
    validate_mmseqs_version,
    validate_recorded_exact_mmseqs_version,
)
from homology_cluster_benchmark.uniref import UniRefIndex

from tests.helpers import FIXTURES, fixture_config


class MMseqsTests(unittest.TestCase):
    def test_command_output_is_streamed_and_preserved_in_stage_log(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            command = CommandSpec(
                "fixture",
                (
                    sys.executable,
                    "-c",
                    "import sys; print('progress-one', flush=True); "
                    "print('progress-two', file=sys.stderr)",
                ),
            )
            streamed = io.StringIO()

            with contextlib.redirect_stdout(streamed):
                execute_commands((command,), root / "logs")

            expected = "progress-one\nprogress-two\n"
            self.assertEqual(streamed.getvalue(), expected)
            self.assertEqual((root / "logs" / "mmseqs_fixture.log").read_text(), expected)

    def test_failed_command_output_is_streamed_and_preserved(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            command = CommandSpec(
                "fixture-failure",
                (
                    sys.executable,
                    "-c",
                    "import sys; print('last-progress'); sys.exit(9)",
                ),
            )
            streamed = io.StringIO()

            with contextlib.redirect_stdout(streamed), self.assertRaisesRegex(
                RuntimeError, "exit code 9"
            ):
                execute_commands((command,), root / "logs")

            self.assertEqual(streamed.getvalue(), "last-progress\n")
            self.assertEqual(
                (root / "logs" / "mmseqs_fixture-failure.log").read_text(),
                "last-progress\n",
            )

    def test_all_six_commands_encode_locked_identity_and_coverage(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for identity in SUPPORTED_IDENTITIES:
                config = fixture_config(root / "out", root / "temp", identity=identity)
                commands = build_mmseqs_commands(config, FIXTURES / "uniref90.fasta", root / "work")
                cluster = commands[1].argv
                self.assertEqual(cluster[cluster.index("--min-seq-id") + 1], f"{identity:.2f}")
                self.assertEqual(cluster[cluster.index("-c") + 1], "0.8")
                self.assertEqual(cluster[cluster.index("--cov-mode") + 1], "0")
                self.assertEqual(cluster[cluster.index("--alignment-mode") + 1], "3")
                self.assertEqual(cluster[cluster.index("--seq-id-mode") + 1], "0")
                self.assertEqual(cluster[cluster.index("--cluster-reassign") + 1], "1")
                self.assertEqual(cluster[cluster.index("-s") + 1], "7.5")
                self.assertEqual(cluster[cluster.index("-e") + 1], "1e-4")
                createdb = commands[0].argv
                self.assertEqual(createdb[createdb.index("--shuffle") + 1], "0")

    def test_mmseqs_release_with_reassign_fix_is_required(self):
        self.assertEqual(validate_mmseqs_version("MMseqs2 Version: 15-6f452"), 15)
        self.assertEqual(validate_mmseqs_version("12-113e3"), 12)
        with self.assertRaisesRegex(ValueError, "too old"):
            validate_mmseqs_version("11-e1a1c")
        with self.assertRaisesRegex(ValueError, "Could not parse"):
            validate_mmseqs_version("unknown-build")

    def test_exact_runtime_version_match_mismatch_unparseable_and_nonzero(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)

            def executable(name: str, output: str, status: int = 0) -> Path:
                path = root / name
                path.write_text(
                    "#!/usr/bin/env bash\n"
                    f"printf '%s\\n' '{output}'\n"
                    f"exit {status}\n"
                )
                path.chmod(0o755)
                return path

            matching = resolve_mmseqs_runtime(str(executable("matching", "15-6f452")))
            self.assertEqual(validate_exact_mmseqs_version("15-6f452", matching), "15-6f452")
            self.assertIsNotNone(matching.executable_sha256)
            with self.assertRaisesRegex(ValueError, "exact version mismatch"):
                validate_exact_mmseqs_version("14-7e284", matching)

            unparseable = resolve_mmseqs_runtime(str(executable("unparseable", "local-build")))
            with self.assertRaisesRegex(ValueError, "unparseable"):
                validate_exact_mmseqs_version("15-6f452", unparseable)

            nonzero = resolve_mmseqs_runtime(str(executable("nonzero", "15-6f452", 9)))
            with self.assertRaisesRegex(ValueError, "exited 9"):
                validate_exact_mmseqs_version("15-6f452", nonzero)

            release_commit = "8cc5ce367b5638c4306c2d7cfc652dd099a4643f"
            commit_runtime = resolve_mmseqs_runtime(
                str(executable("release-commit", release_commit))
            )
            self.assertEqual(
                validate_exact_mmseqs_version("18-8cc5c", commit_runtime),
                release_commit,
            )
            with self.assertRaisesRegex(ValueError, "exact version mismatch"):
                validate_exact_mmseqs_version("18-dead0", commit_runtime)

        self.assertEqual(
            validate_recorded_exact_mmseqs_version("15-6f452", "15-6f452"),
            "15-6f452",
        )
        self.assertEqual(
            validate_recorded_exact_mmseqs_version(
                "18-8cc5c", "8cc5ce367b5638c4306c2d7cfc652dd099a4643f"
            ),
            "8cc5ce367b5638c4306c2d7cfc652dd099a4643f",
        )
        with self.assertRaisesRegex(ValueError, "exact version mismatch"):
            validate_recorded_exact_mmseqs_version("15-6f452", "14-7e284")
        with self.assertRaisesRegex(ValueError, "exactly one version identity"):
            validate_recorded_exact_mmseqs_version(
                "15-6f452", "MMseqs2 Version: 15-6f452"
            )

    def test_unsupported_identity_is_rejected(self):
        for value in (40, 50, 70, 90, 0.4):
            with self.assertRaisesRegex(ValueError, "Unsupported identity"):
                parse_identity(value)

    def test_cluster_tsv_is_complete_and_unique(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            uniref = UniRefIndex.build(FIXTURES / "uniref90.fasta", root / "u.sqlite")
            clusters = ClusterIndex.build(FIXTURES / "clusters.tsv", uniref, root / "c.sqlite")
            self.assertEqual(clusters.member_count(), 7)
            self.assertEqual(clusters.cluster_count(), 6)
            self.assertEqual(clusters.cluster_for("UniRef90_U1B"), "UniRef90_U1")

    def test_missing_member_and_duplicate_member_fail(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            uniref = UniRefIndex.build(FIXTURES / "uniref90.fasta", root / "u.sqlite")
            missing = root / "missing.tsv"
            lines = [line for line in (FIXTURES / "clusters.tsv").read_text().splitlines() if line]
            missing.write_text("\n".join(lines[:-1]) + "\n")
            with self.assertRaisesRegex(ValueError, "missing_members=1"):
                ClusterIndex.build(missing, uniref, root / "missing.sqlite")
            duplicate = root / "duplicate.tsv"
            duplicate.write_text((FIXTURES / "clusters.tsv").read_text() + lines[0] + "\n")
            with self.assertRaisesRegex(ValueError, "more than once"):
                ClusterIndex.build(duplicate, uniref, root / "duplicate.sqlite")

    @unittest.skipUnless(shutil.which("mmseqs"), "MMseqs2 is not installed locally")
    def test_real_mmseqs_version_smoke(self):
        result = subprocess.run(["mmseqs", "version"], capture_output=True, text=True, check=False)
        self.assertEqual(result.returncode, 0)
        self.assertTrue(result.stdout.strip())


if __name__ == "__main__":
    unittest.main()
