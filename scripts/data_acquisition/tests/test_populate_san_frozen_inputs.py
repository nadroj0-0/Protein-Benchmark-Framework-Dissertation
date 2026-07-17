from __future__ import annotations

import functools
import hashlib
import http.server
import os
from pathlib import Path
import subprocess
import tarfile
import tempfile
import threading
import unittest


REPOSITORY = Path(__file__).resolve().parents[3]
SCRIPT = REPOSITORY / "scripts" / "data_acquisition" / "populate_san_frozen_inputs.sh"
PRODUCTION_SPEC = REPOSITORY / "scripts" / "data_acquisition" / "san_frozen_inputs.tsv"


class QuietHandler(http.server.SimpleHTTPRequestHandler):
    requests = 0

    def do_GET(self) -> None:  # noqa: N802 - inherited HTTP method name
        type(self).requests += 1
        super().do_GET()

    def log_message(self, _format: str, *_args: object) -> None:
        return


class SanAcquisitionTest(unittest.TestCase):
    def run_script(
        self, root: Path, spec: Path, *arguments: str, check: bool = True
    ) -> subprocess.CompletedProcess[str]:
        environment = os.environ.copy()
        environment["SAN_INPUT_SPEC"] = str(spec)
        return subprocess.run(
            ["bash", str(SCRIPT), "--root", str(root), "--reserve-gb", "0", *arguments],
            cwd=REPOSITORY,
            env=environment,
            check=check,
            text=True,
            capture_output=True,
        )

    def test_download_publish_rerun_and_verify(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            workspace = Path(temporary)
            served = workspace / "served"
            served.mkdir()
            archive = served / "mmseqs-linux-avx2.tar.gz"
            executable = workspace / "mmseqs"
            executable.write_text("synthetic mmseqs\n", encoding="ascii")
            with tarfile.open(archive, "w:gz") as handle:
                handle.add(executable, arcname="mmseqs/bin/mmseqs")

            handler = functools.partial(QuietHandler, directory=str(served))
            QuietHandler.requests = 0
            server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)
            thread = threading.Thread(target=server.serve_forever, daemon=True)
            thread.start()
            try:
                digest = hashlib.sha256(archive.read_bytes()).hexdigest()
                spec = workspace / "spec.tsv"
                spec.write_text(
                    "# profiles\trole\trelease\trelative_path\turl\texpected_bytes\t"
                    "checksum_algorithm\texpected_checksum\tvalidator\n"
                    f"tools\tmmseqs2\ttest\ttools/mmseqs2/test/archive.tar.gz\t"
                    f"http://127.0.0.1:{server.server_port}/{archive.name}\t"
                    f"{archive.stat().st_size}\tsha256\t{digest}\tmmseqs-archive\n",
                    encoding="ascii",
                )
                root = workspace / "store"

                first = self.run_script(root, spec, "--profile", "tools")
                self.assertIn("downloaded: 1", first.stdout)
                destination = root / "tools" / "mmseqs2" / "test" / "archive.tar.gz"
                self.assertEqual(destination.read_bytes(), archive.read_bytes())
                self.assertTrue(Path(f"{destination}.sha256").is_file())
                self.assertTrue(Path(f"{destination}.provenance.tsv").is_file())
                self.assertTrue((root / "manifests" / "frozen_input_catalog.tsv").is_file())
                self.assertEqual(QuietHandler.requests, 1)

                server.shutdown()
                thread.join(timeout=5)
                second = self.run_script(root, spec, "--profile", "tools")
                self.assertIn("downloaded: 0", second.stdout)
                self.assertIn("skipped:    1", second.stdout)

                verified = self.run_script(
                    root, spec, "--profile", "tools", "--verify-only"
                )
                self.assertIn("full checks: 1", verified.stdout)

                destination.write_bytes(destination.read_bytes() + b"corrupt")
                failed = self.run_script(
                    root, spec, "--profile", "tools", "--verify-only", check=False
                )
                self.assertNotEqual(failed.returncode, 0)
                self.assertIn("Size mismatch", failed.stderr)
            finally:
                server.shutdown()
                server.server_close()
                thread.join(timeout=5)

    def test_dry_run_does_not_create_root(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            workspace = Path(temporary)
            spec = workspace / "spec.tsv"
            spec.write_text(
                "# profiles\trole\trelease\trelative_path\turl\texpected_bytes\t"
                "checksum_algorithm\texpected_checksum\tvalidator\n"
                "references\tfile\ttest\treferences/file.txt\t"
                "https://example.invalid/file.txt\t10\t-\t-\ttext\n",
                encoding="ascii",
            )
            root = workspace / "not-created"
            result = self.run_script(root, spec, "--profile", "references", "--dry-run")
            self.assertIn("No files were changed", result.stdout)
            self.assertFalse(root.exists())

    def test_goa_234_uses_immutable_archive_and_pinned_sha256(self) -> None:
        specification = PRODUCTION_SPEC.read_text(encoding="utf-8")
        matching = [
            line
            for line in specification.splitlines()
            if "\tgoa_t1\t234\t" in line
        ]
        self.assertEqual(len(matching), 1)
        fields = matching[0].split("\t")
        self.assertEqual(
            fields[4],
            "https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/UNIPROT/"
            "goa_uniprot_all.gaf.234.gz",
        )
        self.assertEqual(fields[5], "11664243116")
        self.assertEqual(fields[6], "sha256")
        self.assertEqual(
            fields[7],
            "f315375b07946a0649142b2f4de2e15e282316989677a04e7a561203186dd2ff",
        )
        self.assertNotIn("goa/current_release_numbers.txt", specification)
        self.assertNotIn("goa_t1_md5", specification)


if __name__ == "__main__":
    unittest.main()
