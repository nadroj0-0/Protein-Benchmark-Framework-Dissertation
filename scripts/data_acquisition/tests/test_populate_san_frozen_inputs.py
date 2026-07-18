from __future__ import annotations

import functools
import gzip
import hashlib
import http.server
import os
from pathlib import Path
import subprocess
import sys
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
        self,
        root: Path,
        spec: Path,
        *arguments: str,
        check: bool = True,
        environment_overrides: dict[str, str] | None = None,
    ) -> subprocess.CompletedProcess[str]:
        environment = os.environ.copy()
        environment["SAN_INPUT_SPEC"] = str(spec)
        if environment_overrides:
            environment.update(environment_overrides)
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
                path_catalog = root / "manifests" / "artifact_paths.tsv"
                self.assertTrue(path_catalog.is_file())
                self.assertEqual(
                    path_catalog.read_text(encoding="ascii"),
                    f"artifact_id\tpath\nmmseqs2\t{destination}\n",
                )
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

    def test_temporal_profile_derives_and_reuses_filtered_trembl(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            workspace = Path(temporary)
            served = workspace / "served"
            served.mkdir()

            target_record = (
                "ID   TARGET_HUMAN\n"
                "AC   P00001;\n"
                "OX   NCBI_TaxID=9606;\n"
                "SQ   SEQUENCE   4 AA;\n"
                "     AAAA\n"
                "//\n"
            ).encode("ascii")
            excluded_record = (
                "ID   EXCLUDED_TEST\n"
                "AC   P00002;\n"
                "OX   NCBI_TaxID=999999;\n"
                "SQ   SEQUENCE   4 AA;\n"
                "     CCCC\n"
                "//\n"
            ).encode("ascii")
            trembl_gzip = served / "uniprot_trembl.dat.gz"
            with trembl_gzip.open("wb") as compressed:
                with gzip.GzipFile(fileobj=compressed, mode="wb", mtime=0) as handle:
                    handle.write(target_record + excluded_record)

            knowledgebase = served / "knowledgebase2025_01.tar.gz"
            with tarfile.open(knowledgebase, "w:gz") as handle:
                handle.add(
                    trembl_gzip,
                    arcname="knowledgebase/complete/uniprot_trembl.dat.gz",
                )

            handler = functools.partial(QuietHandler, directory=str(served))
            QuietHandler.requests = 0
            server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)
            thread = threading.Thread(target=server.serve_forever, daemon=True)
            thread.start()
            try:
                spec = workspace / "spec.tsv"
                rows = [
                    "# profiles\trole\trelease\trelative_path\turl\texpected_bytes\t"
                    "checksum_algorithm\texpected_checksum\tvalidator",
                    "\t".join(
                        [
                            "temporal",
                            "uniprot_knowledgebase_t0",
                            "2025_01",
                            "frozen_inputs/uniprot/2025_01/knowledgebase2025_01.tar.gz",
                            f"http://127.0.0.1:{server.server_port}/{knowledgebase.name}",
                            str(knowledgebase.stat().st_size),
                            "sha256",
                            hashlib.sha256(knowledgebase.read_bytes()).hexdigest(),
                            "tar-gzip",
                        ]
                    ),
                    "\t".join(
                        [
                            "temporal",
                            "uniprot_trembl_t1",
                            "2026_02",
                            "frozen_inputs/uniprot/2026_02/uniprot_trembl.dat.gz",
                            f"http://127.0.0.1:{server.server_port}/{trembl_gzip.name}",
                            str(trembl_gzip.stat().st_size),
                            "sha256",
                            hashlib.sha256(trembl_gzip.read_bytes()).hexdigest(),
                            "uniprot-dat-gzip",
                        ]
                    ),
                ]
                spec.write_text("\n".join(rows) + "\n", encoding="ascii")
                root = workspace / "store"
                python_wrapper = workspace / "python-with-bind-check.sh"
                python_wrapper.write_text(
                    "#!/bin/sh\n"
                    "case \"${SINGULARITY_BINDPATH:-}\" in\n"
                    "  *\"${EXPECTED_BIND}\"*) ;;\n"
                    "  *) echo \"missing expected bind: ${EXPECTED_BIND}\" >&2; exit 97 ;;\n"
                    "esac\n"
                    f'exec "{sys.executable}" "$@"\n',
                    encoding="ascii",
                )
                python_wrapper.chmod(0o755)
                run_environment = {
                    "PYTHON_BIN": str(python_wrapper),
                    "EXPECTED_BIND": f"{root.resolve()}:{root.resolve()}",
                }

                first = self.run_script(
                    root,
                    spec,
                    "--profile",
                    "temporal",
                    environment_overrides=run_environment,
                )
                self.assertIn("derived:    2 created, 0 skipped", first.stdout)
                derived_paths = [
                    root
                    / "derived_inputs/uniprot/cafa3_target_taxa/2025_01/"
                    "uniprot_trembl_cafa3_targets.dat.gz",
                    root
                    / "derived_inputs/uniprot/cafa3_target_taxa/2026_02/"
                    "uniprot_trembl_cafa3_targets.dat.gz",
                ]
                for path in derived_paths:
                    with gzip.open(path, "rb") as handle:
                        self.assertEqual(handle.read(), target_record)
                    self.assertTrue(Path(f"{path}.sha256").is_file())
                    self.assertTrue(Path(f"{path}.provenance.tsv").is_file())
                    self.assertTrue(Path(f"{path}.derivation.tsv").is_file())

                path_catalog = (
                    root / "manifests" / "artifact_paths.tsv"
                ).read_text(encoding="ascii")
                self.assertIn("uniprot_trembl_cafa3_targets_t0\t", path_catalog)
                self.assertIn("uniprot_trembl_cafa3_targets_t1\t", path_catalog)
                self.assertEqual(QuietHandler.requests, 2)

                server.shutdown()
                thread.join(timeout=5)
                second = self.run_script(
                    root,
                    spec,
                    "--profile",
                    "temporal",
                    environment_overrides=run_environment,
                )
                self.assertIn("downloaded: 0", second.stdout)
                self.assertIn("derived:    0 created, 2 skipped", second.stdout)

                derived_paths[1].unlink()
                rebuilt = self.run_script(
                    root,
                    spec,
                    "--profile",
                    "temporal",
                    environment_overrides=run_environment,
                )
                self.assertIn("derived:    1 created, 1 skipped", rebuilt.stdout)
                with gzip.open(derived_paths[1], "rb") as handle:
                    self.assertEqual(handle.read(), target_record)
            finally:
                server.shutdown()
                server.server_close()
                thread.join(timeout=5)

    def test_csv_validator_accepts_both_canonical_protein_headers(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            workspace = Path(temporary)
            served = workspace / "served"
            served.mkdir()
            csv_files = {
                "singular.csv": b"protein,sequences,GO:0000001\nP1,AAAA,1\n",
                "plural.csv": b"proteins,sequences,GO:0000001\nP2,CCCC,1\n",
            }
            for name, content in csv_files.items():
                (served / name).write_bytes(content)

            handler = functools.partial(QuietHandler, directory=str(served))
            server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)
            thread = threading.Thread(target=server.serve_forever, daemon=True)
            thread.start()
            try:
                spec = workspace / "spec.tsv"
                lines = [
                    "# profiles\trole\trelease\trelative_path\turl\texpected_bytes\t"
                    "checksum_algorithm\texpected_checksum\tvalidator"
                ]
                for name, content in csv_files.items():
                    lines.append(
                        f"references\t{name}\ttest\treferences/{name}\t"
                        f"http://127.0.0.1:{server.server_port}/{name}\t"
                        f"{len(content)}\tmd5\t{hashlib.md5(content).hexdigest()}\tcsv"
                    )
                spec.write_text("\n".join(lines) + "\n", encoding="ascii")

                root = workspace / "store"
                result = self.run_script(root, spec, "--profile", "references")
                self.assertIn("downloaded: 2", result.stdout)
                for name, content in csv_files.items():
                    self.assertEqual((root / "references" / name).read_bytes(), content)
            finally:
                server.shutdown()
                server.server_close()
                thread.join(timeout=5)

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

    def test_homology_profile_defines_idempotent_common_preprocessing_cache(self) -> None:
        script = SCRIPT.read_text(encoding="utf-8")
        self.assertIn(
            'HOMOLOGY_CACHE_ROLE="homology_common_preprocessing_2026_02"', script
        )
        self.assertIn("process_homology_derived_inputs()", script)
        self.assertIn("homology_cluster_benchmark.common_cache", script)
        self.assertIn('verify "${verify_args[@]}"', script)
        self.assertIn('common_cache build "${build_args[@]}"', script)
        self.assertIn("--replace-existing", script)
        self.assertIn("process_homology_derived_inputs", script)
        self.assertIn("$HOMOLOGY_CACHE_RELATIVE/$CACHE_MARKER", script)


if __name__ == "__main__":
    unittest.main()
