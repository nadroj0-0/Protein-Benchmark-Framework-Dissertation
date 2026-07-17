from __future__ import annotations

import os
from pathlib import Path
import subprocess
import tempfile
import unittest


REPOSITORY = Path(__file__).resolve().parents[3]
HELPER = REPOSITORY / "scripts" / "artifact_catalog.sh"


class ArtifactCatalogTest(unittest.TestCase):
    def run_bash(
        self,
        command: str,
        *,
        environment: dict[str, str] | None = None,
        check: bool = True,
    ) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            ["bash", "-c", command],
            cwd=REPOSITORY,
            env={**os.environ, **(environment or {})},
            check=check,
            text=True,
            capture_output=True,
        )

    def test_comments_before_header_are_valid(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            artifact = root / "artifact.dat"
            artifact.write_text("fixture\n", encoding="ascii")
            catalog = root / "catalog.tsv"
            catalog.write_text(
                "# machine-local paths\n"
                "artifact_id\tpath\n"
                f"fixture\t{artifact}\n",
                encoding="ascii",
            )
            result = self.run_bash(
                f'source "{HELPER}"; '
                f'artifact_catalog_configure "{REPOSITORY}" "{catalog}"; '
                'resolve_artifact_path fixture ""'
            )
            self.assertEqual(result.stdout.strip(), str(artifact))

    def test_explicit_path_wins_over_catalog(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            explicit = root / "explicit.dat"
            catalogued = root / "catalogued.dat"
            explicit.write_text("explicit\n", encoding="ascii")
            catalogued.write_text("catalogued\n", encoding="ascii")
            catalog = root / "catalog.tsv"
            catalog.write_text(
                f"artifact_id\tpath\nfixture\t{catalogued}\n", encoding="ascii"
            )
            result = self.run_bash(
                f'source "{HELPER}"; '
                f'artifact_catalog_configure "{REPOSITORY}" "{catalog}"; '
                f'resolve_artifact_path fixture "{explicit}"'
            )
            self.assertEqual(result.stdout.strip(), str(explicit))

    def test_missing_explicit_path_falls_back_to_catalog(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            catalogued = root / "catalogued.dat"
            catalogued.write_text("catalogued\n", encoding="ascii")
            catalog = root / "catalog.tsv"
            catalog.write_text(
                f"artifact_id\tpath\nfixture\t{catalogued}\n", encoding="ascii"
            )
            missing = root / "missing.dat"
            result = self.run_bash(
                f'source "{HELPER}"; '
                f'artifact_catalog_configure "{REPOSITORY}" "{catalog}"; '
                f'resolve_artifact_path fixture "{missing}"'
            )
            self.assertEqual(result.stdout.strip(), str(catalogued))
            self.assertIn("explicit path for fixture is missing", result.stderr)

    def test_unknown_artifact_returns_nonzero_for_download_fallback(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            catalog = root / "catalog.tsv"
            catalog.write_text("artifact_id\tpath\n", encoding="ascii")
            result = self.run_bash(
                f'source "{HELPER}"; '
                f'artifact_catalog_configure "{REPOSITORY}" "{catalog}"; '
                'resolve_artifact_path unknown ""',
                check=False,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertEqual(result.stdout, "")

    def test_duplicate_and_relative_rows_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            catalog = Path(temporary) / "catalog.tsv"
            catalog.write_text(
                "artifact_id\tpath\n"
                "duplicate\trelative/path\n"
                "duplicate\t/absolute/path\n",
                encoding="ascii",
            )
            result = self.run_bash(
                f'source "{HELPER}"; artifact_catalog_validate "{catalog}"',
                check=False,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("not absolute", result.stderr)
            self.assertIn("Duplicate artifact ID", result.stderr)


if __name__ == "__main__":
    unittest.main()
