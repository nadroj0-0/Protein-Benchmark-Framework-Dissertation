from __future__ import annotations

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "audit_pfp_working_copy.py"


def git(repo: Path, *args: str) -> None:
    subprocess.run(["git", "-C", str(repo), *args], check=True, capture_output=True)


class AuditPfpWorkingCopyTests(unittest.TestCase):
    def test_separates_public_changes_from_local_only_files(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            public = root / "public-source"
            local = root / "local-working"
            public.mkdir()
            git(public, "init", "-b", "main")
            git(public, "config", "user.name", "Test")
            git(public, "config", "user.email", "test@example.com")
            (public / "README.md").write_text("public\n", encoding="utf-8")
            (public / "script.py").write_text("print('public')\n", encoding="utf-8")
            git(public, "add", ".")
            git(public, "commit", "-m", "public")

            subprocess.run(
                ["git", "clone", "--quiet", str(public), str(local)], check=True
            )
            (local / "script.py").write_text("print('local')\n", encoding="utf-8")
            (local / "README.md").unlink()
            (local / "notes.tmp").write_text("private notes\n", encoding="utf-8")

            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    str(local),
                    "--public-repository",
                    str(public),
                    "--skip-environment",
                ],
                check=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )

            self.assertEqual(result.returncode, 1, result.stderr)
            self.assertIn("The local working copy differs", result.stdout)
            self.assertIn("`script.py`", result.stdout)
            self.assertIn("`README.md`", result.stdout)
            self.assertIn("`notes.tmp`", result.stdout)
            self.assertNotIn("private notes", result.stdout)
            self.assertEqual((local / "script.py").read_text(), "print('local')\n")

    def test_identical_clone_returns_zero(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            public = root / "public-source"
            local = root / "local-working"
            public.mkdir()
            git(public, "init", "-b", "main")
            git(public, "config", "user.name", "Test")
            git(public, "config", "user.email", "test@example.com")
            (public / "script.py").write_text("print('same')\n", encoding="utf-8")
            git(public, "add", ".")
            git(public, "commit", "-m", "public")
            subprocess.run(
                ["git", "clone", "--quiet", str(public), str(local)], check=True
            )

            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    str(local),
                    "--public-repository",
                    str(public),
                    "--skip-environment",
                ],
                check=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("byte-identical", result.stdout)


if __name__ == "__main__":
    unittest.main()
