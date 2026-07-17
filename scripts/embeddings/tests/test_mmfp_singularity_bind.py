from __future__ import annotations

import subprocess
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
COMMON = REPO_ROOT / "scripts/reproduction_common.sh"


class MmfpSingularityBindTest(unittest.TestCase):
    def test_adds_existing_paths_once_without_overwriting_existing_binds(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            first = root / "first"
            second = root / "second"
            first.mkdir()
            second.mkdir()
            command = (
                f'source "{COMMON}"; '
                f'add_mmfp_singularity_bind "{first}"; '
                f'add_mmfp_singularity_bind "{first}"; '
                f'add_mmfp_singularity_bind "{second}"; '
                'printf "%s\\n" "$SINGULARITY_BINDPATH"'
            )
            result = subprocess.run(
                ["bash", "-c", command],
                check=True,
                capture_output=True,
                text=True,
            )
            self.assertEqual(
                result.stdout.strip(),
                f"{first.resolve()}:{first.resolve()},"
                f"{second.resolve()}:{second.resolve()}",
            )

    def test_rejects_a_missing_bind_source(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            missing = Path(temporary) / "missing"
            result = subprocess.run(
                [
                    "bash",
                    "-c",
                    f'source "{COMMON}"; add_mmfp_singularity_bind "{missing}"',
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("Cannot bind missing directory", result.stderr)


if __name__ == "__main__":
    unittest.main()
