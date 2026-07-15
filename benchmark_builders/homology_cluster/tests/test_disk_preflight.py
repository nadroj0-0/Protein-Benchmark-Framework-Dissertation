from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
import tempfile
import unittest
from unittest.mock import patch

from homology_cluster_benchmark.models import ResolvedInput
from homology_cluster_benchmark.pipeline import _disk_preflight

from tests.helpers import fixture_config


def _input(name: str, path: Path, size: int) -> ResolvedInput:
    return ResolvedInput(
        name=name,
        resolved_path=path,
        source_url=None,
        release="fixture",
        size_bytes=size,
        sha256="a" * 64,
        expected_sha256=None,
        acquisition="local",
    )


class DiskPreflightTests(unittest.TestCase):
    def _inputs(self, root: Path) -> dict[str, ResolvedInput]:
        return {
            "uniref90_fasta": _input("uniref90_fasta", root / "uniref", 100),
            "idmapping": _input("idmapping", root / "mapping", 40),
            "uniprot_sprot_sequences": _input(
                "uniprot_sprot_sequences", root / "sprot", 30
            ),
            "goa": _input("goa", root / "goa", 50),
            "go_obo": _input("go_obo", root / "obo", 10),
        }

    def _config(self, root: Path, diagnostic: bool):
        return fixture_config(
            root / "output",
            root / "temp",
            diagnostic_pilot=diagnostic,
            scratch_safety_multiplier=1,
            mmseqs_work_multiplier=1,
            publication_safety_multiplier=1,
        )

    def test_diagnostic_pilot_records_but_does_not_enforce_speculative_estimate(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            work = root / "work"
            work.mkdir()
            usage = SimpleNamespace(total=1000, used=990, free=10)
            with patch(
                "homology_cluster_benchmark.pipeline.shutil.disk_usage",
                return_value=usage,
            ):
                report = _disk_preflight(
                    self._config(root, True), work, self._inputs(root)
                )
            self.assertFalse(report["estimate_enforced"])
            self.assertTrue(report["estimate_exceeds_available_space"])

    def test_non_pilot_still_enforces_configured_estimate(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            work = root / "work"
            work.mkdir()
            usage = SimpleNamespace(total=1000, used=990, free=10)
            with patch(
                "homology_cluster_benchmark.pipeline.shutil.disk_usage",
                return_value=usage,
            ):
                with self.assertRaisesRegex(OSError, "Scratch preflight failed"):
                    _disk_preflight(
                        self._config(root, False), work, self._inputs(root)
                    )


if __name__ == "__main__":
    unittest.main()
