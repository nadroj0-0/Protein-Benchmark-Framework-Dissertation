from __future__ import annotations

import json
from pathlib import Path
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from homology_cluster_benchmark.external_clusters import (
    load_external_cluster_provenance,
    validate_external_cluster_counts,
)
from homology_cluster_benchmark.inputs import sha256_file

from tests.helpers import FIXTURES


class ExternalClusterProvenanceTests(unittest.TestCase):
    def _payload(self, assignments: Path) -> dict[str, object]:
        return {
            "schema_name": "homology-external-cluster-assignments",
            "schema_version": 1,
            "artifact": {
                "sha256": sha256_file(assignments),
                "size_bytes": assignments.stat().st_size,
                "members": 6,
                "clusters": 6,
            },
            "producer": {"name": "Daniel Buchan"},
            "input": {
                "uniref_level": 50,
                "release": "2026_02",
                "expected_fasta_sha256": sha256_file(FIXTURES / "uniref50.fasta"),
                "expected_records": 6,
            },
            "method": {
                "identity_percent": 30,
                "coverage": 0.8,
                "coverage_mode": 0,
                "cluster_mode": 0,
                "sensitivity": 4,
                "evalue": "MMseqs2 default",
                "createdb_shuffle": "MMseqs2 default",
                "cluster_reassignment": "MMseqs2 default",
            },
            "usage_policy": {
                "lineage": "supervisor-generated",
                "do_not_merge_with_framework_generated_cluster_cache": True,
            },
        }

    def test_valid_external_artifact_is_bound_to_frozen_uniref(self):
        assignments = FIXTURES / "clusters_uniref50.tsv"
        with tempfile.TemporaryDirectory() as tmp:
            provenance = Path(tmp) / "provenance.json"
            provenance.write_text(json.dumps(self._payload(assignments)))
            payload = load_external_cluster_provenance(
                provenance,
                assignments,
                identity=0.30,
                coverage=0.8,
                cov_mode=0,
                cluster_mode=0,
                sensitivity=4,
                uniref_level=50,
                uniref_release="2026_02",
                uniref_sha256=sha256_file(FIXTURES / "uniref50.fasta"),
                uniref_records=6,
            )
            validate_external_cluster_counts(payload, members=6, clusters=6)

    def test_method_or_count_mismatch_fails_loudly(self):
        assignments = FIXTURES / "clusters_uniref50.tsv"
        with tempfile.TemporaryDirectory() as tmp:
            provenance = Path(tmp) / "provenance.json"
            payload = self._payload(assignments)
            payload["method"]["sensitivity"] = 7.5  # type: ignore[index]
            provenance.write_text(json.dumps(payload))
            with self.assertRaisesRegex(ValueError, "sensitivity"):
                load_external_cluster_provenance(
                    provenance,
                    assignments,
                    identity=0.30,
                    coverage=0.8,
                    cov_mode=0,
                    cluster_mode=0,
                    sensitivity=4,
                    uniref_level=50,
                    uniref_release="2026_02",
                    uniref_sha256=sha256_file(FIXTURES / "uniref50.fasta"),
                    uniref_records=6,
                )

            payload = self._payload(assignments)
            provenance.write_text(json.dumps(payload))
            with self.assertRaisesRegex(ValueError, "counts"):
                validate_external_cluster_counts(payload, members=5, clusters=6)

    def test_assignment_hash_mismatch_fails_loudly(self):
        assignments = FIXTURES / "clusters_uniref50.tsv"
        with tempfile.TemporaryDirectory() as tmp:
            provenance = Path(tmp) / "provenance.json"
            payload = self._payload(assignments)
            payload["artifact"]["sha256"] = "0" * 64  # type: ignore[index]
            provenance.write_text(json.dumps(payload))
            with self.assertRaisesRegex(ValueError, "SHA-256 mismatch"):
                load_external_cluster_provenance(
                    provenance,
                    assignments,
                    identity=0.30,
                    coverage=0.8,
                    cov_mode=0,
                    cluster_mode=0,
                    sensitivity=4,
                    uniref_level=50,
                    uniref_release="2026_02",
                    uniref_sha256=sha256_file(FIXTURES / "uniref50.fasta"),
                    uniref_records=6,
                )


if __name__ == "__main__":
    unittest.main()
