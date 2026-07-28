from __future__ import annotations

import csv
import gzip
import io
import hashlib
import json
import subprocess
import sys
import tarfile
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np


SCRIPT = Path(__file__).resolve().parents[1] / "generate_unambiguous_ppi_delta.py"


def digest(path: Path) -> str:
    value = hashlib.sha256()
    value.update(path.read_bytes())
    return value.hexdigest()


class UnambiguousPpiDeltaTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.details = self.root / "details.tsv.gz"
        self.summary = self.root / "summary.json"
        self.status = self.root / "pair_status.tsv"
        self.h5 = self.root / "string.h5"
        self.work = self.root / "work"
        self.output = self.root / "output"

        with self.status.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["protein_id", "modality", "state"])
            for protein_id in ("P1", "P2", "P3", "P4"):
                writer.writerow(
                    [
                        protein_id,
                        "ppi",
                        "accepted" if protein_id == "P1" else "needs_retry",
                    ]
                )
                writer.writerow([protein_id, "sequence", "accepted"])

        rows = [
            ("P1", "True", "False", "9606.S1", "UniProt_AC"),
            (
                "P2",
                "True",
                "False",
                "9606.S2",
                "Ensembl_UniProt;UniProt_AC",
            ),
            ("P3", "False", "True", "", ""),
            ("P4", "False", "False", "", ""),
            ("OUTSIDE", "True", "False", "9606.S3", "Ensembl_UniProt"),
        ]
        with gzip.open(self.details, "wt", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(
                [
                    "protein_id",
                    "policy",
                    "covered",
                    "resident_string_id_count",
                    "ambiguous",
                    "selected_string_id",
                    "selected_source",
                    "resident_string_ids",
                ]
            )
            for protein_id, covered, ambiguous, string_id, source in rows:
                writer.writerow(
                    [
                        protein_id,
                        "widened_unambiguous",
                        covered,
                        "1" if covered else "0",
                        ambiguous,
                        string_id,
                        source,
                        string_id,
                    ]
                )

        with h5py.File(self.h5, "w") as handle:
            species = handle.create_group("species").create_group("9606")
            species.create_dataset(
                "proteins",
                data=np.asarray([b"9606.S1", b"9606.S2", b"9606.S3"]),
            )
            species.create_dataset(
                "embeddings",
                data=np.asarray(
                    [
                        np.full(512, 1, dtype=np.float32),
                        np.full(512, 2, dtype=np.float32),
                        np.full(512, 3, dtype=np.float32),
                    ]
                ),
            )

        self.summary.write_text(
            json.dumps(
                {
                    "schema_name": "string-alias-policy-coverage-audit",
                    "policies": {
                        "widened_unambiguous": {
                            "description": "fixture",
                            "source_tokens": [
                                "Ensembl_HGNC_uniprot_ids",
                                "Ensembl_UniProt",
                                "Ensembl_flybase_gene_id",
                                "Ensembl_gene",
                                "UniProt_AC",
                                "UniProt_DR_FlyBase",
                                "UniProt_ID",
                            ],
                        }
                    },
                    "input_files": {"string_h5": {"sha256": digest(self.h5)}},
                }
            )
            + "\n",
            encoding="utf-8",
        )

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def command(
        self, *, delta_count: int = 1, final_count: int | None = None
    ) -> list[str]:
        if final_count is None:
            final_count = 1 + delta_count
        return [
            sys.executable,
            str(SCRIPT),
            "--policy-details",
            str(self.details),
            "--audit-summary",
            str(self.summary),
            "--base-pair-status",
            str(self.status),
            "--string-h5",
            str(self.h5),
            "--work-dir",
            str(self.work),
            "--output-root",
            str(self.output),
            "--expected-policy-details-sha256",
            digest(self.details),
            "--expected-audit-summary-sha256",
            digest(self.summary),
            "--expected-base-pair-status-sha256",
            digest(self.status),
            "--expected-string-h5-sha256",
            digest(self.h5),
            "--expected-target-count",
            "4",
            "--expected-base-count",
            "1",
            "--expected-delta-count",
            str(delta_count),
            "--expected-final-count",
            str(final_count),
            "--protein-chunk-size",
            "2",
            "--embedding-batch-size",
            "1",
        ]

    def rewrite_policy_row(self, protein_id: str, **updates: str) -> None:
        with gzip.open(self.details, "rt", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fieldnames = list(reader.fieldnames or [])
            rows = list(reader)
        for row in rows:
            if row["protein_id"] == protein_id:
                row.update(updates)
                break
        else:
            self.fail(f"Missing policy fixture row: {protein_id}")
        with gzip.open(self.details, "wt", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(rows)

    def test_publishes_only_missing_unambiguous_vectors(self) -> None:
        result = subprocess.run(
            self.command(), check=False, capture_output=True, text=True
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        marker = json.loads((self.output / "DELTA_COMPLETE.json").read_text())
        self.assertEqual(marker["base_accepted_count"], 1)
        self.assertEqual(marker["delta_count"], 1)
        self.assertEqual(marker["ambiguous_rejected_count"], 1)
        self.assertEqual(marker["base_overlap_count"], 0)
        self.assertTrue(marker["roundtrip_validated"])
        with tarfile.open(self.output / "ppi_delta.tar.gz", "r:gz") as archive:
            self.assertEqual(archive.getnames(), ["ppi/P2.npy"])
            extracted = archive.extractfile("ppi/P2.npy")
            self.assertIsNotNone(extracted)
            vector = np.load(io.BytesIO(extracted.read()), allow_pickle=False)
        np.testing.assert_array_equal(vector, np.full(512, 2, dtype=np.float32))

    def test_count_mismatch_publishes_nothing(self) -> None:
        result = subprocess.run(
            self.command(delta_count=2, final_count=2),
            check=False,
            capture_output=True,
            text=True,
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Delta count differs", result.stderr)
        self.assertFalse(self.output.exists())

    def test_changed_policy_details_are_rejected(self) -> None:
        command = self.command()
        self.details.write_bytes(self.details.read_bytes() + b"changed")
        result = subprocess.run(command, check=False, capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("policy details SHA-256 mismatch", result.stderr)
        self.assertFalse(self.output.exists())

    def test_widened_policy_must_retain_every_baseline_vector(self) -> None:
        self.rewrite_policy_row(
            "P1",
            covered="False",
            resident_string_id_count="0",
            selected_string_id="",
            selected_source="",
            resident_string_ids="",
        )
        result = subprocess.run(
            self.command(), check=False, capture_output=True, text=True
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Widened policy loses 1 accepted baseline", result.stderr)
        self.assertFalse(self.output.exists())

    def test_selected_source_must_belong_to_frozen_policy(self) -> None:
        self.rewrite_policy_row(
            "P2", selected_source="Ensembl_UniProt;Unreviewed_Alias"
        )
        result = subprocess.run(
            self.command(), check=False, capture_output=True, text=True
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("outside the frozen policy", result.stderr)
        self.assertFalse(self.output.exists())


if __name__ == "__main__":
    unittest.main()
