from __future__ import annotations

import csv
import hashlib
import importlib.util
import io
import json
import sys
import tarfile
import tempfile
import unittest
from pathlib import Path

import numpy as np


SCRIPT = Path(__file__).resolve().parents[1] / "bind_embedding_archive_evidence.py"
WRAPPER = SCRIPT.parents[2] / "hpc_jobs" / "active" / "hpc_bind_embedding_archive_evidence.sh"
FINALIZER = SCRIPT.parents[2] / "hpc_jobs" / "active" / "hpc_homology_embedding_finalize.sh"
SPEC = importlib.util.spec_from_file_location("bind_embedding_archive_evidence", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class ArchiveEvidenceTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.benchmark = self.root / "benchmark"
        self.benchmark.mkdir()
        for aspect in ("bp", "cc", "mf"):
            for split, protein_id, sequence in (
                ("training", "TRAIN", "AAAA"),
                ("validation", "VALID", "CCCC"),
                ("test", "TEST", "GGGG"),
            ):
                path = self.benchmark / f"{aspect}-{split}.csv"
                with path.open("w", encoding="utf-8", newline="") as handle:
                    writer = csv.writer(handle)
                    writer.writerow(("proteins", "sequences", "GO:0000001"))
                    writer.writerow((protein_id, sequence, 1))
        self.config = self.root / "config.json"
        self.config.write_text(
            json.dumps(
                {
                    "modalities": {
                        "sequence": {"directory": "prott5", "dimension": 1024},
                        "text": {"directory": "exp_text_embeddings_temporal", "dimension": 768},
                        "structure": {"directory": "IF1", "dimension": 512},
                        "ppi": {"directory": "ppi", "dimension": 512},
                    }
                }
            )
        )

    def tearDown(self) -> None:
        self.temporary.cleanup()

    @staticmethod
    def array_bytes(dimension: int) -> bytes:
        output = io.BytesIO()
        np.save(output, np.ones((dimension,), dtype=np.float32), allow_pickle=False)
        return output.getvalue()

    def write_archive(self, *, bad_shape: bool = False, duplicate: bool = False) -> Path:
        path = self.root / "cache.tar.gz"
        directories = {
            "sequence": ("prott5", 1024),
            "text": ("exp_text_embeddings_temporal", 768),
            "structure": ("IF1", 512),
            "ppi": ("ppi", 512),
        }
        with tarfile.open(path, "w:gz") as archive:
            for protein_id in ("TRAIN", "VALID", "TEST"):
                for modality, (directory, dimension) in directories.items():
                    if modality == "ppi" and protein_id != "TRAIN":
                        continue
                    malformed = bad_shape and protein_id == "TEST" and modality == "text"
                    actual_dimension = 7 if malformed else dimension
                    data = self.array_bytes(actual_dimension)
                    names = [f"data/embedding_cache/{directory}/{protein_id}.npy"]
                    if duplicate and protein_id == "TRAIN" and modality == "ppi":
                        names.append(f"duplicate/{directory}/{protein_id}.npy")
                    for name in names:
                        info = tarfile.TarInfo(name)
                        info.size = len(data)
                        archive.addfile(info, io.BytesIO(data))
        return path

    def args(self, archive: Path, output: Path):
        return type(
            "Args",
            (),
            {
                "benchmark_dir": self.benchmark,
                "benchmark_id": "fixture",
                "archive": archive,
                "archive_sha256": hashlib.sha256(archive.read_bytes()).hexdigest(),
                "config": self.config,
                "framework_commit": "f" * 40,
                "pfp_commit": MODULE.PFP_COMMIT,
                "output_dir": output,
            },
        )()

    def test_binds_complete_and_missing_pairs_to_exact_csvs(self) -> None:
        archive = self.write_archive()
        output = self.root / "evidence"
        summary = MODULE.publish(self.args(archive, output))
        self.assertEqual(summary["target_count"], 3)
        self.assertEqual(summary["coverage"]["sequence"]["accepted"], 3)
        self.assertEqual(summary["coverage"]["ppi"]["accepted"], 1)
        contract = json.loads((output / "contract.json").read_text())
        recorded = contract.pop("contract_sha256")
        self.assertEqual(recorded, MODULE.canonical_sha256(contract))
        self.assertEqual(len(contract["benchmark_csvs"]), 9)
        with (output / "pair_status.tsv").open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 12)
        self.assertEqual(
            {(row["protein_id"], row["state"]) for row in rows if row["modality"] == "ppi"},
            {("TRAIN", "accepted"), ("VALID", "needs_retry"), ("TEST", "needs_retry")},
        )
        self.assertTrue((output / "RUN_COMPLETE.json").is_file())

    def test_conflicting_sequences_fail(self) -> None:
        path = self.benchmark / "mf-test.csv"
        path.write_text("proteins,sequences,GO:0000001\nTEST,TTTT,1\n")
        with self.assertRaisesRegex(ValueError, "Conflicting sequences"):
            MODULE.publish(self.args(self.write_archive(), self.root / "evidence"))

    def test_malformed_present_array_fails(self) -> None:
        with self.assertRaisesRegex(ValueError, "Wrong array shape"):
            MODULE.publish(self.args(self.write_archive(bad_shape=True), self.root / "evidence"))

    def test_duplicate_archive_pair_fails(self) -> None:
        with self.assertRaisesRegex(ValueError, "repeats target/modality pair"):
            MODULE.publish(self.args(self.write_archive(duplicate=True), self.root / "evidence"))

    def test_wrapper_retries_slow_san_visibility(self) -> None:
        text = WRAPPER.read_text(encoding="utf-8")
        self.assertIn("set -Eeuo pipefail", text)
        self.assertIn("wait_for_file", text)
        self.assertIn("add_mmfp_singularity_bind /SAN/bioinf/bmpfp", text)

        finalizer = FINALIZER.read_text(encoding="utf-8")
        self.assertIn("stage_ledger_with_retry", finalizer)
        self.assertIn('cp -a "$LEDGER_DIR/." "$STAGED_LEDGER/"', finalizer)
        self.assertIn('--ledger-dir "$STAGED_LEDGER"', finalizer)


if __name__ == "__main__":
    unittest.main()
