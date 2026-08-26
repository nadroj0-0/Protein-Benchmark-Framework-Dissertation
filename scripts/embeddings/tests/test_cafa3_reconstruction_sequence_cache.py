import csv
import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "scripts/embeddings/prepare_cafa3_reconstruction_sequence_cache.py"


class Cafa3ReconstructionSequenceCacheTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.benchmark = self.root / "benchmark"
        self.cache = self.root / "cache"
        self.benchmark.mkdir()
        for directory in ("prott5", "exp_text_embeddings_temporal", "IF1", "ppi"):
            (self.cache / directory).mkdir(parents=True)
        self.proteins = {"A": "AAAA", "B": "BBBB", "C": "CCCC"}
        for aspect in ("bp", "cc", "mf"):
            for split, protein_id in zip(
                ("training", "validation", "test"), self.proteins, strict=True
            ):
                with (self.benchmark / f"{aspect}-{split}.csv").open(
                    "w", newline="", encoding="utf-8"
                ) as handle:
                    writer = csv.writer(handle)
                    writer.writerow(("proteins", "sequences", "GO:0000001"))
                    writer.writerow((protein_id, self.proteins[protein_id], 1))
        self.ledger = self.root / "output_checksums.sha256"
        self.ledger.write_text(
            "".join(
                f"{self.sha(path)}  /old/run/{path.name}\n"
                for path in sorted(self.benchmark.glob("*.csv"))
            ),
            encoding="utf-8",
        )
        source_sequences = {"A": "AAAA", "X": "BBBB"}
        self.source_targets = self.root / "targets.tsv"
        self.source_targets.write_text(
            "protein_id\tsequence_sha256\n"
            + "".join(
                f"{protein_id}\t{hashlib.sha256(sequence.encode()).hexdigest()}\n"
                for protein_id, sequence in source_sequences.items()
            ),
            encoding="utf-8",
        )
        np.save(self.cache / "prott5/A.npy", np.ones(1024, dtype=np.float32))
        np.save(self.cache / "prott5/X.npy", np.full(1024, 2, dtype=np.float32))
        for protein_id in ("A", "B", "C"):
            np.save(
                self.cache / "exp_text_embeddings_temporal" / f"{protein_id}.npy",
                np.ones(768, dtype=np.float32),
            )
            np.save(
                self.cache / "IF1" / f"{protein_id}.npy",
                np.ones(512, dtype=np.float32),
            )
            np.save(
                self.cache / "ppi" / f"{protein_id}.npy",
                np.ones(512, dtype=np.float32),
            )

    def tearDown(self) -> None:
        self.temporary.cleanup()

    @staticmethod
    def sha(path: Path) -> str:
        return hashlib.sha256(path.read_bytes()).hexdigest()

    def command(self, *arguments: str) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [sys.executable, str(SCRIPT), *arguments],
            check=False,
            capture_output=True,
            text=True,
        )

    def test_prepare_and_finalize_exact_reuse_boundary(self) -> None:
        actions = self.root / "actions.tsv"
        missing = self.root / "missing.fasta"
        target_manifest = self.root / "target_manifest.tsv"
        prepared = self.root / "prepared.json"
        result = self.command(
            "prepare",
            "--benchmark-dir",
            str(self.benchmark),
            "--benchmark-checksums",
            str(self.ledger),
            "--source-targets",
            str(self.source_targets),
            "--expected-source-targets-sha256",
            self.sha(self.source_targets),
            "--cache-root",
            str(self.cache),
            "--missing-fasta",
            str(missing),
            "--target-manifest",
            str(target_manifest),
            "--actions-tsv",
            str(actions),
            "--report",
            str(prepared),
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        payload = json.loads(prepared.read_text(encoding="utf-8"))
        self.assertEqual(
            payload["actions"],
            {"direct-reuse": 1, "exact-sequence-reuse": 1, "generate": 1},
        )
        self.assertEqual(
            (self.cache / "prott5/B.npy").read_bytes(),
            (self.cache / "prott5/X.npy").read_bytes(),
        )
        self.assertEqual(missing.read_text(encoding="utf-8"), ">C\nCCCC\n")

        np.save(self.cache / "prott5/C.npy", np.full(1024, 3, dtype=np.float32))
        validated = self.root / "validated.json"
        result = self.command(
            "finalize",
            "--benchmark-dir",
            str(self.benchmark),
            "--benchmark-checksums",
            str(self.ledger),
            "--cache-root",
            str(self.cache),
            "--actions-tsv",
            str(actions),
            "--report",
            str(validated),
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(json.loads(validated.read_text())["status"], "passed")

    def test_prepare_rejects_source_id_with_changed_sequence(self) -> None:
        source = self.source_targets.read_text(encoding="utf-8").replace(
            hashlib.sha256(b"AAAA").hexdigest(),
            hashlib.sha256(b"DIFFERENT").hexdigest(),
        )
        self.source_targets.write_text(source, encoding="utf-8")
        result = self.command(
            "prepare",
            "--benchmark-dir",
            str(self.benchmark),
            "--benchmark-checksums",
            str(self.ledger),
            "--source-targets",
            str(self.source_targets),
            "--expected-source-targets-sha256",
            self.sha(self.source_targets),
            "--cache-root",
            str(self.cache),
            "--missing-fasta",
            str(self.root / "missing.fasta"),
            "--target-manifest",
            str(self.root / "target_manifest.tsv"),
            "--actions-tsv",
            str(self.root / "actions.tsv"),
            "--report",
            str(self.root / "prepared.json"),
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("different sequence", result.stderr)


if __name__ == "__main__":
    unittest.main()
