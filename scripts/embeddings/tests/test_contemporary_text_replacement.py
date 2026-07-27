from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
EMBEDDING_SCRIPTS = HERE.parent
REPO_ROOT = HERE.parents[2]
WRAPPER = REPO_ROOT / "hpc_jobs/active/hpc_contemporary_text_replacement_finalize.sh"
INITIALIZER = EMBEDDING_SCRIPTS / "initialize_contemporary_embedding_state.sh"
sys.path.insert(0, str(EMBEDDING_SCRIPTS))


def load_module(name: str, path: Path):
    specification = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(specification)
    assert specification.loader is not None
    specification.loader.exec_module(module)
    return module


replacement = load_module(
    "compose_contemporary_text_replacement",
    EMBEDDING_SCRIPTS / "compose_contemporary_text_replacement.py",
)
archive_manager = load_module(
    "manage_embedding_archive_for_replacement_test",
    EMBEDDING_SCRIPTS / "manage_embedding_archive.py",
)


class ContemporaryTextReplacementTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.config = self.root / "config.json"
        self.policy = self.root / "policy.json"
        self.plan = self.root / "plan"
        self.plan.mkdir()
        self.input_acquisition = self.root / "input_acquisition.tsv"
        self.input_acquisition.write_text("artifact\tsha256\nsource\tvalue\n")
        self.config.write_text(
            json.dumps(
                {
                    "schema_version": 1,
                    "modalities": {
                        "sequence": {"directory": "prott5", "dimension": 3},
                        "text": {
                            "directory": "exp_text_embeddings_temporal",
                            "dimension": 2,
                        },
                        "structure": {"directory": "IF1", "dimension": 4},
                        "ppi": {"directory": "ppi", "dimension": 5},
                    }
                }
            )
        )
        self.policy.write_text(
            json.dumps(
                {
                    "schema_version": 1,
                    "modalities": {
                        "sequence": {
                            "cache_directory": "prott5",
                            "dimension": 3,
                            "min_accepted_fraction": 1.0,
                        },
                        "text": {
                            "cache_directory": "exp_text_embeddings_temporal",
                            "dimension": 2,
                            "min_accepted_fraction": 0.5,
                        },
                        "structure": {
                            "cache_directory": "IF1",
                            "dimension": 4,
                            "min_accepted_fraction": 0.5,
                        },
                        "ppi": {
                            "cache_directory": "ppi",
                            "dimension": 5,
                            "min_accepted_fraction": 0.5,
                        },
                    }
                }
            )
        )
        targets = (("P1", "AAAA"), ("P2", "BBBB"))
        header = "protein_id\tsequence\tsequence_sha256\n"
        rows = "".join(
            f"{protein_id}\t{sequence}\t"
            f"{hashlib.sha256(sequence.encode()).hexdigest()}\n"
            for protein_id, sequence in targets
        )
        (self.plan / "reuse_proteins.tsv").write_text(header + rows)
        (self.plan / "regenerate_proteins.tsv").write_text(header)

        self.base_root = self.root / "base_final"
        self.base_root.mkdir()
        base_cache = self.root / "base_cache"
        self._create_cache_dirs(base_cache)
        for protein_id, offset in (("P1", 1), ("P2", 2)):
            self._save(base_cache / "prott5" / f"{protein_id}.npy", 3, 10 + offset)
            self._save(base_cache / "IF1" / f"{protein_id}.npy", 4, 20 + offset)
            self._save(base_cache / "ppi" / f"{protein_id}.npy", 5, 30 + offset)
            self._save(
                base_cache / "exp_text_embeddings_temporal" / f"{protein_id}.npy",
                2,
                40 + offset,
            )
        base_archive = self.base_root / "contemporary_embedding_cache.tar.gz"
        base_report = archive_manager.create_archive(
            base_cache, base_archive, self.config
        )
        (self.base_root / "FINAL_CACHE_COMPLETE.json").write_text(
            json.dumps(
                {
                    "complete": True,
                    "validated": True,
                    "archive_name": base_archive.name,
                    "archive_sha256": base_report["archive_sha256"],
                }
            )
        )

        self.text_run = self.root / "text_run"
        (self.text_run / "artifacts").mkdir(parents=True)
        replacement_cache = self.root / "replacement_cache"
        (replacement_cache / "exp_text_embeddings_temporal").mkdir(parents=True)
        self._save(
            replacement_cache / "exp_text_embeddings_temporal/P1.npy", 2, 99
        )
        text_archive = self.text_run / "artifacts/corrected_text.tar.gz"
        text_report = archive_manager.create_archive_from_directories(
            replacement_cache,
            text_archive,
            ("exp_text_embeddings_temporal",),
        )
        (self.text_run / "TEXT_GENERATION_COMPLETE.json").write_text(
            json.dumps(
                {
                    "complete": True,
                    "mode": "full-text-generation-only",
                    "requested_cutoff": "2025-03-08",
                    "effective_cutoff": "2025-03-08",
                    "old_text_carried_forward": False,
                    "hydration_performed": False,
                    "state_modified": False,
                    "target_count": 2,
                    "text_available": 1,
                    "archive": "artifacts/corrected_text.tar.gz",
                    "archive_sha256": text_report["archive_sha256"],
                }
            )
        )

    def tearDown(self) -> None:
        self.temporary.cleanup()

    @staticmethod
    def _create_cache_dirs(root: Path) -> None:
        for name in ("prott5", "exp_text_embeddings_temporal", "IF1", "ppi"):
            (root / name).mkdir(parents=True)

    @staticmethod
    def _save(path: Path, dimension: int, value: int) -> None:
        np.save(path, np.full(dimension, value, dtype=np.float32))

    def _arguments(self, output: Path) -> argparse.Namespace:
        return argparse.Namespace(
            base_final_root=self.base_root,
            replacement_run_root=self.text_run,
            plan_dir=self.plan,
            policy=self.policy,
            config=self.config,
            input_acquisition=self.input_acquisition,
            expected_cutoff="2025-03-08",
            variant_name="text-cutoff-2025-03-08__ppi-paper-faithful",
            work_dir=self.root / f"work-{output.name}",
            output_root=output,
            report=None,
        )

    def test_replaces_text_without_changing_non_text(self) -> None:
        output = self.root / "composed"
        result = replacement.compose(self._arguments(output))
        self.assertTrue(result["old_text_carried_forward"] is False)
        self.assertEqual(result["removed_old_text_count"], 2)
        self.assertEqual(
            result["combined_counts"],
            {"ppi": 2, "sequence": 2, "structure": 2, "text": 1},
        )

        extracted = self.root / "extracted"
        archive_manager.extract_archive(
            output / "archive/contemporary_embedding_cache.tar.gz",
            extracted,
            self.config,
        )
        np.testing.assert_array_equal(
            np.load(extracted / "exp_text_embeddings_temporal/P1.npy"),
            np.full(2, 99, dtype=np.float32),
        )
        self.assertFalse((extracted / "exp_text_embeddings_temporal/P2.npy").exists())
        np.testing.assert_array_equal(
            np.load(extracted / "prott5/P2.npy"),
            np.full(3, 12, dtype=np.float32),
        )
        with gzip.open(
            output / "reports/assembly/embedding_assembly.tsv.gz", "rt"
        ) as handle:
            assembly = handle.read()
        self.assertIn("P1\ttext\tavailable\t2", assembly)
        self.assertIn("P2\ttext\tmissing\t2", assembly)
        self.assertTrue((output / "COMPOSITION_COMPLETE.json").is_file())

    def test_rejects_replacement_archive_hash_mismatch(self) -> None:
        marker_path = self.text_run / "TEXT_GENERATION_COMPLETE.json"
        marker = json.loads(marker_path.read_text())
        marker["archive_sha256"] = "0" * 64
        marker_path.write_text(json.dumps(marker))
        output = self.root / "should-not-publish"
        with self.assertRaisesRegex(ValueError, "SHA-256 mismatch"):
            replacement.compose(self._arguments(output))
        self.assertFalse(output.exists())


class ContemporaryTextReplacementWorkflowTests(unittest.TestCase):
    def test_wrapper_keeps_old_final_cache_immutable_and_uses_fresh_roots(self) -> None:
        source = WRAPPER.read_text(encoding="utf-8")
        self.assertIn("compose_contemporary_text_replacement.py", source)
        self.assertIn('BASELINE_ROOT="$VARIANT_ROOT/source_baseline"', source)
        self.assertIn('STATE_ROOT="$VARIANT_ROOT/retry_state"', source)
        self.assertIn('FINAL_ROOT="$VARIANT_ROOT/finalized_pfp_cache"', source)
        self.assertIn("--retire-source-embeddings", source)
        self.assertNotIn('rm -rf -- "$BASE_ROOT"', source)

    def test_initializer_binds_composition_and_corrected_text_provenance(self) -> None:
        source = INITIALIZER.read_text(encoding="utf-8")
        self.assertIn("replacement-composition=", source)
        self.assertIn("corrected-text-generation=", source)


if __name__ == "__main__":
    unittest.main()
