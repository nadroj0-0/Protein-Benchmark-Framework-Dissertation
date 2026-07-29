from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.util
import json
import sys
import tarfile
import tempfile
import unittest
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
EMBEDDING_SCRIPTS = HERE.parent
REPO_ROOT = HERE.parents[2]
WRAPPER = REPO_ROOT / "hpc_jobs/active/hpc_contemporary_widened_ppi_finalize.sh"
INITIALIZER = EMBEDDING_SCRIPTS / "initialize_contemporary_embedding_state.sh"
sys.path.insert(0, str(EMBEDDING_SCRIPTS))


def load_module(name: str, path: Path):
    specification = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(specification)
    assert specification.loader is not None
    specification.loader.exec_module(module)
    return module


composer = load_module(
    "compose_contemporary_ppi_delta",
    EMBEDDING_SCRIPTS / "compose_contemporary_ppi_delta.py",
)
archive_manager = load_module(
    "manage_embedding_archive_for_ppi_composition_test",
    EMBEDDING_SCRIPTS / "manage_embedding_archive.py",
)


class ContemporaryPpiDeltaCompositionTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.config = self.root / "config.json"
        self.policy = self.root / "policy.json"
        self.plan = self.root / "plan"
        self.plan.mkdir()
        self.input_acquisition = self.root / "input_acquisition.tsv"
        self.input_acquisition.write_text("artifact\tsha256\nsource\tvalue\n")
        dimensions = {"sequence": 3, "text": 2, "structure": 4, "ppi": 5}
        directories = {
            "sequence": "prott5",
            "text": "exp_text_embeddings_temporal",
            "structure": "IF1",
            "ppi": "ppi",
        }
        self.config.write_text(
            json.dumps(
                {
                    "schema_version": 1,
                    "modalities": {
                        modality: {
                            "directory": directories[modality],
                            "dimension": dimensions[modality],
                        }
                        for modality in directories
                    },
                }
            )
        )
        self.policy.write_text(
            json.dumps(
                {
                    "schema_version": 1,
                    "modalities": {
                        modality: {
                            "cache_directory": directories[modality],
                            "dimension": dimensions[modality],
                            "min_accepted_fraction": 0.0,
                        }
                        for modality in directories
                    },
                }
            )
        )
        targets = (("P1", "AAAA"), ("P2", "BBBB"), ("P3", "CCCC"))
        header = "protein_id\tsequence\tsequence_sha256\n"
        rows = "".join(
            f"{protein_id}\t{sequence}\t{hashlib.sha256(sequence.encode()).hexdigest()}\n"
            for protein_id, sequence in targets
        )
        (self.plan / "reuse_proteins.tsv").write_text(header + rows)
        (self.plan / "regenerate_proteins.tsv").write_text(header)

        self.base_root = self.root / "base_final"
        self.base_root.mkdir()
        base_cache = self.root / "base_cache"
        for directory in directories.values():
            (base_cache / directory).mkdir(parents=True)
        for index, protein_id in enumerate(("P1", "P2", "P3"), start=1):
            self._save(base_cache / f"prott5/{protein_id}.npy", 3, 10 + index)
            self._save(base_cache / f"exp_text_embeddings_temporal/{protein_id}.npy", 2, 20 + index)
            self._save(base_cache / f"IF1/{protein_id}.npy", 4, 30 + index)
        self._save(base_cache / "ppi/P1.npy", 5, 41)
        base_archive = self.base_root / "contemporary_embedding_cache.tar.gz"
        base_report = archive_manager.create_archive(base_cache, base_archive, self.config)
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
        self.base_sha256 = archive_manager.sha256_file(base_archive)

        self.delta_root = self.root / "delta"
        self.delta_root.mkdir()
        delta_arrays = self.root / "delta_arrays"
        delta_arrays.mkdir()
        hashes: dict[str, str] = {}
        for protein_id, value in (("P2", 52), ("P3", 53)):
            path = delta_arrays / f"{protein_id}.npy"
            self._save(path, 5, value)
            hashes[protein_id] = archive_manager.sha256_file(path)
        delta_archive = self.delta_root / "ppi_delta.tar.gz"
        with tarfile.open(delta_archive, "w:gz") as handle:
            for path in sorted(delta_arrays.glob("*.npy")):
                handle.add(path, arcname=f"ppi/{path.name}")
        mapping = self.delta_root / "mapping.tsv.gz"
        with gzip.open(mapping, "wt", encoding="utf-8", newline="") as handle:
            handle.write("protein_id\tstring_id\tselected_source\tnpy_sha256\n")
            for protein_id in sorted(hashes):
                handle.write(
                    f"{protein_id}\t9606.{protein_id}\tEnsembl_UniProt\t{hashes[protein_id]}\n"
                )
        (self.delta_root / "DELTA_COMPLETE.json").write_text(
            json.dumps(
                {
                    "schema_name": "validated-unambiguous-ppi-delta",
                    "schema_version": 1,
                    "complete": True,
                    "policy": "widened_unambiguous",
                    "roundtrip_validated": True,
                    "target_count": 3,
                    "base_accepted_count": 1,
                    "delta_count": 2,
                    "final_union_count": 3,
                    "base_overlap_count": 0,
                    "archive": delta_archive.name,
                    "archive_sha256": archive_manager.sha256_file(delta_archive),
                    "mapping_report": mapping.name,
                    "mapping_report_sha256": archive_manager.sha256_file(mapping),
                }
            )
        )

    def tearDown(self) -> None:
        self.temporary.cleanup()

    @staticmethod
    def _save(path: Path, dimension: int, value: int) -> None:
        np.save(path, np.full(dimension, value, dtype=np.float32))

    def _arguments(self, output: Path) -> argparse.Namespace:
        return argparse.Namespace(
            base_final_root=self.base_root,
            delta_root=self.delta_root,
            plan_dir=self.plan,
            policy=self.policy,
            config=self.config,
            input_acquisition=self.input_acquisition,
            variant_name="text-cutoff-2025-03-08__ppi-widened-unambiguous",
            work_dir=self.root / f"work-{output.name}",
            output_root=output,
            expected_target_count=3,
            expected_base_ppi_count=1,
            expected_delta_count=2,
            expected_final_ppi_count=3,
            report=None,
        )

    def test_adds_delta_without_replacing_or_changing_other_modalities(self) -> None:
        output = self.root / "composed"
        result = composer.compose(self._arguments(output))
        self.assertEqual(result["delta_count"], 2)
        self.assertEqual(result["replacement_count"], 0)
        self.assertEqual(
            result["combined_counts"],
            {"ppi": 3, "sequence": 3, "structure": 3, "text": 3},
        )
        self.assertEqual(
            archive_manager.sha256_file(
                self.base_root / "contemporary_embedding_cache.tar.gz"
            ),
            self.base_sha256,
        )

        extracted = self.root / "extracted"
        archive_manager.extract_archive(
            output / "archive/contemporary_embedding_cache.tar.gz",
            extracted,
            self.config,
        )
        np.testing.assert_array_equal(
            np.load(extracted / "ppi/P1.npy"), np.full(5, 41, dtype=np.float32)
        )
        np.testing.assert_array_equal(
            np.load(extracted / "ppi/P2.npy"), np.full(5, 52, dtype=np.float32)
        )
        np.testing.assert_array_equal(
            np.load(extracted / "exp_text_embeddings_temporal/P3.npy"),
            np.full(2, 23, dtype=np.float32),
        )

    def test_rejects_delta_archive_hash_mismatch(self) -> None:
        marker_path = self.delta_root / "DELTA_COMPLETE.json"
        marker = json.loads(marker_path.read_text())
        marker["archive_sha256"] = "0" * 64
        marker_path.write_text(json.dumps(marker))
        output = self.root / "should-not-publish"
        with self.assertRaisesRegex(ValueError, "SHA-256 mismatch"):
            composer.compose(self._arguments(output))
        self.assertFalse(output.exists())


class ContemporaryPpiDeltaWorkflowTests(unittest.TestCase):
    def test_wrapper_publishes_separately_after_validation(self) -> None:
        source = WRAPPER.read_text(encoding="utf-8")
        self.assertIn("text-cutoff-2025-03-08__ppi-widened-unambiguous", source)
        self.assertIn("compose_contemporary_ppi_delta.py", source)
        self.assertIn('SCRATCH_VARIANT="$WORK/variant"', source)
        self.assertIn('mv "$PUBLICATION_STAGING" "$VARIANT_ROOT"', source)
        self.assertIn("Paper-faithful archive changed during hydration", source)
        self.assertIn("--expected-delta-count \"$EXPECTED_DELTA_COUNT\"", source)
        self.assertNotIn('rm -rf -- "$BASE_ROOT"', source)

    def test_initializer_binds_base_and_delta_provenance(self) -> None:
        source = INITIALIZER.read_text(encoding="utf-8")
        self.assertIn("base-final-contract=", source)
        self.assertIn("base-variant-completion=", source)
        self.assertIn("ppi-delta-completion=", source)
        self.assertIn("ppi-delta-mapping=", source)


if __name__ == "__main__":
    unittest.main()
