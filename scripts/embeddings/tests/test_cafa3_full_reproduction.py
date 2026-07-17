from __future__ import annotations

import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import scipy.sparse as sparse


REPO_ROOT = Path(__file__).resolve().parents[3]


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


PREFLIGHT = load_module(
    "prepare_cafa3_embedding_preflight_test",
    REPO_ROOT / "scripts/embeddings/prepare_cafa3_embedding_preflight.py",
)
REPORT = load_module(
    "build_cafa3_full_reproduction_report_test",
    REPO_ROOT / "scripts/diagnostics/build_cafa3_full_reproduction_report.py",
)


class PreflightRestoreTest(unittest.TestCase):
    def test_preflight_restore_is_byte_exact(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            data = root / "data"
            backup = root / "backup"
            data.mkdir()
            for aspect in PREFLIGHT.ASPECTS:
                for split in PREFLIGHT.SPLITS:
                    ids = np.asarray(
                        [f"{aspect}_{split}_A", f"{aspect}_{split}_B"], dtype=object
                    )
                    np.save(data / f"{aspect}_{split}_names.npy", ids)
                    sparse.save_npz(
                        data / f"{aspect}_{split}_labels.npz",
                        sparse.csr_matrix([[1, 0], [0, 1]], dtype=np.float32),
                    )
                    sequences = {str(ids[0]): "ACDE", str(ids[1]): "FGHI"}
                    (data / f"{aspect}_{split}_sequences.json").write_text(
                        json.dumps(sequences), encoding="utf-8"
                    )

            PREFLIGHT.write_fasta(data)

            before = {
                path.name: PREFLIGHT.sha256_file(path)
                for path in PREFLIGHT.backup_paths(data)
            }
            created = PREFLIGHT.create_preflight(data, backup, 1)
            self.assertEqual(created["preflight_unique_proteins"], 9)
            for aspect in PREFLIGHT.ASPECTS:
                for split in PREFLIGHT.SPLITS:
                    names = np.load(
                        data / f"{aspect}_{split}_names.npy", allow_pickle=True
                    )
                    self.assertEqual(len(names), 1)

            restored = PREFLIGHT.restore_full(data, backup)
            self.assertTrue(restored["restored"])
            after = {
                path.name: PREFLIGHT.sha256_file(path)
                for path in PREFLIGHT.backup_paths(data)
            }
            self.assertEqual(before, after)
            self.assertEqual(restored["restored_unique_proteins"], 18)


class ExplicitEmbeddingComparisonTest(unittest.TestCase):
    def test_explicit_cache_roots_report_exact_and_different(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            generated = root / "generated"
            published = root / "published"
            directories = {
                "prott5": 4,
                "exp_text_embeddings_temporal": 3,
                "IF1": 2,
                "ppi": 2,
            }
            for directory, dimension in directories.items():
                (generated / directory).mkdir(parents=True)
                (published / directory).mkdir(parents=True)
                array = np.arange(dimension, dtype=np.float32)
                np.save(generated / directory / "P1.npy", array)
                np.save(published / directory / "P1.npy", array)
            np.save(generated / "prott5/P2.npy", np.ones(4, dtype=np.float32))
            np.save(published / "prott5/P2.npy", np.zeros(4, dtype=np.float32))

            out_csv = root / "comparison.csv"
            out_json = root / "summary.json"
            command = [
                sys.executable,
                str(REPO_ROOT / "scripts/diagnostics/compare_embeddings.py"),
                "--generated-cache-root",
                str(generated),
                "--published-cache-root",
                str(published),
                "--out-csv",
                str(out_csv),
                "--out-json",
                str(out_json),
            ]
            subprocess.run(command, check=True, capture_output=True, text=True)
            summary = json.loads(out_json.read_text(encoding="utf-8"))
            self.assertEqual(summary["prott5"]["generated_count"], 2)
            self.assertEqual(summary["prott5"]["published_count"], 2)
            self.assertEqual(summary["prott5"]["statuses"]["exact_match"], 1)
            self.assertEqual(summary["prott5"]["statuses"]["different"], 1)
            self.assertEqual(
                summary["text_temporal"]["statuses"]["exact_match"], 1
            )


class ReportFormattingTest(unittest.TestCase):
    def test_report_retains_nonzero_evaluation_as_observation(self) -> None:
        report = {
            "complete": True,
            "provenance": {
                "pfp_commit": "abc",
                "framework_commit": "def",
                "text_cutoff_date": "2016-02-17",
            },
            "published_cache_discarded": True,
            "embedding_comparison": {
                "prott5": {
                    "generated_count": 1,
                    "published_count": 1,
                    "common_count": 1,
                    "statuses": {"exact_match": 1},
                }
            },
            "evaluation": {
                "results": [
                    {
                        "aspect": aspect,
                        "expected_fmax": values["fmax"],
                        "actual_fmax": values["fmax"] - 0.003,
                        "delta_fmax": -0.003,
                        "expected_wfmax": values["wfmax"],
                        "actual_wfmax": values["wfmax"],
                        "delta_wfmax": 0.0,
                        "passed": False,
                    }
                    for aspect, values in REPORT.EXPECTED.items()
                ]
            },
            "evaluation_exit_status": 1,
            "training": {
                aspect: {
                    "best_epoch": 3,
                    "total_epochs": 8,
                    "best_val_fmax": 0.5,
                    "test_fmax": 0.5,
                    "cafa_fmax": 0.5,
                    "cafa_wfmax": 0.4,
                }
                for aspect in REPORT.ASPECTS
            },
        }
        markdown = REPORT.build_markdown(report)
        self.assertIn("PFP evaluation exit status: `1`", markdown)
        self.assertIn("non-zero status is retained", markdown)
        self.assertIn("Published cache discarded after comparison: `true`", markdown)


class WorkflowContractTest(unittest.TestCase):
    def test_hpc_binds_embedding_state_before_mmfp_python_starts(self) -> None:
        path = REPO_ROOT / "hpc_jobs/active/hpc_cafa3_full_from_scratch_reproduction.sh"
        source = path.read_text(encoding="utf-8")

        bind = 'add_mmfp_singularity_bind "$(dirname "$EMBEDDING_STATE_ROOT")"'
        self.assertIn(bind, source)
        self.assertLess(source.index(bind), source.index("activate_or_create_mmfp_env"))

    def test_full_workflow_avoids_git_dash_c_for_morecambe(self) -> None:
        paths = [
            REPO_ROOT
            / "scripts/reproduction/run_cafa3_full_from_scratch_reproduction.sh",
            REPO_ROOT / "scripts/embeddings/generate_embeddings_dependencies.sh",
            REPO_ROOT
            / "hpc_jobs/active/hpc_cafa3_full_from_scratch_reproduction.sh",
        ]
        for path in paths:
            with self.subTest(path=path):
                source = path.read_text(encoding="utf-8")
                self.assertNotIn("git -C", source)
                self.assertIn("git_in_dir", source)

    def test_published_cache_is_compared_then_deleted_before_training(self) -> None:
        path = REPO_ROOT / "scripts/reproduction/run_cafa3_full_from_scratch_reproduction.sh"
        source = path.read_text(encoding="utf-8")
        compare_at = source.index("compare_embeddings.py")
        discard_at = source.index('rm -rf "$PUBLISHED_ROOT" "$ARCHIVE_STAGE"')
        train_at = source.index('"$PYTHON_BIN" train.py')
        self.assertLess(compare_at, discard_at)
        self.assertLess(discard_at, train_at)
        training_block = source[train_at : source.index("==> [12/13]")]
        self.assertNotIn("PUBLISHED_CACHE", training_block)

    def test_hpc_cleanup_is_unconditional_and_scratch_scoped(self) -> None:
        path = REPO_ROOT / "hpc_jobs/active/hpc_cafa3_full_from_scratch_reproduction.sh"
        source = path.read_text(encoding="utf-8")
        self.assertIn("trap cleanup EXIT", source)
        self.assertIn('WORK="/scratch0/cafa3_full_reproduction_${JOB_TOKEN}"', source)
        self.assertIn('rm -rf "$WORK"', source)
        self.assertIn('[[ -f "$staging/WORKFLOW_COMPLETE.json" ]]', source)

    def test_partial_generation_is_persisted_before_training(self) -> None:
        path = REPO_ROOT / "scripts/reproduction/run_cafa3_full_from_scratch_reproduction.sh"
        source = path.read_text(encoding="utf-8")
        merge_at = source.index("merge_command=(")
        incomplete_at = source.index("publish_incomplete_generation", merge_at)
        train_at = source.index('"$PYTHON_BIN" train.py')
        self.assertLess(merge_at, incomplete_at)
        self.assertLess(incomplete_at, train_at)
        self.assertIn("--embedding-mode initial|resume", source)
        self.assertIn("state_gate_passed", source)

    def test_retry_wrapper_is_modality_specific_and_always_cleans_scratch(self) -> None:
        workflow = (
            REPO_ROOT / "scripts/reproduction/run_cafa3_embedding_retry.sh"
        ).read_text(encoding="utf-8")
        wrapper = (
            REPO_ROOT / "hpc_jobs/active/hpc_cafa3_embedding_retry.sh"
        ).read_text(encoding="utf-8")
        self.assertIn("--requested-pairs \"$REQUESTED\"", workflow)
        self.assertIn("verify_embedding_subset_equivalence.py", workflow)
        self.assertIn("--allowed-extra-pairs \"$CONTROLS\"", workflow)
        self.assertIn("trap cleanup EXIT", wrapper)
        self.assertIn('WORK="/scratch0/cafa3_embedding_retry_${JOB_TOKEN}"', wrapper)
        self.assertIn('rm -rf "$WORK"', wrapper)

    def test_bounded_alphafold_mode_is_opt_in_and_pfp_source_is_unchanged(self) -> None:
        path = REPO_ROOT / "scripts/embeddings/generate_embeddings_structure.sh"
        source = path.read_text(encoding="utf-8")
        self.assertIn('ALPHAFOLD_ACQUISITION_MODE="${ALPHAFOLD_ACQUISITION_MODE:-pfp}"', source)
        self.assertIn("prefetch_alphafold_structures.py", source)
        self.assertIn("python scripts/check_alphafold_coverage.py", source)


if __name__ == "__main__":
    unittest.main()
