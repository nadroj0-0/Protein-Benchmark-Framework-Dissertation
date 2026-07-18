#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import json
import shutil
import subprocess
import sys
import tarfile
import tempfile
import unittest
from pathlib import Path

import numpy as np
import torch


FRAMEWORK_ROOT = Path(__file__).resolve().parents[3]
MODEL_EXECUTION = FRAMEWORK_ROOT / "scripts" / "model_execution"
TEMPORAL_CONFIG = FRAMEWORK_ROOT / "configs" / "pfp_benchmark_run.temporal.json"
HOMOLOGY_CONFIG = FRAMEWORK_ROOT / "configs" / "pfp_benchmark_run.homology.json"
CAFA3_CONFIG = FRAMEWORK_ROOT / "configs" / "pfp_benchmark_run.cafa3.json"

sys.path.insert(0, str(MODEL_EXECUTION))
from prepare_pfp_benchmark import (  # noqa: E402
    bind_reference_archive,
    compare_prepared,
    require_production_homology_publication,
)


FAKE_PREPARER = r'''#!/usr/bin/env python3
import argparse, csv, json
from pathlib import Path
import numpy as np
import scipy.sparse as ssp

p = argparse.ArgumentParser()
p.add_argument("--cafa3-dir", required=True)
p.add_argument("--output-dir", required=True)
a = p.parse_args()
source, output = Path(a.cafa3_dir), Path(a.output_dir)
output.mkdir(parents=True, exist_ok=True)
for aspect, prefix in (("BPO", "bp"), ("CCO", "cc"), ("MFO", "mf")):
    terms = None
    for csv_split, pfp_split in (("training", "train"), ("validation", "valid"), ("test", "test")):
        with (source / f"{prefix}-{csv_split}.csv").open(newline="") as handle:
            rows = list(csv.reader(handle))
        header, values = rows[0], rows[1:]
        if terms is None:
            terms = header[2:]
        lookup = [header.index(term) for term in terms]
        names = np.array([row[0] for row in values], dtype=object)
        labels = np.array([[float(row[index]) for index in lookup] for row in values], dtype=np.float32)
        np.save(output / f"{aspect}_{pfp_split}_names.npy", names)
        ssp.save_npz(output / f"{aspect}_{pfp_split}_labels.npz", ssp.csr_matrix(labels))
        (output / f"{aspect}_{pfp_split}_sequences.json").write_text(
            json.dumps({row[0]: row[1] for row in values})
        )
    (output / f"{aspect}_go_terms.json").write_text(json.dumps(terms))
    (output / f"{aspect}_info.json").write_text(json.dumps({"aspect": aspect}))
'''


def run(
    command: list[str], expected: int = 0, contains: str | None = None
) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(command, text=True, capture_output=True, check=False)
    if result.returncode != expected:
        raise AssertionError(
            f"Expected status {expected}, got {result.returncode}\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
        )
    if contains is not None and contains not in result.stdout + result.stderr:
        raise AssertionError(
            f"Expected output to contain {contains!r}\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
        )
    return result


class ModelExecutionTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.benchmark = self.root / "benchmark"
        self.benchmark.mkdir()
        self.pfp = self.root / "PFP"
        (self.pfp / "scripts").mkdir(parents=True)
        (self.pfp / "scripts" / "prepare_cafa3_data.py").write_text(FAKE_PREPARER)
        (self.pfp / "train.py").write_text("# fake training entrypoint\n")
        self.obo = self.root / "go.obo"
        self.obo.write_text(
            "format-version: 1.2\n\n"
            "[Term]\nid: GO:0000001\nnamespace: biological_process\n\n"
            "[Term]\nid: GO:0000002\nnamespace: cellular_component\n\n"
            "[Term]\nid: GO:0000003\nnamespace: molecular_function\n"
        )
        terms = {"bp": "GO:0000001", "cc": "GO:0000002", "mf": "GO:0000003"}
        split_values = {
            "training": ("TRAIN", "AAAA", "1"),
            "validation": ("VALID", "CCCC", "1"),
            "test": ("TEST", "GGGG", "1"),
        }
        for prefix, term in terms.items():
            for split, row in split_values.items():
                with (self.benchmark / f"{prefix}-{split}.csv").open("w", newline="") as handle:
                    writer = csv.writer(handle)
                    writer.writerow(("proteins", "sequences", term))
                    writer.writerow(row)

    def tearDown(self) -> None:
        self.temp.cleanup()

    def prepare(
        self, config: Path = TEMPORAL_CONFIG, expected: int = 0, contains: str | None = None
    ) -> Path:
        data = self.root / f"data-{len(list(self.root.glob('data-*')))}"
        report = self.root / f"prepare-{data.name}.json"
        self.last_preparation_report = report
        command = [
                sys.executable,
                str(MODEL_EXECUTION / "prepare_pfp_benchmark.py"),
                "--benchmark-dir",
                str(self.benchmark),
                "--data-dir",
                str(data),
                "--obo-file",
                str(self.obo),
                "--pfp-root",
                str(self.pfp),
                "--config",
                str(config),
                "--report",
                str(report),
                "--log-dir",
                str(self.root / f"logs-{data.name}"),
            ]
        if config == HOMOLOGY_CONFIG:
            validation = self.root / "validation_report.json"
            manifest = self.root / "output_manifest.json"
            marker = self.root / "RUN_COMPLETE.json"
            validation.write_text('{"valid": true}\n')
            files = []
            for path in sorted(self.benchmark.glob("*.csv")):
                files.append(
                    {
                        "path": path.name,
                        "size_bytes": path.stat().st_size,
                        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                    }
                )
            files.append(
                {
                    "path": validation.name,
                    "size_bytes": validation.stat().st_size,
                    "sha256": hashlib.sha256(validation.read_bytes()).hexdigest(),
                }
            )
            manifest.write_text(
                json.dumps(
                    {"schema_version": 1, "payload_file_count": len(files), "files": files}
                )
                + "\n"
            )
            digest = hashlib.sha256(manifest.read_bytes()).hexdigest()
            marker.write_text(
                json.dumps(
                    {
                        "complete": True,
                        "manifest_sha256": digest,
                        "scientific_fingerprint": "a" * 64,
                    }
                )
                + "\n"
            )
            for evidence in (validation, manifest, marker):
                command.extend(("--benchmark-evidence", str(evidence)))
        run(command, expected=expected, contains=contains)
        return data

    def test_prepare_materializes_and_verifies_nine_csv_contract(self) -> None:
        data = self.prepare()
        self.assertTrue((data / "proteins.fasta").is_file())
        self.assertEqual(json.loads((data / "BPO_go_terms.json").read_text()), ["GO:0000001"])
        comparison = compare_prepared(data, data)
        self.assertTrue(comparison["passed"])
        self.assertEqual(comparison["artifacts_compared"], 30)
        self.assertEqual(len(comparison["reference_artifact_sha256"]), 30)
        self.assertEqual(len(comparison["reference_fingerprint"]), 64)

        reference = self.root / "published-reference"
        shutil.copytree(data, reference)
        names_path = reference / "BPO_train_names.npy"
        np.save(names_path, np.asarray(["ALTERED"], dtype=object))
        with self.assertRaisesRegex(ValueError, "BPO_train_names.npy"):
            compare_prepared(reference, data)

    def test_reference_archive_binds_extracted_prepared_artifacts(self) -> None:
        data = self.prepare()
        archive = self.root / "fixture_data_splits.tar.gz"
        with tarfile.open(archive, "w:gz") as handle:
            handle.add(data, arcname="data")
        policy = {
            "name": archive.name,
            "size_bytes": archive.stat().st_size,
            "checksum_algorithm": "md5",
            "checksum": hashlib.md5(archive.read_bytes()).hexdigest(),
        }
        binding = bind_reference_archive(data, archive, policy)
        self.assertEqual(binding["member_count"], 30)
        self.assertEqual(binding["checksum"], policy["checksum"])

        np.save(data / "BPO_train_names.npy", np.asarray(["ALTERED"], dtype=object))
        with self.assertRaisesRegex(ValueError, "not extracted"):
            bind_reference_archive(data, archive, policy)

    def test_unexpected_tenth_csv_is_rejected(self) -> None:
        (self.benchmark / "diagnostic.csv").write_text("a,b\n1,2\n")
        self.prepare(expected=1, contains="unexpected CSV")

    def test_inconsistent_go_header_is_rejected(self) -> None:
        path = self.benchmark / "bp-validation.csv"
        path.write_text("proteins,sequences,GO:0000009\nVALID,CCCC,1\n")
        self.prepare(expected=1, contains="GO columns differ")

    def test_temporal_policy_rejects_exact_sequence_leakage(self) -> None:
        for prefix in ("bp", "cc", "mf"):
            path = self.benchmark / f"{prefix}-test.csv"
            with path.open(newline="") as handle:
                rows = list(csv.reader(handle))
            rows[1][1] = "AAAA"
            with path.open("w", newline="") as handle:
                csv.writer(handle).writerows(rows)
        self.prepare(expected=1, contains="exact sequences")

    def test_cafa3_policy_normalizes_legacy_singular_header_without_editing_source(self) -> None:
        path = self.benchmark / "mf-test.csv"
        original = "protein,sequences,GO:0000003\nTEST,GGGG,1\n"
        path.write_text(original)
        data = self.prepare(CAFA3_CONFIG)
        self.assertTrue((data / "MFO_test_names.npy").is_file())
        self.assertEqual(path.read_text(), original)
        report = json.loads(next(self.root.glob("prepare-data-*.json")).read_text())
        self.assertEqual(report["header_compatibility_aliases"][0]["file"], "mf-test.csv")

    def test_homology_policy_retains_zero_training_positive_term(self) -> None:
        path = self.benchmark / "bp-training.csv"
        path.write_text("proteins,sequences,GO:0000001\nTRAIN,AAAA,0\n")
        policy = json.loads(HOMOLOGY_CONFIG.read_text())
        policy["benchmark_contract"].pop("domain_validator")
        policy["benchmark_contract"]["require_benchmark_evidence"] = False
        policy["benchmark_contract"]["required_benchmark_evidence_names"] = []
        fixture_config = self.root / "homology-fixture-config.json"
        fixture_config.write_text(json.dumps(policy) + "\n")
        data = self.prepare(fixture_config)
        self.assertEqual(json.loads((data / "BPO_go_terms.json").read_text()), ["GO:0000001"])
        report = json.loads(next(self.root.glob("prepare-data-*.json")).read_text())
        self.assertEqual(report["csvs"]["bp-training.csv"]["zero_positive_terms"], ["GO:0000001"])

    def test_production_homology_policy_rejects_fabricated_receipts(self) -> None:
        self.prepare(
            HOMOLOGY_CONFIG,
            expected=1,
            contains="Completion marker is missing",
        )

    def test_production_homology_policy_rejects_valid_nonproduction_scope(self) -> None:
        for publication in (
            {
                "benchmark_scope": "fixture-only",
                "production_eligible": False,
                "fixture_mode": True,
            },
            {
                "benchmark_scope": "diagnostic-pilot",
                "production_eligible": False,
                "fixture_mode": False,
            },
        ):
            with self.assertRaisesRegex(ValueError, "dissertation-production"):
                require_production_homology_publication(publication)

        require_production_homology_publication(
            {
                "benchmark_scope": "dissertation-production",
                "production_eligible": True,
                "fixture_mode": False,
            }
        )

    def make_cache(self, data: Path, include_nonsequence: bool = True) -> Path:
        cache = self.root / "cache"
        specifications = {
            "prott5": 1024,
            "exp_text_embeddings_temporal": 768,
            "IF1": 512,
            "ppi": 512,
        }
        ids = {"TRAIN", "VALID", "TEST"}
        for directory, dimension in specifications.items():
            target = cache / directory
            target.mkdir(parents=True)
            if directory != "prott5" and not include_nonsequence:
                continue
            selected = ids if directory == "prott5" else {"TRAIN"}
            for protein_id in selected:
                np.save(target / f"{protein_id}.npy", np.zeros(dimension, dtype=np.float32))
        return cache

    def validate_cache(
        self,
        data: Path,
        cache: Path,
        mode: str,
        expected: int = 0,
        config: Path = TEMPORAL_CONFIG,
        evidence: list[Path] | None = None,
        require_evidence: bool = False,
        ia_dir: Path | None = None,
        aspects: list[str] | None = None,
    ) -> dict:
        report = self.root / f"cache-{mode}.json"
        command = [
                sys.executable,
                str(MODEL_EXECUTION / "validate_pfp_embedding_cache.py"),
                "--data-dir",
                str(data),
                "--cache-root",
                str(cache),
                "--config",
                str(config),
                "--mode",
                mode,
                "--report",
                str(report),
                "--issues-tsv",
                str(self.root / f"issues-{mode}.tsv"),
                "--preparation-report",
                str(self.last_preparation_report),
            ]
        for path in evidence or []:
            command.extend(("--embedding-evidence", str(path)))
        if require_evidence:
            command.append("--require-embedding-evidence")
        if ia_dir:
            command.extend(("--ia-file-dir", str(ia_dir)))
        for aspect in aspects or []:
            command.extend(("--aspect", aspect))
        run(command, expected=expected)
        return json.loads(report.read_text())

    def make_embedding_state_evidence(
        self, cache: Path, wrong_benchmark_hash: bool = False
    ) -> list[Path]:
        targets = {"TEST": "GGGG", "TRAIN": "AAAA", "VALID": "CCCC"}
        target_lines = ["protein_id\tsequence_sha256"]
        for protein_id, sequence in sorted(targets.items()):
            target_lines.append(
                f"{protein_id}\t{hashlib.sha256(sequence.encode()).hexdigest()}"
            )
        target_path = self.root / "targets.tsv"
        target_path.write_text("\n".join(target_lines) + "\n")
        csv_records = []
        for prefix in ("bp", "cc", "mf"):
            for split in ("training", "validation", "test"):
                path = self.benchmark / f"{prefix}-{split}.csv"
                digest = hashlib.sha256(path.read_bytes()).hexdigest()
                if wrong_benchmark_hash and prefix == "bp" and split == "training":
                    digest = "0" * 64
                csv_records.append(
                    {"name": path.name, "sha256": digest, "size_bytes": path.stat().st_size}
                )
        policy = json.loads(
            (FRAMEWORK_ROOT / "configs" / "contemporary_embedding_resume.json").read_text()
        )
        contract = {
            "schema_version": 1,
            "benchmark_id": "fixture",
            "benchmark_csvs": csv_records,
            "targets": {
                "count": len(targets),
                "manifest_sha256": hashlib.sha256(target_path.read_bytes()).hexdigest(),
            },
            "pfp_commit": "1e04fd6d6d3c40458fd41ec1a881ed6e24de768e",
            "framework_commit": "f" * 40,
            "policy": policy,
            "policy_sha256": "1" * 64,
            "environment": None,
            "source_files": [],
            "runtime": {},
        }
        canonical = json.dumps(contract, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
        contract["contract_sha256"] = hashlib.sha256(canonical.encode()).hexdigest()
        contract_path = self.root / "contract.json"
        contract_path.write_text(json.dumps(contract) + "\n")
        coverage = {
            "contract_sha256": contract["contract_sha256"],
            "target_count": 3,
            "coverage": {
                "sequence": {"accepted": 3},
                "text": {"accepted": 1},
                "structure": {"accepted": 1},
                "ppi": {"accepted": 1},
            },
        }
        coverage_path = self.root / "coverage.json"
        coverage_path.write_text(json.dumps(coverage) + "\n")
        pair_status = self.root / "pair_status.tsv"
        lines = [
            "protein_id\tmodality\tstate\tsequence_sha256\tembedding_sha256\t"
            "attempts\tlatest_reason\tlatest_detail"
        ]
        cache_directories = {
            "sequence": "prott5",
            "text": "exp_text_embeddings_temporal",
            "structure": "IF1",
            "ppi": "ppi",
        }
        for protein_id, sequence in sorted(targets.items()):
            sequence_sha = hashlib.sha256(sequence.encode()).hexdigest()
            for modality in ("sequence", "text", "structure", "ppi"):
                accepted = modality == "sequence" or protein_id == "TRAIN"
                embedding_sha = ""
                if accepted:
                    embedding_path = (
                        cache / cache_directories[modality] / f"{protein_id}.npy"
                    )
                    embedding_sha = hashlib.sha256(embedding_path.read_bytes()).hexdigest()
                lines.append(
                    "\t".join(
                        (
                            protein_id,
                            modality,
                            "accepted" if accepted else "needs_retry",
                            sequence_sha,
                            embedding_sha,
                            "0",
                            "" if accepted else "not_attempted",
                            "",
                        )
                    )
                )
        pair_status.write_text("\n".join(lines) + "\n")
        return [coverage_path, contract_path, target_path, pair_status]

    def test_full_cache_allows_missing_nonsequence_and_reports_masks(self) -> None:
        data = self.prepare()
        report = self.validate_cache(data, self.make_cache(data), "full")
        self.assertEqual(report["modalities"]["sequence"]["coverage"], 1.0)
        self.assertEqual(report["modalities"]["ppi"]["valid"], 1)
        self.assertEqual(report["modalities"]["ppi"]["missing"], 2)

    def test_malformed_present_array_fails(self) -> None:
        data = self.prepare()
        cache = self.make_cache(data)
        np.save(cache / "IF1" / "TEST.npy", np.zeros((1, 512), dtype=np.float32))
        self.validate_cache(data, cache, "full", expected=1)

    def test_sequence_only_requires_only_complete_sequence_cache(self) -> None:
        data = self.prepare()
        report = self.validate_cache(data, self.make_cache(data, include_nonsequence=False), "sequence-only")
        self.assertEqual(set(report["modalities"]), {"sequence"})

    def test_aspect_selection_scopes_cache_targets(self) -> None:
        path = self.benchmark / "cc-test.csv"
        path.write_text("proteins,sequences,GO:0000002\nCCONLY,TTTT,1\n")
        data = self.prepare()
        report = self.validate_cache(
            data, self.make_cache(data), "full", aspects=["BPO"]
        )
        self.assertEqual(report["target_count"], 3)

    def test_missing_sequence_fails(self) -> None:
        data = self.prepare()
        cache = self.make_cache(data)
        (cache / "prott5" / "TEST.npy").unlink()
        self.validate_cache(data, cache, "full", expected=1)

    def test_float_array_that_overflows_float32_is_rejected(self) -> None:
        data = self.prepare()
        cache = self.make_cache(data)
        np.save(cache / "prott5" / "TEST.npy", np.full(1024, 1e100, dtype=np.float64))
        report = self.validate_cache(data, cache, "full", expected=1)
        self.assertEqual(
            report["modalities"]["sequence"]["reasons"]["non_finite_after_float32"], 1
        )

    def test_embedding_state_evidence_is_bound_to_csvs_targets_and_cache(self) -> None:
        data = self.prepare()
        cache = self.make_cache(data)
        evidence = self.make_embedding_state_evidence(cache)
        report = self.validate_cache(
            data, cache, "full", evidence=evidence, require_evidence=True
        )
        self.assertTrue(report["provenance_evidence_bound"])

    def test_embedding_state_evidence_rejects_other_benchmark_contract(self) -> None:
        data = self.prepare()
        cache = self.make_cache(data)
        evidence = self.make_embedding_state_evidence(cache, wrong_benchmark_hash=True)
        report = self.validate_cache(
            data, cache, "full", expected=1, evidence=evidence, require_evidence=True
        )
        self.assertFalse(report["provenance_evidence_bound"])

    def test_embedding_state_evidence_rejects_stable_array_substitution(self) -> None:
        data = self.prepare()
        cache = self.make_cache(data)
        evidence = self.make_embedding_state_evidence(cache)
        np.save(cache / "IF1" / "TRAIN.npy", np.ones(512, dtype=np.float32))
        report = self.validate_cache(
            data, cache, "full", expected=1, evidence=evidence, require_evidence=True
        )
        self.assertFalse(report["provenance_evidence_bound"])

    def test_cafa3_policy_requires_complete_precomputed_ia_files(self) -> None:
        data = self.prepare(CAFA3_CONFIG)
        cache = self.make_cache(data)
        self.validate_cache(data, cache, "full", expected=1, config=CAFA3_CONFIG)
        ia_dir = self.root / "ia"
        ia_dir.mkdir()
        for aspect, term in (
            ("BPO", "GO:0000001"),
            ("CCO", "GO:0000002"),
            ("MFO", "GO:0000003"),
        ):
            (ia_dir / f"{aspect}_ia.txt").write_text(f"{term}\t0.0\n")
        report = self.validate_cache(
            data, cache, "full", config=CAFA3_CONFIG, ia_dir=ia_dir
        )
        self.assertEqual(set(report["information_accretion"]), {"BPO", "CCO", "MFO"})

    def test_prepare_only_runner_publishes_completion_marker(self) -> None:
        work = self.root / "runner-work"
        output = self.root / "runner-output"
        run(
            [
                "bash",
                str(MODEL_EXECUTION / "run_pfp_benchmark.sh"),
                "--benchmark-id",
                "fixture",
                "--benchmark-dir",
                str(self.benchmark),
                "--obo-file",
                str(self.obo),
                "--pfp-root",
                str(self.pfp),
                "--work-dir",
                str(work),
                "--output-dir",
                str(output),
                "--config",
                str(TEMPORAL_CONFIG),
                "--execution-mode",
                "prepare-only",
                "--seed",
                "7",
                "--allow-dirty-framework",
                "--allow-unversioned-pfp",
            ]
        )
        self.assertTrue((output / "WORKFLOW_COMPLETE.json").is_file())
        self.assertTrue((output / "output_manifest.json").is_file())
        marker = json.loads((output / "WORKFLOW_COMPLETE.json").read_text())
        self.assertEqual(
            marker["manifest_sha256"],
            hashlib.sha256((output / "output_manifest.json").read_bytes()).hexdigest(),
        )
        run(
            [
                sys.executable,
                str(MODEL_EXECUTION / "manage_output_manifest.py"),
                "verify",
                "--root",
                str(output),
            ]
        )
        run_report = json.loads((output / "reports" / "run_report.json").read_text())
        self.assertEqual(run_report["status"], "passed")
        self.assertEqual(run_report["seed"], 7)

    def test_eval_only_adapter_requires_and_records_cafa_metrics(self) -> None:
        data = self.prepare()
        mmfp = self.pfp / "mmfp"
        mmfp.mkdir()
        (mmfp / "__init__.py").write_text("")
        (mmfp / "dataset.py").write_text(
            "import json, numpy as np, scipy.sparse as ssp, torch\n"
            "from torch.utils.data import Dataset\n"
            "class MultiModalDataset(Dataset):\n"
            "  def __init__(self, data_dir, embedding_dirs, seq_model, aspect, split, normalize='standard', norm_stats=None):\n"
            "    from pathlib import Path\n"
            "    root=Path(data_dir); self.protein_ids=np.load(root/f'{aspect}_{split}_names.npy', allow_pickle=True); self.labels=ssp.load_npz(root/f'{aspect}_{split}_labels.npz').toarray().astype(np.float32); self.norm_stats={} if norm_stats is None else norm_stats\n"
            "  def __len__(self): return len(self.protein_ids)\n"
            "  def __getitem__(self, i):\n"
            "    return {'seq':torch.zeros(1024),'seq_mask':torch.ones(1),'text':torch.zeros(768),'text_mask':torch.zeros(1),'struct':torch.zeros(512),'struct_mask':torch.zeros(1),'ppi':torch.zeros(512),'ppi_mask':torch.zeros(1),'labels':torch.from_numpy(self.labels[i]),'protein_id':str(self.protein_ids[i])}\n"
            "def collate_fn(batch):\n"
            "  return {k: ([x[k] for x in batch] if k=='protein_id' else torch.stack([x[k] for x in batch])) for k in batch[0]}\n"
        )
        (mmfp / "models.py").write_text(
            "import torch\n"
            "class Dummy(torch.nn.Module):\n"
            "  def __init__(self, num_go_terms, **kwargs): super().__init__(); self.bias=torch.nn.Parameter(torch.zeros(num_go_terms))\n"
            "  def forward(self, seq, seq_mask, text, text_mask, struct, struct_mask, ppi, ppi_mask):\n"
            "    batch=seq.shape[0]; return self.bias.expand(batch,-1), torch.zeros((batch,4)), None\n"
            "def create_model(**kwargs): return Dummy(**kwargs)\n"
        )
        (mmfp / "evaluation.py").write_text(
            "def evaluate_with_cafa(**kwargs): return {'fmax':0.6,'wfmax':0.5,'smin':1.2}\n"
        )
        (self.pfp / "train.py").write_text(
            "def set_seed(seed): pass\n"
            "def evaluate(model, loader, criterion, device, compute_weight_stats_detailed=False):\n"
            "  return {'fmax':0.5,'micro_auprc':0.4,'macro_auprc':0.3}\n"
        )
        checkpoint_root = self.root / "checkpoints"
        checkpoint = checkpoint_root / "fusion_comparison" / "prott5" / "BPO" / "gated_bilinear" / "best_model.pt"
        checkpoint.parent.mkdir(parents=True)
        torch.save({"bias": torch.zeros(1)}, checkpoint)
        cache = self.make_cache(data, include_nonsequence=False)
        output = self.root / "evaluation"
        run(
            [
                sys.executable,
                str(MODEL_EXECUTION / "evaluate_pfp_checkpoints.py"),
                "--pfp-root",
                str(self.pfp),
                "--data-dir",
                str(data),
                "--cache-root",
                str(cache),
                "--obo-file",
                str(self.obo),
                "--checkpoint-root",
                str(checkpoint_root),
                "--output-dir",
                str(output),
                "--config",
                str(TEMPORAL_CONFIG),
                "--mode",
                "sequence-only",
                "--aspect",
                "BPO",
            ]
        )
        result = json.loads((output / "BPO" / "results.json").read_text())
        self.assertEqual(result["cafa_fmax"], 0.6)
        self.assertEqual(result["cafa_wfmax"], 0.5)
        self.assertEqual(result["cafa_smin"], 1.2)

    def test_required_reference_comparison_cannot_pass_without_selected_aspect(self) -> None:
        preparation = self.root / "preparation.json"
        preparation.write_text(
            json.dumps({"status": "passed", "benchmark_fingerprint": "fixture"}) + "\n"
        )
        evaluation = self.root / "summary-evaluation" / "BPO"
        evaluation.mkdir(parents=True)
        (evaluation / "results.json").write_text(
            json.dumps(
                {
                    "checkpoint_sha256": "a" * 64,
                    "cafa_fmax": 0.6,
                    "cafa_wfmax": 0.5,
                    "cafa_smin": 1.2,
                }
            )
            + "\n"
        )
        expected = self.root / "expected.json"
        expected.write_text(json.dumps({"metrics": {"CCO": {"cafa_fmax": 0.6}}}) + "\n")
        report = self.root / "run-report.json"
        run(
            [
                sys.executable,
                str(MODEL_EXECUTION / "summarize_pfp_benchmark_run.py"),
                "--benchmark-id",
                "fixture",
                "--execution-mode",
                "eval-only",
                "--modality-mode",
                "full",
                "--seed",
                "42",
                "--preparation-report",
                str(preparation),
                "--evaluation-root",
                str(evaluation.parent),
                "--expected-metrics",
                str(expected),
                "--require-reference-match",
                "--aspect",
                "BPO",
                "--output-json",
                str(report),
                "--output-md",
                str(self.root / "run-report.md"),
            ],
            expected=1,
        )
        self.assertEqual(json.loads(report.read_text())["status"], "failed_reference_match")

    def test_non_finite_mandatory_metric_is_rejected(self) -> None:
        preparation = self.root / "nonfinite-preparation.json"
        preparation.write_text(
            json.dumps({"status": "passed", "benchmark_fingerprint": "fixture"}) + "\n"
        )
        evaluation = self.root / "nonfinite-evaluation" / "BPO"
        evaluation.mkdir(parents=True)
        (evaluation / "results.json").write_text(
            json.dumps(
                {
                    "checkpoint_sha256": "a" * 64,
                    "cafa_fmax": float("nan"),
                    "cafa_wfmax": 0.5,
                    "cafa_smin": 1.2,
                }
            )
            + "\n"
        )
        run(
            [
                sys.executable,
                str(MODEL_EXECUTION / "summarize_pfp_benchmark_run.py"),
                "--benchmark-id",
                "fixture",
                "--execution-mode",
                "eval-only",
                "--modality-mode",
                "full",
                "--seed",
                "42",
                "--preparation-report",
                str(preparation),
                "--evaluation-root",
                str(evaluation.parent),
                "--aspect",
                "BPO",
                "--output-json",
                str(self.root / "nonfinite-report.json"),
                "--output-md",
                str(self.root / "nonfinite-report.md"),
            ],
            expected=1,
            contains="non-finite",
        )


if __name__ == "__main__":
    unittest.main()
