from __future__ import annotations

import csv
import gzip
import hashlib
import importlib.util
import json
import sys
import tarfile
import tempfile
import unittest
from io import BytesIO
from pathlib import Path

import numpy as np


SCRIPT = Path(__file__).resolve().parents[1] / "resolve_embedding_reuse_sources.py"
WRAPPER = SCRIPT.parents[2] / "hpc_jobs" / "active" / "hpc_homology_embedding_source_ledger.sh"
SPEC = importlib.util.spec_from_file_location("embedding_source_resolution", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def digest(path: Path) -> str:
    value = hashlib.sha256()
    value.update(path.read_bytes())
    return value.hexdigest()


def npy_bytes(values: np.ndarray) -> bytes:
    output = BytesIO()
    np.save(output, values, allow_pickle=False)
    return output.getvalue()


def write_archive(path: Path, arrays: dict[tuple[str, str], np.ndarray]) -> str:
    directory = {
        "sequence": "prott5",
        "text": "exp_text_embeddings_temporal",
        "structure": "IF1",
        "ppi": "ppi",
    }
    with tarfile.open(path, "w:gz") as archive:
        for (protein_id, modality), values in sorted(arrays.items()):
            data = npy_bytes(values)
            info = tarfile.TarInfo(
                f"data/embedding_cache/{directory[modality]}/{protein_id}.npy"
            )
            info.size = len(data)
            archive.addfile(info, BytesIO(data))
    return digest(path)


def write_plan(root: Path, reuse: list[dict[str, str]], regenerate: list[dict[str, str]]) -> None:
    columns = [
        "protein_id",
        "sequence",
        "sequence_sha256",
        "action",
        "reason",
        "matching_embedded_benchmarks",
        "embedded_benchmark_memberships",
        "target_memberships",
        "regenerate_modalities",
    ]
    for name, rows in (("reuse_proteins.tsv", reuse), ("regenerate_proteins.tsv", regenerate)):
        with (root / name).open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, columns, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
    files = []
    for name in ("reuse_proteins.tsv", "regenerate_proteins.tsv"):
        path = root / name
        files.append(
            {"path": name, "size_bytes": path.stat().st_size, "sha256": digest(path)}
        )
    manifest = root / "output_manifest.json"
    manifest.write_text(
        json.dumps({"schema_version": 1, "files": files}) + "\n", encoding="utf-8"
    )
    (root / "RUN_COMPLETE.json").write_text(
        json.dumps(
            {
                "complete": True,
                "output_manifest": {
                    "path": "output_manifest.json",
                    "size_bytes": manifest.stat().st_size,
                    "sha256": digest(manifest),
                },
            }
        )
        + "\n",
        encoding="utf-8",
    )


def plan_row(protein_id: str, action: str, benchmarks: list[str]) -> dict[str, str]:
    sequence = "MPEPTIDE" + protein_id[-1]
    return {
        "protein_id": protein_id,
        "sequence": sequence,
        "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
        "action": action,
        "reason": "exact-id-sequence-match" if action == "reuse" else "protein-id-absent",
        "matching_embedded_benchmarks": json.dumps(sorted(benchmarks)),
        "embedded_benchmark_memberships": "[]",
        "target_memberships": '["bp-test.csv"]',
        "regenerate_modalities": "[]" if action == "reuse" else json.dumps(list(MODULE.MODALITIES)),
    }


def all_arrays(protein_id: str, offset: float = 0.0) -> dict[tuple[str, str], np.ndarray]:
    return {
        (protein_id, modality): np.full(
            (MODULE.EXPECTED_DIMENSIONS[modality],), offset + index, dtype=np.float32
        )
        for index, modality in enumerate(MODULE.MODALITIES, start=1)
    }


def read_gzip_tsv(path: Path) -> list[dict[str, str]]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class SourceResolutionTests(unittest.TestCase):
    def test_hpc_wrapper_pins_both_validated_sources_and_preserves_exit_status(self):
        text = WRAPPER.read_text(encoding="utf-8")
        self.assertIn("set -Eeuo pipefail", text)
        self.assertIn(
            "contemporary_paper_faithful=contemporary_hydrated_population=", text
        )
        self.assertIn(
            "cafa3_regenerated_hydrated=cafa3_hydrated_population=", text
        )
        self.assertIn("status=${PIPESTATUS[0]}", text)
        self.assertIn("reuse_ledger/source_resolved", text)

    def test_identical_duplicates_reuse_and_conflict_regenerates_only_one_modality(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            plan = root / "plan"
            plan.mkdir()
            write_plan(plan, [plan_row("P1", "reuse", ["a", "b"])], [])
            first_arrays = all_arrays("P1")
            second_arrays = all_arrays("P1")
            second_arrays[("P1", "text")] = np.full((768,), 99, dtype=np.float32)
            first = root / "first.tar.gz"
            second = root / "second.tar.gz"
            first_sha = write_archive(first, first_arrays)
            second_sha = write_archive(second, second_arrays)
            sources = [
                MODULE.parse_source(f"first=a={first}={first_sha}", 0),
                MODULE.parse_source(f"second=b={second}={second_sha}", 1),
            ]
            output = root / "result"
            MODULE.publish_resolution(plan, output, sources)

            pairs = {row["modality"]: row for row in read_gzip_tsv(output / "resolved_embedding_pairs.tsv.gz")}
            self.assertEqual(pairs["sequence"]["action"], "reuse")
            self.assertEqual(pairs["sequence"]["reason"], "identical-source-arrays")
            self.assertEqual(pairs["sequence"]["selected_source"], "first")
            self.assertEqual(pairs["text"]["action"], "regenerate")
            self.assertEqual(pairs["text"]["reason"], "conflicting-source-arrays")
            self.assertEqual(len(read_gzip_tsv(output / "conflicting_embedding_pairs.tsv.gz")), 1)
            with (output / "regenerate_proteins.tsv").open(
                encoding="utf-8", newline=""
            ) as handle:
                proteins = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(proteins), 1)
            self.assertEqual(json.loads(proteins[0]["regenerate_modalities"]), ["text"])
            self.assertTrue((output / "RUN_COMPLETE.json").is_file())

    def test_missing_invalid_and_coarse_regenerate_are_conservative(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            plan = root / "plan"
            plan.mkdir()
            write_plan(
                plan,
                [plan_row("P1", "reuse", ["a"])],
                [plan_row("P2", "regenerate", [])],
            )
            arrays = all_arrays("P1")
            del arrays[("P1", "ppi")]
            arrays[("P1", "structure")] = np.ones((7,), dtype=np.float32)
            archive = root / "source.tar.gz"
            archive_sha = write_archive(archive, arrays)
            output = root / "result"
            MODULE.publish_resolution(
                plan,
                output,
                [MODULE.parse_source(f"source=a={archive}={archive_sha}", 0)],
            )
            pairs = {(row["protein_id"], row["modality"]): row for row in read_gzip_tsv(output / "resolved_embedding_pairs.tsv.gz")}
            self.assertEqual(pairs[("P1", "ppi")]["reason"], "no-valid-source-array")
            self.assertEqual(pairs[("P1", "structure")]["reason"], "invalid-source-array")
            self.assertTrue(
                all(pairs[("P2", modality)]["action"] == "regenerate" for modality in MODULE.MODALITIES)
            )

    def test_nonmatching_benchmark_source_is_not_a_candidate(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            plan = root / "plan"
            plan.mkdir()
            write_plan(plan, [plan_row("P1", "reuse", ["a"])], [])
            archive = root / "source.tar.gz"
            archive_sha = write_archive(archive, all_arrays("P1"))
            output = root / "result"
            MODULE.publish_resolution(
                plan,
                output,
                [MODULE.parse_source(f"source=b={archive}={archive_sha}", 0)],
            )
            pairs = read_gzip_tsv(output / "resolved_embedding_pairs.tsv.gz")
            self.assertTrue(all(row["reason"] == "no-valid-source-array" for row in pairs))
            self.assertEqual(read_gzip_tsv(output / "source_candidates.tsv.gz"), [])

    def test_archive_hash_mismatch_fails_before_publication(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            plan = root / "plan"
            plan.mkdir()
            write_plan(plan, [plan_row("P1", "reuse", ["a"])], [])
            archive = root / "source.tar.gz"
            write_archive(archive, all_arrays("P1"))
            output = root / "result"
            source = MODULE.SourceSpec("source", "a", archive, "0" * 64, 0)
            with self.assertRaisesRegex(MODULE.ResolutionError, "hash mismatch"):
                MODULE.publish_resolution(plan, output, [source])
            self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
