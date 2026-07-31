from __future__ import annotations

import csv
import json
import subprocess
import sys
import tarfile
import tempfile
import unittest
from io import BytesIO
from pathlib import Path


FRAMEWORK = Path(__file__).parents[3]
SCRIPT = FRAMEWORK / "scripts" / "diagnostics" / "build_cafa3_knowledge_state_census.py"
ROOTS = {"bp": "GO:0008150", "cc": "GO:0005575", "mf": "GO:0003674"}
CHILDREN = {"bp": "GO:1000001", "cc": "GO:2000001", "mf": "GO:3000001"}


def write_csv(path: Path, prefix: str, split: str) -> None:
    root = ROOTS[prefix]
    child = CHILDREN[prefix]
    if split == "test":
        rows = [
            (f"{prefix.upper()}_NK", {root, child}),
            (f"{prefix.upper()}_LK", {root}),
            (f"{prefix.upper()}_TOOFEW", {root, child}),
            (f"{prefix.upper()}_OTHER", {root, child}),
        ]
    else:
        rows = [(f"{prefix.upper()}_{split}", {root, child})]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["proteins", "sequences", root, child])
        for protein_id, positives in rows:
            writer.writerow([protein_id, "AAAA", int(root in positives), int(child in positives)])


def add_tar_bytes(archive: tarfile.TarFile, name: str, value: str) -> None:
    payload = value.encode("utf-8")
    info = tarfile.TarInfo(name)
    info.size = len(payload)
    archive.addfile(info, BytesIO(payload))


def write_archive(path: Path, overlap: bool = False) -> None:
    # Match the organizer artifact: benchmark20171115.tar is uncompressed and
    # contains the benchmark README, not the separate CAFA3 target package.
    with tarfile.open(path, "w:") as archive:
        add_tar_bytes(archive, "benchmark20171115/00README.txt", "official benchmark\n")
        for prefix in ROOTS:
            nk = f"{prefix.upper()}_NK\n"
            lk = f"{prefix.upper()}_LK\n{prefix.upper()}_TOOFEW\n"
            if overlap:
                lk += f"{prefix.upper()}_NK\n"
            add_tar_bytes(
                archive,
                f"benchmark20171115/lists/{prefix}o_all_type1.txt",
                nk,
            )
            add_tar_bytes(
                archive,
                f"benchmark20171115/lists/{prefix}o_all_type2.txt",
                lk,
            )
            add_tar_bytes(
                archive,
                f"benchmark20171115/lists/{prefix}o_all_typex.txt",
                nk + f"{prefix.upper()}_LK\n{prefix.upper()}_TOOFEW\n",
            )
            add_tar_bytes(
                archive,
                f"benchmark20171115/lists/too_few/{prefix}o_FIXTURE_type2.txt",
                f"{prefix.upper()}_TOOFEW\n",
            )
            add_tar_bytes(
                archive,
                f"benchmark20171115/lists/too_few/{prefix}o_EMPTY_type1.txt",
                "",
            )


class Cafa3KnowledgeStateCensusTests(unittest.TestCase):
    def run_fixture(
        self, overlap: bool = False
    ) -> tuple[
        subprocess.CompletedProcess[str],
        Path,
        tempfile.TemporaryDirectory[str],
    ]:
        temporary = tempfile.TemporaryDirectory()
        root = Path(temporary.name)
        csv_dir = root / "csvs"
        csv_dir.mkdir()
        for prefix in ROOTS:
            for split in ("training", "validation", "test"):
                write_csv(csv_dir / f"{prefix}-{split}.csv", prefix, split)
        archive = root / "benchmark20171115.tar"
        write_archive(archive, overlap=overlap)
        output = root / "output"
        completed = subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                "--published-csv-dir",
                str(csv_dir),
                "--official-cafa-archive",
                str(archive),
                "--output-dir",
                str(output),
            ],
            check=False,
            capture_output=True,
            text=True,
        )
        return completed, output, temporary

    def test_classifies_official_type1_type2_and_unclassified_rows(self):
        completed, output, temporary = self.run_fixture()
        self.addCleanup(temporary.cleanup)
        self.assertEqual(completed.returncode, 0, completed.stderr)
        report = json.loads((output / "cafa3_knowledge_state_census.json").read_text())
        for aspect in ("BPO", "CCO", "MFO"):
            states = report["aspects"][aspect]["states"]
            self.assertEqual(states["no_knowledge"]["proteins"], 1)
            self.assertEqual(states["limited_knowledge"]["proteins"], 2)
            self.assertEqual(states["unclassified_by_official_lists"]["proteins"], 1)
            self.assertEqual(states["limited_knowledge"]["root_only"], 1)
            self.assertEqual(states["no_knowledge"]["non_root_observed_truth"], 1)
        self.assertTrue((output / "RUN_COMPLETE.json").is_file())
        self.assertTrue((output / "output_manifest.json").is_file())
        with (output / "cafa3_test_knowledge_states.tsv").open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        too_few = [row for row in rows if row["protein_id"].endswith("_TOOFEW")]
        self.assertEqual(len(too_few), 3)
        self.assertTrue(all(row["official_list_group"] == "too_few" for row in too_few))

    def test_rejects_type1_type2_overlap(self):
        completed, output, temporary = self.run_fixture(overlap=True)
        self.addCleanup(temporary.cleanup)
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("type1/type2 lists overlap", completed.stderr)
        self.assertFalse(output.exists())

    def test_rejects_processed_archive_without_organizer_lists_with_actionable_error(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        root = Path(temporary.name)
        csv_dir = root / "csvs"
        csv_dir.mkdir()
        for prefix in ROOTS:
            for split in ("training", "validation", "test"):
                write_csv(csv_dir / f"{prefix}-{split}.csv", prefix, split)
        archive = root / "data-cafa.tar.gz"
        with tarfile.open(archive, "w:gz") as handle:
            add_tar_bytes(
                handle,
                "benchmark20171115/groundtruth/leafonly_all.txt",
                "BP_NK GO:1000001\n",
            )
        output = root / "output"
        completed = subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                "--published-csv-dir",
                str(csv_dir),
                "--official-cafa-archive",
                str(archive),
                "--output-dir",
                str(output),
            ],
            check=False,
            capture_output=True,
            text=True,
        )
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("benchmark20171115.tar", completed.stderr)
        self.assertIn("not the processed DeepGOPlus", completed.stderr)
        self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
