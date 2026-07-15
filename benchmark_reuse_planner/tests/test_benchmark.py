from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from pfp_benchmark_reuse.benchmark import (  # noqa: E402
    BenchmarkError,
    parse_benchmark,
    verify_input_identities,
)

from helpers import CSV_NAMES, rows_in, write_benchmark  # noqa: E402


class BenchmarkParsingTests(unittest.TestCase):
    def test_missing_required_csv_fails(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            directory = write_benchmark(Path(tmp) / "benchmark")
            (directory / "mf-test.csv").unlink()
            with self.assertRaisesRegex(BenchmarkError, "missing required CSVs: mf-test.csv"):
                parse_benchmark("target", directory)

    def test_identical_occurrences_across_csvs_and_rows_are_deduplicated(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            rows = {
                "bp-training.csv": [("P1", "Aa-*"), ("P1", "Aa-*")],
                "cc-test.csv": [("P1", "Aa-*")],
                "mf-validation.csv": [("P1", "Aa-*")],
            }
            benchmark = parse_benchmark("source", write_benchmark(Path(tmp) / "source", rows))
            self.assertEqual(list(benchmark.proteins), ["P1"])
            self.assertEqual(
                benchmark.proteins["P1"].memberships,
                ("bp-training.csv", "cc-test.csv", "mf-validation.csv"),
            )
            self.assertEqual(benchmark.duplicate_occurrences, 1)

    def test_conflicting_sequences_within_benchmark_fail(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            rows = {
                "bp-training.csv": [("P1", "AAAA")],
                "cc-test.csv": [("P1", "BBBB")],
            }
            directory = write_benchmark(Path(tmp) / "target", rows)
            with self.assertRaisesRegex(BenchmarkError, "conflicting sequences"):
                parse_benchmark("target", directory)

    def test_label_changes_across_csvs_do_not_create_sequence_conflict(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            rows = {
                "bp-training.csv": [("P1", "AAAA", "1")],
                "mf-test.csv": [("P1", "AAAA", "0")],
            }
            benchmark = parse_benchmark("target", write_benchmark(Path(tmp) / "target", rows))
            self.assertEqual(len(benchmark.proteins), 1)

    def test_contradictory_duplicate_rows_within_one_csv_fail(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            rows = rows_in(
                "bp-training.csv",
                ("P1", "AAAA", "1"),
                ("P1", "AAAA", "0"),
            )
            directory = write_benchmark(Path(tmp) / "target", rows)
            with self.assertRaisesRegex(BenchmarkError, "contradictory duplicate rows"):
                parse_benchmark("target", directory)

    def test_malformed_headers_rows_ids_sequences_and_labels_fail(self) -> None:
        cases = {
            "singular header": lambda path: write_benchmark(
                path, headers_by_file={"bp-training.csv": ("protein", "sequences", "GO:1")}
            ),
            "duplicate GO": lambda path: write_benchmark(
                path,
                rows_in("bp-training.csv", ("P1", "AAAA", ("1", "0"))),
                {"bp-training.csv": ("proteins", "sequences", "GO:1", "GO:1")},
            ),
            "bad GO": lambda path: write_benchmark(
                path, headers_by_file={"bp-training.csv": ("proteins", "sequences", "not-go")}
            ),
            "unsafe ID": lambda path: write_benchmark(
                path, rows_in("bp-training.csv", ("bad id", "AAAA"))
            ),
            "empty ID": lambda path: write_benchmark(
                path, rows_in("bp-training.csv", ("", "AAAA"))
            ),
            "empty sequence": lambda path: write_benchmark(
                path, rows_in("bp-training.csv", ("P1", ""))
            ),
            "malformed sequence": lambda path: write_benchmark(
                path, rows_in("bp-training.csv", ("P1", "AA1A"))
            ),
            "nonbinary": lambda path: write_benchmark(
                path, rows_in("bp-training.csv", ("P1", "AAAA", "2"))
            ),
            "row width": lambda path: write_benchmark(
                path,
                rows_in("bp-training.csv", ("P1", "AAAA", "1")),
                {
                    "bp-training.csv": (
                        "proteins",
                        "sequences",
                        "GO:0000001",
                        "GO:0000002",
                    )
                },
            ),
        }
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for index, (name, builder) in enumerate(cases.items()):
                with self.subTest(name=name):
                    directory = builder(root / str(index))
                    with self.assertRaises(BenchmarkError):
                        parse_benchmark("target", directory)

    def test_malformed_csv_quoting_and_empty_overall_benchmark_fail(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            malformed = write_benchmark(root / "malformed")
            (malformed / "bp-training.csv").write_text(
                'proteins,sequences,GO:1\n"P1,AAAA,1\n', encoding="utf-8"
            )
            with self.assertRaises(BenchmarkError):
                parse_benchmark("target", malformed)

            empty = write_benchmark(root / "empty", rows_by_file={})
            with self.assertRaisesRegex(BenchmarkError, "contains no proteins"):
                parse_benchmark("target", empty)

    def test_input_identity_recheck_detects_change(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            directory = write_benchmark(Path(tmp) / "benchmark")
            benchmark = parse_benchmark("target", directory)
            with (directory / CSV_NAMES[0]).open("a", encoding="utf-8") as handle:
                handle.write("P2,BBBB,1\n")
            with self.assertRaisesRegex(BenchmarkError, "changed after planning"):
                verify_input_identities((benchmark,))


if __name__ == "__main__":
    unittest.main()
