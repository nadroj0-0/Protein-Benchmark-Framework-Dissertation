import csv
import sys
import tempfile
import unittest
from pathlib import Path


SRC = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(SRC))

from pfp_embedding_inventory.benchmark import (  # noqa: E402
    BenchmarkError,
    parse_benchmark,
    temporal_text_role,
)

from helpers import make_contract, unique_rows_by_file, write_nine_csvs  # noqa: E402


class BenchmarkParsingTests(unittest.TestCase):
    def test_required_csv_missing_fails_loudly(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(BenchmarkError, "missing required CSVs"):
                parse_benchmark(Path(tmp), make_contract())

    def test_union_and_shared_protein_memberships(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            rows = unique_rows_by_file()
            rows[("cc", "test")].append(rows[("bp", "training")][0])
            write_nine_csvs(directory, rows_by_file=rows)
            benchmark = parse_benchmark(
                directory,
                make_contract("per-ontology-disjoint", "per-ontology-disjoint"),
            )
            self.assertEqual(len(benchmark.proteins), 9)
            shared = benchmark.proteins["P1"]
            self.assertEqual(shared.ontologies, {"BP", "CC"})
            self.assertEqual(shared.splits, {"training", "test"})
            self.assertEqual(len(shared.source_files), 2)
            self.assertEqual(temporal_text_role(shared), "mixed-current-and-test")

    def test_identical_duplicate_row_is_counted_once(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            write_nine_csvs(directory)
            path = directory / "bp-training.csv"
            with path.open("a", newline="") as handle:
                csv.writer(handle).writerow(["P1", "ACDE", "1"])
            benchmark = parse_benchmark(directory, make_contract())
            self.assertEqual(len(benchmark.proteins), 1)
            self.assertEqual(benchmark.duplicate_rows, 1)

    def test_contradictory_duplicate_row_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            write_nine_csvs(directory)
            path = directory / "bp-training.csv"
            with path.open("a", newline="") as handle:
                csv.writer(handle).writerow(["P1", "ACDE", "0"])
            with self.assertRaisesRegex(BenchmarkError, "contradictory duplicate"):
                parse_benchmark(directory, make_contract())

    def test_conflicting_sequence_across_files_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            rows = {(o, s): [("P1", "ACDE", "1")] for o in ("bp", "cc", "mf") for s in ("training", "validation", "test")}
            rows[("mf", "test")] = [("P1", "AAAA", "1")]
            write_nine_csvs(directory, rows_by_file=rows)
            with self.assertRaisesRegex(BenchmarkError, "conflicting sequences"):
                parse_benchmark(directory, make_contract())

    def test_empty_and_malformed_ids_and_sequences_fail(self):
        bad_values = [
            (("", "ACDE", "1"), "protein ID"),
            (("bad/id", "ACDE", "1"), "protein ID"),
            (("P1", "", "1"), "sequence"),
            (("P1", "AC DE", "1"), "sequence"),
        ]
        for bad_row, message in bad_values:
            with self.subTest(row=bad_row), tempfile.TemporaryDirectory() as tmp:
                directory = Path(tmp)
                rows = unique_rows_by_file()
                rows[("bp", "training")] = [bad_row]
                write_nine_csvs(directory, rows_by_file=rows)
                with self.assertRaisesRegex(BenchmarkError, message):
                    parse_benchmark(directory, make_contract())

    def test_train_test_overlap_respects_supplied_contract(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            rows = unique_rows_by_file()
            rows[("bp", "test")] = [rows[("bp", "training")][0]]
            write_nine_csvs(directory, rows_by_file=rows)
            with self.assertRaisesRegex(BenchmarkError, "contract.*violated"):
                parse_benchmark(
                    directory,
                    make_contract("per-ontology-disjoint", "per-ontology-disjoint"),
                )

    def test_custom_regex_cannot_enable_path_traversal(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            rows = unique_rows_by_file()
            rows[("bp", "training")] = [("../escape", "ACDE", "1")]
            write_nine_csvs(directory, rows_by_file=rows)
            contract = make_contract()
            contract = type(contract)(
                id_overlap=contract.id_overlap,
                sequence_overlap=contract.sequence_overlap,
                protein_id_pattern=r"^.*$",
                sequence_pattern=contract.sequence_pattern,
            )
            with self.assertRaisesRegex(BenchmarkError, "malformed protein ID"):
                parse_benchmark(directory, contract)


if __name__ == "__main__":
    unittest.main()
