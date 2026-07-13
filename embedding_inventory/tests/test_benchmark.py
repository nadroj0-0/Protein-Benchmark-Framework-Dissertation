import csv
import sys
import tempfile
import unittest
from pathlib import Path


SRC = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(SRC))

from pfp_embedding_inventory.benchmark import (  # noqa: E402
    BenchmarkError,
    duplicate_row_digest,
    load_aliases,
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

    def test_cross_ontology_target_id_leakage_is_rejected_globally(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            rows = unique_rows_by_file()
            rows[("cc", "test")].append(rows[("bp", "training")][0])
            write_nine_csvs(directory, rows_by_file=rows)
            with self.assertRaisesRegex(BenchmarkError, "global-evaluation.*protein IDs"):
                parse_benchmark(
                    directory,
                    make_contract("global-evaluation-disjoint", "allow"),
                )

    def test_homology_contract_rejects_exact_sequence_under_different_ids(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            rows = unique_rows_by_file()
            rows[("bp", "training")] = [("TRAIN_ID", "ACDEFG", "1")]
            rows[("mf", "test")] = [("TEST_ID", "ACDEFG", "1")]
            write_nine_csvs(directory, rows_by_file=rows)
            with self.assertRaisesRegex(BenchmarkError, "global-evaluation.*exact sequences"):
                parse_benchmark(
                    directory,
                    make_contract("global-evaluation-disjoint", "global-evaluation-disjoint"),
                )

    def test_permissive_source_can_parse_overlap_rejected_for_target(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            rows = unique_rows_by_file()
            rows[("cc", "test")].append(rows[("bp", "training")][0])
            write_nine_csvs(directory, rows_by_file=rows)
            source = parse_benchmark(directory, make_contract("allow", "allow"))
            self.assertIn("P1", source.proteins)
            with self.assertRaises(BenchmarkError):
                parse_benchmark(
                    directory,
                    make_contract("global-evaluation-disjoint", "allow"),
                )

    def test_non_binary_label_is_rejected_incrementally(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            write_nine_csvs(directory)
            path = directory / "bp-training.csv"
            with path.open() as source:
                rows = list(csv.reader(source))
            rows[1][-1] = "2"
            with path.open("w", newline="") as handle:
                csv.writer(handle).writerows(rows)
            with self.assertRaisesRegex(BenchmarkError, "non-binary GO label"):
                parse_benchmark(directory, make_contract())

    def test_duplicate_digest_is_fixed_size_for_wide_labels(self):
        digest = duplicate_row_digest("P1", "ACDE", ("1" for _ in range(10000)))
        self.assertIsInstance(digest, bytes)
        self.assertEqual(len(digest), 32)

    def test_wide_csv_streams_through_parser(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            go_terms = ["GO:%07d" % number for number in range(1, 5001)]
            for ontology in ("bp", "cc", "mf"):
                for split in ("training", "validation", "test"):
                    with (directory / ("%s-%s.csv" % (ontology, split))).open(
                        "w", newline=""
                    ) as handle:
                        writer = csv.writer(handle)
                        writer.writerow(["proteins", "sequences", *go_terms])
                        writer.writerow(["P1", "ACDE", *("0" for _ in go_terms)])
            benchmark = parse_benchmark(directory, make_contract())
            self.assertEqual(len(benchmark.proteins), 1)

    def test_fingerprint_is_order_independent_but_membership_sensitive(self):
        with tempfile.TemporaryDirectory() as first_tmp, tempfile.TemporaryDirectory() as second_tmp:
            first, second = Path(first_tmp), Path(second_tmp)
            rows = {
                (o, s): [("Z", "ACDE", "1"), ("A", "FGHI", "0")]
                for o in ("bp", "cc", "mf") for s in ("training", "validation", "test")
            }
            write_nine_csvs(first, rows_by_file=rows)
            reversed_rows = {key: list(reversed(value)) for key, value in rows.items()}
            write_nine_csvs(second, rows_by_file=reversed_rows)
            self.assertEqual(
                parse_benchmark(first, make_contract()).fingerprint,
                parse_benchmark(second, make_contract()).fingerprint,
            )
            changed = dict(reversed_rows)
            changed[("bp", "training")] = [("Z", "ACDE", "1")]
            write_nine_csvs(second, rows_by_file=changed)
            self.assertNotEqual(
                parse_benchmark(first, make_contract()).fingerprint,
                parse_benchmark(second, make_contract()).fingerprint,
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

    def test_alias_uses_separate_target_and_source_patterns(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "aliases.tsv"
            path.write_text(
                "protein_id\tsource_protein_id\tmodality\tmapping_route\tsource_identity\tmapping_evidence\n"
                "TARGET1\tSOURCE1\tprott5\tcurated\tmodel\tsequence-sha256:abc\n"
            )
            aliases = load_aliases(path, r"^TARGET[0-9]+$", r"^SOURCE[0-9]+$")
            self.assertIn(("TARGET1", "prott5"), aliases)
            with self.assertRaisesRegex(BenchmarkError, "unsafe or malformed"):
                load_aliases(path, r"^TARGET[0-9]+$", r"^UNRELATED[0-9]+$")

    def test_malformed_csv_is_translated_to_benchmark_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            write_nine_csvs(directory)
            (directory / "bp-training.csv").write_text(
                'proteins,sequences,GO:0000001\n"P1,ACDE,1\n'
            )
            with self.assertRaisesRegex(BenchmarkError, "Malformed CSV"):
                parse_benchmark(directory, make_contract())

    def test_control_character_in_id_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            rows = unique_rows_by_file()
            rows[("bp", "training")] = [("P1\x00X", "ACDE", "1")]
            write_nine_csvs(directory, rows_by_file=rows)
            permissive = make_contract()
            permissive = type(permissive)(
                id_overlap=permissive.id_overlap,
                sequence_overlap=permissive.sequence_overlap,
                protein_id_pattern=r"^.*$",
                sequence_pattern=permissive.sequence_pattern,
            )
            with self.assertRaisesRegex(BenchmarkError, "Malformed CSV|malformed protein ID"):
                parse_benchmark(directory, permissive)


if __name__ == "__main__":
    unittest.main()
