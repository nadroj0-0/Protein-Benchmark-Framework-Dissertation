import csv
import gzip
import json
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd


PACKAGE_ROOT = Path(__file__).resolve().parents[1]
SRC = PACKAGE_ROOT / "src"
sys.path.insert(0, str(SRC))

from pfp_benchmark_forensics.analysis import analyze  # noqa: E402
from pfp_benchmark_forensics.config import ConfigError, load_config  # noqa: E402
from pfp_benchmark_forensics.reports import write_reports  # noqa: E402


ROOTS = {
    "bp": ("BPO", "GO:0008150", "GO:1000001", "GO:1000002", "biological_process"),
    "cc": ("CCO", "GO:0005575", "GO:2000001", "GO:2000002", "cellular_component"),
    "mf": ("MFO", "GO:0003674", "GO:3000001", "GO:3000002", "molecular_function"),
}
SPLITS = ("training", "validation", "test")


def write_obo(path: Path) -> None:
    lines = ["format-version: 1.2", ""]
    for _prefix, (_aspect, root, retained, filtered, namespace) in ROOTS.items():
        lines.extend(
            [
                "[Term]",
                f"id: {root}",
                f"name: {namespace} root",
                f"namespace: {namespace}",
                "",
                "[Term]",
                f"id: {retained}",
                "name: retained child",
                f"namespace: {namespace}",
                f"is_a: {root} ! root",
                "",
                "[Term]",
                f"id: {filtered}",
                "name: filtered child",
                f"namespace: {namespace}",
                f"is_a: {root} ! root",
                "",
            ]
        )
    path.write_text("\n".join(lines), encoding="utf-8")


def write_csvs(directory: Path, *, singular_mf: bool = False) -> dict:
    directory.mkdir(parents=True)
    sequences = {}
    for split_index, split in enumerate(SPLITS):
        proteins = [f"P{split_index}A", f"P{split_index}B", f"P{split_index}C"]
        for index, protein in enumerate(proteins):
            sequences[protein] = "M" + "A" * (10 + split_index + index)
        for prefix, (_aspect, root, retained, _filtered, _namespace) in ROOTS.items():
            first = "protein" if singular_mf and prefix == "mf" else "proteins"
            with (directory / f"{prefix}-{split}.csv").open(
                "w", encoding="utf-8", newline=""
            ) as handle:
                writer = csv.writer(handle)
                writer.writerow([first, "sequences", root, retained])
                writer.writerow([proteins[0], sequences[proteins[0]], 1, 0])
                writer.writerow([proteins[1], sequences[proteins[1]], 1, 0])
                writer.writerow([proteins[2], sequences[proteins[2]], 1, 1])
    return sequences


def write_pickles(directory: Path, sequences: dict) -> None:
    directory.mkdir(parents=True)
    files = {
        "training": "train_data_train.pkl",
        "validation": "train_data_valid.pkl",
        "test": "test_data.pkl",
    }
    roots = {item[1] for item in ROOTS.values()}
    filtered = {item[3] for item in ROOTS.values()}
    retained = {item[2] for item in ROOTS.values()}
    for split_index, split in enumerate(SPLITS):
        proteins = [f"P{split_index}A", f"P{split_index}B", f"P{split_index}C"]
        frame = pd.DataFrame(
            {
                "proteins": proteins,
                "sequences": [sequences[protein] for protein in proteins],
                "annotations": [
                    tuple(sorted(roots)),
                    tuple(sorted(roots | filtered)),
                    tuple(sorted(roots | retained)),
                ],
            }
        )
        frame.to_pickle(directory / files[split])


def write_taxonomy(path: Path, sequences: dict) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["protein_id", "taxon_id", "taxon_name"])
        for index, protein in enumerate(sorted(sequences)):
            if protein == "P2C":
                continue
            taxon_id = "9606" if index % 2 == 0 else "10090"
            name = "Homo sapiens" if taxon_id == "9606" else "Mus musculus"
            writer.writerow([protein, taxon_id, name])


def write_taxonomy_rows(path: Path, rows: list) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["protein_id", "taxon_id", "taxon_name"])
        writer.writerows(rows)


def write_uniprot_dat(path: Path, records: list) -> None:
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        for primary, secondary, taxon_id, organism, sequence in records:
            accessions = "; ".join([primary, *secondary]) + ";"
            handle.write(f"ID   {primary}_TEST\n")
            handle.write(f"AC   {accessions}\n")
            handle.write(f"OS   {organism}.\n")
            handle.write(f"OX   NCBI_TaxID={taxon_id};\n")
            if sequence is not None:
                handle.write(
                    f"SQ   SEQUENCE   {len(sequence)} AA;  0 MW;  000000000000000 R;\n"
                )
                handle.write(f"     {sequence}\n")
            handle.write("//\n")


def write_inventory(path: Path, sequences: dict) -> None:
    fields = [
        "protein_id",
        "modality",
        "exists",
        "valid",
        "scientifically_eligible",
        "requested_action",
    ]
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for protein in sorted(sequences):
            for modality in ("prott5", "text", "structure", "ppi"):
                available = modality != "ppi" or protein.endswith("A")
                writer.writerow(
                    {
                        "protein_id": protein,
                        "modality": modality,
                        "exists": str(available).lower(),
                        "valid": str(available).lower(),
                        "scientifically_eligible": str(available).lower(),
                        "requested_action": "reuse" if available else "regenerate",
                    }
                )


def write_categories(path: Path, sequences: dict) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["protein_id", "family_id", "family_name"])
        for protein in sorted(sequences):
            writer.writerow([protein, "IPR0001", "Example family"])


def build_dataset(root: Path, name: str, *, singular_mf: bool = False) -> dict:
    dataset_root = root / name
    csv_dir = dataset_root / "csv"
    pickle_dir = dataset_root / "source"
    obo = dataset_root / "go.obo"
    taxonomy = dataset_root / "taxonomy.tsv"
    inventory = dataset_root / "embedding_inventory.tsv.gz"
    categories = dataset_root / "families.tsv"
    dataset_root.mkdir()
    write_obo(obo)
    sequences = write_csvs(csv_dir, singular_mf=singular_mf)
    write_pickles(pickle_dir, sequences)
    write_taxonomy(taxonomy, sequences)
    write_inventory(inventory, sequences)
    write_categories(categories, sequences)
    return {
        "id": name,
        "benchmark_dir": str(csv_dir),
        "obo_file": str(obo),
        "allow_legacy_singular_protein_header": singular_mf,
        "source_annotations": {
            "type": "pfp-pickle-directory",
            "path": str(pickle_dir),
            "projection_policy": "min_count=50",
        },
        "taxonomy_sources": [
            {
                "type": "tsv",
                "path": str(taxonomy),
                "id_columns": ["protein_id"],
                "taxon_id_column": "taxon_id",
                "taxon_name_column": "taxon_name",
            }
        ],
        "modality_inventory": {
            "type": "embedding-inventory",
            "path": str(inventory),
        },
        "category_sources": [
            {
                "name": "interpro_family",
                "path": str(categories),
                "protein_id_column": "protein_id",
                "category_id_column": "family_id",
                "category_name_column": "family_name",
            }
        ],
    }


class ConfigTests(unittest.TestCase):
    def test_unknown_keys_fail_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "config.json"
            path.write_text(
                json.dumps(
                    {
                        "schema_version": 1,
                        "run_name": "test",
                        "datasets": [],
                        "mystery": True,
                    }
                )
            )
            with self.assertRaisesRegex(ConfigError, "unsupported keys"):
                load_config(path)

    def test_empty_dataset_list_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "config.json"
            path.write_text(
                json.dumps({"schema_version": 1, "run_name": "test", "datasets": []})
            )
            with self.assertRaisesRegex(ConfigError, "non-empty"):
                load_config(path)

    def test_taxonomy_source_names_must_be_unique(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            dataset = build_dataset(root, "dataset")
            source = dataset["taxonomy_sources"][0]
            source["name"] = "duplicate"
            dataset["taxonomy_sources"].append(dict(source))
            path = root / "config.json"
            path.write_text(
                json.dumps(
                    {"schema_version": 1, "run_name": "test", "datasets": [dataset]}
                )
            )
            with self.assertRaisesRegex(ConfigError, "names must be unique"):
                load_config(path)


class AnalysisTests(unittest.TestCase):
    def test_complete_analysis_and_atomic_reports(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            config_path = root / "config.json"
            payload = {
                "schema_version": 1,
                "run_name": "comparison",
                "top_n": 2,
                "datasets": [
                    build_dataset(root, "published", singular_mf=True),
                    build_dataset(root, "contemporary"),
                ],
            }
            config_path.write_text(json.dumps(payload), encoding="utf-8")
            config = load_config(config_path)
            bundle = analyze(config)

            root_rows = [
                row
                for row in bundle.tables["root_only_summary"]
                if row["dataset_id"] == "contemporary"
                and row["aspect"] == "BPO"
                and row["split"] == "test"
            ]
            counts = {row["classification"]: row["proteins"] for row in root_rows}
            self.assertEqual(counts["source_root_only"], 1)
            self.assertEqual(counts["projection_created"], 1)
            self.assertEqual(counts["source_unresolved"], 0)

            taxon_rows = [
                row
                for row in bundle.tables["taxonomy_distribution"]
                if row["dataset_id"] == "contemporary"
                and row["aspect"] == "MFO"
                and row["split"] == "test"
            ]
            self.assertEqual(sum(row["proteins"] for row in taxon_rows), 3)
            self.assertIn("__UNMAPPED__", {row["taxon_id"] for row in taxon_rows})

            ppi = [
                row
                for row in bundle.tables["modality_coverage"]
                if row["dataset_id"] == "contemporary"
                and row["aspect"] == "CCO"
                and row["split"] == "test"
                and row["modality"] == "ppi"
                and row["coverage_state"] == "artifact_valid"
            ]
            self.assertEqual(len(ppi), 1)
            self.assertAlmostEqual(ppi[0]["coverage_fraction"], 1 / 3)
            self.assertTrue(bundle.tables["cross_benchmark_metrics"])
            self.assertTrue(bundle.tables["cross_benchmark_modality"])
            self.assertTrue(bundle.tables["cross_benchmark_taxonomy"])

            output = root / "reports"
            write_reports(output, bundle, config, replace=False)
            self.assertTrue((output / "RUN_COMPLETE.json").is_file())
            self.assertTrue((output / "input_manifest.json").is_file())
            self.assertTrue((output / "output_manifest.json").is_file())
            self.assertTrue((output / "benchmark_forensics.md").is_file())
            self.assertTrue((output / "root_only_provenance.tsv").is_file())
            self.assertTrue((output / "taxonomy_distribution.tsv").is_file())
            self.assertTrue((output / "taxonomy_conflicts.tsv").is_file())
            self.assertTrue((output / "modality_coverage.tsv").is_file())
            manifest = json.loads((output / "output_manifest.json").read_text())
            self.assertTrue(manifest["outputs"])
            self.assertTrue(
                all("/" not in item["path"] for item in manifest["outputs"])
            )
            with self.assertRaises(FileExistsError):
                write_reports(output, bundle, config, replace=False)
            write_reports(output, bundle, config, replace=True)
            self.assertTrue((output / "RUN_COMPLETE.json").is_file())

    def test_source_sequence_mismatch_fails_loudly(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            dataset = build_dataset(root, "dataset")
            frame = pd.read_pickle(
                Path(dataset["source_annotations"]["path"]) / "test_data.pkl"
            )
            frame.at[0, "sequences"] = "MISMATCH"
            frame.to_pickle(
                Path(dataset["source_annotations"]["path"]) / "test_data.pkl"
            )
            config_path = root / "config.json"
            config_path.write_text(
                json.dumps(
                    {"schema_version": 1, "run_name": "test", "datasets": [dataset]}
                )
            )
            with self.assertRaisesRegex(ValueError, "source sequence differs"):
                analyze(load_config(config_path))

    def test_uniprot_dat_taxonomy_adapter(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            dataset = build_dataset(root, "dataset")
            dat = root / "taxonomy.dat.gz"
            with gzip.open(dat, "wt", encoding="utf-8") as handle:
                handle.write(
                    "ID   TEST\nAC   P0A;\nOS   Homo sapiens.\nOX   NCBI_TaxID=9606;\n//\n"
                )
            dataset["taxonomy_sources"] = [{"type": "uniprot-dat", "path": str(dat)}]
            config_path = root / "config.json"
            config_path.write_text(
                json.dumps(
                    {"schema_version": 1, "run_name": "test", "datasets": [dataset]}
                )
            )
            bundle = analyze(load_config(config_path))
            rows = [
                row
                for row in bundle.tables["taxonomy_distribution"]
                if row["aspect"] == "BPO" and row["split"] == "training"
            ]
            human = [row for row in rows if row["taxon_id"] == "9606"]
            self.assertEqual(human[0]["proteins"], 1)

    def test_uniprot_alias_conflict_uses_unique_exact_sequence_match(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            dataset = build_dataset(root, "dataset")
            sequence = pd.read_csv(
                Path(dataset["benchmark_dir"]) / "bp-training.csv"
            ).iloc[0]["sequences"]
            dat = root / "taxonomy.dat.gz"
            records = [
                ("QOLD1", ["P0A"], "11111", "Wrong organism", "MBBBBBBBBBB"),
                ("QOLD2", ["P0A"], "22222", "Matching organism", sequence),
            ]
            write_uniprot_dat(dat, records)
            dataset["taxonomy_sources"] = [
                {
                    "type": "uniprot-dat",
                    "path": str(dat),
                    "name": "historical-uniprot",
                    "priority": 100,
                }
            ]
            config_path = root / "config.json"
            config_path.write_text(
                json.dumps(
                    {"schema_version": 1, "run_name": "test", "datasets": [dataset]}
                )
            )

            bundle = analyze(load_config(config_path))
            membership = {
                row["protein_id"]: row for row in bundle.tables["protein_membership"]
            }
            self.assertEqual(membership["P0A"]["taxon_id"], "22222")
            self.assertEqual(membership["P0A"]["taxonomy_accession_role"], "secondary")
            self.assertTrue(membership["P0A"]["taxonomy_sequence_matches_benchmark"])
            self.assertTrue(membership["P0A"]["taxonomy_conflict_resolved"])
            conflicts = bundle.tables["taxonomy_conflicts"]
            self.assertEqual(len(conflicts), 1)
            self.assertEqual(conflicts[0]["status"], "resolved")
            self.assertEqual(conflicts[0]["alternative_taxon_id"], "11111")

            write_uniprot_dat(dat, list(reversed(records)))
            reversed_bundle = analyze(load_config(config_path))
            reversed_membership = {
                row["protein_id"]: row
                for row in reversed_bundle.tables["protein_membership"]
            }
            self.assertEqual(reversed_membership["P0A"]["taxon_id"], "22222")

    def test_uniprot_alias_ambiguity_is_logged_and_left_unmapped(self):
        cases = {
            "zero-match": [
                ("QOLD1", ["P0A"], "11111", "First organism", "MBBBBBBBBBB"),
                ("QOLD2", ["P0A"], "22222", "Second organism", "MCCCCCCCCCC"),
            ],
            "multiple-match": None,
        }
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            dataset = build_dataset(root, "dataset")
            sequence = pd.read_csv(
                Path(dataset["benchmark_dir"]) / "bp-training.csv"
            ).iloc[0]["sequences"]
            cases["multiple-match"] = [
                ("QOLD1", ["P0A"], "11111", "First organism", sequence),
                ("QOLD2", ["P0A"], "22222", "Second organism", sequence),
            ]
            dat = root / "taxonomy.dat.gz"
            dataset["taxonomy_sources"] = [
                {
                    "type": "uniprot-dat",
                    "path": str(dat),
                    "name": "historical-uniprot",
                    "priority": 100,
                }
            ]
            config_path = root / "config.json"
            for case, records in cases.items():
                with self.subTest(case=case):
                    write_uniprot_dat(dat, records)
                    config_path.write_text(
                        json.dumps(
                            {
                                "schema_version": 1,
                                "run_name": "test",
                                "datasets": [dataset],
                            }
                        )
                    )
                    bundle = analyze(load_config(config_path))
                    membership = {
                        row["protein_id"]: row
                        for row in bundle.tables["protein_membership"]
                    }
                    self.assertFalse(membership["P0A"]["taxonomy_mapped"])
                    self.assertTrue(membership["P0A"]["taxonomy_conflict_unresolved"])
                    conflicts = [
                        row
                        for row in bundle.tables["taxonomy_conflicts"]
                        if row["protein_id"] == "P0A"
                    ]
                    self.assertEqual(len(conflicts), 2)
                    self.assertEqual(
                        {row["status"] for row in conflicts}, {"unresolved"}
                    )
                    self.assertEqual(
                        bundle.summary["datasets"]["dataset"][
                            "taxonomy_unresolved_conflict_proteins"
                        ],
                        1,
                    )

    def test_higher_priority_taxonomy_resolves_and_reports_snapshot_change(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            dataset = build_dataset(root, "dataset")
            old_source = dataset["taxonomy_sources"][0]
            old_source.update({"name": "uniprot-old", "priority": 100})
            new_taxonomy = root / "taxonomy-new.tsv"
            write_taxonomy_rows(
                new_taxonomy,
                [["P0A", "11111", "Updated organism assignment"]],
            )
            new_source = {
                "type": "tsv",
                "path": str(new_taxonomy),
                "name": "uniprot-new",
                "priority": 200,
            }
            dataset["taxonomy_sources"] = [old_source, new_source]
            config_path = root / "config.json"
            config_path.write_text(
                json.dumps(
                    {"schema_version": 1, "run_name": "test", "datasets": [dataset]}
                )
            )

            bundle = analyze(load_config(config_path))
            conflicts = bundle.tables["taxonomy_conflicts"]
            self.assertEqual(len(conflicts), 1)
            self.assertEqual(conflicts[0]["protein_id"], "P0A")
            self.assertEqual(conflicts[0]["selected_taxon_id"], "11111")
            self.assertEqual(conflicts[0]["selected_source_name"], "uniprot-new")
            self.assertEqual(conflicts[0]["alternative_taxon_id"], "9606")
            self.assertEqual(
                bundle.summary["datasets"]["dataset"]["taxonomy_conflict_proteins"],
                1,
            )
            membership = {
                row["protein_id"]: row for row in bundle.tables["protein_membership"]
            }
            self.assertEqual(membership["P0A"]["taxon_id"], "11111")
            self.assertTrue(membership["P0A"]["taxonomy_conflict_resolved"])

            dataset["taxonomy_sources"] = [new_source, old_source]
            config_path.write_text(
                json.dumps(
                    {"schema_version": 1, "run_name": "test", "datasets": [dataset]}
                )
            )
            reversed_bundle = analyze(load_config(config_path))
            self.assertEqual(
                reversed_bundle.tables["taxonomy_conflicts"],
                conflicts,
            )

    def test_equal_priority_taxonomy_conflict_is_logged_and_left_unmapped(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            dataset = build_dataset(root, "dataset")
            first = dataset["taxonomy_sources"][0]
            first.update({"name": "source-a", "priority": 100})
            second_path = root / "taxonomy-second.tsv"
            write_taxonomy_rows(second_path, [["P0A", "11111", "Other organism"]])
            dataset["taxonomy_sources"].append(
                {
                    "type": "tsv",
                    "path": str(second_path),
                    "name": "source-b",
                    "priority": 100,
                }
            )
            config_path = root / "config.json"
            config_path.write_text(
                json.dumps(
                    {"schema_version": 1, "run_name": "test", "datasets": [dataset]}
                )
            )
            bundle = analyze(load_config(config_path))
            membership = {
                row["protein_id"]: row for row in bundle.tables["protein_membership"]
            }
            self.assertFalse(membership["P0A"]["taxonomy_mapped"])
            self.assertTrue(membership["P0A"]["taxonomy_conflict_unresolved"])
            conflicts = [
                row
                for row in bundle.tables["taxonomy_conflicts"]
                if row["protein_id"] == "P0A"
            ]
            self.assertEqual(len(conflicts), 2)
            self.assertEqual({row["status"] for row in conflicts}, {"unresolved"})

    def test_taxonomy_source_priority_must_be_integer(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            dataset = build_dataset(root, "dataset")
            dataset["taxonomy_sources"][0]["priority"] = "latest"
            config_path = root / "config.json"
            config_path.write_text(
                json.dumps(
                    {"schema_version": 1, "run_name": "test", "datasets": [dataset]}
                )
            )
            with self.assertRaisesRegex(ConfigError, "priority must be an integer"):
                load_config(config_path)


if __name__ == "__main__":
    unittest.main()
