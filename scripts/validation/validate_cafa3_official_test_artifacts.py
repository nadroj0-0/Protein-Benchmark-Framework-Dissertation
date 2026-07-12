#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


EXPECTED_CSV_ROWS = {"bp-test.csv": 2392, "cc-test.csv": 1265, "mf-test.csv": 1137}
EXPECTED_TEST_ROWS = 3328


def annotation_map(frame: pd.DataFrame) -> dict[str, frozenset[str]]:
    return {
        str(row.proteins): frozenset(row.annotations)
        for row in frame.itertuples(index=False)
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Hard-gate released CAFA3 test reconstruction against DeepGOPlus/PFP artifacts."
    )
    parser.add_argument("--generated-dir", type=Path, required=True)
    parser.add_argument("--reference-pickle-dir", type=Path, required=True)
    parser.add_argument("--reference-csv-dir", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()

    generated = pd.read_pickle(args.generated_dir / "test_data.pkl")
    reference = pd.read_pickle(args.reference_pickle_dir / "test_data.pkl")
    checks: list[tuple[str, bool, str]] = []

    generated_ids = generated["proteins"].astype(str).tolist()
    reference_ids = reference["proteins"].astype(str).tolist()
    checks.append(("test_data row count", len(generated) == EXPECTED_TEST_ROWS,
                   f"generated={len(generated)} expected={EXPECTED_TEST_ROWS}"))
    checks.append(("test_data unique IDs", len(set(generated_ids)) == len(generated_ids),
                   f"unique={len(set(generated_ids))} rows={len(generated_ids)}"))
    checks.append(("test_data ID set", set(generated_ids) == set(reference_ids),
                   f"generated_only={len(set(generated_ids) - set(reference_ids))} "
                   f"reference_only={len(set(reference_ids) - set(generated_ids))}"))
    checks.append(("test_data row order", generated_ids == reference_ids,
                   f"ordered_rows={len(generated_ids)}"))

    generated_sequences = dict(zip(generated_ids, generated["sequences"], strict=True))
    reference_sequences = dict(zip(reference_ids, reference["sequences"], strict=True))
    shared_ids = set(generated_sequences) & set(reference_sequences)
    sequence_matches = sum(
        generated_sequences[protein_id] == reference_sequences[protein_id]
        for protein_id in shared_ids
    )
    checks.append(("test_data sequences", sequence_matches == len(reference_sequences),
                   f"matching={sequence_matches} reference={len(reference_sequences)}"))

    generated_annotations = annotation_map(generated)
    reference_annotations = annotation_map(reference)
    annotation_matches = sum(
        generated_annotations.get(protein_id) == terms
        for protein_id, terms in reference_annotations.items()
    )
    checks.append(("test_data annotations", annotation_matches == len(reference_annotations),
                   f"matching={annotation_matches} reference={len(reference_annotations)}"))

    for filename, expected_rows in EXPECTED_CSV_ROWS.items():
        generated_csv = pd.read_csv(args.generated_dir / filename)
        reference_csv = pd.read_csv(args.reference_csv_dir / filename)
        generated_csv_ids = set(generated_csv["proteins"].astype(str))
        reference_csv_ids = set(reference_csv["proteins"].astype(str))
        checks.append((f"{filename} row count", len(generated_csv) == expected_rows,
                       f"generated={len(generated_csv)} expected={expected_rows}"))
        checks.append((f"{filename} ID set", generated_csv_ids == reference_csv_ids,
                       f"generated_only={len(generated_csv_ids - reference_csv_ids)} "
                       f"reference_only={len(reference_csv_ids - generated_csv_ids)}"))

    passed = all(ok for _, ok, _ in checks)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# CAFA3 Official Test Artifact Gate",
        "",
        f"Final status: **{'PASS' if passed else 'FAIL'}**",
        "",
        "This gate validates the released `leafonly_all.txt` and target-FASTA path consumed by "
        "DeepGOPlus. It does not claim reconstruction of the unavailable 15-Nov-2017 GOA snapshot.",
        "",
        "|check|status|detail|",
        "|---|---|---|",
    ]
    lines.extend(
        f"|{name}|{'PASS' if ok else 'FAIL'}|{detail}|"
        for name, ok, detail in checks
    )
    args.report.write_text("\n".join(lines) + "\n")
    print(args.report.read_text(), end="")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
