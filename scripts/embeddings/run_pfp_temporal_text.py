#!/usr/bin/env python3
"""Run PFP's temporal text recipe with a configurable historical cutoff."""

from __future__ import annotations

import argparse
import importlib.util
import json
import shutil
import sys
from pathlib import Path


def load_module(path: Path):
    specification = importlib.util.spec_from_file_location("pfp_extract_uniprot_text", path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"Cannot import PFP text extractor: {path}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pfp-root", type=Path, required=True)
    parser.add_argument("--cafa-assessment-dir", type=Path, required=True)
    parser.add_argument("--cutoff-date", required=True, help="Historical cutoff as YYYY-MM-DD")
    parser.add_argument("--workers", type=int, default=5)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    pfp_root = args.pfp_root.resolve()
    data_dir = pfp_root / "data"
    text_dir = data_dir / "embedding_cache" / "uniprot_text"
    temporal_dir = text_dir / "temporal_recipe"
    current = text_dir / "protein_descriptions.tsv"
    historical = text_dir / "protein_descriptions_historical.tsv"
    current_checkpoint = text_dir / "processed_checkpoint.txt"
    historical_checkpoint = text_dir / "historical_checkpoint.txt"
    historical_raw = text_dir / "historical_raw"
    punct = temporal_dir / "protein_descriptions_historical_punct_v1_test.tsv"
    mixed = temporal_dir / "protein_descriptions_mixed.tsv"

    script = pfp_root / "scripts" / "extract_uniprot_text.py"
    if not script.is_file():
        raise SystemExit(f"Missing PFP text extractor: {script}")
    if not args.cafa_assessment_dir.is_dir():
        raise SystemExit(f"Missing CAFA assessment directory: {args.cafa_assessment_dir}")

    module = load_module(script)
    module.CUTOFF_DATE = args.cutoff_date

    current_status = module.run_current_extraction(
        data_dir=data_dir,
        cafa_assessment_dir=args.cafa_assessment_dir,
        output_file=current,
        checkpoint_file=current_checkpoint,
    )
    historical_success, historical_failed, historical_status = module.extract_historical_text(
        data_dir=data_dir,
        cafa_assessment_dir=args.cafa_assessment_dir,
        output_file=historical,
        checkpoint_file=historical_checkpoint,
        raw_dir=historical_raw,
        splits=["test"],
        workers=args.workers,
    )
    historical_created_empty = False
    if not historical.exists():
        historical.touch()
        historical_created_empty = True
    punct_metadata = module.build_historical_punct_v1_test_tsv(
        historical_tsv=historical,
        output_tsv=punct,
        data_dir=data_dir,
    )
    mixed_metadata = module.build_mixed_temporal_tsv(
        current_tsv=current,
        hist_test_tsv=punct,
        output_tsv=mixed,
        bundle_dir=temporal_dir,
        historical_tsv=historical,
        data_dir=data_dir,
    )

    if not mixed.is_file():
        raise SystemExit(f"PFP temporal recipe did not create: {mixed}")
    current_backup = temporal_dir / "protein_descriptions_current_before_mixed.tsv"
    shutil.copyfile(current, current_backup)
    shutil.copyfile(mixed, current)

    report = {
        "schema_version": 1,
        "pfp_text_script": str(script),
        "historical_cutoff": args.cutoff_date,
        "current_status": current_status,
        "historical_success": historical_success,
        "historical_failed": historical_failed,
        "historical_status_counts": historical_status,
        "historical_created_empty": historical_created_empty,
        "punctuation_recipe": punct_metadata,
        "mixed_recipe": mixed_metadata,
        "embedding_input": str(current),
        "embedding_input_source": str(mixed),
    }
    report_path = temporal_dir / "framework_temporal_text_run.json"
    report_path.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
