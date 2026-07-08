from __future__ import annotations

import argparse
from pathlib import Path

from .builder import build_benchmark
from .config import BuildConfig, CAFA3_FINAL_EXP_CODES, normalise_taxa


def read_taxa_file(path: Path | None) -> list[str]:
    if path is None:
        return []
    values = []
    with open(path, "r") as handle:
        for line in handle:
            line = line.strip()
            if line and not line.startswith("#"):
                values.append(line)
    return values


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build 2025->2026 CAFA-style PFP-compatible benchmark CSVs."
    )
    parser.add_argument("--uniprot-t0", action="append", required=True, type=Path,
                        help="UniProt t0 FASTA or DAT file. Repeat for multiple files.")
    parser.add_argument("--uniprot-t1", action="append", required=True, type=Path,
                        help="UniProt t1 FASTA or DAT file. Repeat for multiple files.")
    parser.add_argument("--goa-t0", required=True, type=Path, help="GOA t0 GAF/GAF.gz file.")
    parser.add_argument("--goa-t1", required=True, type=Path, help="GOA t1 GAF/GAF.gz file.")
    parser.add_argument("--go-obo", required=True, type=Path, help="GO ontology OBO file.")
    parser.add_argument("--output-dir", required=True, type=Path, help="Output directory.")
    parser.add_argument("--target-taxon", action="append", default=[],
                        help="NCBI taxon ID to include. Repeatable. Default: all taxa.")
    parser.add_argument("--target-taxa-file", type=Path,
                        help="Optional file containing one taxon ID per line.")
    parser.add_argument("--evidence-code", action="append", default=[],
                        help="Override evidence code set. Repeatable. Default: final CAFA3 policy.")
    parser.add_argument("--min-count", type=int, default=50,
                        help="DeepGOPlus term frequency threshold. Default: 50.")
    parser.add_argument("--split", type=float, default=0.9,
                        help="DeepGOPlus train/valid split. Default: 0.9.")
    parser.add_argument("--seed", type=int, default=0,
                        help="DeepGOPlus split seed. Default: 0.")
    parser.add_argument("--reviewed-only", action="store_true",
                        help="Keep only Swiss-Prot/reviewed records when identifiable.")
    parser.add_argument("--no-rels", action="store_true",
                        help="Do not include OBO relationship parents. Default follows DeepGOPlus with_rels=True.")
    parser.add_argument("--no-intermediates", action="store_true",
                        help="Do not write DeepGOPlus-style pickle intermediates.")
    parser.add_argument("--max-gaf-records", type=int,
                        help="Smoke-test limiter for GAF records. Do not use for full builds.")
    return parser


def config_from_args(args: argparse.Namespace) -> BuildConfig:
    taxa_values = list(args.target_taxon) + read_taxa_file(args.target_taxa_file)
    evidence = frozenset(args.evidence_code) if args.evidence_code else CAFA3_FINAL_EXP_CODES
    return BuildConfig(
        uniprot_t0=tuple(args.uniprot_t0),
        uniprot_t1=tuple(args.uniprot_t1),
        goa_t0=args.goa_t0,
        goa_t1=args.goa_t1,
        go_obo=args.go_obo,
        output_dir=args.output_dir,
        target_taxa=normalise_taxa(taxa_values),
        evidence_codes=evidence,
        min_count=args.min_count,
        split=args.split,
        seed=args.seed,
        reviewed_only=args.reviewed_only,
        include_rels=not args.no_rels,
        write_intermediates=not args.no_intermediates,
        max_gaf_records=args.max_gaf_records,
    )


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    written = build_benchmark(config_from_args(args))
    print("Wrote:")
    for key in sorted(written):
        print(f"  {key}: {written[key]}")
