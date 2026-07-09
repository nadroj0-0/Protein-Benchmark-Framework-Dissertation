from __future__ import annotations

import argparse
from importlib import resources
from pathlib import Path

from .builder import (
    build_benchmark,
    export_from_deepgoplus_pickles,
    generate_deepgoplus_pickles_from_cafa_files,
)
from .config import BuildConfig, EVIDENCE_POLICIES, normalise_taxa


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


def read_packaged_cafa3_taxa() -> list[str]:
    path = resources.files("cafa_benchmark_builder").joinpath("resources/cafa3_target_taxa.txt")
    values = []
    with path.open("r") as handle:
        for line in handle:
            line = line.strip()
            if line and not line.startswith("#"):
                values.append(line)
    return values


def require_args(args: argparse.Namespace, names: list[str]) -> None:
    missing = [name for name in names if getattr(args, name) in (None, [])]
    if missing:
        formatted = ", ".join("--" + name.replace("_", "-") for name in missing)
        raise SystemExit(f"Missing required arguments for --source-mode {args.source_mode}: {formatted}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build 2025->2026 CAFA-style PFP-compatible benchmark CSVs."
    )
    parser.add_argument("--source-mode", choices=("snapshots", "deepgoplus", "cafa3-files"), default="snapshots",
                        help=("Input mode. snapshots parses UniProt/GOA/GO; deepgoplus exports from released "
                              "pickles; cafa3-files regenerates train_data.pkl/test_data.pkl/terms.pkl from "
                              "official CAFA3/DeepGOPlus files."))
    parser.add_argument("--deepgoplus-dir", type=Path,
                        help="Directory containing train_data_train.pkl, train_data_valid.pkl, test_data.pkl, terms.pkl.")
    parser.add_argument("--train-sequences-file", type=Path,
                        help="CAFA/DeepGOPlus training FASTA for --source-mode cafa3-files.")
    parser.add_argument("--train-annotations-file", type=Path,
                        help="CAFA/DeepGOPlus training annotation TSV for --source-mode cafa3-files.")
    parser.add_argument("--test-sequences-file", type=Path,
                        help="CAFA/DeepGOPlus target FASTA for --source-mode cafa3-files.")
    parser.add_argument("--test-annotations-file", type=Path,
                        help="CAFA/DeepGOPlus test/ground-truth annotation TSV for --source-mode cafa3-files.")
    parser.add_argument("--uniprot-t0", action="append", type=Path,
                        help="UniProt t0 FASTA or DAT file. Repeat for multiple files.")
    parser.add_argument("--uniprot-t1", action="append", type=Path,
                        help="UniProt t1 FASTA or DAT file. Repeat for multiple files.")
    parser.add_argument("--goa-t0", type=Path, help="GOA t0 GAF/GAF.gz file.")
    parser.add_argument("--goa-t1", type=Path, help="GOA t1 GAF/GAF.gz file.")
    parser.add_argument("--go-obo", type=Path, required=True, help="GO ontology OBO file.")
    parser.add_argument("--output-dir", type=Path, required=True, help="Output directory.")
    parser.add_argument("--taxon-policy", choices=("all", "cafa3-targets", "custom"), default="all",
                        help="Taxon scope for snapshot mode. all = official broad training; cafa3-targets = CAFA3 target taxa.")
    parser.add_argument("--target-taxon", action="append", default=[],
                        help="NCBI taxon ID to include. Repeatable. Default: all taxa.")
    parser.add_argument("--target-taxa-file", type=Path,
                        help="Optional file containing one taxon ID per line.")
    parser.add_argument("--evidence-policy", choices=sorted(EVIDENCE_POLICIES), default="cafa3-final",
                        help="Named evidence-code policy. Default: cafa3-final.")
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
    require_args(args, ["uniprot_t0", "uniprot_t1", "goa_t0", "goa_t1"])
    taxa_values = []
    if args.taxon_policy == "cafa3-targets":
        taxa_values.extend(read_packaged_cafa3_taxa())
    taxa_values.extend(args.target_taxon)
    taxa_values.extend(read_taxa_file(args.target_taxa_file))
    if args.taxon_policy == "custom" and not taxa_values:
        raise SystemExit("--taxon-policy custom requires --target-taxon or --target-taxa-file")
    evidence = frozenset(args.evidence_code) if args.evidence_code else EVIDENCE_POLICIES[args.evidence_policy]
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
    if args.source_mode == "deepgoplus":
        require_args(args, ["deepgoplus_dir"])
        written = export_from_deepgoplus_pickles(
            deepgoplus_dir=args.deepgoplus_dir,
            go_obo=args.go_obo,
            output_dir=args.output_dir,
            include_rels=not args.no_rels,
            write_intermediates=not args.no_intermediates,
        )
    elif args.source_mode == "cafa3-files":
        require_args(args, [
            "train_sequences_file",
            "train_annotations_file",
            "test_sequences_file",
            "test_annotations_file",
        ])
        written = generate_deepgoplus_pickles_from_cafa_files(
            go_obo=args.go_obo,
            train_sequences_file=args.train_sequences_file,
            train_annotations_file=args.train_annotations_file,
            test_sequences_file=args.test_sequences_file,
            test_annotations_file=args.test_annotations_file,
            output_dir=args.output_dir,
            min_count=args.min_count,
            include_rels=not args.no_rels,
        )
    else:
        written = build_benchmark(config_from_args(args))
    print("Wrote:")
    for key in sorted(written):
        print(f"  {key}: {written[key]}")
