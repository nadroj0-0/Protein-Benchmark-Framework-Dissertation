import argparse
import sys
from pathlib import Path
from typing import Optional, Sequence

from .benchmark import BenchmarkError, load_aliases, parse_benchmark
from .config import ConfigError, load_config
from .inventory import InventoryError, build_inventory
from .reports import write_reports


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Inventory and plan reuse for PFP-compatible benchmark embeddings."
    )
    parser.add_argument("--benchmark-dir", type=Path, required=True)
    parser.add_argument(
        "--source-benchmark-dir",
        type=Path,
        required=True,
        help="CSV benchmark that defines the protein identities represented by the cache",
    )
    parser.add_argument("--embedding-cache", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--policy",
        choices=("paper-faithful", "maximize-coverage"),
        default="paper-faithful",
    )
    parser.add_argument(
        "--aliases",
        type=Path,
        default=None,
        help="Optional explicit tab-separated alias mapping; fuzzy matching is never performed",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        config = load_config(args.config)
        print("[1/4] Parsing benchmark CSVs", flush=True)
        benchmark = parse_benchmark(args.benchmark_dir, config.benchmark_contract)
        print("[2/4] Parsing source benchmark CSVs", flush=True)
        source_benchmark = parse_benchmark(args.source_benchmark_dir, config.benchmark_contract)
        aliases = (
            load_aliases(args.aliases, config.benchmark_contract.protein_id_pattern)
            if args.aliases
            else {}
        )
        print("[3/4] Validating arrays and planning modality-specific actions", flush=True)
        result = build_inventory(
            benchmark=benchmark,
            source_benchmark=source_benchmark,
            embedding_cache=args.embedding_cache,
            config=config,
            policy=args.policy,
            aliases=aliases,
        )
        print("[4/4] Writing inventory reports and manifests", flush=True)
        summary = write_reports(result, args.output_dir, args.embedding_cache)
    except (BenchmarkError, ConfigError, InventoryError, OSError) as exc:
        print("ERROR: %s" % exc, file=sys.stderr)
        return 2

    global_coverage = summary["coverage"]["global"]
    print(
        "Wrote %d proteins x %d modalities to %s"
        % (summary["population"], 4, args.output_dir.resolve())
    )
    print(
        "Complete physical coverage: %d; at least one: %d"
        % (
            global_coverage["complete_four_modalities_present"]["count"],
            global_coverage["at_least_one_modality"]["count"],
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
