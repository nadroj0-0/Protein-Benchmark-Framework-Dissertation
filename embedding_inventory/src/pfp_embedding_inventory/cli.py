import argparse
import sys
from pathlib import Path
from typing import Optional, Sequence

from .benchmark import BenchmarkError, load_aliases, parse_benchmark
from .config import ConfigError, load_config
from .inventory import InventoryError, build_inventory
from .provenance import (
    HashCache,
    ProvenanceError,
    build_run_provenance,
    compute_cache_catalog,
    verify_artifact_scope,
)
from .reports import ReportError, write_reports


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
    parser.add_argument(
        "--artifact-root",
        type=Path,
        default=None,
        help=(
            "Root containing published archives and the pinned PFP reference; "
            "defaults to two parents above --embedding-cache"
        ),
    )
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
    parser.add_argument(
        "--report-level",
        choices=("compact", "full"),
        default="compact",
        help="compact is storage-safe; full adds gzip-compressed row-level inventory",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        config = load_config(args.config)
        print("[1/7] Parsing target benchmark under its leakage contract", flush=True)
        benchmark = parse_benchmark(
            args.benchmark_dir, config.target_benchmark_contract
        )
        print("[2/7] Parsing source benchmark under its identity contract", flush=True)
        source_benchmark = parse_benchmark(
            args.source_benchmark_dir, config.source_benchmark_contract
        )
        aliases = (
            load_aliases(
                args.aliases,
                config.target_benchmark_contract.protein_id_pattern,
                config.source_benchmark_contract.protein_id_pattern,
            )
            if args.aliases
            else {}
        )
        hash_cache = HashCache()
        print("[3/7] Fingerprinting the configured cache catalog", flush=True)
        catalog = compute_cache_catalog(args.embedding_cache, config, hash_cache)
        artifact_root = (
            args.artifact_root
            if args.artifact_root is not None
            else args.embedding_cache.resolve().parent.parent
        )
        print("[4/7] Verifying published artifact scope", flush=True)
        verification = verify_artifact_scope(
            config,
            benchmark,
            source_benchmark,
            catalog,
            args.embedding_cache,
            artifact_root,
            hash_cache,
        )
        print("[5/7] Validating arrays and planning modality-specific actions", flush=True)
        result = build_inventory(
            benchmark=benchmark,
            source_benchmark=source_benchmark,
            embedding_cache=args.embedding_cache,
            config=config,
            policy=args.policy,
            aliases=aliases,
            artifact_verification=verification,
        )
        print("[6/7] Capturing immutable run provenance", flush=True)
        command = list(sys.argv) if argv is None else ["inventory_embeddings.py", *argv]
        repository = Path(__file__).resolve().parents[3]
        provenance = build_run_provenance(
            command=command,
            repository=repository,
            config_path=args.config,
            alias_path=args.aliases,
            target=benchmark,
            source=source_benchmark,
            embedding_cache=args.embedding_cache,
            artifact_root=artifact_root,
            catalog=catalog,
            verification=verification,
            policy=args.policy,
            report_level=args.report_level,
            runtime_options={
                "benchmark_dir": str(args.benchmark_dir.resolve()),
                "source_benchmark_dir": str(args.source_benchmark_dir.resolve()),
                "embedding_cache": str(args.embedding_cache.resolve()),
                "artifact_root": str(artifact_root.resolve()),
                "output_dir": str(args.output_dir.resolve()),
                "aliases": str(args.aliases.resolve()) if args.aliases else None,
                "cache_catalog_reverified_at_report_publication": True,
            },
            hash_cache=hash_cache,
        )
        print("[7/7] Writing storage-safe reports and manifests", flush=True)
        summary = write_reports(
            result,
            args.output_dir,
            args.embedding_cache,
            report_level=args.report_level,
            provenance=provenance,
            protected_roots=(artifact_root,),
        )
    except (
        BenchmarkError,
        ConfigError,
        InventoryError,
        ProvenanceError,
        ReportError,
        UnicodeError,
        OSError,
    ) as exc:
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
    print(
        "Published cache authenticated: %s"
        % str(summary["artifact_verification"]["verified"]).lower()
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
