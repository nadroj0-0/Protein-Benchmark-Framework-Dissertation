"""Command-line entry point for benchmark forensic analysis."""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

from .analysis import analyze
from .config import ConfigError, load_config
from .reports import write_reports


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Characterize PFP benchmark labels, root-only provenance, taxa, "
            "optional protein categories, and modality coverage."
        )
    )
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--replace",
        action="store_true",
        help="Atomically replace an existing output directory after a successful run.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    started = time.perf_counter()
    try:
        config = load_config(args.config)
        bundle = analyze(config)
        write_reports(args.output_dir, bundle, config, replace=args.replace)
    except (
        ConfigError,
        FileNotFoundError,
        FileExistsError,
        ValueError,
        RuntimeError,
    ) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    elapsed = time.perf_counter() - started
    print(
        json.dumps(
            {
                "status": "complete",
                "run_name": config.run_name,
                "datasets": list(bundle.summary["dataset_order"]),
                "output_dir": str(args.output_dir.expanduser().resolve()),
                "elapsed_seconds": elapsed,
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
