from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Optional, Sequence, Tuple

from .benchmark import parse_benchmark, validate_benchmark_name
from .models import ReusePlannerError
from .planner import build_plan
from .reports import write_reports


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="pfp_benchmark_reuse",
        description="Plan benchmark-level reuse from exact protein ID and sequence matches.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    plan = subparsers.add_parser("plan", help="partition target proteins into reuse/regenerate")
    plan.add_argument(
        "--embedded-benchmark",
        action="append",
        required=True,
        metavar="NAME=PATH",
        help="previously embedded benchmark; repeat for multiple references",
    )
    plan.add_argument(
        "--target-benchmark",
        action="append",
        required=True,
        metavar="NAME=PATH",
        help="the one target benchmark",
    )
    plan.add_argument("--output-dir", required=True, type=Path)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if len(args.target_benchmark) != 1:
            raise ReusePlannerError("Exactly one --target-benchmark is required")
        embedded_specs = [_parse_named_path(value) for value in args.embedded_benchmark]
        target_spec = _parse_named_path(args.target_benchmark[0])
        _validate_unique_names(embedded_specs, target_spec)

        embedded = tuple(
            parse_benchmark(name, path)
            for name, path in sorted(embedded_specs, key=lambda item: item[0])
        )
        target = parse_benchmark(*target_spec)
        plan = build_plan(embedded, target)
        output = write_reports(plan, args.output_dir)
    except (ReusePlannerError, OSError) as exc:
        print("error: %s" % exc, file=sys.stderr)
        return 2

    print(
        "Published %s (reuse=%d, regenerate=%d)"
        % (output, len(plan.reuse_records), len(plan.regenerate_records))
    )
    return 0


def _parse_named_path(value: str) -> Tuple[str, Path]:
    name, separator, raw_path = value.partition("=")
    if not separator or not name or not raw_path:
        raise ReusePlannerError("Benchmark arguments must use non-empty NAME=PATH: %r" % value)
    validate_benchmark_name(name)
    return name, Path(raw_path)


def _validate_unique_names(
    embedded: Sequence[Tuple[str, Path]], target: Tuple[str, Path]
) -> None:
    names = [name for name, _ in embedded] + [target[0]]
    duplicates = sorted({name for name in names if names.count(name) > 1})
    if duplicates:
        raise ReusePlannerError("Benchmark names must be unique: %s" % ", ".join(duplicates))
