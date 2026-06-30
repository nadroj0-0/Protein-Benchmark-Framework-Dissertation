#!/usr/bin/env python3
"""
verify_embeddings.py — completeness + correctness gate for generated embeddings.

Verifies that, for every protein ID in the configured dataset splits, a
corresponding {pid}.npy embedding exists for each configured modality. Beyond
completeness, it validates embedding DIMENSION and FINITENESS, so that partial
or corrupted embedding generation is detected BEFORE training rather than
silently propagating into the model.

The verification logic is benchmark-agnostic. All benchmark-specific values
(aspect/split names, cache directory, expected dimensions, coverage thresholds,
sampling parameters) live in a JSON config file. The config fully specifies the
run; adapting to a new benchmark means writing a new config, not editing code.

Usage:
    python verify_embeddings.py --data-dir data
    python verify_embeddings.py --data-dir data --config configs/cafa3.json --strict

Exit code 0 if all checks pass (or only expected-coverage shortfalls within
configured tolerance); 1 if --strict and any check fails.
"""
import argparse
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_CONFIG_PATH = SCRIPT_DIR / "configs" / "cafa3.json"


def load_config(path: str) -> dict[str, Any]:
    with open(path) as f:
        return json.load(f)


def load_split_ids(
    data_dir: Path, aspects: list[str], splits: list[str]
) -> tuple[set[str], dict[tuple[str, str], list[str]]]:
    """Union of protein IDs across aspects/splits — the coverage denominator."""
    ids: set[str] = set()
    per_split: dict[tuple[str, str], list[str]] = {}
    for a in aspects:
        for s in splits:
            f = data_dir / f"{a}_{s}_names.npy"
            if f.exists():
                names = [str(p) for p in np.load(f, allow_pickle=True)]
                per_split[(a, s)] = names
                ids.update(names)
    return ids, per_split


def resolve_dir(cache: Path, candidate_dirs: list[str]) -> Path:
    """Return the first candidate cache dir that exists and is non-empty.

    Resolving from a list of candidates keeps modality-specific naming out of
    the verification logic. Falls back to the first candidate (reported as
    missing downstream) if none qualify.
    """
    for d in candidate_dirs:
        p = cache / d
        if p.exists() and any(p.iterdir()):
            return p
    return cache / candidate_dirs[0]


def sample_check(
    emb_dir: Path, ids: set[str], expected_dim: int, n: int
) -> tuple[int, list[tuple[str, str]], list[str]]:
    """Load up to n existing embeddings; check dimension and finiteness."""
    bad_dim: list[tuple[str, str]] = []
    bad_finite: list[str] = []
    checked = 0
    for pid in ids:
        f = emb_dir / f"{pid}.npy"
        if not f.exists():
            continue
        try:
            arr = np.load(f)
        except Exception as e:  # noqa: BLE001 — report any load failure
            bad_dim.append((pid, f"unreadable: {e}"))
            continue
        flat = np.squeeze(arr)
        if flat.ndim != 1 or flat.shape[0] != expected_dim:
            bad_dim.append((pid, f"shape={arr.shape}"))
        if not np.all(np.isfinite(flat)):
            bad_finite.append(pid)
        checked += 1
        if checked >= n:
            break
    return checked, bad_dim, bad_finite


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="data")
    ap.add_argument("--config", default=DEFAULT_CONFIG_PATH,
                    help="JSON config fully specifying the verification (default: %(default)s)")
    ap.add_argument("--strict", action="store_true")
    args = ap.parse_args()

    cfg = load_config(args.config)

    required_keys = {
        "aspects",
        "splits",
        "cache_dir",
        "modalities",
        "sample_size",
        "catastrophic_factor",
    }

    missing = required_keys - cfg.keys()
    if missing:
        print(
            f"FAIL: configuration '{args.config}' is missing required keys: "
            f"{', '.join(sorted(missing))}"
        )
        sys.exit(1)

    data_dir = Path(args.data_dir)
    cache = data_dir / cfg["cache_dir"]
    print(f"==> Verifying embeddings under {cache}  (config: {args.config})")
    if not cache.exists():
        print(f"FAIL: {cache} does not exist.")
        sys.exit(1)

    all_ids, per_split = load_split_ids(data_dir, cfg["aspects"], cfg["splits"])
    if not all_ids:
        print(f"FAIL: no *_names.npy splits found in {data_dir} — run prepare_data first.")
        sys.exit(1)
    print(f"==> {len(all_ids)} unique protein IDs across {len(per_split)} aspect/split files\n")

    hard_fail = False
    for mod, spec in cfg["modalities"].items():

        required_modality_keys = {"dirs", "dim", "min_coverage"}

        missing = required_modality_keys - spec.keys()
        if missing:
            print(
                f"FAIL: modality '{mod}' is missing required keys: "
                f"{', '.join(sorted(missing))}"
            )
            hard_fail = True
            continue

        emb_dir = resolve_dir(cache, spec["dirs"])
        dim, min_cov = spec["dim"], spec["min_coverage"]

        print("-" * 60)
        print(f"Modality : {mod}")
        print(f"Cache    : {emb_dir.name}")
        print(f"Dim      : {dim}")
        print(f"Coverage : ≥{min_cov:.0%}")
        print("-" * 60)

        if not emb_dir.exists():
            print("  FAIL: directory missing\n")
            hard_fail = True
            continue

        present = sum(1 for pid in all_ids if (emb_dir / f"{pid}.npy").exists())
        cov = present / len(all_ids)
        n_files = sum(1 for _ in emb_dir.glob("*.npy"))
        checked, bad_dim, bad_finite = sample_check(emb_dir, all_ids, dim, cfg["sample_size"])

        dim_ok, fin_ok = not bad_dim, not bad_finite
        print(f"  files on disk:   {n_files}")
        print(f"  coverage:        {present}/{len(all_ids)} = {cov:.1%}  "
              f"[{'OK' if cov >= min_cov else 'LOW — possible partial generation'}]")
        print(f"  dim check:       {checked} sampled  [{'OK' if dim_ok else 'BAD SHAPES'}]")
        for pid, why in bad_dim[:3]:
            print(f"      {pid}: {why}")
        print(f"  finiteness:      [{'OK' if fin_ok else f'NaN/Inf in {len(bad_finite)} sampled'}]")

        if not dim_ok or not fin_ok:
            hard_fail = True
        if cov < min_cov * cfg["catastrophic_factor"]:
            hard_fail = True
        print()

    print("=" * 60)
    if hard_fail:
        print("RESULT: FAIL — wrong shapes, non-finite values, or coverage far below "
              "expected. Do NOT train on these.")
        sys.exit(1 if args.strict else 0)
    print("RESULT: PASS — all modalities present, correct dim, finite, coverage within bounds.")
    sys.exit(0)


if __name__ == "__main__":
    main()