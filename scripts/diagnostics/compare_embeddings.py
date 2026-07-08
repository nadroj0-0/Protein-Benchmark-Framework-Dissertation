#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np
from tqdm import tqdm


@dataclass
class Row:
    modality: str
    protein_id: str
    generated_file: str
    published_file: str
    status: str
    generated_shape: str
    published_shape: str
    generated_dtype: str
    published_dtype: str
    same_shape: bool
    same_dtype: bool
    sha256_equal: bool
    cosine: float | None
    l2: float | None
    max_abs: float | None
    mean_abs: float | None
    generated_finite: bool
    published_finite: bool
    note: str


MODALITIES = {
    "prott5": {
        "generated": "data/embedding_cache/prott5",
        "published": "published/data/embedding_cache/prott5",
    },
    "text_temporal": {
        "generated": "data/embedding_cache/exp_text_embeddings",
        "published": "published/data/embedding_cache/exp_text_embeddings_temporal",
    },
    "text_current": {
        "generated": "data/embedding_cache/exp_text_embeddings",
        "published": "published/data/embedding_cache/exp_text_embeddings",
    },
    "structure": {
        "generated": "data/embedding_cache/IF1",
        "published": "published/data/embedding_cache/IF1",
    },
    "ppi": {
        "generated": "data/embedding_cache/ppi",
        "published": "published/data/embedding_cache/ppi",
    },
}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def flatten_for_numeric(x: np.ndarray) -> np.ndarray:
    return np.asarray(x, dtype=np.float64).reshape(-1)


def numeric_compare(a: np.ndarray, b: np.ndarray) -> tuple[float | None, float | None, float | None, float | None, str]:
    if a.shape != b.shape:
        return None, None, None, None, "shape mismatch"

    af = flatten_for_numeric(a)
    bf = flatten_for_numeric(b)

    if af.size == 0 or bf.size == 0:
        return None, None, None, None, "empty array"

    if not np.isfinite(af).all() or not np.isfinite(bf).all():
        return None, None, None, None, "non-finite values"

    diff = af - bf
    l2 = float(np.linalg.norm(diff))
    max_abs = float(np.max(np.abs(diff)))
    mean_abs = float(np.mean(np.abs(diff)))

    denom = float(np.linalg.norm(af) * np.linalg.norm(bf))
    cosine = float(np.dot(af, bf) / denom) if denom > 0 else None

    return cosine, l2, max_abs, mean_abs, ""


def compare_one(modality: str, protein_id: str, gen_file: Path | None, pub_file: Path | None) -> Row:
    if gen_file is None:
        return Row(modality, protein_id, "", str(pub_file), "missing_generated", "", "", "", "", False, False, False, None, None, None, None, False, False, "")
    if pub_file is None:
        return Row(modality, protein_id, str(gen_file), "", "missing_published", "", "", "", "", False, False, False, None, None, None, None, False, False, "")

    try:
        g = np.load(gen_file, allow_pickle=False)
        p = np.load(pub_file, allow_pickle=False)
    except Exception as e:
        return Row(modality, protein_id, str(gen_file), str(pub_file), "load_error", "", "", "", "", False, False, False, None, None, None, None, False, False, repr(e))

    g_finite = bool(np.isfinite(g).all())
    p_finite = bool(np.isfinite(p).all())
    same_shape = g.shape == p.shape
    same_dtype = g.dtype == p.dtype

    try:
        hash_equal = sha256_file(gen_file) == sha256_file(pub_file)
    except Exception:
        hash_equal = False

    cosine, l2, max_abs, mean_abs, note = numeric_compare(g, p)

    if not g_finite or not p_finite:
        status = "non_finite"
    elif not same_shape:
        status = "shape_mismatch"
    elif hash_equal:
        status = "exact_match"
    elif cosine is not None and cosine >= 0.999999 and max_abs is not None and max_abs <= 1e-5:
        status = "numeric_match"
    else:
        status = "different"

    return Row(
        modality=modality,
        protein_id=protein_id,
        generated_file=str(gen_file),
        published_file=str(pub_file),
        status=status,
        generated_shape=str(tuple(g.shape)),
        published_shape=str(tuple(p.shape)),
        generated_dtype=str(g.dtype),
        published_dtype=str(p.dtype),
        same_shape=same_shape,
        same_dtype=same_dtype,
        sha256_equal=hash_equal,
        cosine=cosine,
        l2=l2,
        max_abs=max_abs,
        mean_abs=mean_abs,
        generated_finite=g_finite,
        published_finite=p_finite,
        note=note,
    )


def npy_map(directory: Path) -> dict[str, Path]:
    if not directory.exists():
        return {}
    return {p.stem: p for p in directory.glob("*.npy")}


def summarise(rows: list[Row]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for modality in sorted({r.modality for r in rows}):
        rs = [r for r in rows if r.modality == modality]
        cosines = [r.cosine for r in rs if r.cosine is not None and math.isfinite(r.cosine)]
        l2s = [r.l2 for r in rs if r.l2 is not None and math.isfinite(r.l2)]
        max_abs = [r.max_abs for r in rs if r.max_abs is not None and math.isfinite(r.max_abs)]

        statuses: dict[str, int] = {}
        for r in rs:
            statuses[r.status] = statuses.get(r.status, 0) + 1

        out[modality] = {
            "total_union": len(rs),
            "statuses": statuses,
            "same_shape": sum(r.same_shape for r in rs),
            "same_dtype": sum(r.same_dtype for r in rs),
            "sha256_equal": sum(r.sha256_equal for r in rs),
            "generated_finite": sum(r.generated_finite for r in rs),
            "published_finite": sum(r.published_finite for r in rs),
            "cosine_mean": float(np.mean(cosines)) if cosines else None,
            "cosine_min": float(np.min(cosines)) if cosines else None,
            "l2_mean": float(np.mean(l2s)) if l2s else None,
            "l2_max": float(np.max(l2s)) if l2s else None,
            "max_abs_mean": float(np.mean(max_abs)) if max_abs else None,
            "max_abs_max": float(np.max(max_abs)) if max_abs else None,
        }
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", default=".", help="PFP repo root")
    parser.add_argument("--out-csv", default="results/embedding_comparison.csv")
    parser.add_argument("--out-json", default="results/embedding_comparison_summary.json")
    parser.add_argument("--limit", type=int, default=None)
    args = parser.parse_args()

    root = Path(args.root).resolve()
    all_rows: list[Row] = []

    for modality, cfg in MODALITIES.items():
        gen_dir = root / cfg["generated"]
        pub_dir = root / cfg["published"]

        gen = npy_map(gen_dir)
        pub = npy_map(pub_dir)

        if not gen and not pub:
            continue

        ids = sorted(set(gen) | set(pub))
        if args.limit:
            ids = ids[: args.limit]

        print(f"\n==> Comparing {modality}")
        print(f"    generated: {gen_dir} ({len(gen)} files)")
        print(f"    published: {pub_dir} ({len(pub)} files)")
        print(f"    union:     {len(ids)} proteins")

        for pid in tqdm(ids, desc=modality):
            all_rows.append(compare_one(modality, pid, gen.get(pid), pub.get(pid)))

    out_csv = root / args.out_csv
    out_json = root / args.out_json
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    out_json.parent.mkdir(parents=True, exist_ok=True)

    with out_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(asdict(all_rows[0]).keys()) if all_rows else list(Row.__annotations__.keys()))
        writer.writeheader()
        for row in all_rows:
            writer.writerow(asdict(row))

    summary = summarise(all_rows)
    with out_json.open("w") as f:
        json.dump(summary, f, indent=2)

    print("\n==> Wrote:")
    print(f"    {out_csv}")
    print(f"    {out_json}")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()