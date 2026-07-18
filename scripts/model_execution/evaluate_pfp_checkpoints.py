#!/usr/bin/env python3
"""Evaluate existing PFP checkpoints against an arbitrary prepared benchmark."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict

import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import DataLoader

from common import (
    atomic_write_json,
    expected_result_dir,
    load_run_config,
    modality_paths,
    selected_aspects,
    sha256_file,
    validate_mandatory_metrics,
)


def strict_cafa_runner(
    obo_file: str,
    pred_file: str,
    truth_file: str,
    output_dir: str,
    ia_file: str | None = None,
) -> Path:
    """Run cafaeval with the full PFP argument contract and no silent fallback."""
    from cafaeval.evaluation import cafa_eval, write_results

    kwargs: Dict[str, str] = {"norm": "cafa", "prop": "max"}
    if ia_file:
        kwargs["ia"] = str(ia_file)
    results = cafa_eval(str(obo_file), str(pred_file), str(truth_file), **kwargs)
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    write_results(*results, out_dir=str(destination))
    return destination


def load_state(path: Path, device: str) -> Dict[str, torch.Tensor]:
    try:
        return torch.load(path, map_location=device, weights_only=True)
    except TypeError:
        return torch.load(path, map_location=device)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pfp-root", type=Path, required=True)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--cache-root", type=Path, required=True)
    parser.add_argument("--obo-file", type=Path, required=True)
    parser.add_argument("--checkpoint-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--mode", choices=("full", "sequence-only"), required=True)
    parser.add_argument("--aspect", action="append", default=[])
    parser.add_argument("--num-workers", type=int, default=0)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--ia-file-dir", type=Path)
    args = parser.parse_args()

    pfp_root = args.pfp_root.resolve()
    if not (pfp_root / "mmfp" / "dataset.py").is_file():
        raise FileNotFoundError(f"Not a PFP checkout: {pfp_root}")
    if not args.obo_file.is_file():
        raise FileNotFoundError(f"GO OBO file is missing: {args.obo_file}")
    sys.path.insert(0, str(pfp_root))
    from mmfp.dataset import MultiModalDataset, collate_fn  # noqa: E402
    import mmfp.evaluation as pfp_evaluation  # noqa: E402
    from mmfp.models import create_model  # noqa: E402
    from train import evaluate, set_seed  # noqa: E402

    pfp_evaluation.run_cafa_evaluator = strict_cafa_runner
    evaluate_with_cafa = pfp_evaluation.evaluate_with_cafa

    config = load_run_config(args.config)
    training = config["training"]
    require_precomputed_ia = bool(
        config.get("evaluation", {}).get("require_precomputed_ia", False)
    )
    aspects = selected_aspects(args.aspect)
    paths = modality_paths(args.cache_root.resolve(), config)
    empty_dir = args.output_dir / "_empty_modality"
    empty_dir.mkdir(parents=True, exist_ok=True)
    embedding_dirs = {
        "prott5": str(paths["sequence"]),
        "text": str(paths["text"]),
        "struct": str(paths["structure"]),
        "ppi": str(paths["ppi"]),
    }
    if args.mode == "sequence-only":
        embedding_dirs.update(
            {"text": str(empty_dir), "struct": str(empty_dir), "ppi": str(empty_dir)}
        )

    device = "cuda" if torch.cuda.is_available() else "cpu"
    set_seed(args.seed)
    summary: Dict[str, Any] = {
        "schema_version": 1,
        "mode": args.mode,
        "device": device,
        "aspects": {},
    }
    for aspect in aspects:
        checkpoint = expected_result_dir(args.checkpoint_root.resolve(), aspect) / "best_model.pt"
        if not checkpoint.is_file():
            raise FileNotFoundError(f"Checkpoint is missing for {aspect}: {checkpoint}")
        train_dataset = MultiModalDataset(
            str(args.data_dir), embedding_dirs, "prott5", aspect, "train", normalize="standard"
        )
        test_dataset = MultiModalDataset(
            str(args.data_dir),
            embedding_dirs,
            "prott5",
            aspect,
            "test",
            normalize="standard",
            norm_stats=train_dataset.norm_stats,
        )
        test_loader = DataLoader(
            test_dataset,
            batch_size=32,
            shuffle=False,
            collate_fn=collate_fn,
            num_workers=args.num_workers,
            pin_memory=(device == "cuda"),
        )
        go_terms = json.loads(
            (args.data_dir / f"{aspect}_go_terms.json").read_text(encoding="utf-8")
        )
        model = create_model(
            fusion_type=training["fusion_type"],
            seq_dim=1024,
            text_dim=768,
            struct_dim=512,
            ppi_dim=512,
            hidden_dim=int(training["hidden_dim"]),
            num_go_terms=len(go_terms),
            dropout=float(training["dropout"]),
            modality_dropout=float(training["modality_dropout"]),
            use_late_fusion=bool(training["use_late_fusion"]),
            late_output_mode=training["late_output_mode"],
        ).to(device)
        model.load_state_dict(load_state(checkpoint, device), strict=True)
        criterion = nn.BCEWithLogitsLoss().to(device)
        ordinary = evaluate(
            model, test_loader, criterion, device, compute_weight_stats_detailed=True
        )
        ia_file = None
        if args.ia_file_dir:
            candidate = args.ia_file_dir / f"{aspect}_ia.txt"
            if candidate.is_file():
                ia_file = str(candidate)
        if require_precomputed_ia and ia_file is None:
            raise FileNotFoundError(
                f"This run config requires precomputed IA for {aspect}"
            )
        output = args.output_dir / aspect
        cafa = evaluate_with_cafa(
            model=model,
            loader=test_loader,
            device=device,
            protein_ids=test_dataset.protein_ids.tolist(),
            go_terms=go_terms,
            obo_file=str(args.obo_file),
            output_dir=str(output / "cafa_eval"),
            model_name=f"prott5_{training['fusion_type']}",
            train_labels=train_dataset.labels,
            ia_file=ia_file,
        )
        result: Dict[str, Any] = {
            "aspect": aspect,
            "mode": args.mode,
            "checkpoint": str(checkpoint),
            "checkpoint_sha256": sha256_file(checkpoint),
            "num_go_terms": len(go_terms),
            "test_fmax": float(ordinary["fmax"]),
            "test_micro_auprc": float(ordinary["micro_auprc"]),
            "test_macro_auprc": float(ordinary["macro_auprc"]),
            "cafa_evaluator_policy": "strict-ia-norm-cafa-prop-max-no-fallback",
        }
        for key, value in cafa.items():
            if isinstance(value, (int, float, np.integer, np.floating)):
                result[f"cafa_{key}"] = float(value)
        validate_mandatory_metrics(result, f"CAFA evaluation for {aspect}")
        result["ia_file"] = ia_file
        result["ia_file_sha256"] = sha256_file(Path(ia_file)) if ia_file else None
        atomic_write_json(output / "results.json", result)
        summary["aspects"][aspect] = result
    atomic_write_json(args.output_dir / "evaluation_summary.json", summary)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
