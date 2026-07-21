#!/usr/bin/env python3
"""Evaluate existing PFP checkpoints against an arbitrary prepared benchmark."""

from __future__ import annotations

import argparse
import json
import re
import shutil
import sys
import tempfile
from pathlib import Path
from typing import Any, Dict

import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import DataLoader

from common import (
    MODALITY_MODES,
    active_modalities,
    atomic_write_json,
    expected_result_dir,
    file_snapshot,
    load_run_config,
    modality_paths,
    require_unchanged,
    selected_aspects,
    sha256_file,
    validate_mandatory_metrics,
)
from prediction_artifacts import (
    EvaluationArrayCapture,
    publish_prediction_artifacts,
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

    kwargs: Dict[str, Any] = {"norm": "cafa", "prop": "max", "no_orphans": False}
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


def load_passed_report(path: Path, label: str) -> tuple[dict[str, Any], dict[str, Any]]:
    resolved = path.resolve()
    snapshot = file_snapshot(resolved)
    value = json.loads(resolved.read_text(encoding="utf-8"))
    require_unchanged(resolved, snapshot, label)
    if not isinstance(value, dict) or value.get("status") != "passed":
        raise ValueError(f"{label} does not declare status=passed: {resolved}")
    return value, snapshot


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pfp-root", type=Path, required=True)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--cache-root", type=Path, required=True)
    parser.add_argument("--obo-file", type=Path, required=True)
    parser.add_argument("--checkpoint-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--mode", choices=MODALITY_MODES, required=True)
    parser.add_argument("--aspect", action="append", default=[])
    parser.add_argument("--num-workers", type=int, default=0)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--ia-file-dir", type=Path)
    parser.add_argument("--benchmark-id")
    parser.add_argument("--framework-commit")
    parser.add_argument("--pfp-commit")
    parser.add_argument("--preparation-report", type=Path)
    parser.add_argument("--embedding-report", type=Path)
    parser.add_argument(
        "--prediction-artifact-dir",
        type=Path,
        help=(
            "Opt in to publishing the exact prediction/truth arrays and IA files "
            "needed for separate sensitivity analyses."
        ),
    )
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
    prediction_destination = (
        args.prediction_artifact_dir.resolve()
        if args.prediction_artifact_dir is not None
        else None
    )
    prediction_stage = None
    prediction_manifest: Dict[str, Any] | None = None
    provenance_inputs: list[tuple[Path, dict[str, Any], str]] = []
    if prediction_destination is not None:
        if prediction_destination.exists():
            raise ValueError(
                f"Prediction artifact directory already exists: {prediction_destination}"
            )
        prediction_destination.parent.mkdir(parents=True, exist_ok=True)
        required_provenance = {
            "--framework-commit": args.framework_commit,
            "--pfp-commit": args.pfp_commit,
            "--preparation-report": args.preparation_report,
            "--embedding-report": args.embedding_report,
        }
        missing_provenance = [
            option for option, value in required_provenance.items() if value is None
        ]
        if missing_provenance:
            raise ValueError(
                "Prediction capture requires provenance arguments: "
                + ", ".join(missing_provenance)
            )
        for option, commit in (
            ("--framework-commit", args.framework_commit),
            ("--pfp-commit", args.pfp_commit),
        ):
            if not isinstance(commit, str) or re.fullmatch(r"[0-9a-f]{40}", commit) is None:
                raise ValueError(f"{option} must be a 40-character lowercase git commit")
        assert args.preparation_report is not None
        assert args.embedding_report is not None
        preparation, preparation_snapshot = load_passed_report(
            args.preparation_report, "Preparation report"
        )
        embedding, embedding_snapshot = load_passed_report(
            args.embedding_report, "Embedding validation report"
        )
        if embedding.get("mode") != args.mode:
            raise ValueError("Embedding validation mode differs from evaluation mode")
        embedded_preparation = embedding.get("preparation_report", {})
        if embedded_preparation.get("sha256") != preparation_snapshot["sha256"]:
            raise ValueError(
                "Embedding validation report is not bound to the selected preparation report"
            )
        benchmark_fingerprint = preparation.get("benchmark_fingerprint")
        if not isinstance(benchmark_fingerprint, str) or not benchmark_fingerprint:
            raise ValueError("Preparation report lacks a benchmark fingerprint")
        provenance_inputs = [
            (Path(preparation_snapshot["path"]), preparation_snapshot, "Preparation report"),
            (Path(embedding_snapshot["path"]), embedding_snapshot, "Embedding validation report"),
        ]
        prediction_stage = Path(
            tempfile.mkdtemp(
                prefix=f".{prediction_destination.name}.stage-",
                dir=str(prediction_destination.parent),
            )
        )
        copied_reports = {}
        try:
            for name, snapshot in (
                ("preparation_report.json", preparation_snapshot),
                ("embedding_validation_report.json", embedding_snapshot),
            ):
                destination = prediction_stage / name
                shutil.copyfile(snapshot["path"], destination)
                if sha256_file(destination) != snapshot["sha256"]:
                    raise ValueError(f"Copied provenance report differs from source: {name}")
                copied_reports[name] = {
                    "artifact_file": name,
                    "source_path": snapshot["path"],
                    "bytes": snapshot["bytes"],
                    "sha256": snapshot["sha256"],
                }
        except BaseException:
            shutil.rmtree(prediction_stage, ignore_errors=True)
            prediction_stage = None
            raise
        prediction_manifest = {
            "schema_version": 2,
            "status": "complete",
            "benchmark_id": args.benchmark_id or str(config.get("name", "unspecified")),
            "mode": args.mode,
            "seed": args.seed,
            "selected_aspects": aspects,
            "pfp_root": str(pfp_root),
            "data_dir": str(args.data_dir.resolve()),
            "cache_root": str(args.cache_root.resolve()),
            "config": {
                "path": str(args.config.resolve()),
                "sha256": sha256_file(args.config.resolve()),
            },
            "obo": {
                "path": str(args.obo_file.resolve()),
                "sha256": sha256_file(args.obo_file.resolve()),
            },
            "provenance": {
                "framework_commit": args.framework_commit,
                "pfp_commit": args.pfp_commit,
                "benchmark_fingerprint": benchmark_fingerprint,
                "source_csv_sha256": preparation.get("source_csv_sha256", {}),
                "preparation_report": copied_reports["preparation_report.json"],
                "embedding_validation_report": copied_reports[
                    "embedding_validation_report.json"
                ],
            },
            "capture_policy": (
                "observe-standard-pfp-cafa-writers-without-extra-inference; "
                "canonical-results-unchanged"
            ),
            "aspects": {},
        }
    paths = modality_paths(args.cache_root.resolve(), config)
    empty_dir = args.output_dir / "_empty_modality"
    empty_dir.mkdir(parents=True, exist_ok=True)
    embedding_dirs = {
        "prott5": str(paths["sequence"]),
        "text": str(paths["text"]),
        "struct": str(paths["structure"]),
        "ppi": str(paths["ppi"]),
    }
    active = set(active_modalities(args.mode))
    pfp_names = {"text": "text", "structure": "struct", "ppi": "ppi"}
    for modality, pfp_name in pfp_names.items():
        if modality not in active:
            embedding_dirs[pfp_name] = str(empty_dir)

    device = "cuda" if torch.cuda.is_available() else "cpu"
    set_seed(args.seed)
    summary: Dict[str, Any] = {
        "schema_version": 1,
        "mode": args.mode,
        "device": device,
        "aspects": {},
    }
    try:
        for aspect in aspects:
            checkpoint = (
                expected_result_dir(args.checkpoint_root.resolve(), aspect)
                / "best_model.pt"
            )
            if not checkpoint.is_file():
                raise FileNotFoundError(
                    f"Checkpoint is missing for {aspect}: {checkpoint}"
                )
            checkpoint_sha256 = sha256_file(checkpoint)
            train_dataset = MultiModalDataset(
                str(args.data_dir),
                embedding_dirs,
                "prott5",
                aspect,
                "train",
                normalize="standard",
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
                (args.data_dir / f"{aspect}_go_terms.json").read_text(
                    encoding="utf-8"
                )
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
                model,
                test_loader,
                criterion,
                device,
                compute_weight_stats_detailed=True,
            )
            ia_file = None
            ia_file_sha256 = None
            if args.ia_file_dir:
                candidate = args.ia_file_dir / f"{aspect}_ia.txt"
                if candidate.is_file():
                    ia_file = str(candidate)
                    ia_file_sha256 = sha256_file(candidate)
            if require_precomputed_ia and ia_file is None:
                raise FileNotFoundError(
                    f"This run config requires precomputed IA for {aspect}"
                )
            output = args.output_dir / aspect
            evaluation_kwargs = {
                "model": model,
                "loader": test_loader,
                "device": device,
                "protein_ids": test_dataset.protein_ids.tolist(),
                "go_terms": go_terms,
                "obo_file": str(args.obo_file),
                "output_dir": str(output / "cafa_eval"),
                "model_name": f"prott5_{training['fusion_type']}",
                "train_labels": train_dataset.labels,
                "ia_file": ia_file,
            }
            capture = None
            if prediction_stage is not None:
                capture = EvaluationArrayCapture(
                    pfp_evaluation, aspect, prediction_stage
                )
                with capture:
                    cafa = evaluate_with_cafa(**evaluation_kwargs)
            else:
                cafa = evaluate_with_cafa(**evaluation_kwargs)
            if sha256_file(checkpoint) != checkpoint_sha256:
                raise ValueError(f"Checkpoint changed during evaluation for {aspect}")
            if ia_file and sha256_file(Path(ia_file)) != ia_file_sha256:
                raise ValueError(f"Precomputed IA changed during evaluation for {aspect}")

            result: Dict[str, Any] = {
                "aspect": aspect,
                "mode": args.mode,
                "checkpoint": str(checkpoint),
                "checkpoint_sha256": checkpoint_sha256,
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
            result["ia_file_sha256"] = ia_file_sha256
            atomic_write_json(output / "results.json", result)
            summary["aspects"][aspect] = result

            if capture is not None:
                assert prediction_manifest is not None
                prediction_manifest["aspects"][aspect] = capture.persist(
                    expected_protein_ids=test_dataset.protein_ids.tolist(),
                    expected_go_terms=go_terms,
                    checkpoint=checkpoint,
                    expected_checkpoint_sha256=checkpoint_sha256,
                    cafa_metrics=cafa,
                    ia_file=Path(ia_file) if ia_file else None,
                    expected_ia_sha256=ia_file_sha256,
                )

        if prediction_stage is not None:
            assert prediction_destination is not None
            assert prediction_manifest is not None
            for path, snapshot, label in provenance_inputs:
                require_unchanged(path, snapshot, label)
            publish_prediction_artifacts(
                prediction_stage, prediction_destination, prediction_manifest
            )
            prediction_stage = None
            try:
                artifact_reference = prediction_destination.relative_to(
                    args.output_dir.resolve()
                ).as_posix()
            except ValueError:
                artifact_reference = str(prediction_destination)
            summary["prediction_artifact_dir"] = artifact_reference
            summary["prediction_artifact_manifest_sha256"] = sha256_file(
                prediction_destination / "prediction_artifact_manifest.json"
            )

        atomic_write_json(args.output_dir / "evaluation_summary.json", summary)
        print(json.dumps(summary, indent=2, sort_keys=True))
        return 0
    finally:
        if prediction_stage is not None:
            shutil.rmtree(prediction_stage, ignore_errors=True)


if __name__ == "__main__":
    raise SystemExit(main())
