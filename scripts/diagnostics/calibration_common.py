#!/usr/bin/env python3
"""Deterministic calibration helpers for captured PFP prediction arrays."""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass
from typing import Any, Mapping, Sequence

import numpy as np
from scipy.optimize import minimize
from scipy.special import expit

from label_space_common import OboGraph


def sha256_json(value: Any) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _validate_matrices(scores: np.ndarray, truth: np.ndarray) -> None:
    if scores.shape != truth.shape or scores.ndim != 2:
        raise ValueError("Calibration scores and truth must be matching matrices")
    if not np.isfinite(scores).all() or np.any(scores < 0) or np.any(scores > 1):
        raise ValueError("Calibration scores must be finite probabilities")
    if not np.isin(truth, (0, 1)).all():
        raise ValueError("Calibration truth must be binary")


@dataclass(frozen=True)
class CalibrationPolicy:
    score_clip_epsilon: float = 1e-6
    bin_l2: float = 5.0
    term_l2: float = 20.0
    minimum_bin_positives: int = 20
    minimum_bin_negatives: int = 20
    minimum_term_positives: int = 20
    minimum_term_negatives: int = 20
    maximum_iterations: int = 200
    optimizer_tolerance: float = 1e-7
    protein_chunk_size: int = 256

    def validate(self) -> None:
        if not 0 < self.score_clip_epsilon < 0.5:
            raise ValueError("score_clip_epsilon must lie between zero and 0.5")
        if self.bin_l2 <= 0 or self.term_l2 <= 0:
            raise ValueError("Calibration L2 penalties must be positive")
        for name in (
            "minimum_bin_positives",
            "minimum_bin_negatives",
            "minimum_term_positives",
            "minimum_term_negatives",
            "maximum_iterations",
            "protein_chunk_size",
        ):
            if getattr(self, name) < 1:
                raise ValueError(f"{name} must be positive")
        if self.optimizer_tolerance <= 0:
            raise ValueError("optimizer_tolerance must be positive")

    def as_dict(self) -> dict[str, Any]:
        return {
            key: getattr(self, key)
            for key in (
                "score_clip_epsilon",
                "bin_l2",
                "term_l2",
                "minimum_bin_positives",
                "minimum_bin_negatives",
                "minimum_term_positives",
                "minimum_term_negatives",
                "maximum_iterations",
                "optimizer_tolerance",
                "protein_chunk_size",
            )
        }


def _logit_scores(scores: np.ndarray, epsilon: float) -> np.ndarray:
    clipped = np.clip(scores, epsilon, 1.0 - epsilon)
    return np.log(clipped) - np.log1p(-clipped)


def propagate_scores_max(
    scores: np.ndarray,
    go_terms: Sequence[str],
    graph: OboGraph,
    aspect: str,
) -> np.ndarray:
    """Apply CAFA-style max propagation over the is_a + part_of output graph."""
    if scores.ndim != 2 or scores.shape[1] != len(go_terms):
        raise ValueError("Score propagation term dimensions differ")
    if len(set(go_terms)) != len(go_terms):
        raise ValueError("Score propagation GO terms must be unique")
    term_index = {term: index for index, term in enumerate(go_terms)}
    closure = graph.ancestor_closure(aspect)
    propagated = np.asarray(scores, dtype=np.float64).copy()
    original = np.asarray(scores, dtype=np.float64)
    for child, child_index in term_index.items():
        ancestors = closure.get(child)
        if not ancestors:
            raise ValueError(
                f"GO term is absent/disconnected during propagation: {child}"
            )
        for ancestor in ancestors:
            parent_index = term_index.get(ancestor)
            if parent_index is not None and parent_index != child_index:
                np.maximum(
                    propagated[:, parent_index],
                    original[:, child_index],
                    out=propagated[:, parent_index],
                )
    if (
        not np.isfinite(propagated).all()
        or np.any(propagated < 0)
        or np.any(propagated > 1)
    ):
        raise RuntimeError("Post-propagation scores are invalid")
    return propagated


def fit_monotone_hierarchical_calibrator(
    scores: np.ndarray,
    truth: np.ndarray,
    term_ids: Sequence[str],
    term_bins: Sequence[str],
    policy: CalibrationPolicy,
) -> dict[str, Any]:
    """Fit a positive-slope logistic model with regularized bin/term intercepts."""
    policy.validate()
    _validate_matrices(scores, truth)
    protein_count, term_count = scores.shape
    if term_count != len(term_ids) or term_count != len(term_bins):
        raise ValueError("Calibration term metadata does not match matrix columns")
    if len(set(term_ids)) != len(term_ids):
        raise ValueError("Calibration term IDs must be unique")
    if any(not value for value in term_bins):
        raise ValueError("Every calibration term requires a non-empty IA-bin label")
    positives = truth.sum(axis=0, dtype=np.int64)
    negatives = protein_count - positives
    total_positive = int(positives.sum())
    total_events = protein_count * term_count
    total_negative = total_events - total_positive
    if not total_positive or not total_negative:
        result = {
            "family": "uncalibrated_insufficient_support",
            "status": "uncalibrated_insufficient_support",
            "reason": "calibration population does not contain both outcome classes",
            "term_ids": list(term_ids),
            "term_bins": list(term_bins),
            "fallback_by_term": ["uncalibrated"] * term_count,
            "support": {
                "proteins": protein_count,
                "events": total_events,
                "positives": total_positive,
                "negatives": total_negative,
                "term_positives": positives.tolist(),
                "term_negatives": negatives.tolist(),
            },
            "policy": policy.as_dict(),
        }
        result["model_sha256"] = sha256_json(result)
        return result

    bin_labels = sorted(set(term_bins))
    bin_lookup = {label: index for index, label in enumerate(bin_labels)}
    term_bin_index = np.asarray(
        [bin_lookup[value] for value in term_bins], dtype=np.int64
    )
    bin_positives = np.bincount(
        term_bin_index, weights=positives, minlength=len(bin_labels)
    ).astype(np.int64)
    bin_events = np.bincount(
        term_bin_index,
        weights=np.full(term_count, protein_count, dtype=np.int64),
        minlength=len(bin_labels),
    ).astype(np.int64)
    bin_negatives = bin_events - bin_positives
    supported_bins = np.flatnonzero(
        (bin_positives >= policy.minimum_bin_positives)
        & (bin_negatives >= policy.minimum_bin_negatives)
    )
    supported_terms = np.flatnonzero(
        (positives >= policy.minimum_term_positives)
        & (negatives >= policy.minimum_term_negatives)
    )
    bin_parameter = {int(value): index for index, value in enumerate(supported_bins)}
    term_parameter = {int(value): index for index, value in enumerate(supported_terms)}
    bin_slice = slice(2, 2 + len(supported_bins))
    term_slice = slice(bin_slice.stop, bin_slice.stop + len(supported_terms))
    parameter_count = term_slice.stop

    prevalence = total_positive / total_events
    initial = np.zeros(parameter_count, dtype=np.float64)
    initial[0] = math.log(prevalence / (1.0 - prevalence))
    initial[1] = 0.0

    def objective(parameters: np.ndarray) -> tuple[float, np.ndarray]:
        alpha = float(parameters[0])
        slope = math.exp(float(parameters[1]))
        bin_effect_by_term = np.zeros(term_count, dtype=np.float64)
        for bin_index, parameter_index in bin_parameter.items():
            bin_effect_by_term[term_bin_index == bin_index] = parameters[
                bin_slice.start + parameter_index
            ]
        term_effect_by_term = np.zeros(term_count, dtype=np.float64)
        if supported_terms.size:
            term_effect_by_term[supported_terms] = parameters[term_slice]

        loss = 0.0
        gradient = np.zeros_like(parameters)
        residual_by_term = np.zeros(term_count, dtype=np.float64)
        for start in range(0, protein_count, policy.protein_chunk_size):
            stop = min(protein_count, start + policy.protein_chunk_size)
            x = _logit_scores(scores[start:stop], policy.score_clip_epsilon)
            y = truth[start:stop]
            eta = (
                alpha
                + slope * x
                + bin_effect_by_term[np.newaxis, :]
                + term_effect_by_term[np.newaxis, :]
            )
            loss += float(np.sum(np.logaddexp(0.0, eta) - y * eta))
            residual = expit(eta) - y
            gradient[0] += float(residual.sum())
            gradient[1] += float(np.sum(residual * (slope * x)))
            residual_by_term += residual.sum(axis=0)

        if supported_bins.size:
            aggregated = np.bincount(
                term_bin_index,
                weights=residual_by_term,
                minlength=len(bin_labels),
            )
            gradient[bin_slice] = aggregated[supported_bins]
            bin_parameters = parameters[bin_slice]
            loss += 0.5 * policy.bin_l2 * float(np.dot(bin_parameters, bin_parameters))
            gradient[bin_slice] += policy.bin_l2 * bin_parameters
        if supported_terms.size:
            gradient[term_slice] = residual_by_term[supported_terms]
            term_parameters = parameters[term_slice]
            loss += (
                0.5 * policy.term_l2 * float(np.dot(term_parameters, term_parameters))
            )
            gradient[term_slice] += policy.term_l2 * term_parameters
        return loss, gradient

    bounds = [(-30.0, 30.0), (-8.0, 8.0)] + [(-20.0, 20.0)] * (parameter_count - 2)
    fit = minimize(
        objective,
        initial,
        method="L-BFGS-B",
        jac=True,
        bounds=bounds,
        options={
            "maxiter": policy.maximum_iterations,
            "ftol": policy.optimizer_tolerance,
            "gtol": policy.optimizer_tolerance,
            "maxls": 30,
        },
    )
    if not fit.success or not np.isfinite(fit.fun):
        raise RuntimeError(
            f"Calibration optimizer failed: status={fit.status}, message={fit.message}"
        )
    parameters = np.asarray(fit.x, dtype=np.float64)
    slope = math.exp(float(parameters[1]))
    if not math.isfinite(slope) or slope <= 0:
        raise RuntimeError("Calibration slope is not finite and positive")

    bin_effects = {
        bin_labels[bin_index]: float(parameters[bin_slice.start + parameter_index])
        for bin_index, parameter_index in bin_parameter.items()
    }
    term_effects = {
        term_ids[term_index]: float(parameters[term_slice.start + parameter_index])
        for term_index, parameter_index in term_parameter.items()
    }
    fallback_by_term = []
    for term_index, bin_index in enumerate(term_bin_index):
        if term_index in term_parameter:
            fallback_by_term.append("term_shrinkage")
        elif int(bin_index) in bin_parameter:
            fallback_by_term.append("aspect_mode_ia_bin")
        else:
            fallback_by_term.append("aspect_mode_platt")

    result = {
        "family": "monotone_multilevel_logistic",
        "formula": (
            "logit(q)=alpha+positive_slope*logit(score)+"
            "regularized_ia_bin_intercept+regularized_term_intercept"
        ),
        "alpha": float(parameters[0]),
        "positive_slope": slope,
        "bin_effects": bin_effects,
        "term_effects": term_effects,
        "term_ids": list(term_ids),
        "term_bins": list(term_bins),
        "fallback_by_term": fallback_by_term,
        "support": {
            "proteins": protein_count,
            "events": total_events,
            "positives": total_positive,
            "negatives": total_negative,
            "bin_labels": bin_labels,
            "bin_positives": bin_positives.tolist(),
            "bin_negatives": bin_negatives.tolist(),
            "term_positives": positives.tolist(),
            "term_negatives": negatives.tolist(),
            "supported_bins": [bin_labels[index] for index in supported_bins],
            "supported_terms": [term_ids[index] for index in supported_terms],
        },
        "policy": policy.as_dict(),
        "optimizer": {
            "name": "scipy_L-BFGS-B",
            "success": bool(fit.success),
            "status": int(fit.status),
            "message": str(fit.message),
            "iterations": int(fit.nit),
            "function_evaluations": int(fit.nfev),
            "objective": float(fit.fun),
        },
    }
    result["model_sha256"] = sha256_json(result)
    return result


def apply_calibrator(
    scores: np.ndarray,
    term_ids: Sequence[str],
    term_bins: Sequence[str],
    model: Mapping[str, Any],
    *,
    protein_chunk_size: int = 256,
) -> tuple[np.ndarray, list[str]]:
    if scores.ndim != 2 or scores.shape[1] != len(term_ids):
        raise ValueError("Calibration application term dimensions differ")
    if list(term_ids) != list(model["term_ids"]):
        raise ValueError("Calibration application GO-term order differs")
    if list(term_bins) != list(model["term_bins"]):
        raise ValueError("Calibration application IA-bin order differs")
    if protein_chunk_size < 1:
        raise ValueError("protein_chunk_size must be positive")
    if model["family"] == "uncalibrated_insufficient_support":
        return (
            np.full(scores.shape, np.nan, dtype=np.float32),
            list(model["fallback_by_term"]),
        )
    epsilon = float(model["policy"]["score_clip_epsilon"])
    alpha = float(model["alpha"])
    slope = float(model["positive_slope"])
    if not slope > 0:
        raise ValueError("Calibration model slope must be positive")
    bin_effects = dict(model["bin_effects"])
    term_effects = dict(model["term_effects"])
    intercept_by_term = np.asarray(
        [
            float(bin_effects.get(term_bin, 0.0))
            + float(term_effects.get(term_id, 0.0))
            for term_id, term_bin in zip(term_ids, term_bins)
        ],
        dtype=np.float64,
    )
    calibrated = np.empty(scores.shape, dtype=np.float32)
    for start in range(0, scores.shape[0], protein_chunk_size):
        stop = min(scores.shape[0], start + protein_chunk_size)
        x = _logit_scores(scores[start:stop], epsilon)
        calibrated[start:stop] = expit(
            alpha + slope * x + intercept_by_term[np.newaxis, :]
        ).astype(np.float32)
    if (
        not np.isfinite(calibrated).all()
        or np.any(calibrated < 0)
        or np.any(calibrated > 1)
    ):
        raise RuntimeError("Calibrated values are invalid")
    return calibrated, list(model["fallback_by_term"])


def calibration_metrics(
    probabilities: np.ndarray, truth: np.ndarray
) -> dict[str, float]:
    if probabilities.shape != truth.shape:
        raise ValueError("Calibration metric arrays differ")
    clipped = np.clip(probabilities.astype(np.float64), 1e-12, 1.0 - 1e-12)
    target = truth.astype(np.float64, copy=False)
    brier = float(np.mean((clipped - target) ** 2))
    log_loss = float(
        -np.mean(target * np.log(clipped) + (1.0 - target) * np.log1p(-clipped))
    )
    return {
        "brier": brier,
        "binary_log_loss": log_loss,
        "mean_probability": float(clipped.mean()),
        "observed_prevalence": float(target.mean()),
    }


def calibration_intercept_slope(
    probabilities: np.ndarray,
    truth: np.ndarray,
    *,
    protein_chunk_size: int = 256,
) -> dict[str, float | None]:
    if probabilities.shape != truth.shape or probabilities.ndim != 2:
        raise ValueError("Calibration slope arrays differ")
    positives = int(truth.sum())
    events = truth.size
    if not positives or positives == events:
        return {"calibration_intercept": None, "calibration_slope": None}

    def objective(parameters: np.ndarray) -> tuple[float, np.ndarray]:
        loss = 0.0
        gradient = np.zeros(2, dtype=np.float64)
        for start in range(0, probabilities.shape[0], protein_chunk_size):
            stop = min(probabilities.shape[0], start + protein_chunk_size)
            x = _logit_scores(probabilities[start:stop], 1e-12)
            y = truth[start:stop]
            eta = parameters[0] + parameters[1] * x
            loss += float(np.sum(np.logaddexp(0.0, eta) - y * eta))
            residual = expit(eta) - y
            gradient[0] += float(residual.sum())
            gradient[1] += float(np.sum(residual * x))
        return loss, gradient

    result = minimize(
        objective,
        np.asarray([0.0, 1.0]),
        method="L-BFGS-B",
        jac=True,
        bounds=[(-30.0, 30.0), (-20.0, 20.0)],
        options={"maxiter": 100, "ftol": 1e-9, "gtol": 1e-8},
    )
    if not result.success:
        return {"calibration_intercept": None, "calibration_slope": None}
    return {
        "calibration_intercept": float(result.x[0]),
        "calibration_slope": float(result.x[1]),
    }


def reliability_rows(
    probabilities: np.ndarray,
    truth: np.ndarray,
    bin_count: int,
) -> list[dict[str, Any]]:
    if probabilities.shape != truth.shape:
        raise ValueError("Reliability arrays differ")
    if bin_count < 2:
        raise ValueError("Reliability bin count must be at least two")
    flattened = probabilities.astype(np.float64, copy=False).reshape(-1)
    target = truth.reshape(-1)
    order = np.argsort(flattened, kind="stable")
    groups = np.array_split(order, bin_count)
    rows = []
    for index, group in enumerate(groups, start=1):
        if not len(group):
            rows.append(
                {
                    "reliability_bin": index,
                    "events": 0,
                    "mean_probability": None,
                    "observed_fraction": None,
                    "minimum_probability": None,
                    "maximum_probability": None,
                }
            )
            continue
        values = flattened[group]
        labels = target[group]
        rows.append(
            {
                "reliability_bin": index,
                "events": int(len(group)),
                "mean_probability": float(values.mean()),
                "observed_fraction": float(labels.mean()),
                "minimum_probability": float(values.min()),
                "maximum_probability": float(values.max()),
            }
        )
    return rows


def expected_calibration_error(rows: Sequence[Mapping[str, Any]]) -> float:
    total = sum(int(row["events"]) for row in rows)
    if not total:
        return float("nan")
    return float(
        sum(
            int(row["events"])
            * abs(float(row["mean_probability"]) - float(row["observed_fraction"]))
            for row in rows
            if row["events"]
        )
        / total
    )
