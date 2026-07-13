"""Shared filesystem-containment checks for immutable inputs."""

from pathlib import Path


class PathSafetyError(ValueError):
    pass


def resolve_within(root: Path, relative: Path, label: str) -> Path:
    root_resolved = root.resolve()
    candidate = (root_resolved / relative).resolve()
    try:
        candidate.relative_to(root_resolved)
    except ValueError as exc:
        raise PathSafetyError("%s escapes its configured root: %s" % (label, candidate)) from exc
    return candidate


def require_resolved_within(root: Path, candidate: Path, label: str) -> Path:
    root_resolved = root.resolve()
    resolved = candidate.resolve()
    try:
        resolved.relative_to(root_resolved)
    except ValueError as exc:
        raise PathSafetyError("%s escapes its configured root: %s" % (label, resolved)) from exc
    return resolved
