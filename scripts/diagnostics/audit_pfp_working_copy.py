#!/usr/bin/env python3
"""Compare a local PFP working directory with a fresh public clone, read-only."""

from __future__ import annotations

import argparse
import difflib
import hashlib
import importlib.metadata
import os
import platform
import subprocess
import sys
import tempfile
from pathlib import Path


DEFAULT_REPOSITORY = "https://github.com/psipred/PFP.git"
KEY_PACKAGES = (
    "torch",
    "numpy",
    "pandas",
    "scipy",
    "scikit-learn",
    "transformers",
    "biopython",
    "biotite",
    "fair-esm",
    "h5py",
)


def run_git(repo: Path, *args: str, check: bool = True) -> str:
    env = os.environ.copy()
    env["GIT_OPTIONAL_LOCKS"] = "0"
    result = subprocess.run(
        ["git", "-C", str(repo), *args],
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        env=env,
    )
    if check and result.returncode:
        raise RuntimeError(result.stderr.strip() or "git command failed")
    return result.stdout.rstrip("\n")


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def tracked_files(repo: Path) -> set[str]:
    output = run_git(repo, "ls-files", "-z")
    return {item for item in output.split("\0") if item}


def is_text(path: Path) -> bool:
    try:
        path.read_text(encoding="utf-8")
    except (UnicodeDecodeError, OSError):
        return False
    return True


def render_items(items: list[str], empty: str = "None.") -> list[str]:
    return [f"- `{item}`" for item in items] if items else [empty]


def package_versions() -> list[str]:
    rows = []
    for name in KEY_PACKAGES:
        try:
            observed = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            observed = "not installed"
        rows.append(f"| `{name}` | `{observed}` |")
    return rows


def pip_freeze() -> tuple[list[str], str | None]:
    result = subprocess.run(
        [sys.executable, "-m", "pip", "freeze", "--all"],
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if result.returncode:
        return [], result.stderr.strip() or "pip freeze failed"
    return sorted(line for line in result.stdout.splitlines() if line), None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pfp_path", type=Path, help="Zijian's local PFP working directory")
    parser.add_argument("--public-repository", default=DEFAULT_REPOSITORY)
    parser.add_argument("--public-ref", default="main")
    parser.add_argument(
        "--max-diff-lines",
        type=int,
        default=500,
        help="Maximum unified-diff lines per modified text file; 0 means unlimited",
    )
    parser.add_argument(
        "--skip-environment",
        action="store_true",
        help="Omit active-Python package information (mainly useful for tests)",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    local = args.pfp_path.expanduser().resolve()
    if not local.is_dir():
        raise SystemExit(f"PFP working directory does not exist: {local}")
    if not (local / ".git").exists():
        raise SystemExit(f"PFP working directory is not a Git checkout: {local}")
    if args.max_diff_lines < 0:
        raise SystemExit("--max-diff-lines cannot be negative")

    with tempfile.TemporaryDirectory(prefix="pfp-public-audit-") as temp_name:
        public = Path(temp_name) / "public"
        clone = subprocess.run(
            [
                "git",
                "clone",
                "--quiet",
                "--depth",
                "1",
                "--branch",
                args.public_ref,
                args.public_repository,
                str(public),
            ],
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        if clone.returncode:
            raise SystemExit(f"Could not clone public PFP repository: {clone.stderr.strip()}")

        public_files = tracked_files(public)
        local_files = tracked_files(local)
        missing: list[str] = []
        modified: list[str] = []
        identical = 0

        for relative in sorted(public_files):
            public_path = public / relative
            local_path = local / relative
            if not local_path.is_file():
                missing.append(relative)
            elif digest(public_path) != digest(local_path):
                modified.append(relative)
            else:
                identical += 1

        local_tracked_only = sorted(local_files - public_files)
        status_lines = run_git(
            local, "status", "--short", "--untracked-files=normal", check=False
        ).splitlines()
        local_untracked = sorted(
            line[3:] for line in status_lines if line.startswith("?? ")
        )

        public_head = run_git(public, "rev-parse", "HEAD")
        local_head = run_git(local, "rev-parse", "HEAD")
        local_branch = run_git(local, "branch", "--show-current") or "detached HEAD"
        local_remote = run_git(local, "remote", "get-url", "origin", check=False) or "not configured"

        materially_different = bool(missing or modified or local_tracked_only)
        verdict = (
            "The local working copy differs from the fresh public release."
            if materially_different
            else "Every public tracked file is byte-identical in the local working copy."
        )

        report = [
            "# PFP Local Working-Copy Audit",
            "",
            "## Safety",
            "",
            "- The supplied PFP directory was opened read-only.",
            "- Git optional locks were disabled for local inspection.",
            "- The public repository was cloned only into an automatically deleted temporary directory.",
            "- Nothing was uploaded and no local-only file contents were included.",
            "- This report was printed to standard output; the script did not create a report file.",
            "",
            "## Identity",
            "",
            f"- Local path: `{local}`",
            f"- Local branch: `{local_branch}`",
            f"- Local HEAD: `{local_head}`",
            f"- Local origin: `{local_remote}`",
            f"- Public repository: `{args.public_repository}`",
            f"- Public ref: `{args.public_ref}`",
            f"- Public HEAD: `{public_head}`",
            "",
            "## Verdict",
            "",
            verdict,
            "",
            "Local-only files are reported separately and are not automatically treated as files that should have been published.",
            "",
            "## File Summary",
            "",
            "| Category | Count |",
            "|---|---:|",
            f"| Public files identical locally | {identical} |",
            f"| Public files modified locally | {len(modified)} |",
            f"| Public files missing locally | {len(missing)} |",
            f"| Local Git-tracked files absent publicly | {len(local_tracked_only)} |",
            f"| Local untracked entries | {len(local_untracked)} |",
            "",
            "## Modified Public Files",
            "",
            *render_items(modified),
            "",
            "## Missing Public Files",
            "",
            *render_items(missing),
            "",
            "## Local Git-Tracked Files Absent From The Public Release",
            "",
            "These may be private development files; their presence alone does not prove they were required for the paper.",
            "",
            *render_items(local_tracked_only),
            "",
            "## Local Untracked Entries",
            "",
            "Only names are reported. Contents are neither read for comparison nor printed.",
            "",
            *render_items(local_untracked),
        ]

        if modified:
            report.extend(["", "## Unified Diffs For Modified Public Text Files"])
            for relative in modified:
                public_path = public / relative
                local_path = local / relative
                report.extend(["", f"### `{relative}`", ""])
                if not (is_text(public_path) and is_text(local_path)):
                    report.extend(
                        [
                            "Binary or non-UTF-8 file; content omitted.",
                            "",
                            f"- Public SHA-256: `{digest(public_path)}`",
                            f"- Local SHA-256: `{digest(local_path)}`",
                        ]
                    )
                    continue

                diff = list(
                    difflib.unified_diff(
                        public_path.read_text(encoding="utf-8").splitlines(),
                        local_path.read_text(encoding="utf-8").splitlines(),
                        fromfile=f"public/{relative}",
                        tofile=f"local/{relative}",
                        lineterm="",
                    )
                )
                truncated = args.max_diff_lines and len(diff) > args.max_diff_lines
                if truncated:
                    diff = diff[: args.max_diff_lines]
                report.extend(["```diff", *diff, "```"])
                if truncated:
                    report.append(
                        f"Diff truncated after {args.max_diff_lines} lines; rerun with `--max-diff-lines 0` for the complete diff."
                    )

        if not args.skip_environment:
            frozen_packages, freeze_error = pip_freeze()
            report.extend(
                [
                    "",
                    "## Active Python Environment",
                    "",
                    "This describes the environment used to run this audit. It is only evidence of the paper environment if that environment was activated first.",
                    "",
                    f"- Python executable: `{sys.executable}`",
                    f"- Python version: `{platform.python_version()}`",
                    f"- Platform: `{platform.platform()}`",
                    f"- `CONDA_PREFIX`: `{os.environ.get('CONDA_PREFIX', 'not set')}`",
                    f"- `VIRTUAL_ENV`: `{os.environ.get('VIRTUAL_ENV', 'not set')}`",
                    "",
                    "| Package | Version |",
                    "|---|---|",
                    *package_versions(),
                    "",
                    "### Complete `pip freeze --all`",
                    "",
                ]
            )
            if freeze_error:
                report.append(f"Could not capture `pip freeze`: `{freeze_error}`")
            else:
                report.extend(["```text", *frozen_packages, "```"])

        print("\n".join(report))
    return 1 if materially_different else 0


if __name__ == "__main__":
    raise SystemExit(main())
