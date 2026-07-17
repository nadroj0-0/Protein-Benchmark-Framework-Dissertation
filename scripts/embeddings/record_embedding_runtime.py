#!/usr/bin/env python3
"""Record hardware and numerical-runtime provenance for an embedding run."""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import os
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def command_output(command: list[str]) -> dict:
    try:
        result = subprocess.run(command, capture_output=True, text=True, check=False)
    except FileNotFoundError:
        return {"available": False, "exit_status": None, "stdout": "", "stderr": ""}
    return {
        "available": True,
        "exit_status": result.returncode,
        "stdout": result.stdout.strip(),
        "stderr": result.stderr.strip(),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--source-file", action="append", default=[])
    args = parser.parse_args()

    packages = {}
    for name in ("torch", "numpy", "transformers", "fair-esm", "biotite"):
        try:
            packages[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            packages[name] = None

    torch_report = {"imported": False}
    try:
        import torch

        torch_report = {
            "imported": True,
            "version": torch.__version__,
            "cuda_version": torch.version.cuda,
            "cuda_available": torch.cuda.is_available(),
            "cudnn_version": torch.backends.cudnn.version(),
            "deterministic_algorithms": torch.are_deterministic_algorithms_enabled(),
            "cudnn_deterministic": torch.backends.cudnn.deterministic,
            "cudnn_benchmark": torch.backends.cudnn.benchmark,
            "cuda_matmul_allow_tf32": torch.backends.cuda.matmul.allow_tf32,
            "cudnn_allow_tf32": torch.backends.cudnn.allow_tf32,
            "device_count": torch.cuda.device_count(),
            "devices": [
                {
                    "index": index,
                    "name": torch.cuda.get_device_name(index),
                    "capability": list(torch.cuda.get_device_capability(index)),
                }
                for index in range(torch.cuda.device_count())
            ],
        }
    except Exception as error:
        torch_report = {
            "imported": False,
            "error": f"{type(error).__name__}:{error}",
        }

    sources = []
    for specification in args.source_file:
        if "=" not in specification:
            raise SystemExit(f"Expected LABEL=PATH for --source-file: {specification}")
        label, raw_path = specification.split("=", 1)
        path = Path(raw_path).resolve()
        if not label or not path.is_file():
            raise SystemExit(f"Missing runtime source {label}: {path}")
        sources.append(
            {
                "label": label,
                "path": str(path),
                "size_bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        )

    report = {
        "schema_version": 1,
        "recorded_at": datetime.now(timezone.utc).isoformat(),
        "hostname": platform.node(),
        "platform": platform.platform(),
        "python": sys.version,
        "executable": sys.executable,
        "cuda_visible_devices": os.environ.get("CUDA_VISIBLE_DEVICES"),
        "sge": {
            key: os.environ.get(key)
            for key in ("JOB_ID", "JOB_NAME", "QUEUE", "HOSTNAME", "SGE_O_WORKDIR")
        },
        "packages": packages,
        "torch": torch_report,
        "nvidia_smi_query": command_output(
            [
                "nvidia-smi",
                "--query-gpu=name,uuid,driver_version,vbios_version",
                "--format=csv,noheader",
            ]
        ),
        "uname": command_output(["uname", "-a"]),
        "sources": sources,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
