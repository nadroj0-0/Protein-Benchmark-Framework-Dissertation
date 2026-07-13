#!/usr/bin/env python3
"""Repository-local entry point for the PFP embedding inventory package."""

import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
PACKAGE_SRC = REPO_ROOT / "embedding_inventory" / "src"
sys.path.insert(0, str(PACKAGE_SRC))

from pfp_embedding_inventory.cli import main  # noqa: E402


if __name__ == "__main__":
    raise SystemExit(main())
