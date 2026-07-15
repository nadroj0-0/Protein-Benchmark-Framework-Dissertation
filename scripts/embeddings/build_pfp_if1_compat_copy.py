#!/usr/bin/env python3
"""Create a provenance-recorded CUDA-safe copy of PFP's IF1 extractor."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path


IMPORT_OLD = '''    except ImportError as e:
        print("ERROR: Could not import `esm` (fair-esm).")
        print("  Try: pip install fair-esm")
        print(f"  Original error: {e}")
        return
'''
IMPORT_NEW = '''    except ImportError as e:
        raise RuntimeError(
            f"Could not import ESM-IF1 dependencies: {e}"
        ) from e
'''

ENCODER_OLD = '''                # Get encoder output (L x 512)
                rep = esm.inverse_folding.util.get_encoder_output(
                    model, alphabet, coords
                )
'''
ENCODER_NEW = '''                # fair-esm 2.0.0 get_encoder_output() creates CPU tensors.
                # Build the same encoder inputs explicitly on the model device.
                batch_converter = esm.inverse_folding.util.CoordBatchConverter(alphabet)
                batch = [(coords, None, None)]
                (
                    encoder_coords,
                    confidence,
                    _,
                    _,
                    padding_mask,
                ) = batch_converter(batch, device=device)
                encoder_out = model.encoder.forward(
                    encoder_coords,
                    padding_mask,
                    confidence,
                    return_all_hiddens=False,
                )
                rep = encoder_out["encoder_out"][0][1:-1, 0]
'''

SUMMARY_OLD = '''    if failed:
        print("\\nExamples of failed proteins:")
        for pid, err in failed[:10]:
            print(f"  {pid}: {err}")
'''
SUMMARY_NEW = '''    if failed:
        print("\\nExamples of failed proteins:")
        for pid, err in failed[:10]:
            print(f"  {pid}: {err}")

    if pdb_files and len(failed) == len(pdb_files):
        raise RuntimeError(f"ESM-IF1 failed for all {len(pdb_files)} PDB files")
'''


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    return parser.parse_args()


def replace_once(source: str, old: str, new: str, label: str) -> str:
    occurrences = source.count(old)
    if occurrences != 1:
        raise SystemExit(
            f"Expected exactly one validated {label} block, found {occurrences}"
        )
    return source.replace(old, new)


def main() -> int:
    args = parse_args()
    if not args.source.is_file():
        raise SystemExit(f"Missing PFP IF1 extractor: {args.source}")
    if args.output.exists():
        raise SystemExit(f"Refusing to overwrite compatibility copy: {args.output}")

    source_text = args.source.read_text(encoding="utf-8")
    output_text = replace_once(
        source_text, IMPORT_OLD, IMPORT_NEW, "IF1 dependency-import"
    )
    output_text = replace_once(
        output_text, ENCODER_OLD, ENCODER_NEW, "IF1 encoder-output"
    )
    output_text = replace_once(
        output_text, SUMMARY_OLD, SUMMARY_NEW, "IF1 failure-summary"
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(output_text, encoding="utf-8")
    shutil.copymode(args.source, args.output)

    report = {
        "schema_version": 1,
        "source": str(args.source.resolve()),
        "source_sha256": digest(args.source),
        "compatibility_copy": str(args.output.resolve()),
        "compatibility_copy_sha256": digest(args.output),
        "changes": [
            "Raise a non-zero error when ESM-IF1 dependencies cannot import.",
            "Create fair-esm encoder tensors explicitly on the model device.",
            "Raise a non-zero error when every PDB extraction fails.",
        ],
        "scientific_output_change": False,
        "upstream_source_modified": False,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
