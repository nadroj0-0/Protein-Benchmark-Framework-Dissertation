#!/usr/bin/env bash
# Validate the published outputs from one multi-linkage benchmark campaign.

#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -l tmem=2G
#$ -l tscratch=1G
#$ -l h_rt=1:0:0
#$ -pe smp 1
#$ -N ml_close
#$ -notify

set -Eeuo pipefail

required_environment=(
  CAMPAIGN_ROOT
  BENCHMARK_DIR
  EMBEDDING_EVIDENCE_DIR
  MODEL_RUN_DIR
  DIAGNOSTICS_DIR
  PREDICTION_CAPTURE_DIR
  CALIBRATION_DIR
  SUBMISSION_LEDGER
  BENCHMARK_ID
  FRAMEWORK_COMMIT
)
for name in "${required_environment[@]}"; do
  [[ -n "${!name:-}" ]] || {
    echo "Missing required environment variable: $name" >&2
    exit 2
  }
done

[[ "$CAMPAIGN_ROOT" == /SAN/* ]] || {
  echo "CAMPAIGN_ROOT must be on SAN" >&2
  exit 2
}
[[ "$FRAMEWORK_COMMIT" =~ ^[0-9a-f]{40}$ ]] || {
  echo "FRAMEWORK_COMMIT must be one lowercase 40-character commit" >&2
  exit 2
}

PYTHON_BIN="$(command -v python3 || command -v python || true)"
[[ -n "$PYTHON_BIN" ]] || {
  echo "Python 3 is required for campaign validation" >&2
  exit 2
}

"$PYTHON_BIN" - \
  "$CAMPAIGN_ROOT" \
  "$BENCHMARK_DIR" \
  "$EMBEDDING_EVIDENCE_DIR" \
  "$MODEL_RUN_DIR" \
  "$DIAGNOSTICS_DIR" \
  "$PREDICTION_CAPTURE_DIR" \
  "$CALIBRATION_DIR" \
  "$SUBMISSION_LEDGER" \
  "$BENCHMARK_ID" \
  "$FRAMEWORK_COMMIT" <<'PY'
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path
import sys


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


campaign_root = Path(sys.argv[1]).resolve()
benchmark_id = sys.argv[9]
framework_commit = sys.argv[10]
ledger = Path(sys.argv[8]).resolve()
specifications = (
    ("benchmark", Path(sys.argv[2]).resolve(), "RUN_COMPLETE.json"),
    ("embedding_evidence", Path(sys.argv[3]).resolve(), "RUN_COMPLETE.json"),
    ("model", Path(sys.argv[4]).resolve(), "WORKFLOW_COMPLETE.json"),
    ("diagnostics", Path(sys.argv[5]).resolve(), "WORKFLOW_COMPLETE.json"),
    ("prediction_capture", Path(sys.argv[6]).resolve(), "WORKFLOW_COMPLETE.json"),
    ("calibration", Path(sys.argv[7]).resolve(), "WORKFLOW_COMPLETE.json"),
)

if not campaign_root.is_dir():
    raise SystemExit(f"Campaign root is missing: {campaign_root}")
if not ledger.is_file():
    raise SystemExit(f"Submission ledger is missing: {ledger}")

validated = []
for label, root, marker_name in specifications:
    try:
        root.relative_to(campaign_root)
    except ValueError as error:
        raise SystemExit(f"{label} output is outside the campaign root: {root}") from error
    marker = root / marker_name
    manifest = root / "output_manifest.json"
    if not marker.is_file() or not manifest.is_file():
        raise SystemExit(f"{label} output lacks its marker or manifest: {root}")
    marker_payload = json.loads(marker.read_text(encoding="utf-8"))
    json.loads(manifest.read_text(encoding="utf-8"))
    if marker_payload.get("complete") is not True:
        raise SystemExit(f"{label} completion marker does not declare complete=true")
    validated.append({
        "stage": label,
        "root": str(root),
        "marker": marker_name,
        "marker_sha256": sha256(marker),
        "manifest_sha256": sha256(manifest),
    })

completion = campaign_root / "RUN_COMPLETE.json"
if completion.exists():
    raise SystemExit(f"Campaign completion marker already exists: {completion}")
payload = {
    "complete": True,
    "schema_name": "multilinkage-campaign-completion",
    "schema_version": 1,
    "completed_at_utc": datetime.now(timezone.utc).isoformat(),
    "benchmark_id": benchmark_id,
    "framework_commit": framework_commit,
    "submission_ledger": str(ledger),
    "submission_ledger_sha256": sha256(ledger),
    "validated_outputs": validated,
}
temporary = campaign_root / f".RUN_COMPLETE.json.{os.getpid()}.tmp"
temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
os.replace(temporary, completion)
print(f"Validated complete multi-linkage campaign: {campaign_root}")
PY
