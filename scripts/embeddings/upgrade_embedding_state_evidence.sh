#!/usr/bin/env bash
# Add immutable per-array evidence hashes to an existing embedding retry state.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_BIN="${PYTHON_BIN:-$(command -v python3 || command -v python || true)}"
STATE_ROOT=""
OUTPUT_DIR=""
RETRIES_FINISHED=0

usage() {
  cat <<'EOF'
Usage: bash scripts/embeddings/upgrade_embedding_state_evidence.sh \
  --state-root PATH --output-dir PATH --confirm-retries-finished

Run this only after every retry job targeting the state has finished. The state
lock prevents simultaneous merge writes, but it cannot detect a generation job
that has not reached its merge step yet. The command adds hashes and refreshes
reports without rebuilding the contract or changing accepted membership.
EOF
}

die() { echo "ERROR: $*" >&2; exit 2; }
require_value() { [[ $# -ge 2 && -n "$2" ]] || die "$1 requires a value"; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --state-root) require_value "$@"; STATE_ROOT="$2"; shift 2 ;;
    --output-dir) require_value "$@"; OUTPUT_DIR="$2"; shift 2 ;;
    --confirm-retries-finished) RETRIES_FINISHED=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) usage >&2; die "Unknown argument: $1" ;;
  esac
done

[[ -n "$PYTHON_BIN" ]] || die "Python is required"
[[ -d "$STATE_ROOT" ]] || die "State root does not exist: $STATE_ROOT"
[[ -f "$STATE_ROOT/contract.json" ]] || die "State contract is missing"
[[ -n "$OUTPUT_DIR" ]] || die "--output-dir is required"
[[ "$RETRIES_FINISHED" == "1" ]] || \
  die "Refusing evidence upgrade without --confirm-retries-finished"
[[ ! -e "$OUTPUT_DIR" ]] || die "Output directory already exists: $OUTPUT_DIR"

mkdir -p "$OUTPUT_DIR"
"$PYTHON_BIN" "$HERE/manage_resumable_embedding_state.py" \
  upgrade-evidence-hashes \
  --state-root "$STATE_ROOT" \
  --report "$OUTPUT_DIR/evidence_hash_upgrade.json" \
  > "$OUTPUT_DIR/command_output.json"

cp -p \
  "$STATE_ROOT/coverage.json" \
  "$STATE_ROOT/EVIDENCE_HASHES_COMPLETE.json" \
  "$OUTPUT_DIR/"

printf '{"complete":true,"state_root":"%s"}\n' "$STATE_ROOT" \
  > "$OUTPUT_DIR/EVIDENCE_UPGRADE_COMPLETE.json"
echo "Embedding evidence hashes upgraded: $STATE_ROOT"
