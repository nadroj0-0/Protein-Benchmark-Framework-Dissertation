#!/usr/bin/env bash
# Compare a target nine-CSV PFP benchmark with one or more benchmarks whose
# proteins have already been processed for embeddings. Heavy staging stays in
# the caller-provided work directory; only planner reports are published.

set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
FRAMEWORK_ROOT="$(cd "${HERE}/../.." && pwd)"

REQUIRED_CSVS=(
  bp-training.csv bp-validation.csv bp-test.csv
  cc-training.csv cc-validation.csv cc-test.csv
  mf-training.csv mf-validation.csv mf-test.csv
)

TARGET_SPEC=""
EMBEDDED_SPECS=()
WORK_DIR=""
OUTPUT_DIR=""
PYTHON_BIN="${PYTHON_BIN:-python3}"
CAFA3_BASE_URL="${CAFA3_BASE_URL:-https://zenodo.org/records/7409660/files}"

usage() {
  cat <<'EOF'
Usage:
  run_benchmark_reuse_plan.sh \
    --target-benchmark NAME=PATH \
    --work-dir PATH \
    --output-dir PATH \
    [--embedded-benchmark NAME=PATH ...] \
    [--python-bin PATH]

If no --embedded-benchmark is supplied, the canonical CAFA3 CSVs represented
by Zijian's published embedding run are downloaded from Zenodo record 7409660,
authenticated, and used under the name cafa3_zijian.

This workflow compares benchmark CSVs only. It does not inspect, download, or
generate embedding arrays.
EOF
}

die() {
  echo "ERROR: $*" >&2
  exit 2
}

split_named_path() {
  local value="$1"
  local name="${value%%=*}"
  local path="${value#*=}"
  [[ "$value" == *=* && -n "$name" && -n "$path" ]] || \
    die "Benchmark arguments must use non-empty NAME=PATH: $value"
  [[ "$name" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || \
    die "Benchmark names must match [A-Za-z0-9][A-Za-z0-9._-]*: $name"
  printf '%s\t%s\n' "$name" "$path"
}

file_size() {
  local path="$1"
  if stat -c '%s' "$path" >/dev/null 2>&1; then
    stat -c '%s' "$path"
  else
    stat -f '%z' "$path"
  fi
}

download_file() {
  local url="$1"
  local destination="$2"
  local partial="${destination}.part"
  echo "Downloading: $url"
  if command -v wget >/dev/null 2>&1; then
    wget --tries=5 --timeout=60 -c "$url" -O "$partial"
  elif command -v curl >/dev/null 2>&1; then
    curl --fail --location --retry 5 --continue-at - --output "$partial" "$url"
  else
    die "Neither wget nor curl is available"
  fi
  mv "$partial" "$destination"
}

stage_benchmark() {
  local source_dir="$1"
  local destination_dir="$2"
  local role="$3"
  local name
  [[ -d "$source_dir" ]] || die "$role benchmark directory does not exist: $source_dir"
  mkdir -p "$destination_dir"
  for name in "${REQUIRED_CSVS[@]}"; do
    [[ -f "${source_dir}/${name}" ]] || die "$role benchmark is missing $name"
    cp -p "${source_dir}/${name}" "${destination_dir}/${name}"
    printf '%s\t%s\t%s\t%s\t%s\n' \
      "$role" "$name" "${source_dir}/${name}" "${destination_dir}/${name}" \
      "$(file_size "${destination_dir}/${name}")" >> "$ACQUISITION_LOG"
  done
}

authenticate_canonical_cafa3() {
  local source_dir="$1"
  "$PYTHON_BIN" - "$source_dir" "$ACQUISITION_LOG" "$CAFA3_BASE_URL" <<'PY'
import hashlib
import os
import shutil
import sys
from pathlib import Path

source_dir = Path(sys.argv[1])
log_path = Path(sys.argv[2])
origin = sys.argv[3]
expected = {
    "bp-test.csv": "e9a4b239cd47a7ac80975f63e259581e",
    "bp-training.csv": "85c19594547a503956226b9c225efc5d",
    "bp-validation.csv": "c2674223770d6a8cf680dd9335d51ebe",
    "cc-test.csv": "0e5dc8528ca95e8897b10cddaa12a775",
    "cc-training.csv": "074b13dd50fad4a6a4f13e4d8d4105d6",
    "cc-validation.csv": "cdc8ceefcab4fb8c9278dd07c184327f",
    "mf-test.csv": "2735e408dd57f6de29b1538f6b150d68",
    "mf-training.csv": "b31a8f22b5934aef61b76ec3b89296da",
    "mf-validation.csv": "897921ce5df8174672200320926ccc87",
}


def digest(path, replacement=None):
    value = hashlib.md5()
    with path.open("rb") as handle:
        first = handle.readline()
        if replacement is not None:
            prefix, new_prefix = replacement
            if not first.startswith(prefix):
                return ""
            first = new_prefix + first[len(prefix):]
        value.update(first)
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


with log_path.open("a", encoding="utf-8") as log:
    for name, wanted in expected.items():
        path = source_dir / name
        observed = digest(path)
        normalized = False
        if observed != wanted and name.startswith("mf-"):
            observed = digest(path, (b"proteins,", b"protein,"))
            normalized = observed == wanted
        if observed != wanted:
            raise SystemExit(
                "Canonical source CSV authentication failed for %s: %s != %s"
                % (name, observed or "unrecognized-header", wanted)
            )

        with path.open("rb") as handle:
            first = handle.readline()
        if first.startswith(b"protein,"):
            temporary = path.with_suffix(path.suffix + ".normalizing")
            with path.open("rb") as source, temporary.open("wb") as destination:
                header = source.readline()
                destination.write(b"proteins," + header[len(b"protein,"):])
                shutil.copyfileobj(source, destination, 1024 * 1024)
            os.replace(str(temporary), str(path))
            normalized = True

        status = "canonical-verified-normalized" if normalized else "canonical-verified"
        log.write(
            "embedded-default\t%s\t%s/%s\t%s\t%s\n"
            % (name, origin.rstrip("/"), name, path, status)
        )
PY
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --target-benchmark)
      [[ $# -ge 2 ]] || die "--target-benchmark requires NAME=PATH"
      [[ -z "$TARGET_SPEC" ]] || die "Exactly one --target-benchmark is accepted"
      TARGET_SPEC="$2"
      shift 2
      ;;
    --embedded-benchmark)
      [[ $# -ge 2 ]] || die "--embedded-benchmark requires NAME=PATH"
      EMBEDDED_SPECS+=("$2")
      shift 2
      ;;
    --work-dir)
      [[ $# -ge 2 ]] || die "--work-dir requires a path"
      WORK_DIR="$2"
      shift 2
      ;;
    --output-dir)
      [[ $# -ge 2 ]] || die "--output-dir requires a path"
      OUTPUT_DIR="$2"
      shift 2
      ;;
    --python-bin)
      [[ $# -ge 2 ]] || die "--python-bin requires a path"
      PYTHON_BIN="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      die "Unknown argument: $1"
      ;;
  esac
done

[[ -n "$TARGET_SPEC" ]] || die "--target-benchmark is required"
[[ -n "$WORK_DIR" ]] || die "--work-dir is required"
[[ -n "$OUTPUT_DIR" ]] || die "--output-dir is required"
command -v "$PYTHON_BIN" >/dev/null 2>&1 || die "Python executable not found: $PYTHON_BIN"
[[ ! -e "$WORK_DIR" ]] || die "Work directory already exists: $WORK_DIR"
[[ ! -e "$OUTPUT_DIR" ]] || die "Output directory already exists: $OUTPUT_DIR"

mkdir -p "$WORK_DIR" "$OUTPUT_DIR"
WORK_DIR="$(cd "$WORK_DIR" && pwd)"
OUTPUT_DIR="$(cd "$OUTPUT_DIR" && pwd)"
ACQUISITION_LOG="${OUTPUT_DIR}/input_acquisition.tsv"
printf 'role\tfile\tsource\tstaged_path\tstatus_or_size\n' > "$ACQUISITION_LOG"

IFS=$'\t' read -r TARGET_NAME TARGET_DIR < <(split_named_path "$TARGET_SPEC")
TARGET_STAGE="${WORK_DIR}/target_${TARGET_NAME}"
stage_benchmark "$TARGET_DIR" "$TARGET_STAGE" target

if [[ "${#EMBEDDED_SPECS[@]}" -eq 0 ]]; then
  DEFAULT_STAGE="${WORK_DIR}/embedded_cafa3_zijian"
  mkdir -p "$DEFAULT_STAGE"
  for name in "${REQUIRED_CSVS[@]}"; do
    download_file "${CAFA3_BASE_URL}/${name}?download=1" "${DEFAULT_STAGE}/${name}"
  done
  authenticate_canonical_cafa3 "$DEFAULT_STAGE"
  EMBEDDED_SPECS+=("cafa3_zijian=${DEFAULT_STAGE}")
fi

PLANNER_ARGS=(plan)
for specification in "${EMBEDDED_SPECS[@]}"; do
  IFS=$'\t' read -r embedded_name embedded_dir < <(split_named_path "$specification")
  embedded_stage="${WORK_DIR}/embedded_${embedded_name}"
  if [[ "$embedded_dir" == "$embedded_stage" ]]; then
    :
  else
    stage_benchmark "$embedded_dir" "$embedded_stage" "embedded:${embedded_name}"
  fi
  PLANNER_ARGS+=(--embedded-benchmark "${embedded_name}=${embedded_stage}")
done
PLANNER_ARGS+=(
  --target-benchmark "${TARGET_NAME}=${TARGET_STAGE}"
  --output-dir "${OUTPUT_DIR}/plan"
)

echo "==> Running exact ID-and-sequence benchmark reuse planner"
printf 'Command:'
printf ' %q' "$PYTHON_BIN" -m pfp_benchmark_reuse "${PLANNER_ARGS[@]}"
printf '\n\n'

PYTHONPATH="${FRAMEWORK_ROOT}/benchmark_reuse_planner/src${PYTHONPATH:+:${PYTHONPATH}}" \
  "$PYTHON_BIN" -m pfp_benchmark_reuse "${PLANNER_ARGS[@]}"

[[ -f "${OUTPUT_DIR}/plan/RUN_COMPLETE.json" ]] || \
  die "Planner did not publish RUN_COMPLETE.json"
[[ -f "${OUTPUT_DIR}/plan/output_manifest.json" ]] || \
  die "Planner did not publish output_manifest.json"

cat > "${OUTPUT_DIR}/WORKFLOW_COMPLETE.json" <<EOF
{
  "complete": true,
  "planner_output": "plan",
  "target_benchmark": "${TARGET_NAME}",
  "schema_version": 1
}
EOF

echo "==> Reuse plan complete: ${OUTPUT_DIR}/plan"
