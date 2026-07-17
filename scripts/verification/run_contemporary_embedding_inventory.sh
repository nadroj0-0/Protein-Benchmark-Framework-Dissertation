#!/usr/bin/env bash
# Stage a nine-CSV benchmark and Zijian's published cache, then run the
# provenance-aware embedding inventory. Heavy artifacts remain in work space.

set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
FRAMEWORK_ROOT="$(cd "${HERE}/../.." && pwd)"
# shellcheck disable=SC1091
source "$FRAMEWORK_ROOT/scripts/artifact_catalog.sh"

REQUIRED_CSVS=(
  bp-training.csv bp-validation.csv bp-test.csv
  cc-training.csv cc-validation.csv cc-test.csv
  mf-training.csv mf-validation.csv mf-test.csv
)
# The empty sentinel keeps older cluster Bash versions happy under `set -u`.
# TARGET_OVERRIDE_COUNT remains the authoritative number of supplied overrides.
TARGET_OVERRIDES=("")
TARGET_OVERRIDE_COUNT=0

BENCHMARK_DIR=""
SOURCE_BENCHMARK_DIR=""
EMBEDDING_ARCHIVE_DIR=""
PFP_REFERENCE_DIR=""
OUTPUT_DIR=""
WORK_DIR=""
ALIASES_FILE=""
CONFIG_PATH="${FRAMEWORK_ROOT}/configs/embedding_inventory.contemporary.json"
POLICY="maximize-coverage"
REPORT_LEVEL="compact"
PYTHON_BIN="${PYTHON_BIN:-python3}"

CAFA3_BASE_URL="${CAFA3_BASE_URL:-https://zenodo.org/records/7409660/files}"
MMFP_BASE_URL="${MMFP_BASE_URL:-https://zenodo.org/records/19498341/files}"
PFP_REFERENCE_URL="${PFP_REFERENCE_URL:-https://github.com/psipred/PFP.git}"
PFP_REFERENCE_COMMIT="1e04fd6d6d3c40458fd41ec1a881ed6e24de768e"

usage() {
  cat <<'EOF'
Usage:
  run_contemporary_embedding_inventory.sh \
    --benchmark-dir PATH \
    --work-dir PATH \
    --output-dir PATH [options]

Required benchmark input:
  --benchmark-dir PATH          Directory containing the nine PFP CSVs.
  --benchmark-csv NAME=PATH     Override one CSV; repeat as needed. With no
                                --benchmark-dir, all nine overrides are required.

Required run paths:
  --work-dir PATH               New/empty workspace for staged CSVs and archives.
  --output-dir PATH             New result root. Inventory is written beneath it.

Optional inputs:
  --source-benchmark-dir PATH   Canonical CAFA3 source CSVs represented by the
                                published cache. Default: catalogue, then Zenodo 7409660.
  --embedding-archive-dir PATH  Directory containing Zijian's three published
                                embedding tarballs. Default: catalogue, then
                                Zenodo 19498341 into the work directory.
  --pfp-reference-dir PATH      Local PFP Git clone containing the pinned
                                reference commit. Default: clone psipred/PFP.
  --aliases PATH                Explicit inventory alias TSV.
  --config PATH                 Inventory config (default: contemporary config).
  --policy NAME                 paper-faithful or maximize-coverage.
  --report-level NAME           compact or full (default: compact).
  --python-bin PATH             Python executable (default: PYTHON_BIN or python3).
  --artifact-catalog PATH       Existing-artifact path map; explicit directories win.
  -h, --help                    Show this help.

This workflow inventories and plans only. It does not generate embeddings.
EOF
}

die() {
  echo "ERROR: $*" >&2
  exit 2
}

warn() {
  echo "WARNING: $*" >&2
}

is_required_csv() {
  local candidate="$1"
  local name
  for name in "${REQUIRED_CSVS[@]}"; do
    if [[ "$candidate" == "$name" ]]; then
      return 0
    fi
  done
  return 1
}

override_for() {
  local requested="$1"
  local entry name
  for entry in "${TARGET_OVERRIDES[@]}"; do
    [[ -n "$entry" ]] || continue
    name="${entry%%=*}"
    if [[ "$name" == "$requested" ]]; then
      printf '%s\n' "${entry#*=}"
      return 0
    fi
  done
  return 1
}

file_size() {
  local path="$1"
  if stat -c '%s' "$path" >/dev/null 2>&1; then
    stat -c '%s' "$path"
  else
    stat -f '%z' "$path"
  fi
}

sha256_file() {
  local path="$1"
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$path" | awk '{print $1}'
  elif command -v shasum >/dev/null 2>&1; then
    shasum -a 256 "$path" | awk '{print $1}'
  else
    "$PYTHON_BIN" - "$path" <<'PY'
import hashlib
import sys

digest = hashlib.sha256()
with open(sys.argv[1], "rb") as handle:
    for block in iter(lambda: handle.read(1024 * 1024), b""):
        digest.update(block)
print(digest.hexdigest())
PY
  fi
}

download_file() {
  local url="$1"
  local destination="$2"
  local partial="${destination}.part"
  mkdir -p "$(dirname "$destination")"
  echo "Downloading: $url"
  if command -v wget >/dev/null 2>&1; then
    wget --tries=5 --timeout=60 -c "$url" -O "$partial"
  elif command -v curl >/dev/null 2>&1; then
    curl --fail --location --retry 5 --continue-at - \
      --output "$partial" "$url"
  else
    die "Neither wget nor curl is available"
  fi
  mv "$partial" "$destination"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --benchmark-dir)
      [[ $# -ge 2 ]] || die "--benchmark-dir requires a path"
      BENCHMARK_DIR="$2"
      shift 2
      ;;
    --benchmark-csv)
      [[ $# -ge 2 ]] || die "--benchmark-csv requires NAME=PATH"
      TARGET_OVERRIDES+=("$2")
      TARGET_OVERRIDE_COUNT=$((TARGET_OVERRIDE_COUNT + 1))
      shift 2
      ;;
    --source-benchmark-dir)
      [[ $# -ge 2 ]] || die "--source-benchmark-dir requires a path"
      SOURCE_BENCHMARK_DIR="$2"
      shift 2
      ;;
    --embedding-archive-dir)
      [[ $# -ge 2 ]] || die "--embedding-archive-dir requires a path"
      EMBEDDING_ARCHIVE_DIR="$2"
      shift 2
      ;;
    --pfp-reference-dir)
      [[ $# -ge 2 ]] || die "--pfp-reference-dir requires a path"
      PFP_REFERENCE_DIR="$2"
      shift 2
      ;;
    --output-dir)
      [[ $# -ge 2 ]] || die "--output-dir requires a path"
      OUTPUT_DIR="$2"
      shift 2
      ;;
    --work-dir)
      [[ $# -ge 2 ]] || die "--work-dir requires a path"
      WORK_DIR="$2"
      shift 2
      ;;
    --aliases)
      [[ $# -ge 2 ]] || die "--aliases requires a path"
      ALIASES_FILE="$2"
      shift 2
      ;;
    --config)
      [[ $# -ge 2 ]] || die "--config requires a path"
      CONFIG_PATH="$2"
      shift 2
      ;;
    --policy)
      [[ $# -ge 2 ]] || die "--policy requires a value"
      POLICY="$2"
      shift 2
      ;;
    --report-level)
      [[ $# -ge 2 ]] || die "--report-level requires a value"
      REPORT_LEVEL="$2"
      shift 2
      ;;
    --python-bin)
      [[ $# -ge 2 ]] || die "--python-bin requires a path"
      PYTHON_BIN="$2"
      shift 2
      ;;
    --artifact-catalog)
      [[ $# -ge 2 ]] || die "--artifact-catalog requires a path"
      ARTIFACT_CATALOG="$2"
      export ARTIFACT_CATALOG
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
artifact_catalog_configure "$FRAMEWORK_ROOT" "${ARTIFACT_CATALOG:-}"

[[ -n "$WORK_DIR" ]] || die "--work-dir is required"
[[ -n "$OUTPUT_DIR" ]] || die "--output-dir is required"
[[ -n "$BENCHMARK_DIR" || "$TARGET_OVERRIDE_COUNT" -gt 0 ]] || \
  die "Supply --benchmark-dir or all nine --benchmark-csv overrides"
[[ "$POLICY" == "paper-faithful" || "$POLICY" == "maximize-coverage" ]] || \
  die "Unsupported policy: $POLICY"
[[ "$REPORT_LEVEL" == "compact" || "$REPORT_LEVEL" == "full" ]] || \
  die "Unsupported report level: $REPORT_LEVEL"
[[ -f "$CONFIG_PATH" ]] || die "Inventory config does not exist: $CONFIG_PATH"
command -v "$PYTHON_BIN" >/dev/null 2>&1 || die "Python executable not found: $PYTHON_BIN"

for entry in "${TARGET_OVERRIDES[@]}"; do
  [[ -n "$entry" ]] || continue
  [[ "$entry" == *=* ]] || die "Invalid --benchmark-csv value: $entry"
  name="${entry%%=*}"
  path="${entry#*=}"
  is_required_csv "$name" || die "Unknown benchmark CSV role: $name"
  [[ -n "$path" ]] || die "Empty path for benchmark CSV: $name"
  occurrences=0
  for comparison in "${TARGET_OVERRIDES[@]}"; do
    [[ "${comparison%%=*}" == "$name" ]] && occurrences=$((occurrences + 1))
  done
  [[ "$occurrences" -eq 1 ]] || die "Duplicate --benchmark-csv override: $name"
done

mkdir -p "$WORK_DIR"
WORK_DIR="$(cd "$WORK_DIR" && pwd)"
OUTPUT_PARENT="$(dirname "$OUTPUT_DIR")"
mkdir -p "$OUTPUT_PARENT"
OUTPUT_PARENT="$(cd "$OUTPUT_PARENT" && pwd)"
OUTPUT_DIR="${OUTPUT_PARENT}/$(basename "$OUTPUT_DIR")"
[[ ! -e "$OUTPUT_DIR" ]] || die "Output path already exists: $OUTPUT_DIR"

TARGET_STAGE="${WORK_DIR}/target_benchmark"
SOURCE_STAGE="${WORK_DIR}/source_benchmark"
ARCHIVE_STAGE="${WORK_DIR}/published_archives"
PUBLISHED_ROOT="${WORK_DIR}/published_cache"
INVENTORY_OUTPUT="${OUTPUT_DIR}/inventory"
ACQUISITION_LOG="${OUTPUT_DIR}/input_acquisition.tsv"

for staging_path in "$TARGET_STAGE" "$SOURCE_STAGE" "$ARCHIVE_STAGE" "$PUBLISHED_ROOT"; do
  [[ ! -e "$staging_path" ]] || die "Work directory contains stale staged data: $staging_path"
done

mkdir -p "$TARGET_STAGE" "$SOURCE_STAGE" "$ARCHIVE_STAGE" "$PUBLISHED_ROOT" "$OUTPUT_DIR"
printf 'kind\titem\tsource\tstaged_path\tchecksum_algorithm\tchecksum\tstatus\n' \
  > "$ACQUISITION_LOG"

echo "==> [1/7] Stage the target nine-CSV benchmark"
for name in "${REQUIRED_CSVS[@]}"; do
  source_path=""
  source_label=""
  if source_path="$(override_for "$name")"; then
    source_label="custom-file"
  elif [[ -n "$BENCHMARK_DIR" ]]; then
    source_path="${BENCHMARK_DIR}/${name}"
    source_label="benchmark-directory"
  else
    die "No input was supplied for $name"
  fi
  [[ -f "$source_path" ]] || die "Target benchmark CSV does not exist: $source_path"
  cp -p "$source_path" "${TARGET_STAGE}/${name}"
  printf 'target-csv\t%s\t%s\t%s\t-\t-\tstaged-%s-size-%s\n' \
    "$name" "$source_path" "${TARGET_STAGE}/${name}" "$source_label" \
    "$(file_size "${TARGET_STAGE}/${name}")" >> "$ACQUISITION_LOG"
done

echo "==> [2/7] Stage and authenticate the canonical cache-source benchmark"
if [[ -n "$SOURCE_BENCHMARK_DIR" && -d "$SOURCE_BENCHMARK_DIR" ]]; then
  for name in "${REQUIRED_CSVS[@]}"; do
    [[ -f "${SOURCE_BENCHMARK_DIR}/${name}" ]] || \
      die "Source benchmark is missing: $name"
    cp -p "${SOURCE_BENCHMARK_DIR}/${name}" "${SOURCE_STAGE}/${name}"
  done
  SOURCE_ORIGIN="$SOURCE_BENCHMARK_DIR"
else
  if [[ -n "$SOURCE_BENCHMARK_DIR" ]]; then
    warn "Source benchmark directory does not exist; falling back to the artifact catalogue or Zenodo: $SOURCE_BENCHMARK_DIR"
  fi
  for name in "${REQUIRED_CSVS[@]}"; do
    source_path="$(resolve_artifact_path "$(canonical_cafa3_artifact_id "$name")" "" || true)"
    if [[ -n "$source_path" ]]; then
      cp -p "$source_path" "${SOURCE_STAGE}/${name}"
    else
      download_file "${CAFA3_BASE_URL}/${name}?download=1" "${SOURCE_STAGE}/${name}"
    fi
  done
  SOURCE_ORIGIN="${ARTIFACT_CATALOG:-$CAFA3_BASE_URL}"
fi

# Accept either pristine Zenodo MF headers (protein) or the sole normalization
# required by the PFP consumer (proteins). Every other byte must authenticate.
"$PYTHON_BIN" - "$SOURCE_STAGE" "$ACQUISITION_LOG" "$SOURCE_ORIGIN" <<'PY'
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
            "source-csv\t%s\t%s/%s\t%s\tmd5\t%s\t%s\n"
            % (name, origin.rstrip("/"), name, path, wanted, status)
        )
PY

echo "==> [3/7] Stage and authenticate Zijian's published embedding archives"
ARCHIVE_NAMES=(
  mmfp_embeddings_prott5.tar.gz
  mmfp_embeddings_struct_ppi.tar.gz
  mmfp_embeddings_text_temporal.tar.gz
)

archive_sha256() {
  case "$1" in
    mmfp_embeddings_prott5.tar.gz)
      printf '%s\n' '30dd88fc4acbe3bc267bd8d5ae05e4d967fa7c169a6f063f12d2395fb0ffb00f'
      ;;
    mmfp_embeddings_struct_ppi.tar.gz)
      printf '%s\n' '6d7a243f2c5e2149c162698b4e6a5e297731a4fea835d57b9049bb31f4af32de'
      ;;
    mmfp_embeddings_text_temporal.tar.gz)
      printf '%s\n' 'df1bf558fab1c018286a5b389665245917ba951f2c6cd8558546d0b1a3b47e36'
      ;;
    *)
      return 1
      ;;
  esac
}

for name in "${ARCHIVE_NAMES[@]}"; do
  destination="${ARCHIVE_STAGE}/${name}"
  if [[ -n "$EMBEDDING_ARCHIVE_DIR" ]]; then
    source_path="$(resolve_artifact_path "$(zijian_embedding_artifact_id "$name")" "${EMBEDDING_ARCHIVE_DIR}/${name}" || true)"
  else
    source_path="$(resolve_artifact_path "$(zijian_embedding_artifact_id "$name")" "" || true)"
  fi
  if [[ -n "$source_path" ]]; then
    cp -p "$source_path" "$destination"
    source_description="$source_path"
  else
    source_description="${MMFP_BASE_URL}/${name}?download=1"
    download_file "$source_description" "$destination"
  fi
  wanted="$(archive_sha256 "$name")"
  observed="$(sha256_file "$destination")"
  [[ "$observed" == "$wanted" ]] || \
    die "Embedding archive SHA-256 mismatch for $name: $observed != $wanted"
  tar -tzf "$destination" >/dev/null
  printf 'embedding-archive\t%s\t%s\t%s\tsha256\t%s\tauthenticated\n' \
    "$name" "$source_description" "$destination" "$observed" >> "$ACQUISITION_LOG"
done

echo "==> [4/7] Extract the authenticated published cache in work space"
if [[ -n "$PFP_REFERENCE_DIR" ]]; then
  [[ -d "${PFP_REFERENCE_DIR}/.git" ]] || \
    die "PFP reference directory is not a Git clone: $PFP_REFERENCE_DIR"
  git clone --no-checkout "$PFP_REFERENCE_DIR" "$PUBLISHED_ROOT"
  PFP_REFERENCE_ORIGIN="$PFP_REFERENCE_DIR"
else
  git clone --no-checkout "$PFP_REFERENCE_URL" "$PUBLISHED_ROOT"
  PFP_REFERENCE_ORIGIN="$PFP_REFERENCE_URL"
fi
git -C "$PUBLISHED_ROOT" checkout --detach "$PFP_REFERENCE_COMMIT"
[[ "$(git -C "$PUBLISHED_ROOT" rev-parse HEAD)" == "$PFP_REFERENCE_COMMIT" ]] || \
  die "PFP reference checkout did not resolve to the pinned commit"
printf 'pfp-reference\tcommit\t%s\t%s\tgit-commit\t%s\tauthenticated\n' \
  "$PFP_REFERENCE_ORIGIN" "$PUBLISHED_ROOT" "$PFP_REFERENCE_COMMIT" >> "$ACQUISITION_LOG"

for name in "${ARCHIVE_NAMES[@]}"; do
  ln "${ARCHIVE_STAGE}/${name}" "${PUBLISHED_ROOT}/${name}" 2>/dev/null || \
    cp -p "${ARCHIVE_STAGE}/${name}" "${PUBLISHED_ROOT}/${name}"
  tar -xzf "${ARCHIVE_STAGE}/${name}" -C "$PUBLISHED_ROOT"
done

CACHE_ROOT="${PUBLISHED_ROOT}/data/embedding_cache"
[[ -d "$CACHE_ROOT" ]] || die "Archives did not create the expected cache: $CACHE_ROOT"

validate_cache_count() {
  local modality="$1"
  local directory="$2"
  local expected="$3"
  local observed
  [[ -d "${CACHE_ROOT}/${directory}" ]] || die "Missing cache directory: $directory"
  observed="$(find "${CACHE_ROOT}/${directory}" -type f -name '*.npy' -print | wc -l | tr -d ' ')"
  [[ "$observed" == "$expected" ]] || \
    die "Unexpected $modality cache count: $observed != $expected"
  printf 'cache-directory\t%s\t%s\t%s\tfile-count\t%s\tauthenticated-archive-count\n' \
    "$modality" "${CACHE_ROOT}/${directory}" "${CACHE_ROOT}/${directory}" \
    "$observed" >> "$ACQUISITION_LOG"
}

validate_cache_count prott5 prott5 69811
validate_cache_count text exp_text_embeddings_temporal 69517
validate_cache_count structure IF1 67948
validate_cache_count ppi ppi 58294

echo "==> [5/7] Run the provenance-aware embedding inventory"
INVENTORY_COMMAND=(
  "$PYTHON_BIN" "${FRAMEWORK_ROOT}/scripts/verification/inventory_embeddings.py"
  --benchmark-dir "$TARGET_STAGE"
  --source-benchmark-dir "$SOURCE_STAGE"
  --embedding-cache "$CACHE_ROOT"
  --artifact-root "$PUBLISHED_ROOT"
  --config "$CONFIG_PATH"
  --output-dir "$INVENTORY_OUTPUT"
  --policy "$POLICY"
  --report-level "$REPORT_LEVEL"
)

if [[ -n "$ALIASES_FILE" ]]; then
  [[ -f "$ALIASES_FILE" ]] || die "Alias TSV does not exist: $ALIASES_FILE"
  STAGED_ALIASES="${WORK_DIR}/aliases.tsv"
  cp -p "$ALIASES_FILE" "$STAGED_ALIASES"
  INVENTORY_COMMAND+=(--aliases "$STAGED_ALIASES")
fi

printf 'Command:'
printf ' %q' "${INVENTORY_COMMAND[@]}"
printf '\n'
PYTHONDONTWRITEBYTECODE=1 "${INVENTORY_COMMAND[@]}"

echo "==> [6/7] Write an operator-facing next-step summary"
"$PYTHON_BIN" - "$INVENTORY_OUTPUT" "$OUTPUT_DIR" <<'PY'
import csv
import gzip
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

inventory = Path(sys.argv[1])
run_root = Path(sys.argv[2])
summary = json.loads((inventory / "embedding_summary.json").read_text(encoding="utf-8"))
coverage = summary["coverage"]["global"]
physical_dir = run_root / "physical_coverage"
physical_dir.mkdir()
physical = {
    modality: {"valid": set(), "not_valid": set(), "missing": set()}
    for modality in ("prott5", "text", "structure", "ppi")
}
with gzip.open(inventory / "embedding_inventory.tsv.gz", "rt", encoding="utf-8", newline="") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        modality = row["modality"]
        if row["valid"] == "true":
            physical[modality]["valid"].add(row["protein_id"])
        else:
            physical[modality]["not_valid"].add(row["protein_id"])
        if row["factual_status"] == "missing":
            physical[modality]["missing"].add(row["protein_id"])
for modality, states in physical.items():
    for state, protein_ids in states.items():
        (physical_dir / ("%s_%s.txt" % (state, modality))).write_text(
            "".join(protein_id + "\n" for protein_id in sorted(protein_ids)),
            encoding="utf-8",
        )


def metric(modality, key):
    return int(coverage["by_modality"][modality][key]["count"])


lines = [
    "# Contemporary benchmark embedding inventory",
    "",
    "- Benchmark proteins: **%s**" % format(int(summary["population"]), ","),
    "- Policy: `%s`" % summary["policy"],
    "- Inventory completion marker: `inventory/RUN_COMPLETE.json`",
    "",
    "## Global result",
    "",
    "| Modality | Physically present | Valid | Reuse | Regenerate |",
    "|---|---:|---:|---:|---:|",
]
for modality in ("prott5", "text", "structure", "ppi"):
    lines.append(
        "| %s | %s | %s | %s | %s |"
        % (
            modality,
            format(metric(modality, "present"), ","),
            format(metric(modality, "valid"), ","),
            format(metric(modality, "reuse"), ","),
            format(metric(modality, "regenerate"), ","),
        )
    )

lines.extend(
    [
        "",
        "## Files for the next stage",
        "",
        "- `inventory/reuse.tsv`: every positively proven protein/modality reuse action with evidence.",
        "- `inventory/regenerate.tsv`: every remaining protein/modality pair with the reason regeneration is required.",
        "- `inventory/reuse/{prott5,text,structure,ppi}.txt`: actionable reuse ID lists.",
        "- `inventory/regenerate/prott5.fasta`: complete sequences requiring ProtT5 generation.",
        "- `inventory/regenerate/{text,structure,ppi}.txt`: actionable generation ID lists.",
        "- `inventory/regenerate_reasons.tsv`: compact reason counts for the regeneration workload.",
        "- `physical_coverage/valid_{modality}.txt`: IDs with a loadable, finite, correctly shaped published array.",
        "- `physical_coverage/not_valid_{modality}.txt`: IDs with no currently usable published array, including missing or invalid files.",
        "- `physical_coverage/missing_{modality}.txt`: IDs with no same-ID/candidate published array.",
        "",
        "## Interpretation",
        "",
        "The operational plan has exactly two buckets. ProtT5 reuse requires an exact complete-sequence SHA-256. "
        "PPI direct-ID reuse requires a valid published array for the same UniProt accession under the fixed STRING v12 source/extractor identity; cross-ID PPI reuse is unsupported. "
        "Text and structure are regenerated unless their original description or AlphaFold structure inputs can be positively proven identical. "
        "Missing, invalid, ambiguous, unavailable, incompatible, and unknown-provenance cases all receive action `regenerate`; their finer reasons remain in the TSV reports.",
        "",
        "This run did not generate, modify, or persist any embedding arrays.",
        "",
    ]
)
(run_root / "job_summary.md").write_text("\n".join(lines), encoding="utf-8")
(run_root / "WORKFLOW_COMPLETE.json").write_text(
    json.dumps(
        {
            "schema_version": 1,
            "completed_at_utc": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
            "inventory_completion_marker": "inventory/RUN_COMPLETE.json",
            "population": int(summary["population"]),
            "policy": summary["policy"],
        },
        indent=2,
        sort_keys=True,
    )
    + "\n",
    encoding="utf-8",
)
PY

echo "==> [7/7] Completed"
echo "Persistent result root: $OUTPUT_DIR"
echo "Heavy staged/downloaded artifacts remain only under: $WORK_DIR"
echo "Read first: ${OUTPUT_DIR}/job_summary.md"
