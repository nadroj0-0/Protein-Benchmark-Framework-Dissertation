#!/usr/bin/env bash
# Copy read-only progress evidence for active homology tasks on one execution host.

#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -l tmem=1G
#$ -l h_rt=00:20:00
#$ -N homology_snapshot

set -Eeuo pipefail

usage() {
  cat <<'EOF'
Usage: qsub -q cpu.q@HOST hpc_jobs/active/hpc_homology_progress_snapshot.sh \
  --main-job-id JOB_ID --expected-host HOST \
  --task TASK_ID:IDENTITY [--task TASK_ID:IDENTITY ...] \
  [--destination-root PATH]

The job reads active task logs from host-local scratch and publishes one atomic
snapshot under the user's home directory. It never modifies the source jobs.
EOF
}

die() { echo "ERROR: $*" >&2; exit 2; }

MAIN_JOB_ID=""
EXPECTED_HOST=""
DESTINATION_ROOT="${HOME}/homology_cluster_progress"
TASK_SPECS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --main-job-id) MAIN_JOB_ID="$2"; shift 2 ;;
    --expected-host) EXPECTED_HOST="$2"; shift 2 ;;
    --task) TASK_SPECS+=("$2"); shift 2 ;;
    --destination-root) DESTINATION_ROOT="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) usage >&2; die "Unknown argument: $1" ;;
  esac
done

[[ "$MAIN_JOB_ID" =~ ^[1-9][0-9]*$ ]] || die "--main-job-id must be numeric"
[[ "$EXPECTED_HOST" =~ ^[A-Za-z0-9._-]+$ ]] || die "--expected-host is required"
[[ "${#TASK_SPECS[@]}" -gt 0 ]] || die "At least one --task is required"
[[ "$DESTINATION_ROOT" == /* && "$DESTINATION_ROOT" != "/" ]] || \
  die "--destination-root must be an absolute non-root path"

host_short="$(hostname -s)"
expected_short="${EXPECTED_HOST%%.*}"
[[ "$host_short" == "$expected_short" ]] || \
  die "Expected host $expected_short but running on $host_short"

SCRATCH_BASE="${HOMOLOGY_SNAPSHOT_SCRATCH_BASE:-/scratch0}"
[[ -d "$SCRATCH_BASE" ]] || die "Scratch base is unavailable: $SCRATCH_BASE"
stamp="$(date -u +%Y%m%dT%H%M%SZ)"
host_root="$DESTINATION_ROOT/job_${MAIN_JOB_ID}/host_${host_short}"
temporary="$host_root/.snapshot_${stamp}.${JOB_ID:-manual}.tmp"
final="$host_root/snapshot_${stamp}"
[[ ! -e "$temporary" && ! -e "$final" ]] || die "Snapshot destination already exists"
mkdir -p "$temporary"

inventory="$temporary/copied_log_inventory.tsv"
printf 'task_id\tidentity\tsource\tsize_before\tmtime_before\tsize_after\tmtime_after\tstable_during_copy\tsnapshot\tsha256\n' \
  > "$inventory"

file_size() {
  stat -c '%s' "$1" 2>/dev/null || stat -f '%z' "$1"
}

file_mtime() {
  stat -c '%Y' "$1" 2>/dev/null || stat -f '%m' "$1"
}

file_sha256() {
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$1" | awk '{print $1}'
  else
    shasum -a 256 "$1" | awk '{print $1}'
  fi
}

write_tree_inventory() {
  local root="$1" depth="$2" output="$3"
  find "$root" -maxdepth "$depth" -print | while IFS= read -r path; do
    if [[ -d "$path" ]]; then
      kind=d
    elif [[ -f "$path" ]]; then
      kind=f
    elif [[ -L "$path" ]]; then
      kind=l
    else
      kind=o
    fi
    printf '%s\t%s\t%s\t%s\n' \
      "$kind" "$(file_size "$path")" "$(file_mtime "$path")" "$path"
  done | sort > "$output"
}

copy_log() {
  local task_id="$1" identity="$2" source="$3" destination="$4"
  local size_before mtime_before size_after mtime_after stable digest
  size_before="$(file_size "$source")"
  mtime_before="$(file_mtime "$source")"
  mkdir -p "$(dirname "$destination")"
  cp -- "$source" "$destination"
  size_after="$(file_size "$source")"
  mtime_after="$(file_mtime "$source")"
  stable=false
  [[ "$size_before" == "$size_after" && "$mtime_before" == "$mtime_after" ]] && stable=true
  digest="$(file_sha256 "$destination")"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$task_id" "$identity" "$source" "$size_before" "$mtime_before" \
    "$size_after" "$mtime_after" "$stable" \
    "${destination#"$temporary"/}" "$digest" >> "$inventory"
}

printf 'field\tvalue\nmain_job_id\t%s\nsnapshot_job_id\t%s\nhost\t%s\nsnapshot_utc\t%s\n' \
  "$MAIN_JOB_ID" "${JOB_ID:-manual}" "$host_short" "$stamp" \
  > "$temporary/snapshot_metadata.tsv"

found_tasks=0
for task_spec in "${TASK_SPECS[@]}"; do
  [[ "$task_spec" =~ ^([1-6]):(30|25|20|15|10|5)$ ]] || \
    die "Invalid task specification: $task_spec"
  task_id="${BASH_REMATCH[1]}"
  identity="${BASH_REMATCH[2]}"
  task_root="$temporary/task_${task_id}_identity_${identity}"
  mkdir -p "$task_root/logs/mmseqs"

  candidates=(
    "$SCRATCH_BASE/homology_runtime_${MAIN_JOB_ID}_${task_id}_${identity}_runtime-${MAIN_JOB_ID}"
  )
  work=""
  for candidate in "${candidates[@]}"; do
    if [[ -d "$candidate" ]]; then
      work="$candidate"
      break
    fi
  done
  if [[ -z "$work" ]]; then
    printf 'source_status\tmissing\nexpected_path\t%s\n' "${candidates[0]}" \
      > "$task_root/status.tsv"
    continue
  fi
  found_tasks=$((found_tasks + 1))
  printf 'source_status\tfound\nsource_path\t%s\n' "$work" > "$task_root/status.tsv"

  for source in "$work"/artifacts/logs/*; do
    [[ -f "$source" ]] || continue
    copy_log "$task_id" "$identity" "$source" "$task_root/logs/$(basename "$source")"
  done
  for log_dir in "$work"/tmp/homology-*/logs/mmseqs; do
    [[ -d "$log_dir" ]] || continue
    for source in "$log_dir"/*.log; do
      [[ -f "$source" ]] || continue
      copy_log "$task_id" "$identity" "$source" \
        "$task_root/logs/mmseqs/$(basename "$source")"
    done
  done

  if ! ps -u "$USER" -o pid,ppid,etime,pcpu,pmem,rss,vsz,state,cmd \
      | awk -v needle="homology_runtime_${MAIN_JOB_ID}_${task_id}_${identity}" \
          'NR == 1 || index($0, needle)' > "$task_root/process_snapshot.txt"; then
    printf 'Process snapshot unavailable on host %s\n' "$host_short" \
      > "$task_root/process_snapshot.txt"
  fi
  write_tree_inventory "$work" 2 "$task_root/work_tree_depth2.tsv"
  for mmseqs_root in "$work"/tmp/homology-*/mmseqs; do
    [[ -d "$mmseqs_root" ]] || continue
    write_tree_inventory "$mmseqs_root" 3 "$task_root/mmseqs_tree_depth3.tsv"
  done
done

[[ "$found_tasks" -gt 0 ]] || die "No requested active task scratch directories were found"
df -h "$SCRATCH_BASE" > "$temporary/scratch_filesystem.txt"
qstat -j "$MAIN_JOB_ID" > "$temporary/grid_engine_job.txt" 2>&1 || true

mkdir -p "$host_root"
mv "$temporary" "$final"
ln -sfn "$(basename "$final")" "$host_root/latest"
echo "Snapshot complete: $final"
echo "Captured $found_tasks of ${#TASK_SPECS[@]} requested task(s)."
