#!/usr/bin/env bash

set -Eeuo pipefail

die() { echo "ERROR: $*" >&2; exit 2; }

FRAMEWORK_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
cd "$FRAMEWORK_ROOT"

BRANCH="$(git symbolic-ref --quiet --short HEAD || true)"
[[ "$BRANCH" == codex/cafa3-organizer-replay ]] || \
  die "Submit from codex/cafa3-organizer-replay, not $BRANCH"
[[ -z "$(git status --porcelain)" ]] || die "Framework checkout must be clean"
FRAMEWORK_COMMIT="$(git rev-parse HEAD)"
REMOTE_COMMIT="$(git ls-remote origin refs/heads/codex/cafa3-organizer-replay | awk 'NR==1 {print $1}')"
[[ "$REMOTE_COMMIT" == "$FRAMEWORK_COMMIT" ]] || \
  die "The pushed organizer-replay branch does not match local HEAD"

command -v qsub >/dev/null 2>&1 || die "qsub is unavailable"
CAMPAIGN_TAG="${CAMPAIGN_TAG:-organizer-replay-$(date -u +%Y%m%dT%H%M%SZ)}"
[[ "$CAMPAIGN_TAG" =~ ^[A-Za-z0-9._-]+$ && "$CAMPAIGN_TAG" =~ [A-Za-z0-9] ]] || \
  die "Unsafe CAMPAIGN_TAG"
RESULTS_ROOT="${RESULTS_ROOT:-$HOME/cafa3_organizer_replay_reports}"
LOG_ROOT="${LOG_ROOT:-$HOME/pfp_cafa3_organizer_replay_logs/$CAMPAIGN_TAG}"
[[ "$RESULTS_ROOT" == /* && "$RESULTS_ROOT" != / ]] || die "RESULTS_ROOT must be absolute"
[[ "$LOG_ROOT" == /* && "$LOG_ROOT" != / ]] || die "LOG_ROOT must be absolute"
conditions=(
  "assigned-date-proxy__is-a-part-of"
  "assigned-date-proxy__is-a"
  "assigned-date-proxy__all"
  "snapshot-membership__is-a-part-of"
)
for condition in "${conditions[@]}"; do
  [[ ! -e "$RESULTS_ROOT/${CAMPAIGN_TAG}-${condition}" ]] || \
    die "Campaign output already exists for condition: $condition"
  [[ ! -e "$RESULTS_ROOT/${CAMPAIGN_TAG}-${condition}.failed" ]] || \
    die "Failed campaign output already exists for condition: $condition"
done

mkdir -p "$RESULTS_ROOT/submissions" "$LOG_ROOT"
output="$(qsub -terse -t 1-4 -tc 4 -N c3orgrep -o "$LOG_ROOT" \
  -v "FRAMEWORK_COMMIT=$FRAMEWORK_COMMIT,RESULTS_ROOT=$RESULTS_ROOT,RUN_TAG=$CAMPAIGN_TAG" \
  hpc_jobs/active/hpc_cafa3_organizer_replay.sh)"
job_id="${output%%.*}"
[[ "$job_id" =~ ^[0-9]+$ ]] || die "Could not parse qsub output: $output"

ledger="$RESULTS_ROOT/submissions/${CAMPAIGN_TAG}.tsv"
{
  printf 'role\tjob_id\ttask_id\tcondition\trun_tag\tframework_commit\texpected_output\n'
  for index in "${!conditions[@]}"; do
    task_id=$((index + 1))
    condition="${conditions[$index]}"
    printf 'organizer_replay_condition\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$job_id" "$task_id" "$condition" "$CAMPAIGN_TAG" "$FRAMEWORK_COMMIT" \
      "$RESULTS_ROOT/${CAMPAIGN_TAG}-${condition}"
  done
} > "$ledger"

printf 'Organizer replay array : %s (tasks 1-4, max 4 concurrent)\n' "$job_id"
printf 'Campaign prefix        : %s/%s\n' "$RESULTS_ROOT" "$CAMPAIGN_TAG"
printf 'Submission ledger      : %s\n' "$ledger"
