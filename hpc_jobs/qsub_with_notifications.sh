#!/usr/bin/env bash

# Submit through Grid Engine with optional machine-local email notifications.
# The address file is deliberately outside the repository.

set -uo pipefail

EMAIL_FILE="${GRID_ENGINE_NOTIFY_EMAIL_FILE:-${HOME:-}/.grid_engine_notify_email}"
QSUB_BIN="${GRID_ENGINE_QSUB_BIN:-}"

warn_and_submit_normally() {
  local reason="$1"
  shift
  printf 'WARNING: Grid Engine email notifications disabled: %s\n' "$reason" >&2
  exec "$QSUB_BIN" "$@"
}

if [[ -z "$QSUB_BIN" ]]; then
  QSUB_BIN="$(type -P qsub 2>/dev/null || true)"
fi

if [[ -z "$QSUB_BIN" || ! -x "$QSUB_BIN" ]]; then
  printf 'ERROR: qsub executable not found\n' >&2
  exit 127
fi

if [[ -z "$EMAIL_FILE" || ! -f "$EMAIL_FILE" || ! -r "$EMAIL_FILE" ]]; then
  warn_and_submit_normally "notification file is missing or unreadable" "$@"
fi

line_count="$(awk 'END { print NR + 0 }' "$EMAIL_FILE" 2>/dev/null)" ||
  warn_and_submit_normally "notification file could not be read" "$@"

email="$(sed -n '1p' "$EMAIL_FILE" 2>/dev/null)" ||
  warn_and_submit_normally "notification file could not be read" "$@"
email="${email%$'\r'}"

local_part="${email%@*}"
domain_part="${email#*@}"
if [[ "$line_count" != "1" || -z "$local_part" || -z "$domain_part" ||
      "$local_part" == "$email" || "$domain_part" == *@* ||
      "$domain_part" != *.* || "$email" == *[[:space:],]* ]]; then
  warn_and_submit_normally "notification file does not contain exactly one valid email address" "$@"
fi

printf 'Grid Engine email notifications enabled for begin, end, abort/reschedule, and suspend events\n' >&2
exec "$QSUB_BIN" -m beas -M "$email" "$@"
