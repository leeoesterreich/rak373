#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 /path/to/runIchorCNA.R [args...]" >&2
  exit 2
fi

RSCRIPT_PATH="$1"
shift

ATTEMPTS="${ATTEMPTS:-3}"
SLEEP_SECONDS="${SLEEP_SECONDS:-10}"

for attempt in $(seq 1 "$ATTEMPTS"); do
  # runIchorCNA.R writes to stdout/stderr; caller typically redirects both to a log.
  if Rscript "$RSCRIPT_PATH" "$@"; then
    exit 0
  fi

  echo "runIchorCNA.R attempt ${attempt}/${ATTEMPTS} failed; sleeping ${SLEEP_SECONDS}s..." >&2
  sleep "$SLEEP_SECONDS"
done

echo "runIchorCNA.R failed after ${ATTEMPTS} attempts." >&2
exit 1

