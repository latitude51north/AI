#!/usr/bin/env bash
# ABOUTME: Checks that the proxied API is reachable, and that each control on it is working.
# ABOUTME: Run it from the allowed machine with the token, and with --blocked from any other.

set -euo pipefail

usage() {
  {
    echo "Usage: TOKEN=<token> verify.sh <base-url>     from the allowed machine"
    echo "       verify.sh --blocked <base-url>         from any other machine"
    echo "  base-url   where the proxy is published, such as http://192.0.2.10:8081"
  } >&2
  exit 2
}

blocked=false
if [[ "${1:-}" == "--blocked" ]]; then
  blocked=true
  shift
fi
[[ $# -eq 1 ]] || usage

url="${1%/}/api/health"
failed=0

# Prints the HTTP status of a GET, or 000 when no answer came within five seconds.
status() {
  curl -s -m 5 -o /dev/null -w '%{http_code}' "$@" || true
}

check() {
  local label=$1 expected=$2 got=$3
  if [[ "$got" == "$expected" ]]; then
    echo "PASS  $label: $got"
  else
    echo "FAIL  $label: $got, expected $expected"
    failed=1
  fi
}

if $blocked; then
  got=$(status "$url")
  if [[ "$got" == "000" ]]; then
    echo "PASS  from this machine: no answer"
  else
    echo "FAIL  from this machine: $got, expected no answer"
    failed=1
  fi
  exit "$failed"
fi

if [[ -z "${TOKEN:-}" ]]; then
  echo "verify.sh: set TOKEN to the API's bearer token" >&2
  exit 2
fi

# A 401 proves the whole path works, from here through the published port and
# the sidecar to the API, and that the API still refuses a caller with no token.
check "without the token" 401 "$(status "$url")"

# The header is passed on standard input, so the token never appears in the
# process list.
check "with the token" 200 "$(status -H @- "$url" <<<"Authorization: Bearer $TOKEN")"

exit "$failed"
