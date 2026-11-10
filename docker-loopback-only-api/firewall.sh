#!/usr/bin/env bash
# ABOUTME: Limits a port published by Docker to one source address, using the DOCKER-USER chain.
# ABOUTME: Safe to run again: it never stacks duplicate rules and leaves ACCEPT above DROP.

set -euo pipefail

usage() {
  {
    echo "Usage: firewall.sh [--dry-run] <source> [port]"
    echo "  source   the one IPv4 address or network allowed to reach the port"
    echo "  port     the port inside the container (default 8081)"
    echo "  --dry-run  print the rules without changing anything"
  } >&2
  exit 2
}

dry_run=false
if [[ "${1:-}" == "--dry-run" ]]; then
  dry_run=true
  shift
fi
[[ $# -ge 1 && $# -le 2 ]] || usage

source_address=$1
port=${2:-8081}

valid_source() {
  local octet
  [[ $1 =~ ^([0-9]{1,3})\.([0-9]{1,3})\.([0-9]{1,3})\.([0-9]{1,3})(/([0-9]{1,2}))?$ ]] || return 1
  for octet in "${BASH_REMATCH[@]:1:4}"; do
    ((10#$octet <= 255)) || return 1
  done
  [[ -z "${BASH_REMATCH[6]}" ]] || ((10#${BASH_REMATCH[6]} <= 32))
}

valid_port() {
  [[ $1 =~ ^[0-9]{1,5}$ ]] && ((10#$1 >= 1 && 10#$1 <= 65535))
}

if ! valid_source "$source_address"; then
  echo "firewall.sh: the source must be an IPv4 address or network, such as 192.0.2.20" >&2
  exit 2
fi
if ! valid_port "$port"; then
  echo "firewall.sh: the port must be a number from 1 to 65535" >&2
  exit 2
fi

# DOCKER-USER sees a request after its destination has been rewritten, so --dport
# is the port inside the container. The rules match forwarded traffic to that
# port for every container, so choose a port nothing else uses.
drop=(DOCKER-USER -p tcp --dport "$port" -j DROP)
accept=(DOCKER-USER -p tcp --dport "$port" -s "$source_address" -j ACCEPT)

if $dry_run; then
  echo "iptables -I ${drop[*]}"
  echo "iptables -I ${accept[*]}"
  echo "iptables -L DOCKER-USER -n --line-numbers"
  exit 0
fi

# -I inserts at the top of the chain, so DROP goes in first and ACCEPT lands above it.
if ! iptables -C "${drop[@]}" 2>/dev/null; then
  iptables -I "${drop[@]}"
fi

# An ACCEPT left by an earlier run may sit anywhere. Remove it and insert it
# again, so it is back on top. The DROP rule stays in place throughout.
if iptables -C "${accept[@]}" 2>/dev/null; then
  iptables -D "${accept[@]}"
fi
iptables -I "${accept[@]}"

# ACCEPT must be listed above DROP.
iptables -L DOCKER-USER -n --line-numbers
