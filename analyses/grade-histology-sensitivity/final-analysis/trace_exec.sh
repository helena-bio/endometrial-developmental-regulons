#!/usr/bin/env bash
set -u
set -o pipefail

if [[ "$#" -lt 4 ]]; then
  echo "usage: trace_exec.sh TRACE STDOUT STDERR COMMAND [ARGS ...]" >&2
  exit 64
fi

trace_path="$1"
stdout_path="$2"
stderr_path="$3"
shift 3

exec /usr/bin/strace -f -qq -e trace=%file,%network -o "${trace_path}" \
  "$@" > "${stdout_path}" 2> "${stderr_path}"
