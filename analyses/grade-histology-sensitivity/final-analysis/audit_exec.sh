#!/usr/bin/env bash
set -u
set -o pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AUDIT_DIR="${ROOT}/audit"
mkdir -p "${AUDIT_DIR}/commands"

if [[ "${1:-}" == "--" ]]; then
  shift
fi
if [[ "$#" -eq 0 ]]; then
  echo "usage: audit_exec.sh -- command [args ...]" >&2
  exit 64
fi

counter_file="${AUDIT_DIR}/counter"
if [[ -f "${counter_file}" ]]; then
  counter="$(<"${counter_file}")"
else
  counter=0
fi
counter=$((counter + 1))
printf '%s\n' "${counter}" > "${counter_file}"
command_id="$(printf '%04d' "${counter}")"

stdout_path="${AUDIT_DIR}/commands/${command_id}.stdout"
stderr_path="${AUDIT_DIR}/commands/${command_id}.stderr"
meta_path="${AUDIT_DIR}/commands/${command_id}.meta"
start_utc="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
start_epoch_ns="$(date -u +%s%N)"
cwd="$(pwd -P)"
printf -v command_q '%q ' "$@"
command_q="${command_q% }"

{
  printf 'command_id=%s\n' "${command_id}"
  printf 'start_utc=%s\n' "${start_utc}"
  printf 'start_epoch_ns=%s\n' "${start_epoch_ns}"
  printf 'cwd=%q\n' "${cwd}"
  printf 'command=%s\n' "${command_q}"
  printf 'wrapper_sha256=%s\n' "$(sha256sum "${ROOT}/audit_exec.sh" | awk '{print $1}')"
} > "${meta_path}"

set +e
"$@" > "${stdout_path}" 2> "${stderr_path}"
exit_status=$?
set -e
if [[ -s "${stdout_path}" ]]; then
  /usr/bin/cat "${stdout_path}"
fi
if [[ -s "${stderr_path}" ]]; then
  /usr/bin/cat "${stderr_path}" >&2
fi

end_utc="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
end_epoch_ns="$(date -u +%s%N)"
{
  printf 'end_utc=%s\n' "${end_utc}"
  printf 'end_epoch_ns=%s\n' "${end_epoch_ns}"
  printf 'exit_status=%s\n' "${exit_status}"
  printf 'stdout_sha256=%s\n' "$(sha256sum "${stdout_path}" | awk '{print $1}')"
  printf 'stderr_sha256=%s\n' "$(sha256sum "${stderr_path}" | awk '{print $1}')"
} >> "${meta_path}"

printf '%s\t%s\t%s\t%s\t%s\n' \
  "${command_id}" "${start_utc}" "${end_utc}" "${exit_status}" "${command_q}" \
  >> "${AUDIT_DIR}/command_log.tsv"

exit "${exit_status}"
