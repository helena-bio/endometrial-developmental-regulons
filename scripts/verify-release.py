from __future__ import annotations

import csv
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / 'release' / 'file-manifest.tsv'
SELF = 'scripts/verify-release.py'
TEXT_SUFFIXES = {
    '.md', '.txt', '.tsv', '.csv', '.json', '.jsonl',
    '.yaml', '.yml', '.py', '.sh', '.log', '.lock',
}
UNICODE_TEXT_ALLOWLIST = {
    'LICENSES/Apache-2.0.txt',
    'LICENSES/CC-BY-4.0.txt',
}
RUNTIME_PARTS = {
    '.venv', 'venv', 'venv_phase1', '__pycache__',
    '.pytest_cache', '.mypy_cache', 'site-packages',
}
FORBIDDEN_PATH_TERMS = (
    'critic', 'chatgpt', 'claude', 'anthropic', 'openai',
    'codex', 'sophiagpt', 'sofiagpt', 'zeusgpt',
    'prompt', 'conversation',
)
FORBIDDEN_TEXT_TERMS = (
    'chatgpt', 'claude', 'anthropic', 'openai', 'codex',
    'sophiagpt', 'sofiagpt', 'zeusgpt',
)
SECRET_TERMS = (
    'github_pat_', 'ghp_', 'x-access-token',
    'authorization: bearer', 'authorization: token',
    'api_key', 'api-key', 'client_secret', 'client-secret',
)
PRIVATE_PATHS = (
    '/srv/sophia-workspace', '/home/sophia',
    '/srv/zeusgpt-workspace', '/home/zeusgpt',
)


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            value.update(chunk)
    return value.hexdigest()


def main() -> None:
    failures = []
    with MANIFEST.open(encoding='utf-8', newline='') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    for row in rows:
        relative = row['path']
        relative_lower = relative.lower()
        path = ROOT / relative
        if any(term in relative_lower for term in FORBIDDEN_PATH_TERMS):
            failures.append(f'forbidden public path: {relative}')
        if not path.is_file():
            failures.append(f'missing file: {relative}')
            continue
        if path.stat().st_size != int(row['size_bytes']):
            failures.append(f'size mismatch: {relative}')
        if digest(path) != row['sha256']:
            failures.append(f'checksum mismatch: {relative}')
        if path.stat().st_size >= 95 * 1024 * 1024:
            failures.append(f'file exceeds release limit: {relative}')
        if any(part in RUNTIME_PARTS or part.startswith('venv_') for part in path.parts):
            failures.append(f'runtime material: {relative}')
        if path.suffix.lower() in {'.pyc', '.pyo'}:
            failures.append(f'compiled runtime file: {relative}')
        if relative == SELF or path.suffix.lower() not in TEXT_SUFFIXES:
            continue
        text = path.read_text(encoding='utf-8')
        lower = text.lower()
        if any(term in lower for term in FORBIDDEN_TEXT_TERMS):
            failures.append(f'forbidden internal-system term: {relative}')
        if any(term in lower for term in SECRET_TERMS):
            failures.append(f'secret-like content: {relative}')
        if any(term in text for term in PRIVATE_PATHS):
            failures.append(f'private absolute path: {relative}')
        if (
            relative not in UNICODE_TEXT_ALLOWLIST
            and any(ord(char) > 127 for char in text)
        ):
            failures.append(f'non-ASCII text: {relative}')
    if failures:
        for failure in failures:
            print(f'FAIL: {failure}')
        raise SystemExit(1)
    print(f'PASS: verified {len(rows)} release files')


if __name__ == '__main__':
    main()
