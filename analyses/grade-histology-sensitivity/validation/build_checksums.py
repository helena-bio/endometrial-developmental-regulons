#!/usr/bin/env python3
"""Write SHA-256 inventory for all critic files except the inventory itself."""

import hashlib
from pathlib import Path


HERE = Path(__file__).resolve().parent
rows = []
for path in sorted(HERE.iterdir()):
    if not path.is_file() or path.name == "CRITIC_SHA256SUMS.txt":
        continue
    rows.append(f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.name}")
(HERE / "CRITIC_SHA256SUMS.txt").write_text("\n".join(rows) + "\n", encoding="ascii")
print(f"checksums={len(rows)}")
