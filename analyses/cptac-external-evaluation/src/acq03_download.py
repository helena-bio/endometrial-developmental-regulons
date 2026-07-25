#!/usr/bin/env python3
"""
TASK-029 ACQUISITION step 03 -- download the 230 selected STAR-Counts files + verify md5.

Downloads each selected file individually from the GDC open data endpoint
(https://api.gdc.cancer.gov/data/<file_id>), verifies the GDC-pinned md5sum,
and records per-file sha256. Any md5 mismatch is a FATAL provenance condition.

Deterministic. No expression parsed here (raw bytes only). Git HEAD unchanged.
"""
import json
import os
import sys
import hashlib
import time
import urllib.request

BASE = os.path.dirname(
    os.path.dirname(os.path.abspath(__file__))
)
WORK = os.environ.get(
    "CPTAC_WORK_DIR",
    os.path.join(BASE, "work"),
)
INTER = os.environ.get(
    "CPTAC_INTERMEDIATE_DIR",
    os.path.join(WORK, "intermediate"),
)
RAWDIR = os.environ.get(
    "CPTAC_STAR_COUNTS_DIR",
    os.path.join(WORK, "rna-star-counts"),
)
os.makedirs(INTER, exist_ok=True)
os.makedirs(RAWDIR, exist_ok=True)
DATA = "https://api.gdc.cancer.gov/data/"

lg = json.load(open(os.path.join(INTER, "acq02_join_ledger.json")))
pick = lg["pick"]
cases = sorted(pick)


def md5_of(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def sha256_of(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


verify = {}
mismatches = []
for i, c in enumerate(cases):
    p = pick[c]
    fid = p["file_id"]
    dest = os.path.join(RAWDIR, f"{fid}.rna_seq.star_counts.tsv")
    if not (os.path.exists(dest) and os.path.getsize(dest) == p["file_size"]):
        # download (retry up to 4x)
        ok = False
        for attempt in range(4):
            try:
                req = urllib.request.Request(DATA + fid)
                with urllib.request.urlopen(req, timeout=180) as r, open(dest, "wb") as out:
                    while True:
                        chunk = r.read(1 << 20)
                        if not chunk:
                            break
                        out.write(chunk)
                ok = True
                break
            except Exception as e:
                time.sleep(2 * (attempt + 1))
                last = str(e)
        if not ok:
            print("DOWNLOAD FAIL", c, fid, last, flush=True)
            sys.exit(2)
    m = md5_of(dest)
    s = sha256_of(dest)
    match = (m == p["md5sum"])
    if not match:
        mismatches.append({"case": c, "file_id": fid,
                           "expected_md5": p["md5sum"], "actual_md5": m})
    verify[c] = {"file_id": fid, "path": dest, "expected_md5": p["md5sum"],
                 "actual_md5": m, "md5_match": match, "sha256": s,
                 "size": os.path.getsize(dest)}
    if (i + 1) % 40 == 0:
        print(f"  {i+1}/{len(cases)} verified", flush=True)

out = {"n_files": len(verify), "n_md5_mismatch": len(mismatches),
       "mismatches": mismatches, "per_case": verify}
with open(os.path.join(INTER, "acq03_download_verify.json"), "w") as f:
    json.dump(out, f, indent=2)
print("done: %d files, %d md5 mismatches" % (len(verify), len(mismatches)), flush=True)
if mismatches:
    print("FATAL: md5 mismatch(es):", mismatches[:5], flush=True)
    sys.exit(3)
