#!/usr/bin/env python3
"""
TASK-027 STEP 3 -- download EXACTLY the 542 pinned UUIDs and verify md5 + size.

Reads download_manifest.tsv (built + asserted by build_download_manifest.py).
For each UUID: GET https://api.gdc.cancer.gov/data/<uuid> into originals/,
confirm HTTP 200, byte length == pinned size, md5 == pinned md5 (HARD STOP on
mismatch -- never substitute or re-fetch a different file), record sha256.

Robust to transient network errors: retry with exponential backoff, a few
attempts. Files already present with a matching md5 are skipped (idempotent
resume; a present-but-mismatched partial is deleted and re-fetched).

Writes checksums.tsv:
  source_set, file_uuid, file_name, pinned_md5, computed_md5, md5_match,
  sha256, pinned_size, computed_size, size_match

At the end ASSERTs 542/542 md5_match, 542/542 size_match,
sum(computed_size)==2287429373. Any failure -> non-zero exit + BLOCKED message.

NO expression value is read. This is a byte-exact fetch + checksum only.
"""
import csv
import hashlib
import os
import sys
import time
import urllib.request
import urllib.error

BASE = os.environ.get(
    "TCGA_ACQUISITION_WORKDIR",
    os.getcwd(),
)
MANIFEST = os.environ.get(
    "TCGA_DOWNLOAD_MANIFEST",
    os.path.join(BASE, "download_manifest.tsv"),
)
ORIG = os.environ.get(
    "TCGA_STAR_ORIGINALS",
    os.path.join(BASE, "originals"),
)
CHECKSUMS = BASE + "/checksums.tsv"
DATA_URL = "https://api.gdc.cancer.gov/data/"

EXPECT_TOTAL_BYTES = 2287429373
EXPECT_N = 542
MAX_ATTEMPTS = 6
BACKOFF_BASE = 3.0  # seconds; 3, 6, 12, 24, 48


def digest_file(path):
    """Return (md5hex, sha256hex, nbytes) streamed."""
    m = hashlib.md5()
    s = hashlib.sha256()
    n = 0
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            m.update(chunk)
            s.update(chunk)
            n += len(chunk)
    return m.hexdigest(), s.hexdigest(), n


def fetch(uuid, dest):
    """Download one UUID to dest with retry/backoff. Returns (http_status, filename_hint)."""
    url = DATA_URL + uuid
    last_err = None
    for attempt in range(1, MAX_ATTEMPTS + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "task027-acquire/1.0"})
            with urllib.request.urlopen(req, timeout=300) as resp:
                status = resp.getcode()
                if status != 200:
                    raise urllib.error.HTTPError(url, status, "non-200", resp.headers, None)
                # GDC sets Content-Disposition: attachment; filename=<gdc filename>
                cd = resp.headers.get("Content-Disposition", "")
                fname_hint = ""
                if "filename=" in cd:
                    fname_hint = cd.split("filename=", 1)[1].strip().strip('"')
                tmp = dest + ".part"
                with open(tmp, "wb") as out:
                    while True:
                        chunk = resp.read(1 << 20)
                        if not chunk:
                            break
                        out.write(chunk)
                os.replace(tmp, dest)
                return status, fname_hint
        except Exception as e:  # noqa: BLE001 -- transient net errors, retry
            last_err = e
            if attempt < MAX_ATTEMPTS:
                wait = BACKOFF_BASE * (2 ** (attempt - 1))
                print("  retry %d/%d for %s after error: %s (sleep %.0fs)"
                      % (attempt, MAX_ATTEMPTS, uuid, e, wait), flush=True)
                time.sleep(wait)
            else:
                raise RuntimeError("download failed after %d attempts: %s"
                                   % (MAX_ATTEMPTS, last_err))
    raise RuntimeError("unreachable")


def main():
    rows = []
    with open(MANIFEST, newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            rows.append(r)
    assert len(rows) == EXPECT_N, "manifest has %d rows, expected %d" % (len(rows), EXPECT_N)

    results = []
    hard_stops = []
    n_done = 0
    for r in rows:
        uuid = r["file_uuid"]
        pinned_md5 = r["md5sum"]
        pinned_size = int(r["file_size_bytes"])
        source_set = r["source_set"]
        # deterministic on-disk name: <uuid>.tsv (GDC filenames recorded separately)
        dest = os.path.join(ORIG, uuid + "__augmented_star_gene_counts.tsv")

        fname_hint = ""
        # idempotent resume: if present and md5 already matches, skip re-fetch
        need_fetch = True
        if os.path.exists(dest):
            cmd5, csha, cn = digest_file(dest)
            if cmd5 == pinned_md5 and cn == pinned_size:
                need_fetch = False
            else:
                os.remove(dest)  # partial/mismatched -> re-fetch clean

        if need_fetch:
            status, fname_hint = fetch(uuid, dest)
            if status != 200:
                hard_stops.append((uuid, "HTTP %s" % status))
                continue

        cmd5, csha, cn = digest_file(dest)
        md5_match = (cmd5 == pinned_md5)
        size_match = (cn == pinned_size)
        results.append({
            "source_set": source_set,
            "file_uuid": uuid,
            "file_name": fname_hint or (uuid + "__augmented_star_gene_counts.tsv"),
            "pinned_md5": pinned_md5,
            "computed_md5": cmd5,
            "md5_match": str(md5_match).lower(),
            "sha256": csha,
            "pinned_size": pinned_size,
            "computed_size": cn,
            "size_match": str(size_match).lower(),
        })
        if not md5_match:
            hard_stops.append((uuid, "MD5 MISMATCH pinned=%s computed=%s" % (pinned_md5, cmd5)))
        if not size_match:
            hard_stops.append((uuid, "SIZE MISMATCH pinned=%d computed=%d" % (pinned_size, cn)))

        n_done += 1
        if n_done % 50 == 0:
            print("  progress: %d/%d files verified" % (n_done, EXPECT_N), flush=True)

    # write checksums.tsv (whatever we have, even on partial failure -- provenance)
    fields = ["source_set", "file_uuid", "file_name", "pinned_md5", "computed_md5",
              "md5_match", "sha256", "pinned_size", "computed_size", "size_match"]
    with open(CHECKSUMS, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for row in results:
            w.writerow(row)

    n_md5 = sum(1 for r in results if r["md5_match"] == "true")
    n_size = sum(1 for r in results if r["size_match"] == "true")
    total_bytes = sum(r["computed_size"] for r in results)

    print("")
    print("SUMMARY")
    print("  files verified   = %d / %d" % (len(results), EXPECT_N))
    print("  md5_match true   = %d / %d" % (n_md5, EXPECT_N))
    print("  size_match true  = %d / %d" % (n_size, EXPECT_N))
    print("  sum computed_size = %d (expect %d)" % (total_bytes, EXPECT_TOTAL_BYTES))

    if hard_stops:
        print("")
        print("BLOCKED: HARD STOP(S) DETECTED:")
        for uuid, why in hard_stops:
            print("  %s -> %s" % (uuid, why))
        sys.exit(3)

    if len(results) != EXPECT_N or n_md5 != EXPECT_N or n_size != EXPECT_N:
        print("BLOCKED: not all 542 files verified clean.")
        sys.exit(4)
    if total_bytes != EXPECT_TOTAL_BYTES:
        print("BLOCKED: byte total %d != %d" % (total_bytes, EXPECT_TOTAL_BYTES))
        sys.exit(5)

    print("")
    print("ALL 542/542 md5_match true, 542/542 size_match true, byte total EXACT.")
    print("WROTE %s" % CHECKSUMS)


if __name__ == "__main__":
    main()
