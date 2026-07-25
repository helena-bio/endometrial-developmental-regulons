#!/usr/bin/env python3
"""
TASK-027 STEP 1 -- build the ONE pinned download manifest (507 primary + 35 normal).

Metadata-only. Reads the two FROZEN source TSVs (SHA-256 verified by caller before
this runs; re-verified here as a backstop) and emits download_manifest.tsv.

NO expression value is read. NO subtype is re-derived. NO threshold is tuned.
If any assertion fails -> print BLOCKED and exit non-zero; caller must NOT download.

Columns emitted (download_manifest.tsv, tab-separated):
  source_set(primary|normal), case_barcode, aliquot_barcode, sample_type,
  file_uuid, md5sum, file_size_bytes, mapped_4way(primary only, else NA)
"""
import csv
import hashlib
import os
import sys

TASK025 = os.environ.get(
    "TCGA_ACQUISITION_CHECKPOINT",
    "upstream/acquisition-checkpoint",
)
WOULD = os.path.join(
    TASK025,
    "would_download_files.tsv",
)
FULL = os.path.join(
    TASK025,
    "gdc_manifest_full.tsv",
)
OUT = os.environ.get(
    "TCGA_DOWNLOAD_MANIFEST",
    "download_manifest.tsv",
)

EXPECT = {
    WOULD: "0801754019fe0b69d88980c650c4f4d604558ec05afd334b714315356acca8af",
    FULL: "85899f606b64c9f53c0115670eb8997889c47ad833e09c6041a66cab8c5234cb",
}
EXPECT_N_PRIMARY = 507
EXPECT_N_NORMAL = 35
EXPECT_N_TOTAL = 542
EXPECT_TOTAL_BYTES = 2287429373
EXPECT_PRIMARY_BYTES = 2139373015
EXPECT_NORMAL_BYTES = 148056358


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fail(msg):
    print("BLOCKED: " + msg, file=sys.stderr)
    sys.exit(2)


# --- backstop SHA-256 re-verification of the two source TSVs ---
for path, exp in EXPECT.items():
    got = sha256(path)
    if got != exp:
        fail("SHA-256 mismatch for %s: got %s expected %s" % (path, got, exp))
    print("SHA-256 OK  %s  %s" % (exp, path))

rows = []  # ordered manifest rows

# --- PRIMARY: would_download_files.tsv ---
# columns: patient_barcode, mapped_4way, file_uuid, md5sum, file_size_bytes, aliquot_barcode
with open(WOULD, newline="") as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        rows.append({
            "source_set": "primary",
            "case_barcode": row["patient_barcode"].strip(),
            "aliquot_barcode": row["aliquot_barcode"].strip(),
            "sample_type": "Primary Tumor",
            "file_uuid": row["file_uuid"].strip(),
            "md5sum": row["md5sum"].strip(),
            "file_size_bytes": int(row["file_size_bytes"]),
            "mapped_4way": row["mapped_4way"].strip(),
        })
n_primary = len(rows)

# --- NORMAL: gdc_manifest_full.tsv rows with sample_type == 'Solid Tissue Normal' ---
# columns: file_id, file_name, md5sum, file_size, state, n_cases, n_samples,
#          n_aliquots, case_barcode, sample_barcode, sample_type, aliquot_barcode
n_normal = 0
with open(FULL, newline="") as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        if row["sample_type"].strip() != "Solid Tissue Normal":
            continue
        rows.append({
            "source_set": "normal",
            "case_barcode": row["case_barcode"].strip(),
            "aliquot_barcode": row["aliquot_barcode"].strip(),
            "sample_type": "Solid Tissue Normal",
            "file_uuid": row["file_id"].strip(),
            "md5sum": row["md5sum"].strip(),
            "file_size_bytes": int(row["file_size"]),
            "mapped_4way": "NA",
        })
        n_normal += 1

n_total = len(rows)

# --- ASSERTIONS (STEP 1) ---
if n_primary != EXPECT_N_PRIMARY:
    fail("n_primary==%d expected %d" % (n_primary, EXPECT_N_PRIMARY))
if n_normal != EXPECT_N_NORMAL:
    fail("n_normal==%d expected %d" % (n_normal, EXPECT_N_NORMAL))
if n_total != EXPECT_N_TOTAL:
    fail("n_total==%d expected %d" % (n_total, EXPECT_N_TOTAL))

total_bytes = sum(r["file_size_bytes"] for r in rows)
primary_bytes = sum(r["file_size_bytes"] for r in rows if r["source_set"] == "primary")
normal_bytes = sum(r["file_size_bytes"] for r in rows if r["source_set"] == "normal")
if total_bytes != EXPECT_TOTAL_BYTES:
    fail("sum(file_size_bytes)==%d expected %d" % (total_bytes, EXPECT_TOTAL_BYTES))
if primary_bytes != EXPECT_PRIMARY_BYTES:
    fail("primary bytes==%d expected %d" % (primary_bytes, EXPECT_PRIMARY_BYTES))
if normal_bytes != EXPECT_NORMAL_BYTES:
    fail("normal bytes==%d expected %d" % (normal_bytes, EXPECT_NORMAL_BYTES))

# all file_uuid / md5sum non-empty and unique
uuids = [r["file_uuid"] for r in rows]
md5s = [r["md5sum"] for r in rows]
if any(not u for u in uuids):
    fail("empty file_uuid present")
if any(not m for m in md5s):
    fail("empty md5sum present")
if len(set(uuids)) != n_total:
    fail("file_uuid not unique: %d unique of %d" % (len(set(uuids)), n_total))
if len(set(md5s)) != n_total:
    fail("md5sum not unique: %d unique of %d" % (len(set(md5s)), n_total))

# --- write manifest ---
fields = ["source_set", "case_barcode", "aliquot_barcode", "sample_type",
          "file_uuid", "md5sum", "file_size_bytes", "mapped_4way"]
with open(OUT, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
    w.writeheader()
    for row in rows:
        w.writerow(row)

print("ASSERTIONS PASS")
print("  n_primary = %d" % n_primary)
print("  n_normal  = %d" % n_normal)
print("  n_total   = %d" % n_total)
print("  total_bytes   = %d" % total_bytes)
print("  primary_bytes = %d" % primary_bytes)
print("  normal_bytes  = %d" % normal_bytes)
print("  all uuid/md5 non-empty and unique: yes")
print("WROTE %s" % OUT)
