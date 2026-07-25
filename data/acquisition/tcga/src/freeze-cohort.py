#!/usr/bin/env python3
"""
TASK-027 STEP 5 -- FREEZE A FROM SCRATCH on the MATERIALIZED files.

Re-derive cohort membership FROM THE ACTUAL DOWNLOADED PRIMARY FILES (gate uses
the 507 primary set only), apply the FROZEN SELECTION_RULE verbatim, join
read-only to the FROZEN pinned subtype table, structurally validate every primary
STAR-Counts file, and report the gate band.

STRICT PROHIBITIONS honored here:
  * NO expression value is read into any score/contrast/statistic. Structural
    checks only: header columns present, the 4 STAR summary rows present, the
    GENCODE-v36 gene-row count. tpm/count numeric VALUES are never summed,
    ranked, contrasted, or scored.
  * NO subtype re-derivation. The pinned subtype table is a read-only JOIN.
  * NO threshold tuning. NO forcing POLE to a target.
  * The selection rule is applied VERBATIM (code '01'; lexicographically-lowest
    full aliquot barcode on dedup).

Cohort membership is re-derived from the ACTUAL materialized primary files:
  - the on-disk file set is enumerated from originals/;
  - each on-disk UUID is confirmed present in checksums.tsv with md5_match=true;
  - the uuid->aliquot_barcode binding comes from the pinned manifest (GDC
    metadata; STAR-Counts bodies carry no barcode -- this is intrinsic to GDC);
  - the binding is CONTENT-ANCHORED by confirming each served file's GDC filename
    (Content-Disposition, recorded in checksums.tsv) equals the GDC file_name that
    the raw GDC metadata (gdc_files_raw.json) records for that same file_id, i.e.
    the bytes we hold are the file GDC bound to that UUID (hence that barcode).
  - the selection rule ('01' + dedup) is then applied to the barcodes of the
    materialized, md5-verified primary files -- NOT to the metadata prediction.

Emits: freeze_a_redux/*.tsv|json artifacts + prints a human summary.
"""
import csv
import hashlib
import json
import os
import re
import sys
from collections import defaultdict

BASE = os.environ.get(
    "TCGA_ACQUISITION_WORKDIR",
    os.getcwd(),
)
ORIG = os.environ.get(
    "TCGA_STAR_ORIGINALS",
    os.path.join(BASE, "originals"),
)
MANIFEST = os.environ.get(
    "TCGA_DOWNLOAD_MANIFEST",
    os.path.join(BASE, "download_manifest.tsv"),
)
CHECKSUMS = BASE + "/checksums.tsv"
OUTDIR = BASE + "/freeze_a_redux"

SUBTYPE = os.environ.get(
    "TCGA_SUBTYPE_TABLE",
    "upstream/subtype-freeze/"
    "subtype_normalized.tsv",
)
SUBTYPE_SHA_EXPECT = "17a899a4869fccf806859342454bde440027786a48bdd50f1d1eace810b1ef9e"

GDC_RAW = os.environ.get(
    "TCGA_GDC_FILES_JSON",
    "upstream/acquisition-checkpoint/"
    "gdc_files_raw.json",
)
GDC_RAW_SHA_EXPECT = "6254d90e56dec110098a293b00ba69780fe443e52b4e88654340bc91824bc0cc"

# STAR-Counts augmented_star_gene_counts.tsv expected structure (GENCODE v36 / GDC DR32+)
EXPECT_DATA_COLS = ["gene_id", "gene_name", "gene_type", "unstranded",
                    "stranded_first", "stranded_second", "tpm_unstranded",
                    "fpkm_unstranded", "fpkm_uq_unstranded"]
EXPECT_SUMMARY_ROWS = ["N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous"]
# GENCODE v36 gene count in GDC augmented STAR files = 60660 gene rows
EXPECT_GENE_ROWS = 60660


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fail(msg):
    print("BLOCKED: " + msg, file=sys.stderr)
    sys.exit(2)


os.makedirs(OUTDIR, exist_ok=True)

# --- 5a. verify frozen subtype table SHA-256 before use ---
got = sha256(SUBTYPE)
if got != SUBTYPE_SHA_EXPECT:
    fail("subtype_normalized.tsv sha256 %s != expected %s" % (got, SUBTYPE_SHA_EXPECT))
print("subtype_normalized.tsv SHA-256 OK: %s" % got)

# provenance: also verify the raw GDC metadata used for the content-anchor cross-check
got_raw = sha256(GDC_RAW)
if got_raw != GDC_RAW_SHA_EXPECT:
    fail("gdc_files_raw.json sha256 %s != expected %s" % (got_raw, GDC_RAW_SHA_EXPECT))
print("gdc_files_raw.json SHA-256 OK: %s" % got_raw)

# --- load pinned manifest (uuid -> barcode / source_set / md5) ---
man = {}
for r in csv.DictReader(open(MANIFEST, newline=""), delimiter="\t"):
    man[r["file_uuid"]] = r

# --- load checksums (uuid -> md5_match, served GDC file_name) ---
chk = {}
for r in csv.DictReader(open(CHECKSUMS, newline=""), delimiter="\t"):
    chk[r["file_uuid"]] = r

# --- load raw GDC metadata (file_id -> canonical GDC file_name) ---
gdc_fname = {}
for rec in json.load(open(GDC_RAW)):
    gdc_fname[rec["file_id"]] = rec["file_name"]

# --- enumerate the ACTUAL on-disk files ---
ondisk = sorted(os.listdir(ORIG))
disk_uuids = []
for name in ondisk:
    m = re.match(r"^([0-9a-f-]{36})__", name)
    if not m:
        fail("unexpected on-disk filename (no uuid prefix): %s" % name)
    disk_uuids.append(m.group(1))

if len(disk_uuids) != 542:
    fail("on-disk file count %d != 542" % len(disk_uuids))
if len(set(disk_uuids)) != 542:
    fail("on-disk uuids not unique")

# --- confirm each on-disk uuid: in manifest, md5_match true, GDC-filename anchored ---
anchor_fail = []
md5_fail = []
for uuid in disk_uuids:
    if uuid not in man:
        fail("on-disk uuid %s not in pinned manifest" % uuid)
    if uuid not in chk or chk[uuid]["md5_match"] != "true":
        md5_fail.append(uuid)
        continue
    # content-anchor: served GDC filename == GDC-recorded file_name for this file_id
    served = chk[uuid]["file_name"]
    canonical = gdc_fname.get(uuid, "")
    if served != canonical:
        anchor_fail.append((uuid, served, canonical))

if md5_fail:
    fail("md5_match not true for %d on-disk files (e.g. %s)" % (len(md5_fail), md5_fail[:3]))
if anchor_fail:
    fail("GDC-filename anchor mismatch for %d files (e.g. %s)"
         % (len(anchor_fail), anchor_fail[:3]))
print("content-anchor OK: all 542 served GDC filenames == GDC-recorded file_name")

# --- split materialized set by source_set (from pinned manifest binding) ---
primary_uuids = [u for u in disk_uuids if man[u]["source_set"] == "primary"]
normal_uuids = [u for u in disk_uuids if man[u]["source_set"] == "normal"]
print("materialized: %d primary, %d normal" % (len(primary_uuids), len(normal_uuids)))

# ============================================================================
# 5b. Apply the FROZEN SELECTION RULE to the barcodes of the MATERIALIZED
#     PRIMARY files (not the metadata prediction).
# ============================================================================
# STEP 1: keep sample-type code == '01'. TCGA barcode: TCGA-XX-XXXX-NNz...
#   sample-type code = the 'NN' field = block index 3 (0-based) split on '-',
#   first two chars.
def sample_type_code(barcode):
    parts = barcode.split("-")
    if len(parts) < 4:
        return None
    return parts[3][:2]

def patient12(barcode):
    return barcode[:12]

# gather materialized primary aliquot barcodes
primary_records = []
for u in primary_uuids:
    bc = man[u]["aliquot_barcode"]
    primary_records.append((u, bc))

# STEP 1 keep '01'
step1 = [(u, bc) for (u, bc) in primary_records if sample_type_code(bc) == "01"]
excluded_not01 = [(u, bc, sample_type_code(bc)) for (u, bc) in primary_records
                  if sample_type_code(bc) != "01"]

# STEP 2 dedup: one per 12-char patient, keep lexicographically-lowest full aliquot barcode
by_patient = defaultdict(list)
for (u, bc) in step1:
    by_patient[patient12(bc)].append((u, bc))

kept = {}          # patient12 -> (uuid, barcode) kept
dropped_by_collapse = []   # (patient12, dropped_uuid, dropped_barcode, kept_barcode)
ties = []          # patients with >1 qualifying primary file
for pat, lst in by_patient.items():
    lst_sorted = sorted(lst, key=lambda x: x[1])  # lexicographic on full barcode
    kept[pat] = lst_sorted[0]
    if len(lst_sorted) > 1:
        keep_bc = lst_sorted[0][1]
        competing = [bc for (_, bc) in lst_sorted]
        ties.append((pat, keep_bc, competing))
        for (du, dbc) in lst_sorted[1:]:
            dropped_by_collapse.append((pat, du, dbc, keep_bc))

n_patients_multi = len(ties)
n_dropped_collapse = len(dropped_by_collapse)
selected_patients = set(kept.keys())  # 12-char patient barcodes surviving '01'+dedup
print("STEP1 kept '01': %d files (excluded non-01: %d)" % (len(step1), len(excluded_not01)))
print("STEP2 dedup: %d independent patients; %d patients with >1 primary; %d files dropped"
      % (len(selected_patients), n_patients_multi, n_dropped_collapse))

# ============================================================================
# 5c. Join (read-only) to the frozen subtype table.
# ============================================================================
subtype = {}   # patient12 -> mapped_4way
for r in csv.DictReader(open(SUBTYPE, newline=""), delimiter="\t"):
    subtype[r["patient_barcode"].strip()] = r["mapped_4way"].strip()
subtype_patients = set(subtype.keys())

matched = sorted(selected_patients & subtype_patients)
primary_no_subtype = sorted(selected_patients - subtype_patients)
subtype_no_primary = sorted(subtype_patients - selected_patients)

arm_counts = defaultdict(int)
for pat in matched:
    arm_counts[subtype[pat]] += 1

# ============================================================================
# 5d. STRUCTURAL integrity of each of the 507 materialized primary files.
#     (validate the 507 selected primary files -- gate set. Also validate ALL
#     507 materialized primary files regardless of dedup, so the count is on the
#     full primary substrate; here selected == all 507 since no primary was
#     dropped by collapse in the materialized set... but we validate every
#     materialized primary file to be safe.)
#     STRUCTURAL ONLY -- no tpm value scored.
# ============================================================================
def validate_star_file(path):
    """Return dict with structural checks; never scores a value."""
    with open(path, "r") as f:
        first = f.readline().rstrip("\n")
        header = f.readline().rstrip("\n").split("\t")
        # next 4 rows are the STAR summary rows
        summ = []
        for _ in range(4):
            summ.append(f.readline().rstrip("\n").split("\t")[0])
        gene_rows = 0
        for line in f:
            if line.strip():
                gene_rows += 1
    cols_ok = (header == EXPECT_DATA_COLS)
    summ_ok = (summ == EXPECT_SUMMARY_ROWS)
    gene_ok = (gene_rows == EXPECT_GENE_ROWS)
    gencode_ok = first.startswith("# gene-model: GENCODE v36")
    return {
        "gencode_header": gencode_ok,
        "columns_ok": cols_ok,
        "summary_rows_ok": summ_ok,
        "gene_rows": gene_rows,
        "gene_rows_ok": gene_ok,
        "all_ok": cols_ok and summ_ok and gene_ok and gencode_ok,
    }

struct_rows = []
struct_fail = []
for u in primary_uuids:
    path = os.path.join(ORIG, u + "__augmented_star_gene_counts.tsv")
    v = validate_star_file(path)
    v["file_uuid"] = u
    v["aliquot_barcode"] = man[u]["aliquot_barcode"]
    struct_rows.append(v)
    if not v["all_ok"]:
        struct_fail.append((u, v))

n_struct_ok = sum(1 for v in struct_rows if v["all_ok"])
gene_row_values = sorted(set(v["gene_rows"] for v in struct_rows))
print("STRUCTURAL: %d/%d primary files valid STAR/GENCODE-v36; distinct gene-row counts: %s"
      % (n_struct_ok, len(struct_rows), gene_row_values))

# ============================================================================
# 5e. Gate band.
# ============================================================================
MANDATORY = ["POLE", "MMRd", "NSMP", "p53abn"]
arm_vals = {a: arm_counts.get(a, 0) for a in MANDATORY}
smallest = min(arm_vals.values())
if smallest >= 40:
    band = "GREEN"
elif smallest >= 20:
    band = "YELLOW"
else:
    band = "RED"
pole = arm_vals["POLE"]

# POLE >= 20 confirmation AFTER preprocessing: every counted POLE patient must have
# a materialized, md5-verified, GENCODE-v36-valid primary file surviving '01'+dedup+join.
pole_patients = [p for p in matched if subtype[p] == "POLE"]
pole_all_valid = True
pole_detail = []
struct_by_uuid = {v["file_uuid"]: v for v in struct_rows}
for p in pole_patients:
    u, bc = kept[p]
    v = struct_by_uuid.get(u, {})
    ok = (chk[u]["md5_match"] == "true") and v.get("all_ok", False)
    pole_detail.append({"patient": p, "uuid": u, "barcode": bc,
                        "md5_match": chk[u]["md5_match"],
                        "struct_ok": v.get("all_ok", False)})
    if not ok:
        pole_all_valid = False

# ============================================================================
# 5f. The 35 normals -- verify + preserve only. Record case barcodes and note a
#     subtype label if the SAME patient has one (informational only; NOT in gate).
# ============================================================================
normal_records = []
for u in normal_uuids:
    bc = man[u]["aliquot_barcode"]
    pat = patient12(bc)
    normal_records.append({
        "file_uuid": u,
        "case_barcode": man[u]["case_barcode"],
        "aliquot_barcode": bc,
        "sample_type": man[u]["sample_type"],
        "md5_match": chk[u]["md5_match"],
        "same_patient_subtype_label": subtype.get(pat, "NA"),
    })

# ============================================================================
# WRITE ARTIFACTS
# ============================================================================
def write_tsv(path, header, rows):
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(header)
        for r in rows:
            w.writerow(r)

write_tsv(OUTDIR + "/cohort_selected_primary.tsv",
          ["patient_barcode", "kept_uuid", "kept_aliquot_barcode", "mapped_4way"],
          [[p, kept[p][0], kept[p][1], subtype.get(p, "NA")] for p in sorted(selected_patients)])

write_tsv(OUTDIR + "/join_matched.tsv",
          ["patient_barcode", "mapped_4way", "kept_uuid", "kept_aliquot_barcode"],
          [[p, subtype[p], kept[p][0], kept[p][1]] for p in matched])

write_tsv(OUTDIR + "/primary_no_subtype.tsv",
          ["patient_barcode", "kept_uuid", "kept_aliquot_barcode"],
          [[p, kept[p][0], kept[p][1]] for p in primary_no_subtype])

write_tsv(OUTDIR + "/subtype_no_primary.tsv",
          ["patient_barcode", "mapped_4way"],
          [[p, subtype[p]] for p in subtype_no_primary])

write_tsv(OUTDIR + "/dedup_ties.tsv",
          ["patient_barcode", "kept_barcode", "competing_barcodes"],
          [[p, keep, ";".join(comp)] for (p, keep, comp) in ties])

write_tsv(OUTDIR + "/dropped_by_collapse.tsv",
          ["patient_barcode", "dropped_uuid", "dropped_barcode", "kept_barcode"],
          dropped_by_collapse)

write_tsv(OUTDIR + "/excluded_not_primary01.tsv",
          ["file_uuid", "aliquot_barcode", "sample_type_code"],
          excluded_not01)

write_tsv(OUTDIR + "/structural_integrity_primary.tsv",
          ["file_uuid", "aliquot_barcode", "gencode_header", "columns_ok",
           "summary_rows_ok", "gene_rows", "gene_rows_ok", "all_ok"],
          [[v["file_uuid"], v["aliquot_barcode"], v["gencode_header"], v["columns_ok"],
            v["summary_rows_ok"], v["gene_rows"], v["gene_rows_ok"], v["all_ok"]]
           for v in struct_rows])

with open(OUTDIR + "/normals_preserved.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["file_uuid", "case_barcode", "aliquot_barcode",
                                      "sample_type", "md5_match",
                                      "same_patient_subtype_label"], delimiter="\t")
    w.writeheader()
    for r in normal_records:
        w.writerow(r)

summary = {
    "materialized_primary_files": len(primary_uuids),
    "materialized_normal_files": len(normal_uuids),
    "step1_kept_01": len(step1),
    "excluded_not_primary01": len(excluded_not01),
    "independent_primary_patients": len(selected_patients),
    "patients_with_multiple_primary_rna": n_patients_multi,
    "files_dropped_by_collapse": n_dropped_collapse,
    "join_matched": len(matched),
    "join_primary_no_subtype": len(primary_no_subtype),
    "join_subtype_no_primary": len(subtype_no_primary),
    "arm_counts": {a: arm_vals[a] for a in MANDATORY},
    "smallest_mandatory_arm": smallest,
    "pole_count": pole,
    "pole_ge_20": pole >= 20,
    "pole_all_files_md5_and_struct_valid": pole_all_valid,
    "gate_band": band,
    "structural_primary_valid": n_struct_ok,
    "structural_primary_total": len(struct_rows),
    "distinct_gene_row_counts": gene_row_values,
    "normals_verified_preserved": len(normal_records),
    "normals_all_md5_true": all(r["md5_match"] == "true" for r in normal_records),
}
json.dump(summary, open(OUTDIR + "/summary.json", "w"), indent=2, sort_keys=True)
json.dump(pole_detail, open(OUTDIR + "/pole_patient_detail.json", "w"), indent=2)

# ============================================================================
# HUMAN SUMMARY
# ============================================================================
print("")
print("=================== FREEZE A REDUX -- DERIVED (no number forced) ===================")
print("materialized primary files      : %d" % len(primary_uuids))
print("materialized normal files       : %d" % len(normal_uuids))
print("STEP1 keep '01'                 : %d files (excluded non-01: %d)"
      % (len(step1), len(excluded_not01)))
print("independent primary patients    : %d" % len(selected_patients))
print("  patients with >1 primary RNA  : %d (files dropped by collapse: %d)"
      % (n_patients_multi, n_dropped_collapse))
print("join to frozen subtype table:")
print("  matched                       : %d" % len(matched))
print("  primary-no-subtype            : %d" % len(primary_no_subtype))
print("  subtype-no-primary            : %d" % len(subtype_no_primary))
print("per-subtype (matched):")
for a in MANDATORY:
    print("  %-7s : %d" % (a, arm_vals[a]))
print("smallest mandatory arm          : %d" % smallest)
print("POLE count                      : %d   (>=20 ? %s)" % (pole, pole >= 20))
print("POLE all files md5+struct valid : %s" % pole_all_valid)
print("GATE BAND                       : %s" % band)
print("structural primary valid        : %d/%d (gene-row counts seen: %s)"
      % (n_struct_ok, len(struct_rows), gene_row_values))
print("normals verified+preserved      : %d (all md5 true: %s)"
      % (len(normal_records), all(r["md5_match"] == "true" for r in normal_records)))

if struct_fail:
    print("")
    print("BLOCKED: %d primary files FAILED structural validation:" % len(struct_fail))
    for u, v in struct_fail[:10]:
        print("  %s -> %s" % (u, v))
    sys.exit(6)
if not pole_all_valid:
    print("")
    print("BLOCKED: at least one counted POLE patient lacks a valid/verified primary file.")
    sys.exit(7)

print("")
print("WROTE artifacts under %s" % OUTDIR)
