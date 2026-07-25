#!/usr/bin/env python3
"""
TASK-029 ACQUISITION step 02 -- deterministic patient->sample->RNA join + dedup.

From the authoritative GDC STAR-Counts response (acq01), for each CPTAC-3 uterine
case build the set of Primary-Tumour RNA files, then apply the frozen dedup rule k
(one patient = one primary-tumour analytic unit; deterministic tie-break) to select
EXACTLY ONE RNA file per joined case.

Join target = the sealed 230 subtyped cases (Discovery 95 + Confirmatory 135) from
the SHA-sealed join tables. Reports:
  - sample_type / tissue_type census (to separate tumour from normal-adjacent),
  - per-case candidate primary-tumour RNA files,
  - the deterministic single pick per case + tie-break audit,
  - cross-stratum overlap check (0 required),
  - any case with 0 primary-tumour RNA files (fatal attrition) or subtype conflict.

Metadata only. Deterministic. Git HEAD unchanged.
"""
import json
import os
from collections import defaultdict, Counter

BASE = os.path.dirname(
    os.path.dirname(os.path.abspath(__file__))
)
WORK = os.environ.get(
    "CPTAC_WORK_DIR",
    os.path.join(BASE, "work"),
)
SRC = os.environ.get(
    "CPTAC_ACQUISITION_OUTPUT",
    os.path.join(WORK, "acquisition"),
)
JOINT = os.environ.get(
    "CPTAC_DATA_LINKAGE_DIR",
    os.path.join(BASE, "data-linkage"),
)
OUT = os.environ.get(
    "CPTAC_INTERMEDIATE_DIR",
    os.path.join(WORK, "intermediate"),
)
os.makedirs(OUT, exist_ok=True)

resp = json.load(open(os.path.join(SRC, "acq01_gdc_star_counts_raw.json")))
hits = resp["data"]["hits"]

disc = json.load(open(os.path.join(JOINT, "discovery_join.json")))["case_to_subtype_joined"]
conf = json.load(open(os.path.join(JOINT, "confirmatory_join.json")))["case_to_subtype_joined"]
subtype_of = {}
stratum_of = {}
for c, s in disc.items():
    subtype_of[c] = s
    stratum_of[c] = "Discovery"
for c, s in conf.items():
    if c in subtype_of:
        raise SystemExit(f"FATAL: case {c} in BOTH strata rosters")
    subtype_of[c] = s
    stratum_of[c] = "Confirmatory"
needed = set(subtype_of)
assert len(needed) == 230, len(needed)

# ---- flatten hits: one row per (file, case, sample) ----
# a STAR-Counts file is per-aliquot; each hit carries case(s)+sample(s).
rows = []
sample_type_census = Counter()
tissue_type_census = Counter()
for h in hits:
    fid = h["file_id"]; fname = h["file_name"]
    md5 = h.get("md5sum"); size = h.get("file_size"); state = h.get("state")
    access = h.get("access")
    for case in h.get("cases", []):
        cid = case.get("submitter_id")
        for samp in case.get("samples", []):
            st = samp.get("sample_type")
            tt = samp.get("tissue_type")
            ssid = samp.get("submitter_id")
            sample_type_census[st] += 1
            tissue_type_census[tt] += 1
            # aliquot ids
            aliquots = []
            for portion in samp.get("portions", []) or []:
                for analyte in portion.get("analytes", []) or []:
                    for al in analyte.get("aliquots", []) or []:
                        aliquots.append(al.get("submitter_id"))
            rows.append({
                "file_id": fid, "file_name": fname, "md5sum": md5,
                "file_size": size, "state": state, "access": access,
                "case_submitter_id": cid, "sample_submitter_id": ssid,
                "sample_type": st, "tissue_type": tt,
                "aliquot_submitter_ids": aliquots,
            })

print("total (file,case,sample) rows:", len(rows))
print("sample_type census (all uterine STAR-Counts):", dict(sample_type_census))
print("tissue_type census:", dict(tissue_type_census))

# ---- restrict to the 230 needed cases + PRIMARY TUMOUR samples ----
# Primary tumour = sample_type == 'Primary Tumor' AND tissue_type == 'Tumor'.
# Normal-adjacent (sample_type contains 'Normal' or tissue_type=='Normal') EXCLUDED.
def is_primary_tumor(r):
    return (r["sample_type"] == "Primary Tumor") and (r["tissue_type"] == "Tumor")

by_case_primary = defaultdict(list)
by_case_any = defaultdict(list)
excluded_nonprimary = []
for r in rows:
    if r["case_submitter_id"] not in needed:
        continue
    by_case_any[r["case_submitter_id"]].append(r)
    if is_primary_tumor(r):
        by_case_primary[r["case_submitter_id"]].append(r)
    else:
        excluded_nonprimary.append(r)

# ---- per-case: dedup to ONE primary-tumour RNA file ----
# Deterministic tie-break (frozen, outcome-blind, no expression seen):
#   1. must be state=='released' and access=='open';
#   2. prefer the LOWEST sample_submitter_id suffix number (the study-primary /
#      first-collected vial: e.g. -01 before -02; -02 before -12), matching the
#      "keep study-primary aliquot" rule (lower vial = primary tube);
#   3. then LOWEST sample_submitter_id lexicographically;
#   4. then LOWEST file_id lexicographically (final deterministic tie-break).
# All ties + drops are logged.
def suffix_num(ssid):
    # CPTAC sample submitter id: C3L-XXXXX-NN ; NN = vial/sample code
    try:
        return int(ssid.rsplit("-", 1)[-1])
    except Exception:
        return 9999

pick = {}
dedup_audit = {}
zero_primary = []
for c in sorted(needed):
    cands = by_case_primary.get(c, [])
    # dedup identical file_ids (the recon cache had duplicate rows)
    seen = {}
    for r in cands:
        seen[r["file_id"]] = r
    cands = list(seen.values())
    if not cands:
        zero_primary.append(c)
        continue
    ranked = sorted(cands, key=lambda r: (
        0 if (r["state"] == "released" and r["access"] == "open") else 1,
        suffix_num(r["sample_submitter_id"]),
        r["sample_submitter_id"],
        r["file_id"],
    ))
    chosen = ranked[0]
    pick[c] = chosen
    dedup_audit[c] = {
        "n_primary_candidates_distinct_files": len(cands),
        "chosen_file_id": chosen["file_id"],
        "chosen_sample": chosen["sample_submitter_id"],
        "chosen_sample_type": chosen["sample_type"],
        "chosen_tissue_type": chosen["tissue_type"],
        "all_candidate_samples": sorted(r["sample_submitter_id"] for r in cands),
        "all_candidate_file_ids": sorted(r["file_id"] for r in cands),
    }

print("\ncases with a primary-tumour RNA pick:", len(pick), "/ 230")
print("cases with ZERO primary-tumour RNA file:", len(zero_primary), zero_primary)

# ---- cross-stratum overlap check at case / sample / aliquot / file ----
disc_cases = set(disc); conf_cases = set(conf)
case_overlap = disc_cases & conf_cases
disc_files = {pick[c]["file_id"] for c in disc if c in pick}
conf_files = {pick[c]["file_id"] for c in conf if c in pick}
file_overlap = disc_files & conf_files
disc_samples = {pick[c]["sample_submitter_id"] for c in disc if c in pick}
conf_samples = {pick[c]["sample_submitter_id"] for c in conf if c in pick}
sample_overlap = disc_samples & conf_samples
disc_aliq = set()
conf_aliq = set()
for c in disc:
    if c in pick:
        disc_aliq |= set(pick[c]["aliquot_submitter_ids"])
for c in conf:
    if c in pick:
        conf_aliq |= set(pick[c]["aliquot_submitter_ids"])
aliq_overlap = disc_aliq & conf_aliq
print("\ncross-stratum overlap: case=%d sample=%d aliquot=%d file=%d" % (
    len(case_overlap), len(sample_overlap), len(aliq_overlap), len(file_overlap)))

# ---- per-stratum subtype counts after join ----
sub_ct = {"Discovery": Counter(), "Confirmatory": Counter()}
for c in pick:
    sub_ct[stratum_of[c]][subtype_of[c]] += 1
print("\nDiscovery subtype counts (picked):", dict(sub_ct["Discovery"]))
print("Confirmatory subtype counts (picked):", dict(sub_ct["Confirmatory"]))

# ---- multiplicity summary among the 230 ----
mult = Counter(len(dedup_audit[c]["all_candidate_file_ids"]) for c in dedup_audit)
print("\nprimary-tumour candidate-files-per-case distribution:", dict(mult))

# ---- write join ledger ----
join_out = {
    "gdc_data_release": resp.get("data", {}).get("release"),
    "n_needed_cases": 230,
    "n_picked": len(pick),
    "zero_primary_cases": zero_primary,
    "cross_stratum_overlap": {
        "case": sorted(case_overlap), "sample": sorted(sample_overlap),
        "aliquot": sorted(aliq_overlap), "file": sorted(file_overlap)},
    "discovery_subtype_counts": dict(sub_ct["Discovery"]),
    "confirmatory_subtype_counts": dict(sub_ct["Confirmatory"]),
    "candidate_files_per_case_distribution": {str(k): v for k, v in mult.items()},
    "pick": {c: {
        "stratum": stratum_of[c], "subtype": subtype_of[c],
        "file_id": pick[c]["file_id"], "file_name": pick[c]["file_name"],
        "md5sum": pick[c]["md5sum"], "file_size": pick[c]["file_size"],
        "state": pick[c]["state"], "access": pick[c]["access"],
        "sample_submitter_id": pick[c]["sample_submitter_id"],
        "sample_type": pick[c]["sample_type"], "tissue_type": pick[c]["tissue_type"],
        "aliquot_submitter_ids": pick[c]["aliquot_submitter_ids"],
    } for c in sorted(pick)},
    "dedup_audit": dedup_audit,
    "excluded_nonprimary_n": len(excluded_nonprimary),
    "sample_type_census_all": dict(sample_type_census),
    "tissue_type_census_all": dict(tissue_type_census),
}
os.makedirs(OUT, exist_ok=True)
with open(os.path.join(OUT, "acq02_join_ledger.json"), "w") as f:
    json.dump(join_out, f, indent=2)
print("\nwrote", os.path.join(OUT, "acq02_join_ledger.json"))
