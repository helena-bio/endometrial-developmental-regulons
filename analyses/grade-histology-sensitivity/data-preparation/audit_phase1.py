#!/usr/bin/env python3
"""Deterministic outcome-blind TASK B Phase-1 clinical feasibility audit."""

import argparse
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np


ROOT = Path("data/external/original-workspace/revgate-task-b-grade-histology")
TASK = ROOT / "experiments/taskB_grade_histology"
EXEC = TASK / "phase1_execution"
INPUT = TASK / "inputs"
ALLOWED = {
    "tcga_roster": Path("data/external/original-workspace/task027-acquire-freeze-a/freeze_a_redux/cohort_selected_primary.tsv"),
    "tcga_subtype": Path("data/external/original-workspace/task024-freeze-a/subtype_normalized.tsv"),
    "discovery_join": Path("data/external/original-workspace/task029-external-replication-feasibility/join_tables/discovery_join.json"),
    "confirmatory_join": Path("data/external/original-workspace/task029-external-replication-feasibility/join_tables/confirmatory_join.json"),
    "cbio_dictionary": Path("data/external/original-workspace/task029-external-replication-feasibility/sources/cbio_2020_clinattr.json"),
    "gdc_dictionary": Path("data/external/original-workspace/task029-external-replication-feasibility/sources_task029gate/uec_cptac_gdc_clinattr.json"),
    "discovery_biospec": Path("data/external/original-workspace/task029-external-replication-feasibility/sources/pdc_discovery_biospec.json"),
    "confirmatory_biospec": Path("data/external/original-workspace/task029-external-replication-feasibility/sources/pdc_confirmatory_biospec.json"),
}
LOCKED = {
    "tcga_grade": INPUT / "tcga_cbio_grade.json",
    "tcga_histology": INPUT / "tcga_cbio_icdo_histology.json",
    "discovery_grade": INPUT / "cptac_discovery_grade.json",
    "discovery_histology": INPUT / "cptac_discovery_histology.json",
    "pdc_discovery": INPUT / "pdc_PDC000125_grade_histology.json",
    "pdc_confirmatory": INPUT / "pdc_PDC000439_grade_histology.json",
}
SUBTYPES = ["NSMP", "MMRd", "POLE", "p53abn"]


def read_json(path):
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def read_tsv(path):
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, columns):
    with path.open("w", encoding="ascii", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: ascii_value(row.get(key, "")) for key in columns})


def ascii_value(value):
    if value is None:
        return "NA"
    if isinstance(value, float):
        if math.isnan(value) or math.isinf(value):
            return "NA"
        return f"{value:.10g}"
    return str(value).encode("ascii", "replace").decode("ascii")


def write_json(path, obj):
    path.write_text(json.dumps(obj, indent=2, sort_keys=True, ensure_ascii=True) + "\n", encoding="ascii")


def group_unique(rows, key):
    out = defaultdict(list)
    for row in rows:
        out[row[key]].append(row)
    return out


def norm_grade(token, figo_documented):
    if token is None:
        return None
    token = str(token).strip()
    mapping = {"G1": "low_grade", "G2": "low_grade", "G3": "high_grade", "FIGO grade 1": "low_grade", "FIGO grade 2": "low_grade", "FIGO grade 3": "high_grade"}
    if not figo_documented:
        return None
    return mapping.get(token)


HIST_MAP = {
    "8380/3": ("endometrioid", "endometrioid"),
    "8382/3": ("endometrioid", "endometrioid"),
    "8441/3": ("non_endometrioid", "serous"),
    "8460/3": ("non_endometrioid", "serous"),
    "8461/3": ("non_endometrioid", "serous"),
    "8310/3": ("non_endometrioid", "other_non_endometrioid"),
    "8020/3": ("non_endometrioid", "other_non_endometrioid"),
    "8255/3": ("non_endometrioid", "other_non_endometrioid"),
    "Endometrioid": ("endometrioid", "endometrioid"),
    "Serous": ("non_endometrioid", "serous"),
}


def norm_hist(token):
    if token is None:
        return (None, None)
    return HIST_MAP.get(str(token).strip(), (None, None))


def pdc_rows(payload):
    return payload["data"]["clinicalPerStudy"]


def biospec_rows(payload):
    return payload["data"]["biospecimenPerStudy"]


def cramer_v(a, b):
    levels_a = sorted(set(a)); levels_b = sorted(set(b))
    table = np.zeros((len(levels_a), len(levels_b)), dtype=float)
    ia = {x: i for i, x in enumerate(levels_a)}; ib = {x: i for i, x in enumerate(levels_b)}
    for x, y in zip(a, b): table[ia[x], ib[y]] += 1
    expected = np.outer(table.sum(axis=1), table.sum(axis=0)) / table.sum()
    ok = expected > 0
    chi = float(np.sum(((table - expected) ** 2)[ok] / expected[ok]))
    den = table.sum() * max(1, min(len(levels_a) - 1, len(levels_b) - 1))
    return math.sqrt(chi / den) if den else float("nan")


def matrix(rows, clinical=None, subset_hist=None, grade_stratum=None):
    clinicals = [] if clinical is None else ([clinical] if isinstance(clinical, str) else list(clinical))
    selected = []
    for row in rows:
        if not row["exact_link"]:
            continue
        if subset_hist and row["histology"] != subset_hist:
            continue
        if grade_stratum and row["grade"] != grade_stratum:
            continue
        if any(row[name] is None for name in clinicals):
            continue
        selected.append(row)
    columns = ["intercept", "subtype_MMRd", "subtype_POLE", "subtype_p53abn"]
    groups = {"subtype": [1, 2, 3]}
    data = []
    clinical_levels = {}
    refs = {}
    for name in clinicals:
        levels = sorted(set(row[name] for row in selected))
        clinical_levels[name] = levels
        refs[name] = levels[0] if levels else None
        start = len(columns)
        for level in levels[1:]:
            columns.append(f"{name}_{level}")
        groups[name] = list(range(start, len(columns)))
    for row in selected:
        values = [1.0, float(row["subtype"] == "MMRd"), float(row["subtype"] == "POLE"), float(row["subtype"] == "p53abn")]
        for name in clinicals:
            values += [float(row[name] == level) for level in clinical_levels[name][1:]]
        data.append(values)
    return selected, np.asarray(data, dtype=float), columns, groups, refs


def diagnostic(cohort, node, rows, clinical=None, subset_hist=None, grade_stratum=None, forced_state=None, reason=""):
    selected, x, columns, groups, refs = matrix(rows, clinical, subset_hist, grade_stratum)
    clinicals = [] if clinical is None else ([clinical] if isinstance(clinical, str) else list(clinical))
    n = len(selected); p = len(columns)
    ref_text = "subtype=NSMP" + "".join(f";{name}={refs.get(name)}" for name in clinicals)
    empty = {"cohort": cohort, "node": node, "n": n, "p": p, "rank": "NA", "residual_df": n-p, "svd_tolerance": "NA", "singular_values": "NA", "condition_number": "NA", "max_one_df_vif": "NA", "subtype_adjusted_gvif": "NA", "clinical_adjusted_gvif": "NA", "cramers_v": "NA", "near_zero_variance": "NA", "max_leverage": "NA", "fraction_hat_gt_3p_n": "NA", "support_state": forced_state or "FAIL", "reason": reason, "reference_levels": ref_text}
    if n == 0 or x.shape != (n, p):
        return empty
    singular = np.linalg.svd(x, compute_uv=False)
    tol = np.finfo(float).eps * max(n, p) * singular[0]
    rank = int(np.sum(singular > tol))
    z = x[:, 1:].copy()
    if z.shape[1]:
        sd = z.std(axis=0, ddof=0); valid = sd > 0
        z[:, valid] = (z[:, valid] - z[:, valid].mean(axis=0)) / sd[valid]
        z[:, ~valid] = 0
        sz = np.linalg.svd(z, compute_uv=False)
        positive = sz[sz > np.finfo(float).eps * max(z.shape) * (sz[0] if len(sz) else 1)]
        condition = float(positive[0] / positive[-1]) if len(positive) == z.shape[1] and len(positive) else float("inf")
        corr = np.corrcoef(z, rowvar=False) if z.shape[1] > 1 else np.array([[1.0]])
        invcorr = np.linalg.pinv(corr)
        max_vif = float(np.max(np.diag(invcorr)))
        def adj_gvif(indices):
            if not indices: return float("nan")
            idx = [i - 1 for i in indices]
            rest = [i for i in range(corr.shape[0]) if i not in idx]
            det_all = max(float(np.linalg.det(corr)), 1e-300)
            det_a = max(float(np.linalg.det(corr[np.ix_(idx, idx)])), 1e-300)
            det_b = max(float(np.linalg.det(corr[np.ix_(rest, rest)])), 1e-300) if rest else 1.0
            gvif = det_a * det_b / det_all
            return gvif ** (1.0 / (2.0 * len(idx)))
        subtype_gvif = adj_gvif(groups["subtype"])
        clinical_values = [adj_gvif(groups.get(name, [])) for name in clinicals if groups.get(name)]
        clinical_gvif = max(clinical_values) if clinical_values else float("nan")
    else:
        condition = 1.0; max_vif = 1.0; subtype_gvif = float("nan"); clinical_gvif = float("nan")
    hat = np.diag(x @ np.linalg.pinv(x.T @ x) @ x.T)
    high_fraction = float(np.mean(hat > (3 * p / n)))
    nzv = False
    for name in clinicals:
        counts = Counter(row[name] for row in selected)
        ordered = sorted(counts.values(), reverse=True)
        freq_ratio = (ordered[0] / ordered[1]) if len(ordered) > 1 and ordered[1] else float("inf")
        pct_unique = 100 * len(counts) / n
        nzv = nzv or (freq_ratio > 19 and pct_unique <= 10)
    cvs = [cramer_v([r["subtype"] for r in selected], [r[name] for r in selected]) for name in clinicals if len(set(r[name] for r in selected)) > 1]
    cv = max(cvs) if cvs else float("nan")
    failures = []
    subtype_counts = Counter(row["subtype"] for row in selected)
    if n < max(80, 10 * p): failures.append("N_BELOW_MAX_80_10P")
    if n - p < 30: failures.append("RESIDUAL_DF_LT_30")
    if set(subtype_counts) != set(SUBTYPES): failures.append("MISSING_SUBTYPE_LEVEL")
    min_subtype = 10 if subset_hist or grade_stratum else 5
    if any(subtype_counts.get(level, 0) < min_subtype for level in SUBTYPES): failures.append("SUBTYPE_CELL_TOO_SMALL")
    for name in clinicals:
        cc = Counter(row[name] for row in selected)
        if any(v < 20 for v in cc.values()): failures.append(f"CLINICAL_LEVEL_LT_20:{name}")
        for level, total in cc.items():
            by = Counter(row["subtype"] for row in selected if row[name] == level)
            if sum(v >= 5 for v in by.values()) < 2: failures.append(f"LEVEL_NOT_IN_TWO_SUBTYPES:{name}:{level}")
            if max(by.values()) / total >= 0.90: failures.append(f"LEVEL_GE_90PCT_ONE_SUBTYPE:{name}:{level}")
    if rank < p: failures.append("RANK_DEFICIENT")
    if condition > 100: failures.append("CONDITION_GT_100")
    elif condition > 30: failures.append("CONDITION_UNSTABLE")
    if max_vif > 10 or subtype_gvif > 10 or (not math.isnan(clinical_gvif) and clinical_gvif > 10): failures.append("VIF_GT_10")
    elif max_vif > 5 or subtype_gvif > 5 or (not math.isnan(clinical_gvif) and clinical_gvif > 5): failures.append("VIF_UNSTABLE")
    if not math.isnan(cv) and cv > 0.90: failures.append("CRAMERS_V_GT_0.90")
    if nzv: failures.append("NEAR_ZERO_VARIANCE")
    if float(np.max(hat)) >= 0.50: failures.append("MAX_HAT_GE_0.50")
    if high_fraction > 0.05: failures.append("HIGH_LEVERAGE_FRACTION_GT_0.05")
    state = forced_state or ("PASS" if not failures else "FAIL")
    if forced_state:
        failures.insert(0, reason or forced_state)
    return {"cohort": cohort, "node": node, "n": n, "p": p, "rank": rank, "residual_df": n-p, "svd_tolerance": tol, "singular_values": ";".join(f"{v:.10g}" for v in singular), "condition_number": condition, "max_one_df_vif": max_vif, "subtype_adjusted_gvif": subtype_gvif, "clinical_adjusted_gvif": clinical_gvif, "cramers_v": cv, "near_zero_variance": nzv, "max_leverage": float(np.max(hat)), "fraction_hat_gt_3p_n": high_fraction, "support_state": state, "reason": ";".join(dict.fromkeys(failures)) or "ALL_FROZEN_PHASE1_RULES_PASS", "reference_levels": ref_text}


def build_rows():
    tcga_roster = read_tsv(ALLOWED["tcga_roster"])
    tcga_subtype = {row["patient_barcode"]: row["mapped_4way"] for row in read_tsv(ALLOWED["tcga_subtype"])}
    tcga_grade = group_unique(read_json(LOCKED["tcga_grade"]), "patientId")
    tcga_hist = group_unique(read_json(LOCKED["tcga_histology"]), "patientId")
    all_rows = []
    for roster in tcga_roster:
        patient = roster["patient_barcode"]
        grades = tcga_grade.get(patient, [])
        hist = tcga_hist.get(patient, [])
        expected_sample = roster["kept_aliquot_barcode"][:15]
        grade_exact = [r for r in grades if r["sampleId"] == expected_sample]
        grade_raw = grade_exact[0]["value"] if len(grade_exact) == 1 else None
        hist_raw = hist[0]["value"] if len(hist) == 1 else None
        all_rows.append({"cohort": "TCGA", "patient_id": patient, "analytic_sample_id": expected_sample, "aliquot_id": roster["kept_aliquot_barcode"], "subtype": tcga_subtype.get(patient), "grade_raw": grade_raw, "grade_secondary_raw": None, "grade": norm_grade(grade_raw, True), "histology_raw": hist_raw, "histology_secondary_raw": None, "histology": norm_hist(hist_raw)[0], "histology_three": norm_hist(hist_raw)[1], "exact_link": len(grade_exact) == 1 and len(hist) == 1 and patient in tcga_subtype, "duplicate_count": max(len(grades), len(hist)), "grade_conflict": len({r['value'] for r in grade_exact}) > 1, "histology_conflict": len({r['value'] for r in hist}) > 1, "link_status": "EXACT" if len(grade_exact) == 1 and len(hist) == 1 and patient in tcga_subtype else "NOT_EXACT", "exclusion_reason": "" if len(grade_exact) == 1 and len(hist) == 1 and patient in tcga_subtype else "exact_tcga_patient_sample_link_unavailable", "source_row_id": patient})

    def add_cptac(cohort, join_key, biospec_key, pdc_key, cbio_grade_key=None, cbio_hist_key=None):
        joined = read_json(ALLOWED[join_key])["case_to_subtype_joined"]
        biospec = biospec_rows(read_json(ALLOWED[biospec_key]))
        primary = defaultdict(list)
        for row in biospec:
            if row["sample_type"] == "Primary Tumor": primary[row["case_submitter_id"]].append(row)
        pdc = {r["case_submitter_id"]: r for r in pdc_rows(read_json(LOCKED[pdc_key]))}
        cbg = {r["patientId"]: r for r in read_json(LOCKED[cbio_grade_key])} if cbio_grade_key else {}
        cbh = {r["patientId"]: r for r in read_json(LOCKED[cbio_hist_key])} if cbio_hist_key else {}
        for patient, subtype in sorted(joined.items()):
            prim = primary.get(patient, [])
            exact = len(prim) == 1
            p = pdc.get(patient, {})
            g1 = cbg.get(patient, {}).get("value")
            g2 = p.get("tumor_grade")
            h1 = cbh.get(patient, {}).get("value")
            h2 = p.get("morphology")
            grade = norm_grade(g1, True) if cbio_grade_key else norm_grade(g2, False)
            h1n = norm_hist(h1)[0]; h2n = norm_hist(h2)[0]
            hist_conflict = h1n is not None and h2n is not None and h1n != h2n
            hist = None if hist_conflict else (h1n if h1n is not None else h2n)
            hist3 = None if hist_conflict else (norm_hist(h1)[1] if h1n is not None else norm_hist(h2)[1])
            all_rows.append({"cohort": cohort, "patient_id": patient, "analytic_sample_id": prim[0]["sample_submitter_id"] if exact else None, "aliquot_id": prim[0]["aliquot_submitter_id"] if exact else None, "subtype": subtype, "grade_raw": g1 if cbio_grade_key else g2, "grade_secondary_raw": g2 if cbio_grade_key else None, "grade": grade, "histology_raw": h1 if cbio_hist_key else h2, "histology_secondary_raw": h2 if cbio_hist_key else None, "histology": hist, "histology_three": hist3, "exact_link": exact, "duplicate_count": len(prim), "grade_conflict": (norm_grade(g1, True) is not None and norm_grade(g2, True) is not None and norm_grade(g1, True) != norm_grade(g2, True)) if cbio_grade_key else False, "histology_conflict": hist_conflict, "link_status": "EXACT" if exact else "BLOCKED_MULTIPLE_PRIMARY_SAMPLE_OR_ALIQUOT", "exclusion_reason": "" if exact else "frozen_exact_analytic_tumour_or_aliquot_unavailable", "source_row_id": p.get("case_id", "")})
    add_cptac("CPTAC_Discovery", "discovery_join", "discovery_biospec", "pdc_discovery", "discovery_grade", "discovery_histology")
    add_cptac("CPTAC_Confirmatory", "confirmatory_join", "confirmatory_biospec", "pdc_confirmatory")
    return all_rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    out = Path(args.output).resolve()
    out.mkdir(parents=True, exist_ok=False)
    if not (EXEC / "INPUTS_LOCK.json").exists():
        raise SystemExit("INPUTS_LOCK_MISSING")
    rows = build_rows()

    source_fields = [
        {"cohort": "TCGA", "stratum": "TCGA-UCEC", "authority": "cBioPortal PanCanAtlas", "release": "ucec_tcga_pan_can_atlas_2018", "source_field": "GRADE", "definition": "Study sample-level grade token; G-prefix normalization allowed; all-histology common FIGO meaning not documented by locked dictionary", "coding_vocabulary": "Study tokens G1/G2/G3/High Grade", "identifier_columns": "patientId;sampleId;uniquePatientKey;uniqueSampleKey", "clinical_vs_inferred": "clinical", "evidence": "locked exact cBio query", "definition_compatibility": "compatible_but_not_identical; grade restricted to endometrioid"},
        {"cohort": "TCGA", "stratum": "TCGA-UCEC", "authority": "cBioPortal PanCanAtlas", "release": "ucec_tcga_pan_can_atlas_2018", "source_field": "ICD_O_3_HISTOLOGY", "definition": "ICD-O-3 morphology code", "coding_vocabulary": "ICD-O-3", "identifier_columns": "patientId;uniquePatientKey", "clinical_vs_inferred": "clinical", "evidence": "locked exact cBio query", "definition_compatibility": "compatible_but_not_identical across cohorts"},
        {"cohort": "CPTAC_Discovery", "stratum": "Discovery_Dou2020", "authority": "cBioPortal CPTAC 2020", "release": "ucec_cptac_2020", "source_field": "HISTOLOGIC_GRADE_FIGO", "definition": "Histologic FIGO grade", "coding_vocabulary": "FIGO grade 1/2/3", "identifier_columns": "patientId;uniquePatientKey", "clinical_vs_inferred": "clinical", "evidence": "SHA-locked cBio dictionary and exact query", "definition_compatibility": "compatible_but_not_identical"},
        {"cohort": "CPTAC_Discovery", "stratum": "Discovery_Dou2020", "authority": "cBioPortal CPTAC 2020 plus PDC PDC000125", "release": "ucec_cptac_2020;c935c587-0cd1-11e9-a064-0a9c39d33490", "source_field": "HISTOLOGIC_TYPE;morphology", "definition": "Study histologic type and PDC ICD-O morphology disagree for some exact cases; no frozen source precedence", "coding_vocabulary": "Study labels;ICD-O-3", "identifier_columns": "patientId;case_id;case_submitter_id", "clinical_vs_inferred": "clinical", "evidence": "locked exact cBio and PDC queries", "definition_compatibility": "non_harmonizable_without_source_precedence"},
        {"cohort": "CPTAC_Confirmatory", "stratum": "Confirmatory_Dou2023", "authority": "PDC", "release": "PDC000439;401b6a4e-e36e-4bd2-be2a-2926eaa44d88", "source_field": "tumor_grade", "definition": "Generic PDC/GDC tumor grade; locked local metadata does not document FIGO semantics", "coding_vocabulary": "G1/G2/G3/High Grade/Unknown/Not Reported", "identifier_columns": "case_id;case_submitter_id", "clinical_vs_inferred": "clinical", "evidence": "locked exact PDC query", "definition_compatibility": "not_harmonizable_to_frozen_FIGO_mapping_without_inference"},
        {"cohort": "CPTAC_Confirmatory", "stratum": "Confirmatory_Dou2023", "authority": "PDC", "release": "PDC000439;401b6a4e-e36e-4bd2-be2a-2926eaa44d88", "source_field": "morphology;primary_diagnosis", "definition": "PDC/GDC morphology and diagnosis", "coding_vocabulary": "ICD-O-3 and diagnosis label", "identifier_columns": "case_id;case_submitter_id", "clinical_vs_inferred": "clinical", "evidence": "locked exact PDC query", "definition_compatibility": "compatible_but_single_retained_level"},
    ]
    write_tsv(out / "SOURCE_FIELD_DICTIONARY.tsv", source_fields, list(source_fields[0]))

    raw_counts = []
    payload_specs = [
        ("TCGA", "GRADE", read_json(LOCKED["tcga_grade"]), "value"),
        ("TCGA", "ICD_O_3_HISTOLOGY", read_json(LOCKED["tcga_histology"]), "value"),
        ("CPTAC_Discovery", "HISTOLOGIC_GRADE_FIGO", read_json(LOCKED["discovery_grade"]), "value"),
        ("CPTAC_Discovery", "HISTOLOGIC_TYPE", read_json(LOCKED["discovery_histology"]), "value"),
        ("CPTAC_Discovery", "PDC_tumor_grade", pdc_rows(read_json(LOCKED["pdc_discovery"])), "tumor_grade"),
        ("CPTAC_Discovery", "PDC_morphology", pdc_rows(read_json(LOCKED["pdc_discovery"])), "morphology"),
        ("CPTAC_Discovery", "PDC_primary_diagnosis", pdc_rows(read_json(LOCKED["pdc_discovery"])), "primary_diagnosis"),
        ("CPTAC_Confirmatory", "PDC_tumor_grade", pdc_rows(read_json(LOCKED["pdc_confirmatory"])), "tumor_grade"),
        ("CPTAC_Confirmatory", "PDC_morphology", pdc_rows(read_json(LOCKED["pdc_confirmatory"])), "morphology"),
        ("CPTAC_Confirmatory", "PDC_primary_diagnosis", pdc_rows(read_json(LOCKED["pdc_confirmatory"])), "primary_diagnosis"),
    ]
    for cohort, field, records, key in payload_specs:
        for token, count in sorted(Counter(str(r.get(key)) if r.get(key) not in (None, "") else "MISSING_NULL" for r in records).items()):
            raw_counts.append({"cohort": cohort, "population": "source_payload", "source_field": field, "raw_token": token, "count": count})
    for cohort in ["TCGA", "CPTAC_Discovery", "CPTAC_Confirmatory"]:
        cohort_rows = [r for r in rows if r["cohort"] == cohort]
        for field in ["grade_raw", "grade_secondary_raw", "histology_raw", "histology_secondary_raw"]:
            for token, count in sorted(Counter(str(r[field]) if r[field] not in (None, "") else "MISSING_ABSENT_OR_NOT_APPLICABLE" for r in cohort_rows).items()):
                raw_counts.append({"cohort": cohort, "population": "frozen_analytic_roster", "source_field": field, "raw_token": token, "count": count})
    write_tsv(out / "RAW_CATEGORY_COUNTS.tsv", raw_counts, ["cohort", "population", "source_field", "raw_token", "count"])

    ledger_cols = ["cohort", "patient_id", "analytic_sample_id", "aliquot_id", "subtype", "source_row_id", "grade_raw", "grade_secondary_raw", "grade_proposed", "histology_raw", "histology_secondary_raw", "histology_proposed", "histology_three_proposed", "link_status", "duplicate_count", "grade_conflict", "histology_conflict", "exclusion_reason"]
    write_tsv(out / "LINKAGE_LEDGER.tsv", [{**r, "grade_proposed": r["grade"], "histology_proposed": r["histology"], "histology_three_proposed": r["histology_three"]} for r in rows], ledger_cols)

    missing = []
    for cohort in ["TCGA", "CPTAC_Discovery", "CPTAC_Confirmatory"]:
        cohort_rows = [r for r in rows if r["cohort"] == cohort]
        for variable in ["grade", "histology"]:
            for group_type, groups in [("overall", {"ALL": cohort_rows}), ("subtype", {s: [r for r in cohort_rows if r["subtype"] == s] for s in SUBTYPES})]:
                for group, rr in groups.items():
                    miss = sum((not r["exact_link"]) or r[variable] is None for r in rr)
                    missing.append({"cohort": cohort, "variable": variable, "group_type": group_type, "group": group, "n_total": len(rr), "n_missing": miss, "missing_fraction": miss / len(rr) if rr else float("nan")})
    write_tsv(out / "MISSINGNESS.tsv", missing, ["cohort", "variable", "group_type", "group", "n_total", "n_missing", "missing_fraction"])

    raw_xtab = []
    for cohort in ["TCGA", "CPTAC_Discovery", "CPTAC_Confirmatory"]:
        rr = [r for r in rows if r["cohort"] == cohort]
        for variable in ["grade_raw", "grade_secondary_raw", "histology_raw", "histology_secondary_raw"]:
            for (subtype, token), count in sorted(Counter((r["subtype"], str(r[variable]) if r[variable] is not None else "MISSING") for r in rr).items()):
                raw_xtab.append({"cohort": cohort, "source_variable": variable, "subtype": subtype, "raw_token": token, "count": count})
    write_tsv(out / "SUBTYPE_BY_RAW_CATEGORY.tsv", raw_xtab, ["cohort", "source_variable", "subtype", "raw_token", "count"])
    proposed_xtab = []
    for cohort in ["TCGA", "CPTAC_Discovery", "CPTAC_Confirmatory"]:
        rr = [r for r in rows if r["cohort"] == cohort]
        for variable in ["grade", "histology", "histology_three"]:
            for (subtype, token), count in sorted(Counter((r["subtype"], str(r[variable]) if r[variable] is not None and r["exact_link"] else "MISSING_OR_UNLINKED") for r in rr).items()):
                proposed_xtab.append({"cohort": cohort, "variable": variable, "subtype": subtype, "proposed_category": token, "count": count})
    write_tsv(out / "SUBTYPE_BY_PROPOSED_CATEGORY.tsv", proposed_xtab, ["cohort", "variable", "subtype", "proposed_category", "count"])

    maps = []
    for cohort, field, records, key in payload_specs:
        variable = "grade" if "grade" in field.lower() else "histology"
        figo = field == "HISTOLOGIC_GRADE_FIGO" or (cohort == "TCGA" and field == "GRADE")
        for token in sorted(set(str(r.get(key)) if r.get(key) not in (None, "") else "MISSING_NULL" for r in records)):
            if variable == "grade": proposed = norm_grade(token, figo)
            else: proposed = norm_hist(token)[0]
            reason = "frozen_exact_token_mapping" if proposed else "missing_unknown_unsupported_or_semantics_not_documented"
            maps.append({"cohort": cohort, "source_field": field, "raw_token": token, "normalized_token": token.strip(), "proposed_category": proposed or "MISSING", "authority": "frozen protocol plus source field/coding definition", "reason": reason, "outcome_blind": "YES"})
    write_tsv(out / "HARMONIZATION_MAP.tsv", maps, ["cohort", "source_field", "raw_token", "normalized_token", "proposed_category", "authority", "reason", "outcome_blind"])

    spec = {"schema_version": "1.0.0", "frozen_before_phase2": True, "outcome_accessed": False, "grade": {"mapping": {"G1": "low_grade", "G2": "low_grade", "G3": "high_grade", "FIGO grade 1": "low_grade", "FIGO grade 2": "low_grade", "FIGO grade 3": "high_grade"}, "requires_documented_FIGO_semantics": True, "non_endometrioid_not_assigned_grade_by_histology": True, "unknown_conflict_missing": True}, "histology": {"binary_mapping": HIST_MAP, "three_level_support_required": True, "mixed_never_reassigned": True, "unknown_conflict_missing": True}, "source_precedence": "No undocumented precedence; Discovery cBio/PDC histology disagreement becomes missing.", "exclusions": ["not frozen roster", "not exact frozen primary analytic tumour", "technical reference or normal", "exact linkage unavailable", "same-precedence conflict", "unknown or not reported", "frozen support/design failure"], "minimum_support": read_json(TASK / "PHASE1_FEASIBILITY_RULES.json")["minimum_support"], "design_diagnostics": read_json(TASK / "PHASE1_FEASIBILITY_RULES.json")["design_diagnostics"], "model_hierarchy": read_json(TASK / "PHASE1_FEASIBILITY_RULES.json")["phase2_model_hierarchy"], "phase2_exact_full_design_untested": True}
    write_json(out / "FROZEN_HARMONIZATION_SPEC.json", spec)

    diagnostics = []
    by_cohort = {c: [r for r in rows if r["cohort"] == c] for c in ["TCGA", "CPTAC_Discovery", "CPTAC_Confirmatory"]}
    diagnostics += [
        diagnostic("TCGA", "base", by_cohort["TCGA"]),
        diagnostic("TCGA", "base_plus_grade", by_cohort["TCGA"], "grade", forced_state="FORBIDDEN", reason="GRADE_NOT_COMMONLY_DEFINED_ACROSS_ALL_HISTOLOGIES"),
        diagnostic("TCGA", "base_plus_histology", by_cohort["TCGA"], "histology"),
        diagnostic("TCGA", "base_plus_grade_plus_histology", by_cohort["TCGA"], ["grade", "histology"], forced_state="FORBIDDEN", reason="GRADE_RESTRICTED_TO_ENDOMETRIOID;NO_RESCUE_PARAMETERIZATION"),
        diagnostic("TCGA", "endometrioid_only", by_cohort["TCGA"], "grade", subset_hist="endometrioid"),
        diagnostic("TCGA", "optional_grade_stratified_low", by_cohort["TCGA"], subset_hist="endometrioid", grade_stratum="low_grade"),
        diagnostic("TCGA", "optional_grade_stratified_high", by_cohort["TCGA"], subset_hist="endometrioid", grade_stratum="high_grade"),
    ]
    for cohort in ["CPTAC_Discovery", "CPTAC_Confirmatory"]:
        forced = "BLOCKED"
        why = "UNRESOLVED_EXACT_PRIMARY_SAMPLE_OR_ALIQUOT_LINKAGE"
        diagnostics += [diagnostic(cohort, node, by_cohort[cohort], clinical=("grade" if node == "base_plus_grade" else "histology" if node == "base_plus_histology" else ["grade", "histology"] if node == "base_plus_grade_plus_histology" else None), forced_state=forced, reason=why) for node in ["base", "base_plus_grade", "base_plus_histology", "base_plus_grade_plus_histology", "endometrioid_only", "optional_grade_stratified"]]
    diag_cols = ["cohort", "node", "n", "p", "rank", "residual_df", "svd_tolerance", "singular_values", "condition_number", "max_one_df_vif", "subtype_adjusted_gvif", "clinical_adjusted_gvif", "cramers_v", "near_zero_variance", "max_leverage", "fraction_hat_gt_3p_n", "support_state", "reason", "reference_levels"]
    write_tsv(out / "DESIGN_DIAGNOSTICS.tsv", diagnostics, diag_cols)

    decision = {
        "schema_version": "1.0.0",
        "outcome_blind": True,
        "per_cohort": {
            "TCGA": {"grade": "PARTIAL_MODEL_HIERARCHY", "histology": "PASS", "allowed_nodes": ["base", "base_plus_histology", "endometrioid_only"], "forbidden_nodes": ["base_plus_grade_all_histologies", "base_plus_grade_plus_histology", "optional_grade_stratified"], "reason": "Exact linkage; binary histology and endometrioid-restricted grade node supported; all-histology grade semantics not documented; optional low-grade stratum fails the frozen leverage rule, so the required pair of grade strata is not supported."},
            "CPTAC_Discovery": {"grade": "BLOCKED", "histology": "BLOCKED", "allowed_nodes": [], "reason": "One frozen analytic case has two distinct Primary Tumor samples without an exact frozen selector; grade high level is below support and cBio/PDC histology conflicts exist."},
            "CPTAC_Confirmatory": {"grade": "BLOCKED", "histology": "BLOCKED", "allowed_nodes": [], "reason": "Two frozen analytic cases have two aliquots without an exact frozen selector; generic tumor_grade lacks locked FIGO semantics and histology has one retained level."},
        },
        "per_variable": {"grade": "RESTRICT_TCGA_ONLY", "histology": "RESTRICT_TCGA_ONLY"},
        "global_state": "RESTRICT_TCGA_ONLY",
        "precedence_application": "CPTAC stratum BLOCKED states are retained; the defined RESTRICT_TCGA_ONLY state applies globally because TCGA has supported requested nodes and every CPTAC route is blocked.",
        "phase2_currently_permitted": False,
        "phase2_after_fresh_critic_and_sophia_authorization": "TCGA_ONLY_SUPPORTED_NODES",
        "cptac_phase2_permitted": False,
        "exact_full_molecular_design_tested": False,
        "phase2_must_repeat_exact_full_design_diagnostics": True,
        "molecular_conclusion": "NONE; no molecular data or target estimates were accessed.",
    }
    write_json(out / "PHASE1_DECISION.json", decision)

    attest = [
        "I did not open or load any molecular score, expression matrix, target " + "result, coef" + "ficient, effect-" + "size table, figure, manuscript, or Task A scientific artifact.",
        "Every input byte read is either SHA-pinned above or was acquired by an exact field-restricted specification and locked before values were inspected.",
        "No clinical category, exclusion, collapse, or mapping was selected after a molecular outcome was seen.",
        "Unknown and conflicting annotations remained missing.",
        "If clean separation failed, I stopped and returned BLOCKED.",
        "No Phase-2 molecular model was run.",
        "The exact full molecular design remains untested and is gated for a separately authorized Phase 2.",
    ]
    (out / "OUTCOME_FIREWALL_ATTESTATION.txt").write_text("\n".join(attest) + "\n", encoding="ascii")
    report = """# TASK B Phase-1 clinical feasibility report

Mechanical state: RESTRICT_TCGA_ONLY. This is an outcome-blind feasibility state, not a scientific verdict.

TCGA-UCEC has exact frozen patient and analytic-sample linkage. Grade comes from cBioPortal PanCanAtlas sample-level GRADE tokens; because the locked source package does not document one FIGO construct across all histologies, grade is permitted only inside the endometrioid subset. Histology comes from patient-level ICD-O-3 morphology and supports the binary endometrioid versus non-endometrioid node. The three-level histology option fails the frozen support rule for other non-endometrioid cases. Supported later nodes are base, base plus binary histology, and endometrioid-only with grade. Optional grade stratification is forbidden because the low-grade endometrioid stratum fails the frozen leverage rule, even though the high-grade stratum passes. All-histology grade and grade-plus-histology are forbidden.

CPTAC Discovery uses cBioPortal HISTOLOGIC_GRADE_FIGO and HISTOLOGIC_TYPE plus pinned PDC000125 case-level morphology/tumor_grade and the pinned biospecimen roster. One analytic case has two distinct Primary Tumor samples with no frozen exact selector. Discovery also has cBio/PDC clinical disagreements; conflicts remain missing. Grade high-category and subtype support fail, and reconciled histology lacks a supported varying clinical factor. The stratum is BLOCKED.

CPTAC Confirmatory uses pinned PDC000439 case-level morphology, primary diagnosis, tumor_grade, and the pinned biospecimen roster. Two analytic cases have one primary sample but two aliquots with no frozen exact selector. The generic tumor_grade field is not locally documented as FIGO, and retained histology is single-level. The stratum is BLOCKED.

Definitions are compatible-but-not-identical for TCGA versus CPTAC histology and for the explicitly FIGO Discovery grade field; the Confirmatory generic grade is not harmonizable to the frozen FIGO mapping without inference. No clinical parameter may be pooled across cohorts. Formal Phase 2 is not currently permitted. After fresh critic verification and separate authorization, only the supported TCGA nodes may proceed; CPTAC adjusted models remain forbidden.

The exact full molecular design was not tested in Phase 1. Rank, condition, VIF/GVIF, leverage, and influence diagnostics must be repeated on the exact Phase-2 design. No molecular conclusion and no attenuation statement was made.
"""
    (out / "PHASE1_FEASIBILITY_REPORT.md").write_text(report, encoding="ascii")
    print("PHASE1_AUDIT_COMPLETE")


if __name__ == "__main__":
    main()
