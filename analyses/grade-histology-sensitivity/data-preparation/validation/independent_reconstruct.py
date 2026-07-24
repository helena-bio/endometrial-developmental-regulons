#!/usr/bin/env python3
"""Independent clinical/subtype-only reconstruction for TASK B Phase 1.

No producer analysis module is imported.  Inputs are explicit and default-deny.
The producer tables are opened only after the independent rows and diagnostics
have been constructed, for field-level comparison.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import pathlib
from collections import Counter, defaultdict

import numpy as np


ROOT = pathlib.Path("data/external/original-workspace/revgate-task-b-grade-histology")
TASK = ROOT / "experiments/taskB_grade_histology"
EXEC = TASK / "phase1_execution"
OUT = TASK / "phase1_critic"

ALLOW = [
    ("data/external/original-workspace/task027-acquire-freeze-a/freeze_a_redux/cohort_selected_primary.tsv", "e6f9d390730bff4248379b2280b12e2914ec800e2f56516b466f2151b169f3be"),
    ("data/external/original-workspace/task024-freeze-a/subtype_normalized.tsv", "17a899a4869fccf806859342454bde440027786a48bdd50f1d1eace810b1ef9e"),
    ("data/external/original-workspace/task029-external-replication-feasibility/join_tables/discovery_join.json", "315a9d29f1add19fd90db85188762937ac2e49262d64a51813aca5fb5ccb48ad"),
    ("data/external/original-workspace/task029-external-replication-feasibility/join_tables/confirmatory_join.json", "95b9f754970e7f0691faddf6e350aca6cfe465c5c0c98600a4ed34ebbfe43743"),
    ("data/external/original-workspace/task029-external-replication-feasibility/sources/cbio_2020_clinattr.json", "25161910befe62fea99ba5d754e0e43f8d2d80413538b8bfcac0c52aa5e58e6e"),
    ("data/external/original-workspace/task029-external-replication-feasibility/sources_task029gate/uec_cptac_gdc_clinattr.json", "96794b79b83a663a9caf5dc51dce3d71a29983040a2b4125f972a100c7627628"),
    ("data/external/original-workspace/task029-external-replication-feasibility/sources/pdc_clin_PDC000125.json", "61d05e1298312fbdf98e23b027b93e9b0310612baf1fc9fb6fe088913bc08291"),
    ("data/external/original-workspace/task029-external-replication-feasibility/sources/pdc_clin_PDC000439.json", "6ff6418b4f201d7f2db3b1c69617d3435a01b4877cf37a06f4c56a4fb7ba4012"),
    ("data/external/original-workspace/task029-external-replication-feasibility/sources/pdc_discovery_biospec.json", "c98f51f1ec9916252a210c40fae93adfe2ae42e998decdac5adbf82c8bd3b74e"),
    ("data/external/original-workspace/task029-external-replication-feasibility/sources/pdc_confirmatory_biospec.json", "bc159c1d0500cce635b7631a069c814a9f3c802f947dc3c9bf828573e0b7aae9"),
    ("data/external/original-workspace/task029-external-replication-feasibility/sources_task029gate/gdc_case_mapping.json", "4e5fd1736d2b5a17f6cd66929953a6893f4858a0a480480faa26ccafd2e798c8"),
    ("data/external/original-workspace/task029-external-replication-feasibility/sources_task029gate/gdc_rna_primary_by_case.json", "87514a428a9cf4810621ff7c46f5c6696f78bda2e3fecf66e7fe9920cbe9fbda"),
]

ACQUIRED = {
    "tcga_cbio_grade.json": "4af93f837dda71a50b6e5010fb63404b851f29ed4b671f611351c1d6af3c3edc",
    "tcga_cbio_icdo_histology.json": "fce9c7b47ef10a85876cc92c8b1eee330246db91a4c20b0403f9b9ce0f4a89b9",
    "cptac_discovery_grade.json": "287b4a02df4f41c3ff64727215bf56a003b6d048e7d20e06af1e359d9d5075d9",
    "cptac_discovery_histology.json": "4133634004e840605a9e1fc9c3b4a1cafeeb8dd182fcf4a5dc9f5cfef7e17393",
    "pdc_PDC000125_grade_histology.json": "8143943f3084a93ee1e49a62f892c8ac02b16efd74ee7a5e42df8109157e442e",
    "pdc_PDC000439_grade_histology.json": "de0dc78e2705070d8b21d0fb3db8c2a777a246ef49398036872fc0c6abebe61d",
}

PROTO = {
    "PHASE1_OUTCOME_BLIND_PROTOCOL.md": "09c6950e997b9c2534180d34c86e246ba749943b07a4f6bed62f5cbcea5ee2b9",
    "PHASE1_FEASIBILITY_RULES.json": "aa1eaea70991d9e915f499243a11964e01b7242c2f190cac783661b896e1eb42",
    "ALLOWED_INPUTS_AND_OUTCOME_FIREWALL.json": "40e11800d4d5f56388ca7563886215d07e219dc6d8f440b7069e0d3be772615b",
    "RESEARCH_ACCESS_LOG.tsv": "bd754088fe6d8528b5f5ddbdb0246a956834d06e75fc8b5cebd6ffa20934352b",
}

TCGA_HIST = {
    "8020/3": ("non_endometrioid", "other_non_endometrioid"),
    "8255/3": ("non_endometrioid", "other_non_endometrioid"),
    "8310/3": ("non_endometrioid", "other_non_endometrioid"),
    "8380/3": ("endometrioid", "endometrioid"),
    "8382/3": ("endometrioid", "endometrioid"),
    "8441/3": ("non_endometrioid", "serous"),
    "8460/3": ("non_endometrioid", "serous"),
    "8461/3": ("non_endometrioid", "serous"),
}
GRADE = {"G1": "low_grade", "G2": "low_grade", "G3": "high_grade", "FIGO grade 1": "low_grade", "FIGO grade 2": "low_grade", "FIGO grade 3": "high_grade"}
SUBTYPES = ["NSMP", "MMRd", "POLE", "p53abn"]


def sha_path(path: pathlib.Path | str) -> str:
    return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()


def load_json(path: pathlib.Path | str):
    with pathlib.Path(path).open(encoding="utf-8") as handle:
        return json.load(handle)


def load_tsv(path: pathlib.Path | str) -> list[dict[str, str]]:
    with pathlib.Path(path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: pathlib.Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="ascii") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def na(value):
    return "NA" if value is None else value


def pdc_records(name: str) -> list[dict]:
    obj = load_json(TASK / "inputs" / name)
    assert set(obj) == {"data"} and set(obj["data"]) == {"clinicalPerStudy"}
    records = obj["data"]["clinicalPerStudy"]
    assert all(set(r) == {"case_id", "case_submitter_id", "morphology", "primary_diagnosis", "tumor_grade"} for r in records)
    return records


def exact_one_primary(biospec_path: str, roster: dict[str, str]) -> tuple[dict[str, list[dict]], list[str]]:
    obj = load_json(biospec_path)
    assert set(obj) == {"data"} and set(obj["data"]) == {"biospecimenPerStudy"}
    expected = {"aliquot_id", "aliquot_submitter_id", "case_id", "case_submitter_id", "disease_type", "sample_id", "sample_submitter_id", "sample_type"}
    assert all(set(r) == expected for r in obj["data"]["biospecimenPerStudy"])
    grouped = defaultdict(list)
    for r in obj["data"]["biospecimenPerStudy"]:
        if r["sample_type"] == "Primary Tumor" and r["case_submitter_id"] in roster:
            grouped[r["case_submitter_id"]].append(r)
    ambiguous = []
    for case in roster:
        samples = {r["sample_submitter_id"] for r in grouped[case]}
        aliquots = {r["aliquot_submitter_id"] for r in grouped[case]}
        if len(samples) != 1 or len(aliquots) != 1:
            ambiguous.append(case)
    return grouped, sorted(ambiguous)


def build_rows() -> list[dict]:
    # TCGA: only the permitted identifier/subtype columns are retained.
    tcga_roster = load_tsv(ALLOW[0][0])
    tcga_sub = {r["patient_barcode"]: r["mapped_4way"] for r in load_tsv(ALLOW[1][0])}
    assert len(tcga_roster) == 507 and len(tcga_sub) == 507
    assert {r["patient_barcode"]: r["mapped_4way"] for r in tcga_roster} == tcga_sub
    tg = load_json(TASK / "inputs/tcga_cbio_grade.json")
    th = load_json(TASK / "inputs/tcga_cbio_icdo_histology.json")
    assert all(set(r) == {"clinicalAttributeId", "patientId", "sampleId", "studyId", "uniquePatientKey", "uniqueSampleKey", "value"} and r["clinicalAttributeId"] == "GRADE" for r in tg)
    assert all(set(r) == {"clinicalAttributeId", "patientId", "studyId", "uniquePatientKey", "value"} and r["clinicalAttributeId"] == "ICD_O_3_HISTOLOGY" for r in th)
    tg_map = {r["sampleId"]: r["value"] for r in tg}
    th_map = {r["patientId"]: r["value"] for r in th}
    rows = []
    for r in tcga_roster:
        patient = r["patient_barcode"]
        sample = r["kept_aliquot_barcode"][:15]
        raw_g, raw_h = tg_map.get(sample), th_map.get(patient)
        hb, h3 = TCGA_HIST.get(raw_h, (None, None))
        rows.append({"cohort": "TCGA", "patient_id": patient, "analytic_sample_id": sample, "aliquot_id": r["kept_aliquot_barcode"], "subtype": r["mapped_4way"], "source_row_id": patient, "grade_raw": raw_g, "grade_secondary_raw": None, "grade": GRADE.get(raw_g), "histology_raw": raw_h, "histology_secondary_raw": None, "histology": hb, "histology_three": h3, "link_status": "EXACT", "duplicate_count": 1, "grade_conflict": False, "histology_conflict": False, "exclusion_reason": ""})

    # CPTAC Discovery: cBio FIGO is the grade source; cBio and PDC histology
    # are same-precedence clinical sources and disagreement remains missing.
    dj = load_json(ALLOW[2][0])
    dsub = dj["case_to_subtype_joined"]
    assert len(dsub) == 95 and set(dsub.values()) == set(SUBTYPES)
    dg = {r["patientId"]: r["value"] for r in load_json(TASK / "inputs/cptac_discovery_grade.json")}
    dh = {r["patientId"]: r["value"] for r in load_json(TASK / "inputs/cptac_discovery_histology.json")}
    dp = {r["case_submitter_id"]: r for r in pdc_records("pdc_PDC000125_grade_histology.json")}
    dspec, damb = exact_one_primary(ALLOW[8][0], dsub)
    assert len(damb) == 1
    for patient, subtype in dsub.items():
        pr = dp[patient]
        primary = dspec[patient]
        raw_g, sec_g = dg.get(patient), pr["tumor_grade"]
        raw_h, sec_h = dh.get(patient), pr["morphology"]
        # Compare only where both encode an ordinal grade; PDC remains secondary
        # because its generic semantics are not used for harmonization.
        gnorm = {"FIGO grade 1": "G1", "FIGO grade 2": "G2", "FIGO grade 3": "G3"}.get(raw_g)
        binary = {"G1": "low_grade", "G2": "low_grade", "G3": "high_grade"}
        gconf = gnorm is not None and sec_g in binary and binary[gnorm] != binary[sec_g]
        pdc_hist = "endometrioid" if sec_h == "8380/3" else None
        cbio_hist = "endometrioid" if raw_h == "Endometrioid" else "non_endometrioid" if raw_h == "Serous" else None
        hconf = cbio_hist is not None and pdc_hist is not None and cbio_hist != pdc_hist
        exact = patient not in damb
        one = primary[0] if exact else None
        rows.append({"cohort": "CPTAC_Discovery", "patient_id": patient, "analytic_sample_id": one["sample_submitter_id"] if one else None, "aliquot_id": one["aliquot_submitter_id"] if one else None, "subtype": subtype, "source_row_id": pr["case_id"], "grade_raw": raw_g, "grade_secondary_raw": sec_g, "grade": GRADE.get(raw_g), "histology_raw": raw_h, "histology_secondary_raw": sec_h, "histology": None if hconf else cbio_hist, "histology_three": None if hconf else ("endometrioid" if raw_h == "Endometrioid" else "serous" if raw_h == "Serous" else None), "link_status": "EXACT" if exact else "BLOCKED_MULTIPLE_PRIMARY_SAMPLE_OR_ALIQUOT", "duplicate_count": len(primary), "grade_conflict": gconf, "histology_conflict": hconf, "exclusion_reason": "" if exact else "frozen_exact_analytic_tumour_or_aliquot_unavailable"})

    # Confirmatory: PDC generic tumour_grade is descriptive but cannot be mapped
    # to FIGO.  Morphology is endometrioid; two cases have no exact aliquot.
    cj = load_json(ALLOW[3][0])
    csub = cj["case_to_subtype_joined"]
    assert len(csub) == 135 and set(csub.values()) == set(SUBTYPES)
    cp = {r["case_submitter_id"]: r for r in pdc_records("pdc_PDC000439_grade_histology.json")}
    cspec, camb = exact_one_primary(ALLOW[9][0], csub)
    assert len(camb) == 2
    for patient, subtype in csub.items():
        pr = cp[patient]
        primary = cspec[patient]
        exact = patient not in camb
        one = primary[0] if exact else None
        hb = "endometrioid" if pr["morphology"] == "8380/3" else None
        rows.append({"cohort": "CPTAC_Confirmatory", "patient_id": patient, "analytic_sample_id": one["sample_submitter_id"] if one else None, "aliquot_id": one["aliquot_submitter_id"] if one else None, "subtype": subtype, "source_row_id": pr["case_id"], "grade_raw": pr["tumor_grade"], "grade_secondary_raw": None, "grade": None, "histology_raw": pr["morphology"], "histology_secondary_raw": None, "histology": hb, "histology_three": hb, "link_status": "EXACT" if exact else "BLOCKED_MULTIPLE_PRIMARY_SAMPLE_OR_ALIQUOT", "duplicate_count": len(primary), "grade_conflict": False, "histology_conflict": False, "exclusion_reason": "" if exact else "frozen_exact_analytic_tumour_or_aliquot_unavailable"})
    return rows


def matrix(rows: list[dict], clinical: str | None = None) -> tuple[np.ndarray, list[str], dict[str, list[int]]]:
    columns = [np.ones(len(rows), dtype=float)]
    names = ["intercept"]
    terms = {"subtype": []}
    for subtype in ["MMRd", "POLE", "p53abn"]:
        terms["subtype"].append(len(columns))
        columns.append(np.array([r["subtype"] == subtype for r in rows], float))
        names.append(f"subtype:{subtype}")
    if clinical in {"grade", "both"}:
        terms["grade"] = [len(columns)]
        columns.append(np.array([r["grade"] == "low_grade" for r in rows], float))
        names.append("grade:low_grade")
    if clinical in {"histology", "both"} and len({r["histology"] for r in rows}) > 1:
        terms["histology"] = [len(columns)]
        columns.append(np.array([r["histology"] == "non_endometrioid" for r in rows], float))
        names.append("histology:non_endometrioid")
    return np.column_stack(columns), names, terms


def cramers_v(rows: list[dict], clinical: str) -> float | None:
    levels = sorted({r[clinical] for r in rows})
    if len(levels) < 2:
        return None
    table = np.array([[sum(r["subtype"] == s and r[clinical] == c for r in rows) for c in levels] for s in SUBTYPES], float)
    expected = table.sum(1)[:, None] * table.sum(0)[None, :] / table.sum()
    chi = np.sum((table - expected) ** 2 / expected)
    return math.sqrt(chi / (table.sum() * min(table.shape[0] - 1, table.shape[1] - 1)))


def diagnostic(cohort: str, node: str, rows: list[dict], clinical: str | None, forced: str | None = None, correct_protocol_node: bool = True) -> dict:
    work = [r for r in rows if r["link_status"] == "EXACT"]
    if node == "endometrioid_only" and correct_protocol_node:
        work = [r for r in work if r["histology"] == "endometrioid" and r["grade"] is not None]
        clinical = "grade"
    elif node.startswith("optional_grade_stratified_"):
        target = node.rsplit("_", 1)[-1] + "_grade"
        work = [r for r in work if r["histology"] == "endometrioid" and r["grade"] == target]
        clinical = None
    elif clinical == "both":
        work = [r for r in work if r["grade"] is not None and r["histology"] is not None]
    elif clinical:
        work = [r for r in work if r[clinical] is not None]
    n = len(work)
    if n == 0:
        return {"cohort": cohort, "node": node, "n": 0, "p": 4, "rank": None, "residual_df": -4, "condition_number": None, "max_one_df_vif": None, "subtype_adjusted_gvif": None, "clinical_adjusted_gvif": None, "cramers_v": None, "near_zero_variance": None, "max_leverage": None, "fraction_hat_gt_3p_n": None, "support_state": forced or "FAIL"}
    X, names, terms = matrix(work, clinical)
    n, p = X.shape
    singular = np.linalg.svd(X, compute_uv=False)
    rank = int(np.linalg.matrix_rank(X, tol=np.finfo(float).eps * max(n, p) * singular[0]))
    z = (X[:, 1:] - X[:, 1:].mean(0)) / X[:, 1:].std(0, ddof=0)
    condition = float(np.linalg.cond(z)) if z.shape[1] else 1.0
    corr = np.corrcoef(X[:, 1:], rowvar=False)
    if np.ndim(corr) == 0:
        corr = np.array([[1.0]])
    corr_singular = not np.all(np.isfinite(corr)) or np.linalg.matrix_rank(corr) < corr.shape[0]
    vifs = np.full(corr.shape[0], np.inf) if corr_singular else np.diag(np.linalg.inv(corr))
    subtype_idx = [i - 1 for i in terms["subtype"]]
    other_idx = [i for i in range(corr.shape[0]) if i not in subtype_idx]
    if other_idx and not corr_singular:
        gvif = np.linalg.det(corr[np.ix_(subtype_idx, subtype_idx)]) * np.linalg.det(corr[np.ix_(other_idx, other_idx)]) / np.linalg.det(corr)
        subtype_gvif = float(gvif ** (1 / (2 * len(subtype_idx))))
        # Each clinical factor is a distinct one-df term; report the maximum
        # adjusted GVIF across those terms, while subtype is one multi-df term.
        clinical_gvif = float(max(math.sqrt(vifs[i]) for i in other_idx))
    elif not other_idx:
        subtype_gvif, clinical_gvif = 1.0, None
    else:
        subtype_gvif, clinical_gvif = math.inf, math.inf
    hat = np.sum(X * (X @ np.linalg.pinv(X.T @ X)), axis=1)
    clinical_fields = ["grade", "histology"] if clinical == "both" else [clinical] if clinical else []
    clinical_levels = {f: sorted({r[f] for r in work}) for f in clinical_fields}
    nzv = any(len(levels) < 2 or (max(Counter(r[f] for r in work).values()) / min(Counter(r[f] for r in work).values()) > 19 and 100 * len(levels) / n <= 10) for f, levels in clinical_levels.items())
    cv_values = [cramers_v(work, f) for f in clinical_fields if len(clinical_levels[f]) > 1]
    state = forced or "PASS"
    return {"cohort": cohort, "node": node, "n": n, "p": p, "rank": rank, "residual_df": n - p, "svd_tolerance": float(np.finfo(float).eps * max(n, p) * singular[0]), "singular_values": ";".join(f"{x:.10g}" for x in singular), "condition_number": condition, "max_one_df_vif": float(max(vifs)), "subtype_adjusted_gvif": subtype_gvif, "clinical_adjusted_gvif": clinical_gvif, "cramers_v": max(cv_values) if cv_values else None, "near_zero_variance": nzv, "max_leverage": float(max(hat)), "fraction_hat_gt_3p_n": float(np.mean(hat > 3 * p / n)), "support_state": state}


def diagnostics(rows: list[dict]) -> tuple[list[dict], list[dict]]:
    by = {c: [r for r in rows if r["cohort"] == c] for c in ["TCGA", "CPTAC_Discovery", "CPTAC_Confirmatory"]}
    correct = []
    correct += [diagnostic("TCGA", "base", by["TCGA"], None), diagnostic("TCGA", "base_plus_grade", by["TCGA"], "grade", "FORBIDDEN"), diagnostic("TCGA", "base_plus_histology", by["TCGA"], "histology"), diagnostic("TCGA", "base_plus_grade_plus_histology", by["TCGA"], "both", "FORBIDDEN"), diagnostic("TCGA", "endometrioid_only", by["TCGA"], None), diagnostic("TCGA", "optional_grade_stratified_low", by["TCGA"], None, "FAIL"), diagnostic("TCGA", "optional_grade_stratified_high", by["TCGA"], None)]
    for cohort in ["CPTAC_Discovery", "CPTAC_Confirmatory"]:
        correct += [diagnostic(cohort, "base", by[cohort], None, "BLOCKED"), diagnostic(cohort, "base_plus_grade", by[cohort], "grade", "BLOCKED"), diagnostic(cohort, "base_plus_histology", by[cohort], "histology", "BLOCKED"), diagnostic(cohort, "base_plus_grade_plus_histology", by[cohort], "both", "BLOCKED"), diagnostic(cohort, "endometrioid_only", by[cohort], None, "BLOCKED"), diagnostic(cohort, "optional_grade_stratified_low", by[cohort], None, "BLOCKED"), diagnostic(cohort, "optional_grade_stratified_high", by[cohort], None, "BLOCKED")]
    # Producer-labelled diagnostics recreated independently to compare numeric rows.
    labelled = list(correct[:7])
    for cohort in ["CPTAC_Discovery", "CPTAC_Confirmatory"]:
        labelled += [diagnostic(cohort, "base", by[cohort], None, "BLOCKED"), diagnostic(cohort, "base_plus_grade", by[cohort], "grade", "BLOCKED"), diagnostic(cohort, "base_plus_histology", by[cohort], "histology", "BLOCKED"), diagnostic(cohort, "base_plus_grade_plus_histology", by[cohort], "both", "BLOCKED"), diagnostic(cohort, "endometrioid_only", by[cohort], None, "BLOCKED", False), diagnostic(cohort, "optional_grade_stratified", by[cohort], None, "BLOCKED", False)]
    return correct, labelled


def aggregate_summary(rows: list[dict], correct_diag: list[dict]) -> dict:
    by = {c: [r for r in rows if r["cohort"] == c] for c in ["TCGA", "CPTAC_Discovery", "CPTAC_Confirmatory"]}
    def counts(rs, field):
        return dict(sorted(Counter(na(r[field]) for r in rs).items()))
    tcga = by["TCGA"]
    disc = by["CPTAC_Discovery"]
    conf = by["CPTAC_Confirmatory"]
    dgrade_disagree = sum(r["grade_raw"] is not None and r["grade_secondary_raw"] in {"G1", "G2", "G3"} and {"FIGO grade 1": "G1", "FIGO grade 2": "G2", "FIGO grade 3": "G3"}.get(r["grade_raw"]) != r["grade_secondary_raw"] for r in disc)
    boundary = sum(r["grade_conflict"] for r in disc)
    return {
        "TCGA": {"n": len(tcga), "grade_raw": counts(tcga, "grade_raw"), "grade_proposed": counts(tcga, "grade"), "histology_proposed": counts(tcga, "histology")},
        "CPTAC_Discovery": {"n": len(disc), "ambiguous_linkage": sum(r["link_status"] != "EXACT" for r in disc), "grade_disagreements_any_ordinal": dgrade_disagree, "grade_binary_boundary_conflicts": boundary, "histology_conflicts": sum(r["histology_conflict"] for r in disc), "grade_proposed": counts(disc, "grade"), "histology_proposed": counts(disc, "histology"), "subtype_counts": counts(disc, "subtype")},
        "CPTAC_Confirmatory": {"n": len(conf), "ambiguous_linkage": sum(r["link_status"] != "EXACT" for r in conf), "grade_proposed": counts(conf, "grade"), "histology_proposed": counts(conf, "histology"), "subtype_counts": counts(conf, "subtype")},
        "diagnostic_nodes": [{k: r.get(k) for k in ["cohort", "node", "n", "p", "rank", "residual_df", "condition_number", "max_leverage", "fraction_hat_gt_3p_n", "support_state"]} for r in correct_diag],
        "mechanical_state": {"TCGA_candidates": ["base", "base_plus_binary_histology", "endometrioid_only_plus_grade"], "TCGA_forbidden": ["all_histology_plus_grade", "all_histology_plus_grade_plus_histology", "optional_grade_stratified"], "CPTAC_adjustment": "BLOCKED", "global": "RESTRICT_TCGA_ONLY", "phase2_currently_permitted": False},
    }


def reconstructed_tables(rows: list[dict]) -> dict[str, list[dict]]:
    raw = []
    def add_counts(cohort, population, field, values, missing_token=None):
        tokens = [missing_token if v is None else v for v in values]
        for token, count in sorted(Counter(tokens).items()):
            raw.append({"cohort": cohort, "population": population, "source_field": field, "raw_token": token, "count": count})
    # Immutable source-payload counts.
    for field, filename in [("GRADE", "tcga_cbio_grade.json"), ("ICD_O_3_HISTOLOGY", "tcga_cbio_icdo_histology.json")]:
        add_counts("TCGA", "source_payload", field, [r["value"] for r in load_json(TASK / "inputs" / filename)])
    for field, filename in [("HISTOLOGIC_GRADE_FIGO", "cptac_discovery_grade.json"), ("HISTOLOGIC_TYPE", "cptac_discovery_histology.json")]:
        add_counts("CPTAC_Discovery", "source_payload", field, [r["value"] for r in load_json(TASK / "inputs" / filename)])
    for cohort, filename in [("CPTAC_Discovery", "pdc_PDC000125_grade_histology.json"), ("CPTAC_Confirmatory", "pdc_PDC000439_grade_histology.json")]:
        recs = pdc_records(filename)
        for source, key in [("PDC_tumor_grade", "tumor_grade"), ("PDC_morphology", "morphology"), ("PDC_primary_diagnosis", "primary_diagnosis")]:
            add_counts(cohort, "source_payload", source, [r[key] for r in recs])
    by = {c: [r for r in rows if r["cohort"] == c] for c in ["TCGA", "CPTAC_Discovery", "CPTAC_Confirmatory"]}
    for cohort, cohort_rows in by.items():
        for field in ["grade_raw", "grade_secondary_raw", "histology_raw", "histology_secondary_raw"]:
            add_counts(cohort, "frozen_analytic_roster", field, [r[field] for r in cohort_rows], "MISSING_ABSENT_OR_NOT_APPLICABLE")
    miss = []
    proposed = []
    raw_cross = []
    for cohort, cohort_rows in by.items():
        for variable in ["grade", "histology"]:
            for group_type, groups in [("overall", {"ALL": cohort_rows}), ("subtype", {s: [r for r in cohort_rows if r["subtype"] == s] for s in ["NSMP", "MMRd", "POLE", "p53abn"]})]:
                for group, rs in groups.items():
                    n_missing = sum(r[variable] is None or r["link_status"] != "EXACT" for r in rs)
                    miss.append({"cohort": cohort, "variable": variable, "group_type": group_type, "group": group, "n_total": len(rs), "n_missing": n_missing, "missing_fraction": n_missing / len(rs)})
        for variable in ["grade", "histology", "histology_three"]:
            c = Counter((r["subtype"], "MISSING_OR_UNLINKED" if r[variable] is None or r["link_status"] != "EXACT" else r[variable]) for r in cohort_rows)
            for (subtype, category), count in sorted(c.items()): proposed.append({"cohort": cohort, "variable": variable, "subtype": subtype, "proposed_category": category, "count": count})
        for variable in ["grade_raw", "grade_secondary_raw", "histology_raw", "histology_secondary_raw"]:
            c = Counter((r["subtype"], "MISSING" if r[variable] is None else r[variable]) for r in cohort_rows)
            for (subtype, token), count in sorted(c.items()): raw_cross.append({"cohort": cohort, "source_variable": variable, "subtype": subtype, "raw_token": token, "count": count})
    return {"RAW_CATEGORY_COUNTS.tsv": raw, "MISSINGNESS.tsv": miss, "SUBTYPE_BY_PROPOSED_CATEGORY.tsv": proposed, "SUBTYPE_BY_RAW_CATEGORY.tsv": raw_cross}


def compare_reconstructed_tables(tables: dict[str, list[dict]]) -> list[dict]:
    specs = {
        "RAW_CATEGORY_COUNTS.tsv": (["cohort", "population", "source_field", "raw_token"], "count"),
        "MISSINGNESS.tsv": (["cohort", "variable", "group_type", "group"], "n_total,n_missing,missing_fraction"),
        "SUBTYPE_BY_PROPOSED_CATEGORY.tsv": (["cohort", "variable", "subtype", "proposed_category"], "count"),
        "SUBTYPE_BY_RAW_CATEGORY.tsv": (["cohort", "source_variable", "subtype", "raw_token"], "count"),
    }
    out = []
    for name, independent in tables.items():
        keys, value_fields = specs[name]
        vals = value_fields.split(",")
        p_rows = load_tsv(EXEC / "run1" / name)
        pi = {tuple(r[k] for k in keys): r for r in p_rows}
        ii = {tuple(str(r[k]) for k in keys): r for r in independent}
        for key in sorted(set(pi) | set(ii)):
            for field in vals:
                iv = ii.get(key, {}).get(field, "MISSING_ROW")
                pv = pi.get(key, {}).get(field, "MISSING_ROW")
                try:
                    diff = abs(float(iv) - float(pv)); match = diff <= 5e-9 * max(1.0, abs(float(iv)))
                except (ValueError, TypeError):
                    diff = None; match = str(iv) == str(pv)
                out.append({"artifact": name, "row_key": "|".join(key), "field": field, "independent": str(iv), "producer": str(pv), "match": str(match).upper(), "absolute_difference": "NA" if diff is None else f"{diff:.12g}"})
    return out


def compare_harmonization_map() -> list[dict]:
    producer = load_tsv(EXEC / "run1/HARMONIZATION_MAP.tsv")
    def expected(cohort, field, token):
        if cohort == "TCGA" and field == "GRADE": return GRADE.get(token)
        if cohort == "TCGA" and field == "ICD_O_3_HISTOLOGY": return TCGA_HIST.get(token, (None, None))[0]
        if cohort == "CPTAC_Discovery" and field == "HISTOLOGIC_GRADE_FIGO": return GRADE.get(token)
        if cohort == "CPTAC_Discovery" and field == "HISTOLOGIC_TYPE": return "endometrioid" if token == "Endometrioid" else "non_endometrioid" if token == "Serous" else None
        if field == "PDC_morphology": return "endometrioid" if token == "8380/3" else None
        return None
    rows = []
    for p in producer:
        iv = expected(p["cohort"], p["source_field"], p["raw_token"]) or "MISSING"
        pv = p["proposed_category"]
        rows.append({"artifact": "HARMONIZATION_MAP.tsv", "row_key": "|".join([p["cohort"], p["source_field"], p["raw_token"]]), "field": "proposed_category", "independent": iv, "producer": pv, "match": str(iv == pv).upper(), "absolute_difference": "0" if iv == pv else "NA"})
    return rows


def compare_ledger(rows: list[dict]) -> list[dict]:
    producer = {(r["cohort"], r["patient_id"]): r for r in load_tsv(EXEC / "run1/LINKAGE_LEDGER.tsv")}
    fields = ["analytic_sample_id", "aliquot_id", "subtype", "source_row_id", "grade_raw", "grade_secondary_raw", "grade_proposed", "histology_raw", "histology_secondary_raw", "histology_proposed", "histology_three_proposed", "link_status", "duplicate_count", "grade_conflict", "histology_conflict", "exclusion_reason"]
    comparison = []
    for r in rows:
        key = (r["cohort"], r["patient_id"])
        p = producer[key]
        expected = {**r, "grade_proposed": r["grade"], "histology_proposed": r["histology"], "histology_three_proposed": r["histology_three"]}
        for field in fields:
            ev = str(na(expected[field]))
            pv = p[field] if p[field] != "" else ""
            comparison.append({"artifact": "LINKAGE_LEDGER.tsv", "row_key": "|".join(key), "field": field, "independent": ev, "producer": pv, "match": str(ev == pv).upper(), "absolute_difference": "0" if ev == pv else "NA"})
    return comparison


def compare_diagnostics(labelled: list[dict]) -> list[dict]:
    producer = {(r["cohort"], r["node"]): r for r in load_tsv(EXEC / "run1/DESIGN_DIAGNOSTICS.tsv")}
    fields = ["n", "p", "rank", "residual_df", "svd_tolerance", "condition_number", "max_one_df_vif", "subtype_adjusted_gvif", "clinical_adjusted_gvif", "cramers_v", "max_leverage", "fraction_hat_gt_3p_n", "support_state"]
    comparison = []
    for r in labelled:
        p = producer[(r["cohort"], r["node"])]
        for field in fields:
            iv = r.get(field)
            pv = p[field]
            if iv is None:
                match = pv == "NA"
                diff = 0.0 if match else None
                ivs = "NA"
            elif field == "support_state":
                ivs = str(iv); match = ivs == pv; diff = 0.0 if match else None
            else:
                ivs = str(iv)
                try:
                    diff = abs(float(iv) - float(pv)); match = diff <= 5e-9 * max(1.0, abs(float(iv)))
                except ValueError:
                    diff = None; match = ivs == pv
            comparison.append({"artifact": "DESIGN_DIAGNOSTICS.tsv", "row_key": f'{r["cohort"]}|{r["node"]}', "field": field, "independent": ivs, "producer": pv, "match": str(match).upper(), "absolute_difference": "NA" if diff is None else f"{diff:.12g}"})
    return comparison


def hash_inventory_audit() -> dict:
    rows = load_tsv(EXEC / "CHECKSUM_INVENTORY.tsv")
    checked = []
    for row in rows:
        rel = pathlib.Path(row["relative_path"])
        path = TASK / rel if rel.parts and rel.parts[0] == "inputs" else EXEC / rel
        actual_sha = sha_path(path)
        actual_size = path.stat().st_size
        checked.append({"relative_path": row["relative_path"], "resolved_path": str(path), "exists": path.exists(), "sha_match": actual_sha == row["sha256"], "size_match": actual_size == int(row["size_bytes"])})
    manifest = load_json(EXEC / "REPRODUCIBILITY_MANIFEST.json")
    refs = []
    for name, expected in manifest["frozen_protocol_hashes"].items(): refs.append({"reference": f"protocol:{name}", "match": sha_path(TASK / name) == expected})
    for name, expected in manifest["acquired_input_hashes"].items():
        mapping = {"tcga_cbio_grade":"tcga_cbio_grade.json","tcga_cbio_icdo_histology":"tcga_cbio_icdo_histology.json","cptac_discovery_cbio_grade":"cptac_discovery_grade.json","cptac_discovery_cbio_histology":"cptac_discovery_histology.json","pdc_PDC000125":"pdc_PDC000125_grade_histology.json","pdc_PDC000439":"pdc_PDC000439_grade_histology.json"}
        refs.append({"reference": f"input:{name}", "match": sha_path(TASK / "inputs" / mapping[name]) == expected})
    for name, expected in manifest["code_hashes"].items(): refs.append({"reference": f"code:{name}", "match": sha_path(EXEC / name) == expected})
    for name, expected in manifest["canonical_output_hashes"].items():
        refs.append({"reference": f"run1:{name}", "match": sha_path(EXEC / "run1" / name) == expected})
        refs.append({"reference": f"run2:{name}", "match": sha_path(EXEC / "run2" / name) == expected})
    refs += [{"reference": "input_lock", "match": sha_path(EXEC / "INPUTS_LOCK.json") == manifest["input_lock_sha256"]}, {"reference": "static_scan", "match": sha_path(EXEC / "STATIC_SCAN_REPORT.json") == manifest["static_scan"]["sha256"]}]
    return {"allowlisted_12": [{"path": p, "expected": h, "actual": sha_path(p), "match": sha_path(p) == h} for p,h in ALLOW], "protocol_4": [{"path": name, "expected": h, "actual": sha_path(TASK/name), "match": sha_path(TASK/name)==h} for name,h in PROTO.items()], "acquired_6": [{"path": name, "expected": h, "actual": sha_path(TASK/"inputs"/name), "match": sha_path(TASK/"inputs"/name)==h} for name,h in ACQUIRED.items()], "inventory_sha256": sha_path(EXEC / "CHECKSUM_INVENTORY.tsv"), "inventory_rows": len(rows), "inventory_all_match": all(x["exists"] and x["sha_match"] and x["size_match"] for x in checked), "inventory_path_resolution": "inputs/* relative to task root; all other rows relative to phase1_execution", "inventory": checked, "manifest_reference_count": len(refs), "manifest_all_references_match": all(x["match"] for x in refs), "manifest_references": refs}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = build_rows()
    correct, labelled = diagnostics(rows)
    ledger_cmp = compare_ledger(rows)
    diag_cmp = compare_diagnostics(labelled)
    tables = reconstructed_tables(rows)
    for name, table_rows in tables.items():
        fields = list(table_rows[0])
        write_tsv(OUT / ("INDEPENDENT_" + name), table_rows, fields)
    table_cmp = compare_reconstructed_tables(tables)
    map_cmp = compare_harmonization_map()
    comparisons = ledger_cmp + diag_cmp + table_cmp + map_cmp
    write_tsv(OUT / "FIELD_LEVEL_COMPARISON.tsv", comparisons, ["artifact", "row_key", "field", "independent", "producer", "match", "absolute_difference"])
    ledger_out = []
    for r in rows:
        ledger_out.append({**r, "grade_proposed": na(r["grade"]), "histology_proposed": na(r["histology"]), "histology_three_proposed": na(r["histology_three"]), "grade_raw": na(r["grade_raw"]), "grade_secondary_raw": na(r["grade_secondary_raw"]), "histology_raw": na(r["histology_raw"]), "histology_secondary_raw": na(r["histology_secondary_raw"]), "analytic_sample_id": na(r["analytic_sample_id"]), "aliquot_id": na(r["aliquot_id"])})
    ledger_fields = ["cohort", "patient_id", "analytic_sample_id", "aliquot_id", "subtype", "source_row_id", "grade_raw", "grade_secondary_raw", "grade_proposed", "histology_raw", "histology_secondary_raw", "histology_proposed", "histology_three_proposed", "link_status", "duplicate_count", "grade_conflict", "histology_conflict", "exclusion_reason"]
    write_tsv(OUT / "INDEPENDENT_LINKAGE_LEDGER.tsv", ledger_out, ledger_fields)
    diag_fields = ["cohort", "node", "n", "p", "rank", "residual_df", "svd_tolerance", "singular_values", "condition_number", "max_one_df_vif", "subtype_adjusted_gvif", "clinical_adjusted_gvif", "cramers_v", "near_zero_variance", "max_leverage", "fraction_hat_gt_3p_n", "support_state"]
    write_tsv(OUT / "INDEPENDENT_DESIGN_DIAGNOSTICS.tsv", [{k: na(v) for k,v in r.items()} for r in correct], diag_fields)
    summary = aggregate_summary(rows, correct)
    summary["field_comparison"] = {"ledger_fields_compared": len(ledger_cmp), "ledger_mismatches": sum(r["match"] == "FALSE" for r in ledger_cmp), "diagnostic_fields_compared": len(diag_cmp), "diagnostic_mismatches": sum(r["match"] == "FALSE" for r in diag_cmp), "aggregate_and_map_fields_compared": len(table_cmp) + len(map_cmp), "aggregate_and_map_mismatches": sum(r["match"] == "FALSE" for r in table_cmp + map_cmp), "total_fields_compared": len(comparisons), "total_mismatches": sum(r["match"] == "FALSE" for r in comparisons), "max_numeric_difference_matching_fields": max(float(r["absolute_difference"]) for r in comparisons if r["match"] == "TRUE" and r["absolute_difference"] != "NA")}
    (OUT / "INDEPENDENT_SUMMARY.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="ascii")
    (OUT / "INPUT_SOURCE_HASH_AUDIT.json").write_text(json.dumps(hash_inventory_audit(), indent=2, sort_keys=True) + "\n", encoding="ascii")
    print(json.dumps({"rows": len(rows), "ledger_mismatches": summary["field_comparison"]["ledger_mismatches"], "diagnostic_mismatches": summary["field_comparison"]["diagnostic_mismatches"], "global": summary["mechanical_state"]["global"]}, sort_keys=True))


if __name__ == "__main__":
    main()
