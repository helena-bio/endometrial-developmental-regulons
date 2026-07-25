#!/usr/bin/env python3
"""
TASK-029 RESULT FREEZE -- assemble the frozen result tables from primary + sensitivity JSON.
Writes: primary per-target table (TSV), per-stratum effects table (TSV), meta table (TSV),
sample-flow ledger (TSV), model-diagnostics (JSON), consolidated verdict (JSON/MD).
No new computation; pure assembly from sealed result JSON. Git HEAD unchanged.
"""
import os
import json
import numpy as np

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
RESULTS = os.environ.get(
    "CPTAC_RESULTS_DIR",
    os.path.join(BASE, "results"),
)

pr = json.load(open(os.path.join(RESULTS, "primary_results.json")))
sn = json.load(open(os.path.join(RESULTS, "sensitivity_results.json")))
lg = json.load(open(os.path.join(INTER, "acq02_join_ledger.json")))
verify = json.load(open(os.path.join(INTER, "acq03_download_verify.json")))

C1_T = ["LHX1", "PAX8"]; C2_T = ["GATA2", "HOXA9", "SOX9", "WT1"]
ORDER = C2_T + C1_T

# ---- primary per-target table ----
with open(os.path.join(RESULTS, "table_primary_per_target.tsv"), "w") as f:
    f.write("target\tfamily\tcontrast\tpred_direction\td_meta\tmeta_ci_lo\tmeta_ci_hi\t"
            "meta_ci_excl0\tBH_q_of_6\tdirection_match\tabs_d_ge_0.5\tstratum_consistent\t"
            "opposite_direction_veto\tstatus\n")
    for tf in ORDER:
        m = pr["meta"][tf]; v = pr["per_target_verdict"][tf]
        preddir = "p53abn_enrichment(d>0)" if v["pred_sign"] == 1 else "POLE+MMRd_depletion(d<0)"
        f.write(f"{tf}\t{v['family']}\t{v['contrast']}\t{preddir}\t{m['d_meta']:.4f}\t"
                f"{m['ci_lo']:.4f}\t{m['ci_hi']:.4f}\t{m['meta_ci_excludes_0']}\t"
                f"{m['BH_q_of_6']:.6f}\t{v['direction_match']}\t{v['abs_d_ge_0.5']}\t"
                f"{v['stratum_consistent']}\t{m['opposite_direction_veto']}\t{v['status']}\n")

# ---- per-stratum effects table ----
with open(os.path.join(RESULTS, "table_per_stratum_effects.tsv"), "w") as f:
    f.write("target\tcontrast\tstratum\td\tboot_se\tboot_ci_lo\tboot_ci_hi\tn_boot_used\tdirection\n")
    for tf in ORDER:
        m = pr["meta"][tf]
        for strat in ["Discovery", "Confirmatory"]:
            s = m["strata"][strat]
            direction = "neg" if s["d"] < 0 else "pos"
            f.write(f"{tf}\t{m['contrast']}\t{strat}\t{s['d']:.4f}\t{s['se']:.4f}\t"
                    f"{s['ci_lo']:.4f}\t{s['ci_hi']:.4f}\t{s['n_boot_used']}\t{direction}\n")

# ---- meta table ----
with open(os.path.join(RESULTS, "table_meta_fixedeffect.tsv"), "w") as f:
    f.write("target\tcontrast\td_disc\td_conf\td_meta\tse_meta\tz\tp_two_sided\t"
            "meta_ci_lo\tmeta_ci_hi\topposite_direction_veto\tBH_q_of_6\n")
    for tf in ORDER:
        m = pr["meta"][tf]
        f.write(f"{tf}\t{m['contrast']}\t{m['strata']['Discovery']['d']:.4f}\t"
                f"{m['strata']['Confirmatory']['d']:.4f}\t{m['d_meta']:.4f}\t{m['se_meta']:.4f}\t"
                f"{m['z']:.4f}\t{m['p_two_sided']:.3e}\t{m['ci_lo']:.4f}\t{m['ci_hi']:.4f}\t"
                f"{m['opposite_direction_veto']}\t{m['BH_q_of_6']:.6f}\n")

# ---- M1 (F2 sensitivity-only) table ----
with open(os.path.join(RESULTS, "table_M1_F2_sensitivity.tsv"), "w") as f:
    f.write("stratum\tn\tsigma_resid\tcontrast\td\tboot_ci_lo\tboot_ci_hi\tboot_se\n")
    for strat in ["Discovery", "Confirmatory"]:
        m1 = pr["strata"][strat]["M1_model"]
        for c in ["C1", "C2", "C3"]:
            cc = m1["contrasts"][c]
            f.write(f"{strat}\t{m1['n']}\t{m1['sigma_resid']:.4f}\t{c}\t{cc['d']:.4f}\t"
                    f"{cc['d_ci_lo']:.4f}\t{cc['d_ci_hi']:.4f}\t{cc['d_boot_se']:.4f}\n")

# ---- sample-flow ledger ----
with open(os.path.join(RESULTS, "sample_flow_ledger.tsv"), "w") as f:
    f.write("stratum\tsubtype\tcase_submitter_id\tchosen_sample\tfile_id\tmd5\tsha256\t"
            "n_primary_candidate_files\n")
    for c in sorted(lg["pick"]):
        p = lg["pick"][c]
        aud = lg["dedup_audit"][c]
        f.write(f"{p['stratum']}\t{p['subtype']}\t{c}\t{p['sample_submitter_id']}\t"
                f"{p['file_id']}\t{p['md5sum']}\t{verify['per_case'][c]['sha256']}\t"
                f"{aud['n_primary_candidates_distinct_files']}\n")

# ---- sensitivity (a) table ----
sa = sn["sensitivity_a_estimate_purity"]
with open(os.path.join(RESULTS, "table_sensitivity_estimate_purity.tsv"), "w") as f:
    f.write("target\td_meta\tmeta_ci_lo\tmeta_ci_hi\tmeta_ci_excl0\tBH_q_of_6\t"
            "opposite_direction_veto\n")
    for tf in ORDER:
        m = sa["meta"][tf]
        f.write(f"{tf}\t{m['d_meta']:.4f}\t{m['ci_lo']:.4f}\t{m['ci_hi']:.4f}\t{m['ci_excl0']}\t"
                f"{m['BH_q_of_6']:.6f}\t{m['opposite_veto']}\n")

# ---- model diagnostics ----
diag = {
    "primary_model": pr["model"],
    "per_stratum_sigma_resid": {
        s: {"M1": pr["strata"][s]["M1_model"]["sigma_resid"],
            **{tf: pr["strata"][s]["m3_model"][tf]["sigma_resid"] for tf in ORDER}}
        for s in ["Discovery", "Confirmatory"]},
    "per_stratum_omnibus_p_F": {
        s: {"M1": pr["strata"][s]["M1_model"]["omnibus_p_F"],
            **{tf: pr["strata"][s]["m3_model"][tf]["omnibus_p_F"] for tf in ORDER}}
        for s in ["Discovery", "Confirmatory"]},
    "estimate_purity_collinearity": sa["collinearity"],
    "estimate_purity_model_instability_VIF_gt10": sa["model_instability_VIF_gt10"],
    "m3_edges_used_per_stratum": {s: pr["strata"][s]["m3_edges_used"]
                                  for s in ["Discovery", "Confirmatory"]},
    "allzero_genes_per_stratum": {s: pr["strata"][s]["allzero_n"]
                                  for s in ["Discovery", "Confirmatory"]},
}


def clean(o):
    if isinstance(o, dict):
        return {k: clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(x) for x in o]
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, float) and not np.isfinite(o):
        return None
    return o


json.dump(clean(diag), open(os.path.join(RESULTS, "model_diagnostics.json"), "w"), indent=2)
print("frozen tables written to", RESULTS)
