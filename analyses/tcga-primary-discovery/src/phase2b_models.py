"""
TASK-028 Freeze B v3 -- PHASE 2b: LOCKED PRIMARY MODELS + SENSITIVITIES + D_t.
Spec: B2.4/B2.5/B2.7/B2.8/B2.12. Deterministic seeds (B2.11). No tuning. ASCII only.

CONFIGS (each = score-source x sample-mask x covariate-set):
  PRIMARY: M1 ssGSEA + 20 M3 signed-VIPER regulons, model on CPE complete-case n=506,
           covariates = [M4, purity_CPE, composition].
  SENS no-purity: same scores, n=507, covariates = [M4, composition].
  SENS ABSOLUTE:  same scores, n=502, covariates = [M4, purity_ABSOLUTE, composition].
  SENS compact-M1: compact-M1 ssGSEA, n=506 CPE model.
  SENS signed-ssGSEA: signed-ssGSEA regulon scores, n=506 CPE model (M3 only).
  SENS S1: S1 regulon scores, n=506 CPE model (M3 only).
  SENS S2: S2 regulon scores, n=506 CPE model (M3 only).

FAMILIES (B2.7): F1 = omnibus tests (M1 + 20 regulons) = 21 modules PRIMARY config.
  F2 = module x {C1,C2,C3} for F1-gated modules. F3 = D_t. F4 = module x C1 equivalence.
  BH q<=0.05 confirmatory. C1-C3 ALWAYS reported (effect + CI) regardless of gate.

H1 floor |d|>=0.5. H3 equivalence CONFIRMATORY for C1 only (TOST/CI-in-SESOI, SESOI |d|<0.30).
D_t (H4): residualization A ~ z(TFexpr) + M4 + purity + composition (NO subtype), then the
  same full-model contrast on residual R. Declared iff |D|>=0.5 AND boot CI excludes 0 AND
  perm p passes F3 AND direction condition (opposite-sign non-trivial TF-expr contrast).
"""
import os
import sys
import json
import numpy as np
import pandas as pd
from datetime import datetime, timezone

sys.path.insert(0, os.path.dirname(__file__))
from common_v3 import (MASTER_SEED, GIT_HEAD, substep_seed, INTER, RESULTS_V3, SUBTYPES)
import model_engine as ME

N_BOOT = 2000
N_PERM = 2000
SESOI = 0.30       # H3 equivalence margin |d|<0.30
FLOOR = 0.50       # H1 biological floor |d|>=0.50


def utcnow():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


print("=== PHASE 2b (models) start", utcnow(), "===", flush=True)

# ---- load scores + covariates ----
sc = np.load(os.path.join(INTER, "scores_v3.npz"), allow_pickle=True)
cov = pd.read_csv(os.path.join(INTER, "covariates_v3.tsv"), sep="\t")
patient_order = list(sc["patient_order"])
assert list(cov["patient_barcode"]) == patient_order, "covariate/score patient order mismatch"

subtype = cov["subtype"].to_numpy()
M4 = cov["M4_prolif"].to_numpy(float)
CPE = cov["purity_CPE"].to_numpy(float)
ABS = cov["purity_ABSOLUTE"].to_numpy(float)
COMP = cov["composition"].to_numpy(float)
mask506 = cov["cpe_complete_case"].to_numpy() == 1
mask507 = np.ones(507, dtype=bool)
mask502 = cov["absolute_complete_case"].to_numpy() == 1
CORE_TFS = list(sc["CORE_TFS"])

cov_dict = {"M4": M4, "purity_CPE": CPE, "purity_ABSOLUTE": ABS, "composition": COMP}

# ---- config table ----
# each entry: config_name -> (cov_names, sample_mask, dict module_key -> score_array)
def m3_scores(prefix):
    return {f"M3_{tf}": sc[f"{prefix}__{tf}"] for tf in CORE_TFS}

CONFIGS = {}
# PRIMARY (M1 + M3 primary) on n506 CPE
CONFIGS["PRIMARY"] = {
    "cov_names": ["M4", "purity_CPE", "composition"], "mask": mask506,
    "scores": {"M1": sc["M1"], **m3_scores("M3primary")},
    "label": "primary (CPE n=506; M1 ssGSEA + signed-VIPER regulons)"}
# no-purity n507
CONFIGS["SENS_nopurity"] = {
    "cov_names": ["M4", "composition"], "mask": mask507,
    "scores": {"M1": sc["M1"], **m3_scores("M3primary")},
    "label": "sensitivity no-purity (n=507; drop b_purity)"}
# ABSOLUTE n502
CONFIGS["SENS_absolute"] = {
    "cov_names": ["M4", "purity_ABSOLUTE", "composition"], "mask": mask502,
    "scores": {"M1": sc["M1"], **m3_scores("M3primary")},
    "label": "sensitivity ABSOLUTE-purity (n=502)"}
# compact-M1 (M1 only) n506 CPE
CONFIGS["SENS_compactM1"] = {
    "cov_names": ["M4", "purity_CPE", "composition"], "mask": mask506,
    "scores": {"M1": sc["compactM1"]},
    "label": "sensitivity compact-M1 ssGSEA (CPE n=506)"}
# signed-ssGSEA (M3 only) n506 CPE
CONFIGS["SENS_signedssgsea"] = {
    "cov_names": ["M4", "purity_CPE", "composition"], "mask": mask506,
    "scores": m3_scores("M3signedssgsea"),
    "label": "sensitivity signed-ssGSEA regulons (CPE n=506)"}
# S1 (M3 only) n506 CPE
CONFIGS["SENS_S1"] = {
    "cov_names": ["M4", "purity_CPE", "composition"], "mask": mask506,
    "scores": m3_scores("M3S1"),
    "label": "sensitivity M3 S1 flag-split (all 79 both-flagged; CPE n=506)"}
# S2 (M3 only) n506 CPE
CONFIGS["SENS_S2"] = {
    "cov_names": ["M4", "purity_CPE", "composition"], "mask": mask506,
    "scores": m3_scores("M3S2"),
    "label": "sensitivity M3 S2 exclude-both-flagged (CPE n=506)"}


def run_config(cfg_name, cfg, do_boot_perm):
    """Fit all modules in a config; return dict module_key -> result.
    do_boot_perm: if True, run bootstrap CIs + permutation p (expensive)."""
    cov_names = cfg["cov_names"]; mask = cfg["mask"]
    res = {}
    for mkey, y in cfg["scores"].items():
        if np.sum(~np.isfinite(np.asarray(y)[mask])) > 0:
            res[mkey] = {"error": "nonfinite_score_in_masked_sample",
                         "n_nonfinite": int(np.sum(~np.isfinite(np.asarray(y)[mask])))}
            continue
        fit = ME.fit_module(y, subtype, cov_dict, cov_names, mask)
        if do_boot_perm:
            bseed = substep_seed(f"boot__{cfg_name}__{mkey}")
            pseed = substep_seed(f"perm__{cfg_name}__{mkey}")
            fit["bootstrap"] = ME.bootstrap_contrasts(
                y, subtype, cov_dict, cov_names, mask, bseed, n_boot=N_BOOT)
            fit["permutation"] = ME.permutation_omnibus_and_contrasts(
                y, subtype, cov_dict, cov_names, mask, pseed, n_perm=N_PERM)
        res[mkey] = fit
    return res


# -----------------------------------------------------------------
# Run PRIMARY (full boot+perm) first -> it drives F1/F2 families.
# -----------------------------------------------------------------
print("[run] PRIMARY config (boot+perm)", utcnow(), flush=True)
primary = run_config("PRIMARY", CONFIGS["PRIMARY"], do_boot_perm=True)

# F1 family: omnibus p across the 21 primary modules (M1 + 20 regulons)
f1_p = {mkey: primary[mkey]["omnibus"]["p_F"] for mkey in primary if "omnibus" in primary[mkey]}
f1_q = ME.bh_fdr(f1_p)
for mkey in primary:
    if "omnibus" in primary[mkey]:
        primary[mkey]["omnibus"]["BH_q_F1"] = f1_q[mkey]
        primary[mkey]["omnibus"]["F1_gate_pass"] = bool(f1_q[mkey] <= 0.05)

# F2 family: module x contrast tests (ONLY for F1-gated modules), using permutation p
# per-contrast. Confirmatory q<=0.05.
f2_p = {}
for mkey in primary:
    if primary[mkey].get("omnibus", {}).get("F1_gate_pass"):
        for c in ["C1", "C2", "C3"]:
            f2_p[f"{mkey}::{c}"] = primary[mkey]["permutation"]["perm_p_contrast"][c]
f2_q = ME.bh_fdr(f2_p)
for key, q in f2_q.items():
    mkey, c = key.split("::")
    primary[mkey]["contrasts"][c]["BH_q_F2"] = q
    primary[mkey]["contrasts"][c]["perm_p"] = f2_p[key]

# F4 family: module x C1 equivalence (TOST) for F1-gated modules.
# TOST on d for C1: two one-sided tests against +/-SESOI using bootstrap CI-in-SESOI rule.
# Clean-null (CAT5) = the 90% CI on d (equiv. two-sided 95% for each one-sided) lies ENTIRELY
# within (-SESOI, +SESOI). We use the bootstrap 95% CI on d and require it inside the SESOI.
def equivalence_C1(fit):
    d = fit["contrasts"]["C1"]["d"]
    ci = fit.get("bootstrap", {}).get("C1", {})
    lo = ci.get("d_ci_lo"); hi = ci.get("d_ci_hi")
    if lo is None or hi is None or not np.isfinite(lo) or not np.isfinite(hi):
        return {"tost_pass": None, "d": d, "d_ci_lo": lo, "d_ci_hi": hi}
    inside = (lo > -SESOI) and (hi < SESOI)
    return {"tost_pass": bool(inside), "d": d, "d_ci_lo": lo, "d_ci_hi": hi,
            "SESOI": SESOI}

for mkey in primary:
    if primary[mkey].get("omnibus", {}).get("F1_gate_pass"):
        primary[mkey]["equivalence_C1"] = equivalence_C1(primary[mkey])

# -----------------------------------------------------------------
# B3 category assignment per (module, contrast) -- CONFIRMATORY gate (B2.12).
# CAT1 enrichment / CAT2 depletion require: F1 gate pass AND |d|>=0.5 AND F2 q<=0.05.
# CAT5 clean-null via H3 equivalence (C1 only). Else DESCRIPTIVE/EXPLORATORY or inconclusive.
# -----------------------------------------------------------------
def categorize(fit, mkey):
    out = {}
    gate = fit.get("omnibus", {}).get("F1_gate_pass", False)
    for c in ["C1", "C2", "C3"]:
        cc = fit["contrasts"][c]
        d = cc["d"]
        q = cc.get("BH_q_F2")
        boot = fit.get("bootstrap", {}).get(c, {})
        d_lo = boot.get("d_ci_lo"); d_hi = boot.get("d_ci_hi")
        ci_excl0 = (d_lo is not None and d_hi is not None and
                    np.isfinite(d_lo) and np.isfinite(d_hi) and (d_lo > 0 or d_hi < 0))
        cat = "DESCRIPTIVE"
        reason = []
        if gate and q is not None and q <= 0.05 and abs(d) >= FLOOR and ci_excl0:
            cat = "CAT1_ENRICHMENT" if d > 0 else "CAT2_DEPLETION"
        else:
            # equivalence clean-null only C1, confirmatory-admissible
            if c == "C1":
                eq = fit.get("equivalence_C1", {})
                if eq.get("tost_pass"):
                    cat = "CAT5_CLEAN_NULL"
                    reason.append("H3 equivalence (C1) CI-in-SESOI")
            if cat == "DESCRIPTIVE":
                if not gate:
                    reason.append("omnibus F1 gate FAIL -> descriptive/exploratory")
                else:
                    if q is None or q > 0.05:
                        reason.append(f"F2 q>{0.05}")
                    if abs(d) < FLOOR:
                        reason.append(f"|d|<{FLOOR}")
                    if not ci_excl0:
                        reason.append("bootstrap CI includes 0")
        out[c] = {"category": cat, "d": d, "q_F2": q,
                  "ci_excludes_0": bool(ci_excl0), "reason": "; ".join(reason)}
    return out

for mkey in primary:
    if "contrasts" in primary[mkey]:
        primary[mkey]["b3_category"] = categorize(primary[mkey], mkey)

# -----------------------------------------------------------------
# D_t (H4) TF-MODULE DISCORDANCE -- for each of the 20 regulons, PRIMARY config (n=506 CPE).
# A = signed regulon activity z-standardized across the fitted sample; residualize on
# z(TFexpr) + M4 + purity + composition (NO subtype). D = full-model contrast on residual R.
# Direction: e_X = same contrast on covariate-adjusted subtype means of z(TFexpr).
# -----------------------------------------------------------------
def zscore(v, mask):
    x = np.asarray(v, float)[mask]
    mu = x.mean(); sd = x.std(ddof=1)
    return (x - mu) / sd if sd > 0 else np.zeros_like(x)


def run_dt(tf):
    mask = mask506
    idx = np.where(mask)[0]
    A_full = np.asarray(sc[f"M3primary__{tf}"], float)
    A = zscore(A_full, mask)                          # z(A) across the 506
    zTF = zscore(sc[f"TFexprLog__{tf}"], mask)        # z(TFexpr) across the 506
    subt = subtype[idx]
    cov_local = {"zTF": zTF, "M4": M4[idx], "purity_CPE": CPE[idx], "composition": COMP[idx]}
    # residualization model: A ~ intercept + zTF + M4 + purity + composition (NO subtype)
    Xr = np.column_stack([np.ones(len(idx)), cov_local["zTF"], cov_local["M4"],
                          cov_local["purity_CPE"], cov_local["composition"]])
    br, resid, sr, xi, dof, s2 = ME.fit_ols(A, Xr)
    R = resid                                         # residual in SD(A) units already (A is z)
    # D = full-model contrast on residual: fit R ~ subtype + covariates, EMMs, contrast/sigma
    cov_R = {"M4": M4, "purity_CPE": CPE, "composition": COMP}
    # need R aligned to full-index; embed into 507 with NaN off-mask
    Rfull = np.full(507, np.nan); Rfull[idx] = R
    dt_fit = ME.fit_module(Rfull, subtype, cov_R, ["M4", "purity_CPE", "composition"], mask)
    # bootstrap + permutation on D (F3)
    bseed = substep_seed(f"dt_boot__{tf}")
    pseed = substep_seed(f"dt_perm__{tf}")
    dt_boot = ME.bootstrap_contrasts(Rfull, subtype, cov_R,
                                     ["M4", "purity_CPE", "composition"], mask, bseed, n_boot=N_BOOT)
    dt_perm = ME.permutation_omnibus_and_contrasts(Rfull, subtype, cov_R,
                                     ["M4", "purity_CPE", "composition"], mask, pseed, n_perm=N_PERM)
    # e_X: contrast on covariate-adjusted subtype means of z(TFexpr)
    zTFfull = np.full(507, np.nan); zTFfull[idx] = zTF
    ex_fit = ME.fit_module(zTFfull, subtype, cov_R, ["M4", "purity_CPE", "composition"], mask)
    exseed = substep_seed(f"dt_ex_boot__{tf}")
    ex_boot = ME.bootstrap_contrasts(zTFfull, subtype, cov_R,
                                     ["M4", "purity_CPE", "composition"], mask, exseed, n_boot=N_BOOT)
    out = {"tf": tf, "n": int(mask.sum()), "residual_model_R2_zTF_beta": float(br[1]),
           "sigma_resid_R": dt_fit["sigma_resid"], "D": {}, "e_X": {}, "declared": {}}
    for c in ["C1", "C2", "C3"]:
        D = dt_fit["contrasts"][c]["d"]
        D_lo = dt_boot[c]["d_ci_lo"]; D_hi = dt_boot[c]["d_ci_hi"]
        D_ci_excl0 = (np.isfinite(D_lo) and np.isfinite(D_hi) and (D_lo > 0 or D_hi < 0))
        pperm = dt_perm["perm_p_contrast"][c]
        eX = ex_fit["contrasts"][c]["d"]
        eX_lo = ex_boot[c]["d_ci_lo"]; eX_hi = ex_boot[c]["d_ci_hi"]
        eX_ci_excl0 = (np.isfinite(eX_lo) and np.isfinite(eX_hi) and (eX_lo > 0 or eX_hi < 0))
        out["D"][c] = {"D": D, "ci_lo": D_lo, "ci_hi": D_hi, "ci_excludes_0": bool(D_ci_excl0),
                       "perm_p": pperm}
        out["e_X"][c] = {"e_X": eX, "ci_lo": eX_lo, "ci_hi": eX_hi,
                         "ci_excludes_0": bool(eX_ci_excl0)}
    return out


print("[run] D_t H4 for 20 regulons", utcnow(), flush=True)
dt_results = {}
for tf in CORE_TFS:
    dt_results[tf] = run_dt(tf)
    print("  D_t", tf, "done", utcnow(), flush=True)

# F3 family: BH across all D_t contrast permutation p (20 TF x 3 contrasts = 60)
f3_p = {}
for tf in CORE_TFS:
    for c in ["C1", "C2", "C3"]:
        f3_p[f"{tf}::{c}"] = dt_results[tf]["D"][c]["perm_p"]
f3_q = ME.bh_fdr(f3_p)
for tf in CORE_TFS:
    for c in ["C1", "C2", "C3"]:
        q = f3_q[f"{tf}::{c}"]
        d = dt_results[tf]["D"][c]
        d["BH_q_F3"] = q
        # DECLARE discordance iff ALL FOUR: |D|>=0.5; boot CI excl 0; perm passes F3 (q<=0.05);
        # direction: |e_X|>=0.5 AND e_X CI excl 0 AND sign(D) OPPOSITE sign(e_X).
        eX = dt_results[tf]["e_X"][c]
        cond_i = abs(d["D"]) >= FLOOR
        cond_ii = d["ci_excludes_0"]
        cond_iii = (q is not None and q <= 0.05)
        cond_iv = (abs(eX["e_X"]) >= FLOOR and eX["ci_excludes_0"] and
                   np.sign(d["D"]) == -np.sign(eX["e_X"]) and eX["e_X"] != 0 and d["D"] != 0)
        if cond_i and cond_ii and cond_iii and cond_iv:
            verdict = "DISCORDANCE_DECLARED"
        elif cond_i and cond_ii and cond_iii and not cond_iv:
            verdict = "regulon activity unexplained by TF abundance"
        else:
            verdict = "no discordance (criteria not met)"
        dt_results[tf]["declared"][c] = {
            "verdict": verdict,
            "cond_i_absD_ge_0.5": bool(cond_i),
            "cond_ii_bootCI_excl0": bool(cond_ii),
            "cond_iii_permF3_q_le_0.05": bool(cond_iii),
            "cond_iv_direction": bool(cond_iv),
            "BH_q_F3": q}

# -----------------------------------------------------------------
# Run remaining SENSITIVITY configs (boot+perm too, for direction/mag/category concordance).
# -----------------------------------------------------------------
sens_results = {}
for cfg_name in ["SENS_nopurity", "SENS_absolute", "SENS_compactM1",
                 "SENS_signedssgsea", "SENS_S1", "SENS_S2"]:
    print(f"[run] {cfg_name} (boot+perm)", utcnow(), flush=True)
    r = run_config(cfg_name, CONFIGS[cfg_name], do_boot_perm=True)
    # F1 omnibus BH within the config's own family
    p = {mk: r[mk]["omnibus"]["p_F"] for mk in r if "omnibus" in r[mk]}
    q = ME.bh_fdr(p)
    for mk in r:
        if "omnibus" in r[mk]:
            r[mk]["omnibus"]["BH_q_F1"] = q[mk]
            r[mk]["omnibus"]["F1_gate_pass"] = bool(q[mk] <= 0.05)
    # F2 contrast BH within config
    p2 = {}
    for mk in r:
        if r[mk].get("omnibus", {}).get("F1_gate_pass"):
            for c in ["C1", "C2", "C3"]:
                p2[f"{mk}::{c}"] = r[mk]["permutation"]["perm_p_contrast"][c]
    q2 = ME.bh_fdr(p2)
    for key, qq in q2.items():
        mk, c = key.split("::")
        r[mk]["contrasts"][c]["BH_q_F2"] = qq
        r[mk]["contrasts"][c]["perm_p"] = p2[key]
    for mk in r:
        if "contrasts" in r[mk]:
            if r[mk].get("omnibus", {}).get("F1_gate_pass"):
                r[mk]["equivalence_C1"] = equivalence_C1(r[mk])
            r[mk]["b3_category"] = categorize(r[mk], mk)
    sens_results[cfg_name] = r

# -----------------------------------------------------------------
# Persist everything
# -----------------------------------------------------------------
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

allout = {
    "phase": "2b", "version": "v3", "git_HEAD": GIT_HEAD, "master_seed": MASTER_SEED,
    "n_boot": N_BOOT, "n_perm": N_PERM, "SESOI": SESOI, "H1_floor": FLOOR,
    "config_labels": {k: v["label"] for k, v in CONFIGS.items()},
    "primary": primary,
    "dt": dt_results,
    "sensitivities": sens_results,
    "families": {
        "F1_primary_omnibus_p": f1_p, "F1_primary_omnibus_q": f1_q,
        "F2_primary_contrast_p": f2_p, "F2_primary_contrast_q": f2_q,
        "F3_dt_perm_p": f3_p, "F3_dt_perm_q": f3_q},
    "timestamp_utc": utcnow(),
}
with open(os.path.join(RESULTS_V3, "phase2b_models.json"), "w") as f:
    json.dump(clean(allout), f, indent=2)
print("=== PHASE 2b done", utcnow(), "===", flush=True)
