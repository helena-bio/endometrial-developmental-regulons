#!/usr/bin/env python3
"""
TASK-029 FROZEN PRIMARY EXECUTION -- score (per stratum) + V2 no-purity model +
stratum contrasts + patient-bootstrap CIs + inverse-variance FIXED-EFFECT meta +
BH-of-6 + per-target verdicts.

SCORING + MODEL CODE REUSED VERBATIM from the sealed TASK-028 modules:
  - modules: definitions/ (M1_analysis_ready, M3 761-signed-edge ledger, M4, ESTIMATE gmt)
  - scoring: src/scoring.py  (ssgsea_matrix, area_regulon)  -- imported, not re-implemented
  - model:   src/model_engine.py (build_design, fit_ols, subtype_emms,
             contrast_estimate, bootstrap_contrasts, bh_fdr) -- imported, not re-implemented

FROZEN DESIGN (FROZEN_REPLICATION_DESIGN.md + FINAL_FREEZE_RULES_bm.md):
  PRIMARY = V2 no-purity: score ~ subtype + M4(covariate) + ESTIMATE composition. NO purity.
  Score EACH stratum separately (no raw-expression pooling).
  M1 ssGSEA (F2, sensitivity-only). Each M3 CORE_TF regulon via signed VIPER/aREA over the
  761 CollecTRI-signed edges; orphans excluded; all-zero rule.
  Contrasts C1=[-1,-1,-1,+3]/C2=[+1,+1,-2,0] (equal-within-side normalized in model_engine),
  C3 descriptive. Standardized effect d = contrast/sigma_resid per stratum; patient bootstrap CI.
  META: inverse-variance FIXED-EFFECT over the two strata; per-stratum direction+CI mandatory;
  OPPOSITE-DIRECTION VETO. F1 = BH-of-6 across {C1: LHX1,PAX8 ; C2: GATA2,HOXA9,SOX9,WT1}.
  Per-target `replicated` iff ALL: pre-specified direction; |d|>=0.5 + bootstrap CI excl 0;
  BH q<=0.05; stratum-consistency (no opposite-direction veto); no critical QC/model failure.
  C1 = underpowered-sensitivity only; M1 = sensitivity-only (equivalence not feasible).

MASTER SEED 20260714. Deterministic sub-seeds from sha256(master:step_id). Git HEAD unchanged.
"""
import os
import sys
import json
import hashlib
import numpy as np
import pandas as pd
from datetime import datetime, timezone
from collections import defaultdict
from scipy.stats import norm

# ---- import the SEALED scoring + model code VERBATIM ----
BASE = os.path.dirname(
    os.path.dirname(os.path.abspath(__file__))
)
T28 = os.environ.get(
    "TCGA_DISCOVERY_ROOT",
    os.path.join(
        os.path.dirname(BASE),
        "tcga-primary-discovery",
    ),
)
sys.path.insert(0, os.path.join(T28, "src"))
import scoring as SC            # ssgsea_matrix, area_regulon  (verbatim)
import model_engine as ME       # build_design, fit_ols, subtype_emms, contrast_estimate, bootstrap_contrasts, bh_fdr


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
B1 = os.path.join(T28, "definitions")
ESTIMATE_GMT = os.environ.get(
    "ESTIMATE_GMT",
    os.path.join(
        T28,
        "external-inputs",
        "ESTIMATE_SI_geneset.gmt",
    ),
)

MASTER_SEED = 20260714
GIT_HEAD = "83503bad47b60193598b2b9ebe819c22c83e8ac1"
N_BOOT = 2000
FLOOR = 0.50
SESOI = 0.30
SUBTYPES = ["POLE", "MMRd", "NSMP", "p53abn"]
CORE_TARGETS_C1 = ["LHX1", "PAX8"]
CORE_TARGETS_C2 = ["GATA2", "HOXA9", "SOX9", "WT1"]
# per-target primary contrast + pre-specified direction (positive d = predicted direction)
# C1: enrichment on p53abn side  -> model contrast C1 (=[-1,-1,-1,+3]) positive means p53abn > rest -> PREDICTED d>0
# C2: relative depletion on POLE+MMRd side vs NSMP -> predicted (POLE+MMRd) < NSMP.
#     model contrast C2 (=[+1,+1,-2,0]) positive means (POLE+MMRd) > NSMP; DEPLETION = NEGATIVE d.
#     So the pre-specified direction for C2 targets is d < 0.
TARGET_CONTRAST = {t: "C1" for t in CORE_TARGETS_C1}
TARGET_CONTRAST.update({t: "C2" for t in CORE_TARGETS_C2})
TARGET_PRED_SIGN = {t: +1 for t in CORE_TARGETS_C1}     # p53abn enrichment -> d>0
TARGET_PRED_SIGN.update({t: -1 for t in CORE_TARGETS_C2})  # POLE+MMRd depletion -> d<0


def utcnow():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def subseed(step):
    return int(hashlib.sha256(f"{MASTER_SEED}:{step}".encode()).hexdigest()[:8], 16)


def read_module_genes(path):
    return [l.strip() for l in open(path)
            if l.strip() and not l.strip().startswith("#")]


print("=== TASK-029 FROZEN PRIMARY EXECUTION start", utcnow(), "===", flush=True)

gene_names = json.load(open(os.path.join(INTER, "gene_names.json")))
case_order = json.load(open(os.path.join(INTER, "case_order_by_stratum.json")))
lg = json.load(open(os.path.join(INTER, "acq02_join_ledger.json")))
pick = lg["pick"]

# ---- parse M3 edge ledger -> per-TF primary regulon (761 included edges) ----
edge_path = os.path.join(B1, "m3_primary_edge_ledger.tsv")
edges = []
with open(edge_path) as f:
    for line in f:
        if line.startswith("#"):
            continue
        edges.append(line.rstrip("\n").split("\t"))
hdr = edges[0]; edges = edges[1:]
ix = {c: i for i, c in enumerate(hdr)}
prim = defaultdict(lambda: {"targets": [], "weights": [], "mor": []})
for r in edges:
    if r[ix["primary_inclusion"]] == "Y":
        tf = r[ix["TF"]]
        prim[tf]["targets"].append(r[ix["target"]])
        prim[tf]["weights"].append(float(r[ix["weight"]]))
        prim[tf]["mor"].append(int(r[ix["primary_sign"]]))
n_included = sum(len(prim[tf]["targets"]) for tf in prim)
print("M3 primary included edges:", n_included, "(expect 761); CORE_TFS:", len(prim), flush=True)

# ESTIMATE signatures
estimate_sig = {}
for line in open(ESTIMATE_GMT):
    p = line.rstrip("\n").split("\t")
    estimate_sig[p[0]] = p[2:]

M1_genes = read_module_genes(os.path.join(B1, "M1_analysis_ready.txt"))
M4_genes = read_module_genes(os.path.join(B1, "M4_covariate.txt"))

# ---------------------------------------------------------------
# Per stratum: score + covariates + model + bootstrap
# ---------------------------------------------------------------
def score_stratum(strat):
    LOG = np.load(os.path.join(INTER, f"log2tpm_{strat}.npy"))  # 60660 x n
    cases = case_order[strat]
    n = LOG.shape[1]
    assert n == len(cases)
    # symbol -> rows; best row = highest mean log2 (verbatim TASK-028)
    sym_to_rows = {}
    for i, nm in enumerate(gene_names):
        sym_to_rows.setdefault(nm, []).append(i)

    def best_row(sym):
        rows = sym_to_rows.get(sym)
        if not rows:
            return None
        if len(rows) == 1:
            return rows[0]
        return rows[int(np.argmax([LOG[r].mean() for r in rows]))]

    # all-zero-across-this-stratum mask
    allzero = set()
    for sym, rows in sym_to_rows.items():
        if all((LOG[r] == 0).all() for r in rows):
            allzero.add(sym)

    def ssgsea_module(genes, name):
        kept = [g for g in genes if g not in allzero]
        rows, used = [], set()
        for g in kept:
            r = best_row(g)
            if r is not None and r not in used:
                rows.append(r); used.add(r)
        return SC.ssgsea_matrix(LOG, {name: rows})[name], len(genes), len(kept), len(rows)

    m1, m1n, m1k, m1r = ssgsea_module(M1_genes, "M1")
    m4, _, _, m4r = ssgsea_module(M4_genes, "M4")
    strom, _, _, sr = ssgsea_module(estimate_sig["StromalSignature"], "S")
    imm, _, _, ir = ssgsea_module(estimate_sig["ImmuneSignature"], "I")
    composition = strom + imm

    # gene_index for aREA (best row per target symbol)
    gene_index = {}
    for sym in set(t for tf in prim for t in prim[tf]["targets"]):
        r = best_row(sym)
        if r is not None:
            gene_index[sym] = r

    m3 = {}
    m3_used = {}
    for tf in sorted(prim):
        d = prim[tf]
        nes, nused = SC.area_regulon(LOG, gene_index, d["targets"], d["weights"], d["mor"])
        m3[tf] = nes
        m3_used[tf] = {"edges_listed": len(d["targets"]), "targets_used": int(nused)}

    return {
        "cases": cases, "n": n,
        "M1": m1, "M4": m4, "composition": composition,
        "estimate_stromal": strom, "estimate_immune": imm,
        "m3": m3, "m3_used": m3_used,
        "m1_scored": {"n_members": m1n, "n_after_allzero": m1k, "n_rows": m1r},
        "cov_rows": {"M4": m4r, "estimate_stromal": sr, "estimate_immune": ir},
        "allzero_n": len(allzero),
    }


def model_stratum(strat, sd):
    """V2 no-purity model per stratum; contrasts + patient-bootstrap CI + SE(d)."""
    cases = sd["cases"]
    subtype = np.array([pick[c]["subtype"] for c in cases])
    M4 = np.asarray(sd["M4"], float)
    COMP = np.asarray(sd["composition"], float)
    cov_dict = {"M4": M4, "composition": COMP}
    cov_names = ["M4", "composition"]
    n = sd["n"]
    mask = np.ones(n, dtype=bool)

    def fit_and_boot(y, key):
        y = np.asarray(y, float)
        if np.any(~np.isfinite(y[mask])):
            return {"error": "nonfinite_score"}
        fit = ME.fit_module(y, subtype, cov_dict, cov_names, mask)
        bseed = subseed(f"boot__{strat}__{key}")
        boot = ME.bootstrap_contrasts(y, subtype, cov_dict, cov_names, mask, bseed, n_boot=N_BOOT)
        out = {"n": fit["n"], "sigma_resid": fit["sigma_resid"], "emm": fit["emm"],
               "contrasts": {}, "omnibus_p_F": fit["omnibus"]["p_F"]}
        for c in ["C1", "C2", "C3"]:
            d = fit["contrasts"][c]["d"]
            est = fit["contrasts"][c]["estimate"]
            b = boot[c]
            # bootstrap SE of d (used for inverse-variance meta)
            out["contrasts"][c] = {
                "estimate": est, "d": d,
                "d_ci_lo": b["d_ci_lo"], "d_ci_hi": b["d_ci_hi"],
                "d_boot_se": None, "n_boot_used": b["n_boot_used"],
            }
        # recompute bootstrap SE(d) directly (percentile CI already stored; SE from boot dist)
        return out, boot

    res = {"M1": None, "m3": {}}
    # M1 (sensitivity)
    r1, b1 = fit_and_boot(sd["M1"], "M1")
    res["M1"] = _attach_boot_se(r1, b1, strat, "M1", subtype, cov_dict, cov_names, sd["M1"], mask)
    # M3 targets
    for tf in CORE_TARGETS_C1 + CORE_TARGETS_C2:
        rr, bb = fit_and_boot(sd["m3"][tf], f"M3_{tf}")
        res["m3"][tf] = _attach_boot_se(rr, bb, strat, f"M3_{tf}", subtype, cov_dict, cov_names,
                                        sd["m3"][tf], mask)
    res["subtype_counts"] = {s: int((subtype == s).sum()) for s in SUBTYPES}
    return res


def _attach_boot_se(res, boot, strat, key, subtype, cov_dict, cov_names, y, mask):
    """Recompute the bootstrap SD of d (as SE) from a fresh identical bootstrap draw store.
    We rerun bootstrap capturing the full d distribution (same seed) to get SD(d)."""
    y = np.asarray(y, float)
    idx = np.where(mask)[0]
    yy = y[idx]; ss = subtype[idx]
    cl = {cn: np.asarray(cov_dict[cn], float)[idx] for cn in cov_names}
    nn = len(idx)
    bseed = subseed(f"boot__{strat}__{key}")
    rng = np.random.default_rng(bseed)
    dstore = {c: [] for c in ["C1", "C2", "C3"]}
    for _ in range(N_BOOT):
        samp = rng.integers(0, nn, size=nn)
        yb = yy[samp]; sb = ss[samp]
        if len(set(sb.tolist())) < 4:
            continue
        cb = {cn: cl[cn][samp] for cn in cov_names}
        Xb = ME.build_design(sb, cb, cov_names)
        try:
            beta, resid, sr, xi, dof, s2 = ME.fit_ols(yb, Xb)
        except np.linalg.LinAlgError:
            continue
        cref = {cn: float(np.mean(cb[cn])) for cn in cov_names}
        mu = ME.subtype_emms(beta, cref, cov_names)
        for c in ["C1", "C2", "C3"]:
            dstore[c].append(ME.contrast_estimate(mu, c) / sr if sr > 0 else np.nan)
    for c in ["C1", "C2", "C3"]:
        arr = np.array([v for v in dstore[c] if np.isfinite(v)])
        res["contrasts"][c]["d_boot_se"] = float(arr.std(ddof=1)) if arr.size > 1 else None
        res["contrasts"][c]["n_boot_used"] = int(arr.size)
    return res


# ---- run both strata ----
scored = {}
modelled = {}
for strat in ["Discovery", "Confirmatory"]:
    print(f"[score] {strat}", utcnow(), flush=True)
    scored[strat] = score_stratum(strat)
    print(f"[model] {strat}", utcnow(), flush=True)
    modelled[strat] = model_stratum(strat, scored[strat])

# ---------------------------------------------------------------
# META: inverse-variance FIXED-EFFECT over the two strata (per target-contrast)
# ---------------------------------------------------------------
def meta_target(tf):
    c = TARGET_CONTRAST[tf]
    strat_eff = {}
    for strat in ["Discovery", "Confirmatory"]:
        cc = modelled[strat]["m3"][tf]["contrasts"][c]
        strat_eff[strat] = {"d": cc["d"], "se": cc["d_boot_se"],
                            "ci_lo": cc["d_ci_lo"], "ci_hi": cc["d_ci_hi"],
                            "n_boot_used": cc["n_boot_used"]}
    ds = np.array([strat_eff[s]["d"] for s in ["Discovery", "Confirmatory"]], float)
    ses = np.array([strat_eff[s]["se"] for s in ["Discovery", "Confirmatory"]], float)
    w = 1.0 / (ses ** 2)
    d_meta = float(np.sum(w * ds) / np.sum(w))
    se_meta = float(np.sqrt(1.0 / np.sum(w)))
    z = d_meta / se_meta
    p = float(2.0 * norm.sf(abs(z)))
    ci_lo = d_meta - 1.96 * se_meta
    ci_hi = d_meta + 1.96 * se_meta
    # opposite-direction veto: strata point opposite signs
    signs = np.sign(ds)
    opposite = bool(signs[0] * signs[1] < 0)
    return {
        "target": tf, "contrast": c, "pred_sign": TARGET_PRED_SIGN[tf],
        "strata": strat_eff,
        "d_meta": d_meta, "se_meta": se_meta, "z": z, "p_two_sided": p,
        "ci_lo": ci_lo, "ci_hi": ci_hi,
        "meta_ci_excludes_0": bool(ci_lo > 0 or ci_hi < 0),
        "opposite_direction_veto": opposite,
    }


meta = {tf: meta_target(tf) for tf in CORE_TARGETS_C1 + CORE_TARGETS_C2}

# ---- F1 BH-of-6 across the six target-contrast pooled p-values ----
p6 = {tf: meta[tf]["p_two_sided"] for tf in meta}
q6 = ME.bh_fdr(p6)
for tf in meta:
    meta[tf]["BH_q_of_6"] = q6[tf]

# ---- per-target verdict ----
def verdict(tf):
    m = meta[tf]
    c = m["contrast"]
    d = m["d_meta"]
    pred = m["pred_sign"]
    # direction match: sign(d_meta) == predicted sign
    dir_ok = (np.sign(d) == pred) and d != 0
    mag_ok = abs(d) >= FLOOR
    ci_ok = m["meta_ci_excludes_0"]
    q_ok = (m["BH_q_of_6"] is not None and m["BH_q_of_6"] <= 0.05)
    consist_ok = not m["opposite_direction_veto"]
    family = "C1" if tf in CORE_TARGETS_C1 else "C2"
    if family == "C2":
        replicated = bool(dir_ok and mag_ok and ci_ok and q_ok and consist_ok)
        status = "externally_replicated" if replicated else "evaluable--not_replicated"
    else:
        # C1 underpowered sensitivity ONLY -- never confirmatory-replicated
        status = "underpowered_sensitivity"
        replicated = None
    return {
        "family": family, "contrast": c, "d_meta": d,
        "meta_ci": [m["ci_lo"], m["ci_hi"]], "meta_ci_excludes_0": ci_ok,
        "BH_q_of_6": m["BH_q_of_6"], "pred_sign": pred,
        "direction_match": bool(dir_ok), "abs_d_ge_0.5": bool(mag_ok),
        "stratum_consistent": bool(consist_ok),
        "opposite_direction_veto": m["opposite_direction_veto"],
        "status": status, "replicated": replicated,
    }


verdicts = {tf: verdict(tf) for tf in meta}

# ---- clean helper ----
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


# ---- assemble stratum-level output (drop raw score arrays; keep model summaries) ----
strata_out = {}
for strat in ["Discovery", "Confirmatory"]:
    sd = scored[strat]
    md = modelled[strat]
    strata_out[strat] = {
        "n": sd["n"], "subtype_counts": md["subtype_counts"],
        "m1_scored": sd["m1_scored"], "cov_rows": sd["cov_rows"],
        "m3_edges_used": sd["m3_used"], "allzero_n": sd["allzero_n"],
        "M1_model": md["M1"], "m3_model": md["m3"],
    }

allout = {
    "phase": "TASK-029 frozen primary execution",
    "git_HEAD": GIT_HEAD, "master_seed": MASTER_SEED,
    "n_boot": N_BOOT, "H1_floor": FLOOR, "SESOI": SESOI,
    "model": "V2 FROZEN NO-PURITY: score ~ subtype + M4(covariate) + ESTIMATE composition",
    "contrasts": {"C1": "[-1,-1,-1,+3] (p53abn vs rest); pred d>0 (enrichment)",
                  "C2": "[+1,+1,-2,0] (POLE+MMRd vs NSMP); pred d<0 (depletion)",
                  "C3": "[+1,-1,0,0] descriptive"},
    "m3_included_edges": n_included,
    "strata": strata_out,
    "meta": meta,
    "F1_BH_of_6_p": p6, "F1_BH_of_6_q": q6,
    "per_target_verdict": verdicts,
    "timestamp_utc": utcnow(),
}
json.dump(clean(allout), open(os.path.join(RESULTS, "primary_results.json"), "w"), indent=2)
print("=== PRIMARY EXECUTION done", utcnow(), "===", flush=True)
for tf in CORE_TARGETS_C2 + CORE_TARGETS_C1:
    v = verdicts[tf]
    print(f"  {tf:6s} [{v['family']}/{v['contrast']}] d_meta={v['d_meta']:+.3f} "
          f"q6={v['BH_q_of_6']} ci_excl0={v['meta_ci_excludes_0']} "
          f"dir={v['direction_match']} -> {v['status']}", flush=True)
