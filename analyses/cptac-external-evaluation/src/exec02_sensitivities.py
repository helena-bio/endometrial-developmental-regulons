#!/usr/bin/env python3
"""
TASK-029 FROZEN PRE-FROZEN SENSITIVITIES (labelled; CANNOT override the primary verdict).

(a) ESTIMATE-purity sensitivity: add an ESTIMATE-derived tumour-purity covariate to the V2
    model (score ~ subtype + M4 + ESTIMATE composition + ESTIMATE_purity). ESTIMATE-purity =
    Yoshihara 2013 monotone transform of the ESTIMATE score:
       purity = cos( 0.6049872018 + 0.0001467884 * ESTIMATEScore ).
    Expression-derived, non-identical to Aran-CPE, possibly collinear with composition; per the
    frozen sub-rules this is a DISCLOSED sensitivity ONLY. Both strata (ESTIMATE computable).
    Report collinearity (VIF), model-instability, and re-run C1/C2 meta + BH-of-6.

(b) MMRd/classifier robustness sensitivity: run ONLY IF an alternative TCGA-style MMRd/p53abn
    call is available OUTCOME-BLIND with pre-documented provenance. Availability is assessed
    here (do NOT improvise if unavailable). Provenance: Dou2023 mmc2 Meta_table raw genomic
    fields (MSI_status, MSIsensor_ratio, POLE, CNV_status, TP53) + IHC (MLH1/MSH2/MSH6/PMS2/p53).
    These are RNA-independent, per-case, pre-existing. LIMITATION: they exist only for the
    Confirmatory stratum in the sealed open sources; Discovery has only the cBioPortal
    GENOMICS_SUBTYPE (the primary label). So the alt call is Confirmatory-only -> a single-stratum
    classifier-robustness sensitivity, disclosed as such.

(c) center/TMT-batch PROXY sensitivity: pre-authorized scope only (disclosed proxy, NOT a true
    site correction). Confirmatory Batch/Plex from mmc2; assessed for availability + scope.

Reuses sealed scoring/model code VERBATIM. MASTER SEED 20260714. Git HEAD unchanged.
"""
import os
import sys
import json
import hashlib
import numpy as np
from datetime import datetime, timezone
from scipy.stats import norm

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
import scoring as SC
import model_engine as ME


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
MMC2 = os.environ.get(
    "CPTAC_MMC2_XLSX",
    os.path.join(BASE, "external-inputs", "dou2023_mmc2.xlsx"),
)

MASTER_SEED = 20260714
N_BOOT = 2000
FLOOR = 0.50
SUBTYPES = ["POLE", "MMRd", "NSMP", "p53abn"]
C1_T = ["LHX1", "PAX8"]
C2_T = ["GATA2", "HOXA9", "SOX9", "WT1"]
TARGET_CONTRAST = {t: "C1" for t in C1_T}; TARGET_CONTRAST.update({t: "C2" for t in C2_T})
PRED = {t: +1 for t in C1_T}; PRED.update({t: -1 for t in C2_T})


def utcnow():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def subseed(step):
    return int(hashlib.sha256(f"{MASTER_SEED}:{step}".encode()).hexdigest()[:8], 16)


def read_module_genes(path):
    return [l.strip() for l in open(path) if l.strip() and not l.strip().startswith("#")]


gene_names = json.load(open(os.path.join(INTER, "gene_names.json")))
case_order = json.load(open(os.path.join(INTER, "case_order_by_stratum.json")))
lg = json.load(open(os.path.join(INTER, "acq02_join_ledger.json")))
pick = lg["pick"]

estimate_sig = {}
for line in open(ESTIMATE_GMT):
    p = line.rstrip("\n").split("\t")
    estimate_sig[p[0]] = p[2:]

# edge ledger primary regulons
from collections import defaultdict
edges = [l.rstrip("\n").split("\t") for l in open(os.path.join(B1, "m3_primary_edge_ledger.tsv"))
         if not l.startswith("#")]
hdr = edges[0]; edges = edges[1:]; ixe = {c: i for i, c in enumerate(hdr)}
prim = defaultdict(lambda: {"targets": [], "weights": [], "mor": []})
for r in edges:
    if r[ixe["primary_inclusion"]] == "Y":
        tf = r[ixe["TF"]]
        prim[tf]["targets"].append(r[ixe["target"]])
        prim[tf]["weights"].append(float(r[ixe["weight"]]))
        prim[tf]["mor"].append(int(r[ixe["primary_sign"]]))

M1_genes = read_module_genes(os.path.join(B1, "M1_analysis_ready.txt"))
M4_genes = read_module_genes(os.path.join(B1, "M4_covariate.txt"))


def build_scores(strat):
    LOG = np.load(os.path.join(INTER, f"log2tpm_{strat}.npy"))
    cases = case_order[strat]
    sym_to_rows = {}
    for i, nm in enumerate(gene_names):
        sym_to_rows.setdefault(nm, []).append(i)

    def best_row(sym):
        rows = sym_to_rows.get(sym)
        if not rows:
            return None
        return rows[0] if len(rows) == 1 else rows[int(np.argmax([LOG[r].mean() for r in rows]))]

    allzero = set(sym for sym, rows in sym_to_rows.items()
                  if all((LOG[r] == 0).all() for r in rows))

    def ssg(genes, name):
        rows, used = [], set()
        for g in genes:
            if g in allzero:
                continue
            r = best_row(g)
            if r is not None and r not in used:
                rows.append(r); used.add(r)
        return SC.ssgsea_matrix(LOG, {name: rows})[name]

    m1 = ssg(M1_genes, "M1")
    m4 = ssg(M4_genes, "M4")
    strom = ssg(estimate_sig["StromalSignature"], "S")
    imm = ssg(estimate_sig["ImmuneSignature"], "I")
    comp = strom + imm
    estimate_score = strom + imm  # ESTIMATEScore = stromal + immune (Yoshihara)
    # ESTIMATE-purity (Yoshihara 2013 monotone transform)
    est_purity = np.cos(0.6049872018 + 0.0001467884 * estimate_score)
    gene_index = {}
    for sym in set(t for tf in prim for t in prim[tf]["targets"]):
        r = best_row(sym)
        if r is not None:
            gene_index[sym] = r
    m3 = {}
    for tf in sorted(prim):
        d = prim[tf]
        nes, _ = SC.area_regulon(LOG, gene_index, d["targets"], d["weights"], d["mor"])
        m3[tf] = nes
    subtype = np.array([pick[c]["subtype"] for c in cases])
    return {"cases": cases, "subtype": subtype, "M1": m1, "M4": m4,
            "composition": comp, "estimate_score": estimate_score,
            "estimate_purity": est_purity, "m3": m3, "n": len(cases)}


def d_and_bootse(strat, key, y, subtype, cov_dict, cov_names):
    y = np.asarray(y, float)
    n = len(y); mask = np.ones(n, dtype=bool)
    fit = ME.fit_module(y, subtype, cov_dict, cov_names, mask)
    bseed = subseed(f"sens_boot__{strat}__{key}")
    rng = np.random.default_rng(bseed)
    idx = np.arange(n)
    cl = {cn: np.asarray(cov_dict[cn], float) for cn in cov_names}
    store = {c: [] for c in ["C1", "C2", "C3"]}
    for _ in range(N_BOOT):
        s = rng.integers(0, n, size=n)
        if len(set(subtype[s].tolist())) < 4:
            continue
        cb = {cn: cl[cn][s] for cn in cov_names}
        Xb = ME.build_design(subtype[s], cb, cov_names)
        try:
            beta, resid, sr, xi, dof, s2 = ME.fit_ols(y[s], Xb)
        except np.linalg.LinAlgError:
            continue
        cref = {cn: float(np.mean(cb[cn])) for cn in cov_names}
        mu = ME.subtype_emms(beta, cref, cov_names)
        for c in ["C1", "C2", "C3"]:
            store[c].append(ME.contrast_estimate(mu, c) / sr if sr > 0 else np.nan)
    out = {}
    for c in ["C1", "C2", "C3"]:
        arr = np.array([v for v in store[c] if np.isfinite(v)])
        out[c] = {"d": fit["contrasts"][c]["d"],
                  "d_boot_se": float(arr.std(ddof=1)) if arr.size > 1 else None,
                  "d_ci_lo": float(np.percentile(arr, 2.5)) if arr.size else None,
                  "d_ci_hi": float(np.percentile(arr, 97.5)) if arr.size else None,
                  "n_boot_used": int(arr.size)}
    return out


def meta_and_bh(per_strat_d):
    """per_strat_d: {tf: {'Discovery':{d,se}, 'Confirmatory':{d,se}}} -> meta + BH-of-6."""
    meta = {}
    for tf, sd in per_strat_d.items():
        ds = np.array([sd["Discovery"]["d"], sd["Confirmatory"]["d"]], float)
        ses = np.array([sd["Discovery"]["se"], sd["Confirmatory"]["se"]], float)
        w = 1.0 / ses ** 2
        dm = float(np.sum(w * ds) / np.sum(w)); sem = float(np.sqrt(1.0 / np.sum(w)))
        z = dm / sem; p = float(2.0 * norm.sf(abs(z)))
        meta[tf] = {"d_meta": dm, "se_meta": sem, "z": z, "p": p,
                    "ci_lo": dm - 1.96 * sem, "ci_hi": dm + 1.96 * sem,
                    "ci_excl0": bool((dm - 1.96 * sem) > 0 or (dm + 1.96 * sem) < 0),
                    "opposite_veto": bool(ds[0] * ds[1] < 0),
                    "strata": sd}
    q = ME.bh_fdr({tf: meta[tf]["p"] for tf in meta})
    for tf in meta:
        meta[tf]["BH_q_of_6"] = q[tf]
    return meta


# ================= (a) ESTIMATE-purity sensitivity =================
from statsmodels.stats.outliers_influence import variance_inflation_factor

sens_a = {"label": "ESTIMATE-purity sensitivity (V2 + ESTIMATE-derived purity covariate)",
          "cannot_override_primary": True, "per_strat": {}, "collinearity": {}}
scored = {s: build_scores(s) for s in ["Discovery", "Confirmatory"]}
per_strat_d = {tf: {} for tf in C1_T + C2_T}
for strat in ["Discovery", "Confirmatory"]:
    sc = scored[strat]
    cov_dict = {"M4": sc["M4"], "composition": sc["composition"],
                "estimate_purity": sc["estimate_purity"]}
    cov_names = ["M4", "composition", "estimate_purity"]
    # collinearity: VIF among covariates (+ subtype dummies), and purity~composition correlation
    subt = sc["subtype"]
    dummies = np.column_stack([(subt == s).astype(float) for s in ["MMRd", "NSMP", "p53abn"]])
    X = np.column_stack([np.ones(sc["n"]), dummies, sc["M4"], sc["composition"], sc["estimate_purity"]])
    names = ["intercept", "MMRd", "NSMP", "p53abn", "M4", "composition", "estimate_purity"]
    vifs = {names[j]: float(variance_inflation_factor(X, j)) for j in range(1, X.shape[1])}
    corr = float(np.corrcoef(sc["composition"], sc["estimate_purity"])[0, 1])
    sens_a["collinearity"][strat] = {"VIF": vifs, "corr_composition_purity": corr,
                                     "max_VIF": max(vifs.values())}
    for tf in C1_T + C2_T:
        r = d_and_bootse(strat, f"estpurity_M3_{tf}", sc["m3"][tf], subt, cov_dict, cov_names)
        c = TARGET_CONTRAST[tf]
        per_strat_d[tf][strat] = {"d": r[c]["d"], "se": r[c]["d_boot_se"],
                                  "ci_lo": r[c]["d_ci_lo"], "ci_hi": r[c]["d_ci_hi"]}
sens_a["meta"] = meta_and_bh(per_strat_d)
# instability flag: any covariate VIF>10
sens_a["model_instability_VIF_gt10"] = any(
    sens_a["collinearity"][s]["max_VIF"] > 10 for s in ["Discovery", "Confirmatory"])

# ================= (b) MMRd/classifier robustness sensitivity: AVAILABILITY =================
import openpyxl
wb = openpyxl.load_workbook(MMC2, read_only=True, data_only=True)
ws = wb["Meta_table"]
rows = ws.iter_rows(values_only=True)
mh = next(rows)
midx = {h: i for i, h in enumerate(mh)}
alt_fields = ["MSI_status", "MSIsensor_ratio", "POLE", "CNV_status", "TP53",
              "Ancillary_studies_mlh1", "Ancillary_studies_msh2",
              "Ancillary_studies_msh6", "Ancillary_studies_pms2", "Ancillary_studies_p53",
              "ABSOLUTE_tumor_purity", "Batch", "Plex"]
present = {f: (f in midx) for f in alt_fields}
# tally non-empty ABSOLUTE purity + MSI_status for confirmatory joined cases
conf_cases = set(case_order["Confirmatory"])
absol = 0; msi_present = 0; batch_vals = {}
for row in rows:
    cid = row[midx["Case_id"]]
    if cid not in conf_cases:
        continue
    av = row[midx["ABSOLUTE_tumor_purity"]]
    if av not in (None, "", "NA"):
        absol += 1
    mv = row[midx["MSI_status"]]
    if mv not in (None, "", "NA"):
        msi_present += 1
    bv = row[midx["Batch"]]
    batch_vals[bv] = batch_vals.get(bv, 0) + 1
wb.close()

sens_b = {
    "label": "MMRd/p53abn classifier-robustness sensitivity (alternative TCGA-style call)",
    "cannot_override_primary": True,
    "fields_present_in_confirmatory_mmc2": present,
    "confirmatory_ABSOLUTE_purity_nonempty": absol,
    "confirmatory_MSI_status_nonempty": msi_present,
    "availability_verdict": (
        "AVAILABLE but CONFIRMATORY-STRATUM ONLY. RNA-independent alternative-classifier "
        "components (MSI_status/MSIsensor_ratio/POLE/CNV_status/TP53 + MLH1/MSH2/MSH6/PMS2/p53 "
        "IHC) exist per-case in the Dou2023 mmc2 Meta_table for the Confirmatory stratum. The "
        "Discovery stratum has only the cBioPortal GENOMICS_SUBTYPE (the primary label) in the "
        "sealed open sources -- no independent raw MSI/CNV/POLE components. Per the frozen rule "
        "the alternative call is therefore a SINGLE-STRATUM (Confirmatory-only) classifier-"
        "robustness sensitivity; a two-stratum meta alt-classifier replication is NOT available "
        "and is NOT improvised. This is a DISCLOSED LIMITATION. The primary Dou-2023 genomic "
        "subtype labels remain PRIMARY and define the confirmatory C2 analysis."),
    "not_run_as_full_meta_reason": "no RNA-independent alt classifier for Discovery in sealed sources",
}

# ================= (c) center / TMT-batch PROXY sensitivity: SCOPE =================
sens_c = {
    "label": "center / TMT-batch PROXY sensitivity (disclosed proxy, not a site correction)",
    "cannot_override_primary": True,
    "confirmatory_batch_distribution": {str(k): v for k, v in batch_vals.items()},
    "scope_verdict": (
        "PROXY ONLY, pre-authorized scope. Confirmatory Batch(b1-b4)/Plex are TMT PROTEOMIC "
        "processing batches, NOT the RNA-seq batch and NOT a true recruiting tissue_source_site "
        "(unpopulated in GDC/PDC). Accrual-family C3L/C3N appears in BOTH strata (not fully "
        "confounded with the Discovery-vs-Confirmatory split). A DISCLOSED batch/center proxy "
        "sensitivity is admissible in scope, but it is NOT a true-site correction and CANNOT "
        "change the primary verdict. A true tissue_source_site correction remains a LIMITATION "
        "(not improvised)."),
}

out = {"phase": "TASK-029 pre-frozen sensitivities", "master_seed": MASTER_SEED,
       "git_HEAD": "83503bad47b60193598b2b9ebe819c22c83e8ac1",
       "sensitivity_a_estimate_purity": sens_a,
       "sensitivity_b_mmrd_classifier": sens_b,
       "sensitivity_c_center_batch_proxy": sens_c,
       "timestamp_utc": utcnow()}


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


json.dump(clean(out), open(os.path.join(RESULTS, "sensitivity_results.json"), "w"), indent=2)
print("=== SENSITIVITIES done", utcnow(), "===", flush=True)
print("(a) ESTIMATE-purity meta (cannot override primary):", flush=True)
for tf in C2_T + C1_T:
    m = sens_a["meta"][tf]
    print(f"    {tf:6s} d_meta={m['d_meta']:+.3f} ci_excl0={m['ci_excl0']} q6={m['BH_q_of_6']:.4f} "
          f"veto={m['opposite_veto']}", flush=True)
print("    max VIF:", {s: round(sens_a['collinearity'][s]['max_VIF'], 2) for s in ['Discovery','Confirmatory']},
      "instability_VIF>10:", sens_a["model_instability_VIF_gt10"], flush=True)
print("(b) alt-classifier availability:", sens_b["availability_verdict"][:60], "...", flush=True)
print("(c) batch proxy scope: confirmatory-only proxy, disclosed", flush=True)
