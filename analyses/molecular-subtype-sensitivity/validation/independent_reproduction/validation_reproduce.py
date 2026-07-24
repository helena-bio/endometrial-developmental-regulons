#!/usr/bin/env python3
"""Independent TASK-A reconstruction.

This critic implementation does not import producer code. It reconstructs CPTAC
scores from frozen expression/regulon bytes, builds every design matrix directly,
and independently implements OLS, bootstrap, permutation, BH-18, taxonomy, and
fixed-effect synthesis.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import platform
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("PYTHONHASHSEED", "0")

import numpy as np
import pandas as pd
import scipy
from scipy import stats

ROOT = Path("data/external/original-workspace/revgate-task-a-perclass")
TASK = ROOT / "experiments/taskA_perclass_c2"
OUT = TASK / "critic/independent_reproduction"
T28 = Path("data/external/original-workspace/task028-freeze-b-draft")
T29 = Path("data/external/original-workspace/task029-external-replication-feasibility")
T30 = Path("data/external/original-workspace/revgate-tcga-no-purity-verify")
ORIGINAL = Path("data/external/original-workspace/revgate")
TARGETS = ["GATA2", "SOX9", "HOXA9", "WT1", "PAX8", "LHX1"]
MODELS = ["TCGA_PRIMARY_CPE_N506", "TCGA_NOPURITY_N507",
          "CPTAC_DISCOVERY_N95", "CPTAC_CONFIRMATORY_N135"]
META = "CPTAC_FIXED_EFFECT_META"
SUBTYPES = ["POLE", "MMRd", "NSMP", "p53abn"]
CONTRASTS = {
    "POLE_vs_NSMP": np.array([1., 0., -1., 0.]),
    "MMRd_vs_NSMP": np.array([0., 1., -1., 0.]),
    "POLE_vs_MMRd": np.array([1., -1., 0., 0.]),
    "FROZEN_C2": np.array([.5, .5, -1., 0.]),
}
MASTER = 20260722
N_BOOT = N_PERM = 2000


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def seed(step: str) -> int:
    return int(hashlib.sha256(f"{MASTER}:{step}".encode()).hexdigest()[:8], 16)


def write_tsv(name: str, rows) -> None:
    df = rows if isinstance(rows, pd.DataFrame) else pd.DataFrame(rows)
    df.to_csv(OUT / name, sep="\t", index=False, na_rep="NA", lineterminator="\n",
              quoting=csv.QUOTE_MINIMAL)


def write_json(name: str, obj) -> None:
    def clean(x):
        if isinstance(x, dict): return {str(k): clean(v) for k, v in x.items()}
        if isinstance(x, (list, tuple)): return [clean(v) for v in x]
        if isinstance(x, np.integer): return int(x)
        if isinstance(x, np.floating): x = float(x)
        if isinstance(x, np.bool_): return bool(x)
        if isinstance(x, float) and not math.isfinite(x): return None
        return x
    (OUT / name).write_text(json.dumps(clean(obj), indent=2, sort_keys=True,
                                      allow_nan=False) + "\n", encoding="ascii")


def command(args, cwd: Path) -> bytes:
    return subprocess.run(args, cwd=cwd, check=True, capture_output=True).stdout


def git_snapshot(repo: Path) -> dict:
    status = command(["git", "status", "--porcelain=v1"], repo)
    diff = command(["git", "diff", "--binary", "HEAD"], repo)
    staged = command(["git", "diff", "--cached", "--binary", "HEAD"], repo)
    return {
        "path": str(repo),
        "head": command(["git", "rev-parse", "HEAD"], repo).decode().strip(),
        "branch": command(["git", "branch", "--show-current"], repo).decode().strip(),
        "status_sha256": hashlib.sha256(status).hexdigest(),
        "unstaged_binary_diff_sha256": hashlib.sha256(diff).hexdigest(),
        "staged_binary_diff_sha256": hashlib.sha256(staged).hexdigest(),
    }


def ssgsea(expr: np.ndarray, sets: dict[str, list[int]]) -> dict[str, np.ndarray]:
    """Independent literal implementation of frozen ssGSEA alpha=0.25."""
    ng, ns = expr.shape
    result = {k: np.zeros(ns) for k in sets}
    masks = {}
    for key, idx in sets.items():
        mask = np.zeros(ng, dtype=bool); mask[idx] = True; masks[key] = mask
    for j in range(ns):
        order = np.argsort(expr[:, j], kind="mergesort")
        rank = np.empty(ng); rank[order] = np.arange(1., ng + 1.)
        desc = order[::-1]
        rank_weight = rank[desc] ** .25
        for key, mask in masks.items():
            hit = mask[desc]
            up = np.where(hit, rank_weight / rank_weight[hit].sum(), 0.)
            down = np.where(hit, 0., 1. / (ng - hit.sum()))
            result[key][j] = np.cumsum(up - down).sum()
    for key, val in result.items():
        result[key] = val / (val.max() - val.min())
    return result


def cptac_data(stratum: str, regulons: dict) -> dict:
    inter = T29 / "execution/intermediate"
    expr = np.load(inter / f"log2tpm_{stratum}.npy")
    genes = json.loads((inter / "gene_names.json").read_text())
    cases = json.loads((inter / "case_order_by_stratum.json").read_text())[stratum]
    pick = json.loads((inter / "acq02_join_ledger.json").read_text())["pick"]
    symbol_rows = defaultdict(list)
    for i, g in enumerate(genes): symbol_rows[g].append(i)
    means = expr.mean(1)
    best = {g: rows[int(np.argmax(means[rows]))] for g, rows in symbol_rows.items()}
    allzero = {g for g, rows in symbol_rows.items() if all(np.all(expr[r] == 0) for r in rows)}
    modules = {"M4": [x for x in (T28 / "sealed_v3/B1/M4_covariate.txt").read_text().splitlines()
                       if x and not x.startswith("#")]}
    for line in (T28 / "experimenter_final/sources/ESTIMATE_SI_geneset.gmt").read_text().splitlines():
        f = line.split("\t")
        if f[0] == "StromalSignature": modules["stromal"] = f[2:]
        if f[0] == "ImmuneSignature": modules["immune"] = f[2:]
    sets = {}
    for key, symbols in modules.items():
        rows, seen = [], set()
        for g in symbols:
            r = best.get(g)
            if g not in allzero and r is not None and r not in seen:
                rows.append(r); seen.add(r)
        sets[key] = rows
    cov = ssgsea(expr, sets)

    # Independent aREA reconstruction over the frozen 60,660-gene universe.
    mu = expr.mean(1, keepdims=True)
    sd = expr.std(1, ddof=1, keepdims=True)
    z = (expr - mu) / sd
    finite = np.all(np.isfinite(z), axis=1)
    pos = -np.ones(expr.shape[0], dtype=int)
    pos[finite] = np.arange(finite.sum())
    z = z[finite]
    prepared = {}
    for tf, reg in regulons.items():
        rows, weights, signs = [], [], []
        for g, w, s in zip(reg["target"], reg["weight"], reg["sign"]):
            r = best.get(g)
            if r is not None and finite[r]:
                rows.append(pos[r]); weights.append(float(w)); signs.append(float(s))
        raww = np.asarray(weights) / max(weights)
        prepared[tf] = (np.asarray(rows), raww / raww.sum(), np.asarray(signs),
                        math.sqrt(float(raww @ raww)))
    scores = {t: np.zeros(expr.shape[1]) for t in TARGETS}
    n = z.shape[0]
    for j in range(expr.shape[1]):
        v = z[:, j]
        two = stats.norm.ppf(stats.rankdata(v, method="average") / (n + 1.))
        one = stats.norm.ppf((stats.rankdata(np.abs(v), method="average") / (n + 1.) + 1.) / 2.)
        for tf, (rows, wt, sign, normer) in prepared.items():
            directional = float(np.sum(wt * sign * two[rows]))
            magnitude = float(np.sum(wt * one[rows]))
            scores[tf][j] = (abs(directional) + max(magnitude, 0.)) * np.sign(directional) * normer
    return {
        "subtype": np.asarray([pick[c]["subtype"] for c in cases]),
        "cov": {"M4": cov["M4"], "composition": cov["stromal"] + cov["immune"]},
        "cov_names": ["M4", "composition"], "scores": scores,
        "cases": cases, "targets_used": {t: len(prepared[t][0]) for t in TARGETS},
    }


def load_data() -> dict:
    inter = T28 / "execution/intermediate"
    npz = np.load(inter / "scores_v3.npz", allow_pickle=True)
    covdf = pd.read_csv(inter / "covariates_v3.tsv", sep="\t")
    assert list(npz["patient_order"]) == list(covdf.patient_barcode)
    subtype = covdf.subtype.to_numpy()
    score = {t: np.asarray(npz[f"M3primary__{t}"], float) for t in TARGETS}
    cov = {"M4": covdf.M4_prolif.to_numpy(float),
           "Aran_CPE": covdf.purity_CPE.to_numpy(float),
           "composition": covdf.composition.to_numpy(float)}
    keep = covdf.cpe_complete_case.to_numpy() == 1
    data = {
        MODELS[0]: {"subtype": subtype[keep], "cov": {k: v[keep] for k, v in cov.items()},
                    "cov_names": ["M4", "Aran_CPE", "composition"],
                    "scores": {k: v[keep] for k, v in score.items()}},
        MODELS[1]: {"subtype": subtype, "cov": {"M4": cov["M4"], "composition": cov["composition"]},
                    "cov_names": ["M4", "composition"], "scores": score},
    }
    edge = pd.read_csv(T28 / "sealed_v3/B1/m3_primary_edge_ledger.tsv", sep="\t", comment="#")
    edge = edge[edge.primary_inclusion == "Y"]
    regulons = {}
    for t in TARGETS:
        e = edge[edge.TF == t]
        regulons[t] = {"target": e.target.tolist(), "weight": e.weight.to_numpy(float),
                       "sign": e.primary_sign.to_numpy(float)}
    data[MODELS[2]] = cptac_data("Discovery", regulons)
    data[MODELS[3]] = cptac_data("Confirmatory", regulons)
    return data


def matrix(subtype, cov, names):
    return np.column_stack([np.ones(len(subtype)),
                            subtype == "MMRd", subtype == "NSMP", subtype == "p53abn",
                            *[cov[x] for x in names]]).astype(float)


def group_rows(cov, names):
    means = [float(np.mean(cov[x])) for x in names]
    return np.asarray([[1., s == "MMRd", s == "NSMP", s == "p53abn", *means]
                       for s in SUBTYPES], float)


def solve(y, subtype, cov, names, uncertainty=False):
    X = matrix(subtype, cov, names); n, p = X.shape
    beta, _, rank, _ = np.linalg.lstsq(X, y, rcond=None)
    if rank != p: raise np.linalg.LinAlgError("rank deficient")
    resid = y - X @ beta; sigma = math.sqrt(float(resid @ resid) / (n - p))
    E = group_rows(cov, names)
    vals = {key: {"b": float(w @ E @ beta), "d": float(w @ E @ beta / sigma)}
            for key, w in CONTRASTS.items()}
    if uncertainty:
        vb = sigma * sigma * np.linalg.inv(X.T @ X); tc = stats.t.ppf(.975, n - p)
        for key, w in CONTRASTS.items():
            l = w @ E; se = math.sqrt(float(l @ vb @ l)); b = vals[key]["b"]
            vals[key].update(se=se, t_lo=b-tc*se, t_hi=b+tc*se,
                             wald_p=float(2*stats.t.sf(abs(b/se), n-p)))
    return sigma, n-p, vals


def analyze(y, subtype, cov, names, model, target):
    sigma, df, obs = solve(y, subtype, cov, names, True)
    rng = np.random.default_rng(seed(f"boot__{model}__{target}"))
    store = {c: {"b": [], "d": []} for c in CONTRASTS}; fail = 0
    max_b = max_d = 0.
    for _ in range(N_BOOT):
        ix = rng.integers(0, len(y), len(y)); st = subtype[ix]
        if len(set(st)) != 4: fail += 1; continue
        try:
            _, _, val = solve(y[ix], st, {k: v[ix] for k, v in cov.items()}, names)
        except (np.linalg.LinAlgError, ValueError):
            fail += 1; continue
        for c in CONTRASTS:
            store[c]["b"].append(val[c]["b"]); store[c]["d"].append(val[c]["d"])
        b = [val[c]["b"] for c in ["POLE_vs_NSMP", "MMRd_vs_NSMP", "FROZEN_C2"]]
        d = [val[c]["d"] for c in ["POLE_vs_NSMP", "MMRd_vs_NSMP", "FROZEN_C2"]]
        max_b = max(max_b, abs(b[2]-.5*b[0]-.5*b[1])/max(1.,abs(b[2]),.5*abs(b[0])+.5*abs(b[1])))
        max_d = max(max_d, abs(d[2]-.5*d[0]-.5*d[1])/max(1.,abs(d[2]),.5*abs(d[0])+.5*abs(d[1])))
    for c in CONTRASTS:
        ba, da = np.asarray(store[c]["b"]), np.asarray(store[c]["d"])
        obs[c].update(b_boot_lo=float(np.percentile(ba,2.5)), b_boot_hi=float(np.percentile(ba,97.5)),
                      d_boot_lo=float(np.percentile(da,2.5)), d_boot_hi=float(np.percentile(da,97.5)),
                      d_boot_se=float(np.std(da,ddof=1)))
    rng = np.random.default_rng(seed(f"perm__{model}__{target}"))
    ge = {c: 1 for c in list(CONTRASTS)[:3]}
    for _ in range(N_PERM):
        _, _, val = solve(y, rng.permutation(subtype), cov, names)
        for c in ge:
            ge[c] += abs(val[c]["d"]) >= abs(obs[c]["d"]) - 1e-12
    for c in ge: obs[c]["perm_p"] = ge[c] / (N_PERM + 1.)
    return {"sigma": sigma, "df": df, "values": obs, "boot_used": N_BOOT-fail,
            "boot_failed": fail, "max_boot_b_scaled": max_b, "max_boot_d_scaled": max_d}


def taxonomy(v, meta=False):
    pn, mn, pm = [v[x] for x in list(CONTRASTS)[:3]]; prefix = "CPTAC_META_" if meta else ""
    a, b = pn["value"], mn["value"]
    if a == 0 or b == 0: return prefix+"ZERO_BOUNDARY_UNRESOLVED"
    ex = lambda x: x["lo"] > 0 or x["hi"] < 0
    if a*b < 0:
        return prefix+("OPPOSITE_DIRECTION_HETEROGENEITY_SUPPORTED" if ex(pm) and ex(pn) and ex(mn)
                       else "OPPOSITE_DIRECTION_POINT_HETEROGENEITY_UNRESOLVED")
    if ex(pm) and abs(pm["d_pm"]) >= .5:
        if abs(a) > abs(b): return prefix+"SAME_DIRECTION_MATERIALLY_DIFFERENT_POLE_DOMINANT"
        if abs(b) > abs(a): return prefix+"SAME_DIRECTION_MATERIALLY_DIFFERENT_MMRD_DOMINANT"
        return prefix+"SAME_DIRECTION_MATERIALLY_DIFFERENT_TIE_ERROR"
    if ex(pm): return prefix+"SAME_DIRECTION_DISTINGUISHABLE_BELOW_MATERIALITY_FLOOR"
    if ex(pn) and ex(mn): return prefix+"SAME_DIRECTION_COMPATIBLE_BOTH_RESOLVED"
    return prefix+"SAME_DIRECTION_POINT_ESTIMATES_UNRESOLVED"


def bh18(p):
    a=np.asarray(p,float); order=np.argsort(a,kind="mergesort"); ranked=a[order]
    qrank=np.minimum.accumulate((ranked*18/np.arange(1,19))[::-1])[::-1]
    q=np.empty(18); q[order]=np.minimum(qrank,1.); return q


def compare_frames(critic: pd.DataFrame, producer: pd.DataFrame, keys, fields, source):
    p = producer.set_index(keys); rows=[]
    for _, cr in critic.iterrows():
        key=tuple(cr[k] for k in keys); key=key[0] if len(key)==1 else key
        pr=p.loc[key]
        for field in fields:
            cv, pv=cr[field], pr[field]
            if pd.isna(cv) and pd.isna(pv): diff=rel=0.; ok=True
            elif isinstance(cv,(str,bool,np.bool_)) or isinstance(pv,(str,bool,np.bool_)):
                diff=rel=np.nan; ok=str(cv)==str(pv)
            else:
                diff=abs(float(cv)-float(pv)); rel=diff/max(1e-300,abs(float(pv))); ok=diff <= 1e-11*max(1.,abs(float(pv)))
            rows.append({"source":source, **{k:cr[k] for k in keys}, "field":field,
                         "critic":cv,"producer":pv,"abs_diff":diff,"rel_diff":rel,"match":ok})
    return rows


def main():
    spec=json.loads((TASK/"FROZEN_POSTHOC_SPEC.json").read_text())
    start={"original":git_snapshot(ORIGINAL),"task030":git_snapshot(T30)}
    data=load_data(); counts=[]
    expected={m:spec["models"][m]["expected_counts"] for m in MODELS}
    for m in MODELS:
        got={s:int(np.sum(data[m]["subtype"]==s)) for s in SUBTYPES}; got["total"]=len(data[m]["subtype"])
        assert got==expected[m], (m,got,expected[m]); counts.append({"model":m,**got})
    fits={}
    for m in MODELS:
        for t in TARGETS:
            print("critic fit",m,t,flush=True)
            d=data[m]; fits[m,t]=analyze(d["scores"][t],d["subtype"],d["cov"],d["cov_names"],m,t)
    rows=[]; identity=[]; states={}
    for m in MODELS:
        for t in TARGETS:
            f=fits[m,t]; v=f["values"]
            state=taxonomy({c:{"value":v[c]["b"],"lo":v[c]["t_lo"],"hi":v[c]["t_hi"],
                                   "d_pm":v["POLE_vs_MMRd"]["d"]} for c in list(CONTRASTS)[:3]})
            states[m,t]=state
            bc=v["FROZEN_C2"]["b"]; dc=v["FROZEN_C2"]["d"]
            bre=.5*v["POLE_vs_NSMP"]["b"]+.5*v["MMRd_vs_NSMP"]["b"]
            dre=.5*v["POLE_vs_NSMP"]["d"]+.5*v["MMRd_vs_NSMP"]["d"]
            identity.append({"model":m,"target":t,"b_abs_error":abs(bc-bre),"d_abs_error":abs(dc-dre),
                             "bootstrap_max_b_scaled_error":f["max_boot_b_scaled"],
                             "bootstrap_max_d_scaled_error":f["max_boot_d_scaled"]})
            cnt={s:int(np.sum(data[m]["subtype"]==s)) for s in SUBTYPES}
            for c in list(CONTRASTS)[:3]:
                x=v[c]; rows.append({"model":m,"target":t,"contrast":c,"analytic_n":len(data[m]["subtype"]),
                    "residual_df":f["df"],**{f"n_{s}":cnt[s] for s in SUBTYPES},"coefficient":x["b"],
                    "coefficient_SE":x["se"],"residual_SD":f["sigma"],"d":x["d"],
                    "coefficient_t_CI_lo":x["t_lo"],"coefficient_t_CI_hi":x["t_hi"],
                    "coefficient_boot_CI_lo":x["b_boot_lo"],"coefficient_boot_CI_hi":x["b_boot_hi"],
                    "d_boot_CI_lo":x["d_boot_lo"],"d_boot_CI_hi":x["d_boot_hi"],"d_boot_SE":x["d_boot_se"],
                    "perm_p_raw":x["perm_p"],"wald_t_p_diagnostic":x["wald_p"],"n_boot_used":f["boot_used"],
                    "n_boot_failed":f["boot_failed"],"interpretation_state":state})
    df=pd.DataFrame(rows)
    # CPTAC FE synthesis, including C2 for discrepancy audit.
    meta=[]; meta_lookup={}
    for t in TARGETS:
        for c in CONTRASTS:
            a=fits[MODELS[2],t]["values"][c]; b=fits[MODELS[3],t]["values"][c]
            w=np.asarray([1/a["d_boot_se"]**2,1/b["d_boot_se"]**2]); ds=np.asarray([a["d"],b["d"]])
            est=float(w@ds/w.sum()); se=math.sqrt(1/w.sum())
            r={"target":t,"contrast":c,"d_FE":est,"SE_FE":se,"CI_lo":est-1.96*se,"CI_hi":est+1.96*se,
               "p_raw":float(2*stats.norm.sf(abs(est/se))),"discovery_coefficient":a["b"],
               "discovery_coefficient_SE":a["se"],"discovery_residual_SD":fits[MODELS[2],t]["sigma"],
               "discovery_d":a["d"],"discovery_d_boot_CI_lo":a["d_boot_lo"],"discovery_d_boot_CI_hi":a["d_boot_hi"],
               "discovery_boot_SE_d":a["d_boot_se"],"confirmatory_coefficient":b["b"],
               "confirmatory_coefficient_SE":b["se"],"confirmatory_residual_SD":fits[MODELS[3],t]["sigma"],
               "confirmatory_d":b["d"],"confirmatory_d_boot_CI_lo":b["d_boot_lo"],"confirmatory_d_boot_CI_hi":b["d_boot_hi"],
               "confirmatory_boot_SE_d":b["d_boot_se"],"weight_discovery":w[0],"weight_confirmatory":w[1],
               "weight_fraction_discovery":w[0]/w.sum(),"weight_fraction_confirmatory":w[1]/w.sum(),
               "opposite_direction_flag":bool(ds[0]*ds[1]<0)}
            meta.append(r); meta_lookup[t,c]=r
    for t in TARGETS:
        vals={c:{"value":meta_lookup[t,c]["d_FE"],"lo":meta_lookup[t,c]["CI_lo"],"hi":meta_lookup[t,c]["CI_hi"],
                 "d_pm":meta_lookup[t,"POLE_vs_MMRd"]["d_FE"]} for c in list(CONTRASTS)[:3]}
        state=taxonomy(vals,True); states[META,t]=state
        disc=meta_lookup[t,"FROZEN_C2"]["d_FE"]-.5*meta_lookup[t,"POLE_vs_NSMP"]["d_FE"]-.5*meta_lookup[t,"MMRd_vs_NSMP"]["d_FE"]
        for c in CONTRASTS: meta_lookup[t,c]["interpretation_state"]=state; meta_lookup[t,c]["meta_C2_direct_minus_reconstruction"]=disc
        for c in list(CONTRASTS)[:3]:
            x=meta_lookup[t,c]; df.loc[len(df)]={"model":META,"target":t,"contrast":c,"analytic_n":230,"residual_df":np.nan,
                "n_POLE":13,"n_MMRd":72,"n_NSMP":109,"n_p53abn":36,"coefficient":np.nan,"coefficient_SE":np.nan,
                "residual_SD":np.nan,"d":x["d_FE"],"coefficient_t_CI_lo":np.nan,"coefficient_t_CI_hi":np.nan,
                "coefficient_boot_CI_lo":np.nan,"coefficient_boot_CI_hi":np.nan,"d_boot_CI_lo":x["CI_lo"],
                "d_boot_CI_hi":x["CI_hi"],"d_boot_SE":x["SE_FE"],"perm_p_raw":x["p_raw"],
                "wald_t_p_diagnostic":np.nan,
                "n_boot_used":min(fits[MODELS[2],t]["boot_used"],fits[MODELS[3],t]["boot_used"]),
                "n_boot_failed":fits[MODELS[2],t]["boot_failed"]+fits[MODELS[3],t]["boot_failed"],
                "interpretation_state":state}
    for m in MODELS+[META]:
        ix=df.model==m; df.loc[ix,"BH18_q_descriptive"]=bh18(df.loc[ix,"perm_p_raw"].to_numpy())

    prod=pd.read_csv(TASK/"results/PER_CLASS_CONTRASTS_LONG.tsv",sep="\t",na_values=["NA"])
    numfields=["analytic_n","residual_df","n_POLE","n_MMRd","n_NSMP","n_p53abn","coefficient","coefficient_SE",
        "residual_SD","d","coefficient_t_CI_lo","coefficient_t_CI_hi","coefficient_boot_CI_lo","coefficient_boot_CI_hi",
        "d_boot_CI_lo","d_boot_CI_hi","d_boot_SE","perm_p_raw","BH18_q_descriptive","wald_t_p_diagnostic",
        "n_boot_used","n_boot_failed","interpretation_state"]
    comparisons=compare_frames(df,prod,["model","target","contrast"],numfields,"PER_CLASS_CONTRASTS_LONG.tsv")
    metadf=pd.DataFrame(meta)
    prodmeta=pd.read_csv(TASK/"results/CPTAC_FIXED_EFFECT_META.tsv",sep="\t")
    mf=[c for c in metadf.columns if c in prodmeta.columns and c not in ["target","contrast"]]
    comparisons+=compare_frames(metadf,prodmeta,["target","contrast"],mf,"CPTAC_FIXED_EFFECT_META.tsv")

    # Audit reporting-only correction by comparing every common structured field.
    old=pd.read_csv(TASK/"run1/results/PER_CLASS_CONTRASTS_LONG.tsv",sep="\t",na_values=["NA"])
    new=pd.read_csv(TASK/"run1_definitive/results/PER_CLASS_CONTRASTS_LONG.tsv",sep="\t",na_values=["NA"])
    common=[c for c in old.columns if c in new.columns]
    scientific_common=[c for c in common if c not in ["meta_discovery_coefficient","meta_discovery_coefficient_SE",
        "meta_discovery_residual_SD","meta_confirmatory_coefficient","meta_confirmatory_coefficient_SE","meta_confirmatory_residual_SD"]]
    correction_equal=old[scientific_common].equals(new[scientific_common])
    added=[c for c in new.columns if c not in old.columns]
    runpair=[]
    for rel in sorted(p.relative_to(TASK/"run1_definitive") for p in (TASK/"run1_definitive").rglob("*") if p.is_file()):
        q=TASK/"run2_definitive"/rel; runpair.append({"relative_path":str(rel),"run1_sha256":sha(TASK/"run1_definitive"/rel),
            "run2_sha256":sha(q),"byte_identical":sha(TASK/"run1_definitive"/rel)==sha(q)})

    end={"original":git_snapshot(ORIGINAL),"task030":git_snapshot(T30)}
    input_rows=pd.read_csv(TASK/"INPUT_CHECKSUMS.tsv",sep="\t")
    input_recheck=[]
    for _,r in input_rows.iterrows():
        got=sha(Path(r.path)); input_recheck.append({"path":r.path,"expected":r.expected_sha256,"observed":got,"match":got==r.expected_sha256})
    cdf=pd.DataFrame(comparisons)
    write_tsv("INDEPENDENT_RESULTS.tsv",df)
    write_tsv("INDEPENDENT_META.tsv",metadf)
    write_tsv("INDEPENDENT_IDENTITIES.tsv",identity)
    write_tsv("FULL_MACHINE_COMPARISON.tsv",cdf)
    write_json("FULL_MACHINE_COMPARISON.json",comparisons)
    write_tsv("INPUT_RECHECK.tsv",input_recheck)
    write_tsv("DEFINITIVE_RUN_PAIR.tsv",runpair)
    summary={"counts":counts,"targets_used_cptac":{m:data[m]["targets_used"] for m in MODELS[2:]},
      "comparison_fields":len(cdf),"comparison_matches":int(cdf.match.sum()),"comparison_failures":int((~cdf.match).sum()),
      "max_abs_diff":float(cdf.abs_diff.max(skipna=True)),"max_rel_diff":float(cdf.rel_diff.max(skipna=True)),
      "identity_max":{k:max(r[k] for r in identity) for k in ["b_abs_error","d_abs_error","bootstrap_max_b_scaled_error","bootstrap_max_d_scaled_error"]},
      "schema_correction":{"scientific_common_fields_equal":correction_equal,"new_columns":added},
      "definitive_runs":{"files":len(runpair),"all_byte_identical":all(r["byte_identical"] for r in runpair)},
      "input_recheck":{"rows":len(input_recheck),"all_match":all(r["match"] for r in input_recheck)},
      "preservation":{"start":start,"end":end,"unchanged":start==end},
      "environment":{"python":sys.version,"numpy":np.__version__,"pandas":pd.__version__,"scipy":scipy.__version__,
                     "platform":platform.platform(),"OPENBLAS_NUM_THREADS":os.environ.get("OPENBLAS_NUM_THREADS"),
                     "OMP_NUM_THREADS":os.environ.get("OMP_NUM_THREADS"),"PYTHONHASHSEED":os.environ.get("PYTHONHASHSEED")}}
    write_json("SUMMARY.json",summary)
    print(json.dumps(summary,indent=2,default=str))


if __name__ == "__main__": main()
