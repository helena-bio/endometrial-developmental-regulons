"""
TASK-028 Freeze B v3 -- model engine (LOCKED; spec-exact B2.4/B2.5/B2.8).
OLS with subtype (4-level, 3 df) + M4 + [purity] + composition. EMMs at reference
covariate values. Fixed full-model contrasts C1/C2/C3 with EQUAL-within-side pooling
(unweighted subtype-EMM means). 3-df omnibus (Wald on the subtype factor) + Kruskal-
Wallis companion. Standardized effects d = contrast / sigma_resid. Patient-level
bootstrap CIs. Subtype-label permutation nulls. Deterministic seeds.

No tuning. ASCII only.
"""
import numpy as np
from scipy import stats

SUBTYPES = ["POLE", "MMRd", "NSMP", "p53abn"]
# contrast weights on [POLE, MMRd, NSMP, p53abn]; difference-of-means (equal-within-side)
# C1 = p53abn vs mean(POLE,MMRd,NSMP); C2 = mean(POLE,MMRd) vs NSMP; C3 = POLE vs MMRd
CONTRAST_L = {
    "C1": np.array([-1/3, -1/3, -1/3, 1.0]),
    "C2": np.array([0.5, 0.5, -1.0, 0.0]),
    "C3": np.array([1.0, -1.0, 0.0, 0.0]),
}


def build_design(subtype_vec, cov_dict, cov_names):
    """Full design: intercept + 3 subtype dummies (ref POLE) + covariates (in cov_names order).
    subtype_vec: array of subtype labels (n,). cov_dict: name->array(n,).
    Returns X (n x p), and info to map subtype EMMs.
    Subtype coded with POLE as reference; dummy order [MMRd, NSMP, p53abn]."""
    n = len(subtype_vec)
    dummies = np.zeros((n, 3))
    for j, s in enumerate(["MMRd", "NSMP", "p53abn"]):
        dummies[:, j] = (subtype_vec == s).astype(float)
    cols = [np.ones(n), dummies]
    for cn in cov_names:
        cols.append(np.asarray(cov_dict[cn], dtype=float).reshape(n))
    X = np.column_stack(cols)
    return X


def _cov_ref_values(cov_dict, cov_names, idx):
    """Reference covariate values = mean over the FITTED sample (EMM convention)."""
    return {cn: float(np.mean(np.asarray(cov_dict[cn])[idx])) for cn in cov_names}


def fit_ols(y, X):
    """Plain OLS via lstsq. Returns beta, resid, sigma_resid (sqrt of MSE with p df),
    XtX_inv, dof_resid."""
    n, p = X.shape
    beta, _, rank, _ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    dof = n - p
    sse = float(resid @ resid)
    sigma2 = sse / dof
    XtX_inv = np.linalg.inv(X.T @ X)
    return beta, resid, np.sqrt(sigma2), XtX_inv, dof, sigma2


def subtype_emms(beta, cov_ref, cov_names):
    """Covariate-adjusted subtype marginal means at reference covariate values.
    beta layout: [intercept, dMMRd, dNSMP, dp53abn, cov1, cov2, ...]."""
    b0 = beta[0]
    dMMRd, dNSMP, dp53 = beta[1], beta[2], beta[3]
    cov_contrib = sum(beta[4 + k] * cov_ref[cov_names[k]] for k in range(len(cov_names)))
    mu = {
        "POLE": b0 + cov_contrib,
        "MMRd": b0 + dMMRd + cov_contrib,
        "NSMP": b0 + dNSMP + cov_contrib,
        "p53abn": b0 + dp53 + cov_contrib,
    }
    return mu


def contrast_estimate(mu, cname):
    L = CONTRAST_L[cname]
    return float(sum(L[i] * mu[SUBTYPES[i]] for i in range(4)))


def omnibus_wald(beta, XtX_inv, sigma2, dof):
    """3-df Wald test on the subtype factor (dummy coefficients 1,2,3 jointly = 0).
    F stat and p. Also return the chi2/LRT-style via F*3."""
    idx = [1, 2, 3]
    R = np.zeros((3, len(beta)))
    for r, j in enumerate(idx):
        R[r, j] = 1.0
    Rb = R @ beta
    M = R @ (sigma2 * XtX_inv) @ R.T
    Minv = np.linalg.inv(M)
    wald_chi2 = float(Rb @ Minv @ Rb)          # ~ chi2(3) asymptotically
    F = wald_chi2 / 3.0                          # F(3, dof)
    p_F = float(stats.f.sf(F, 3, dof))
    return {"wald_chi2_3df": wald_chi2, "F": F, "df1": 3, "df2": int(dof), "p_F": p_F}


def kruskal_companion(y, subtype_vec):
    groups = [y[subtype_vec == s] for s in SUBTYPES]
    H, p = stats.kruskal(*groups)
    return {"H": float(H), "p": float(p)}


def fit_module(y_all, subtype_all, cov_dict_all, cov_names, sample_mask):
    """Fit one module on the complete-case subset given by sample_mask.
    Returns a result dict: EMMs, contrasts (estimate + d), sigma_resid, omnibus, KW,
    and the arrays needed for bootstrap/permutation."""
    idx = np.where(sample_mask)[0]
    y = np.asarray(y_all, dtype=float)[idx]
    subt = np.asarray(subtype_all)[idx]
    cov_local = {cn: np.asarray(cov_dict_all[cn], dtype=float)[idx] for cn in cov_names}
    X = build_design(subt, cov_local, cov_names)
    beta, resid, sigma_resid, XtX_inv, dof, sigma2 = fit_ols(y, X)
    cov_ref = _cov_ref_values(cov_dict_all, cov_names, idx)
    mu = subtype_emms(beta, cov_ref, cov_names)
    out = {"n": int(len(idx)), "sigma_resid": sigma_resid,
           "emm": mu, "contrasts": {}, "beta": beta.tolist()}
    for cname in ["C1", "C2", "C3"]:
        est = contrast_estimate(mu, cname)
        d = est / sigma_resid if sigma_resid > 0 else np.nan
        out["contrasts"][cname] = {"estimate": est, "d": d}
    out["omnibus"] = omnibus_wald(beta, XtX_inv, sigma2, dof)
    out["kruskal"] = kruskal_companion(y, subt)
    return out


def bootstrap_contrasts(y_all, subtype_all, cov_dict_all, cov_names, sample_mask,
                        seed, n_boot=2000):
    """Patient-level (case) bootstrap of the contrast estimates and d.
    Resample patients WITH replacement from the fitted (complete-case) sample; refit;
    recompute contrasts + d. Percentile 95% CIs. Deterministic (seed)."""
    idx = np.where(sample_mask)[0]
    y = np.asarray(y_all, dtype=float)[idx]
    subt = np.asarray(subtype_all)[idx]
    cov_local = {cn: np.asarray(cov_dict_all[cn], dtype=float)[idx] for cn in cov_names}
    n = len(idx)
    rng = np.random.default_rng(seed)
    boot = {c: {"est": np.full(n_boot, np.nan), "d": np.full(n_boot, np.nan)}
            for c in ["C1", "C2", "C3"]}
    n_fail = 0
    for b in range(n_boot):
        samp = rng.integers(0, n, size=n)
        ys = y[samp]; ss = subt[samp]
        # require all 4 subtypes present for identifiability
        if len(set(ss.tolist())) < 4:
            n_fail += 1
            continue
        cs = {cn: cov_local[cn][samp] for cn in cov_names}
        Xb = build_design(ss, cs, cov_names)
        try:
            beta, resid, sigma_resid, XtX_inv, dof, sigma2 = fit_ols(ys, Xb)
        except np.linalg.LinAlgError:
            n_fail += 1
            continue
        cov_ref = {cn: float(np.mean(cs[cn])) for cn in cov_names}
        mu = subtype_emms(beta, cov_ref, cov_names)
        for cname in ["C1", "C2", "C3"]:
            est = contrast_estimate(mu, cname)
            boot[cname]["est"][b] = est
            boot[cname]["d"][b] = est / sigma_resid if sigma_resid > 0 else np.nan
    ci = {}
    for cname in ["C1", "C2", "C3"]:
        e = boot[cname]["est"][np.isfinite(boot[cname]["est"])]
        d = boot[cname]["d"][np.isfinite(boot[cname]["d"])]
        ci[cname] = {
            "est_ci_lo": float(np.percentile(e, 2.5)) if e.size else np.nan,
            "est_ci_hi": float(np.percentile(e, 97.5)) if e.size else np.nan,
            "d_ci_lo": float(np.percentile(d, 2.5)) if d.size else np.nan,
            "d_ci_hi": float(np.percentile(d, 97.5)) if d.size else np.nan,
            "n_boot_used": int(e.size),
        }
    ci["_n_boot_failed"] = int(n_fail)
    return ci


def permutation_omnibus_and_contrasts(y_all, subtype_all, cov_dict_all, cov_names,
                                      sample_mask, seed, n_perm=2000):
    """Subtype-label permutation null (B2.4 / F3). Permute subtype labels among the
    fitted patients; refit; record omnibus F and |contrast d|. Two-sided permutation
    p-values for each contrast (|d_perm| >= |d_obs|), and omnibus p (F_perm >= F_obs).
    Deterministic (seed)."""
    idx = np.where(sample_mask)[0]
    y = np.asarray(y_all, dtype=float)[idx]
    subt = np.asarray(subtype_all)[idx]
    cov_local = {cn: np.asarray(cov_dict_all[cn], dtype=float)[idx] for cn in cov_names}
    cov_ref = {cn: float(np.mean(cov_local[cn])) for cn in cov_names}
    # observed
    X = build_design(subt, cov_local, cov_names)
    beta, resid, sigma_resid, XtX_inv, dof, sigma2 = fit_ols(y, X)
    mu = subtype_emms(beta, cov_ref, cov_names)
    obs = {c: abs(contrast_estimate(mu, c) / sigma_resid) for c in ["C1", "C2", "C3"]}
    obs_F = omnibus_wald(beta, XtX_inv, sigma2, dof)["F"]
    rng = np.random.default_rng(seed)
    ge_contrast = {c: 1 for c in ["C1", "C2", "C3"]}  # +1 for observed (add-one)
    ge_F = 1
    for _ in range(n_perm):
        perm = rng.permutation(subt)
        Xp = build_design(perm, cov_local, cov_names)
        bp, rp, sp, xi, dfp, s2p = fit_ols(y, Xp)
        mup = subtype_emms(bp, cov_ref, cov_names)
        for c in ["C1", "C2", "C3"]:
            if sp > 0 and abs(contrast_estimate(mup, c) / sp) >= obs[c] - 1e-12:
                ge_contrast[c] += 1
        Fp = omnibus_wald(bp, xi, s2p, dfp)["F"]
        if Fp >= obs_F - 1e-12:
            ge_F += 1
    denom = n_perm + 1
    return {
        "perm_p_omnibus": ge_F / denom,
        "perm_p_contrast": {c: ge_contrast[c] / denom for c in ["C1", "C2", "C3"]},
        "n_perm": n_perm,
    }


def bh_fdr(pvals):
    """Benjamini-Hochberg q-values. pvals: dict key->p. Returns dict key->q."""
    items = [(k, v) for k, v in pvals.items() if v is not None and np.isfinite(v)]
    m = len(items)
    if m == 0:
        return {k: np.nan for k in pvals}
    items.sort(key=lambda kv: kv[1])
    q = {}
    prev = 1.0
    for rank in range(m - 1, -1, -1):
        k, p = items[rank]
        val = p * m / (rank + 1)
        prev = min(prev, val)
        q[k] = min(prev, 1.0)
    for k in pvals:
        if k not in q:
            q[k] = np.nan
    return q
