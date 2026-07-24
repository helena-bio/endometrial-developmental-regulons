#!/usr/bin/env python3
"""TASK-030 cycle 2: frozen any-single-mapped-gene LOO gate.

Reconstruct the sealed signed aREA regulons, delete every admitted edge for one
mapped target symbol at a time, refit the locked PRIMARY and SENS_nopurity
models, and reconstruct F1/F2 by substituting only the altered module into the
corresponding sealed full family. No target, edge, threshold, orientation,
contrast, model, covariate, family, or category is tuned. ASCII only.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy import stats
from scipy.stats import norm, rankdata


ROOT = Path("data/external/original-workspace")
WORKTREE = ROOT / ".language-model tool/work/revgate-tcga-no-purity-verify"
OUT = WORKTREE / "experiments/task030_verify_current_main/gene_loo"
RESULTS = OUT / "results"
MANUSCRIPT = OUT / "manuscript"
T28 = ROOT / ".language-model tool/work/task028-freeze-b-draft"
T29 = ROOT / ".language-model tool/work/task029-external-replication-feasibility"
T30C1 = ROOT / ".language-model tool/work/task030-tcga-nopurity-audit"
INTER = T28 / "execution/intermediate"
R28 = T28 / "execution/results_v3"
R29 = T29 / "execution/results"
LEDGER = T28 / "sealed_v3/B1/m3_primary_edge_ledger.tsv"
SEED_SCHEME = OUT / "SEED_SCHEME.md"

MASTER_SEED = 20260713
N_BOOT = 2000
N_PERM = 2000
FLOOR = 0.50
Q_CUTOFF = 0.05
SEED_SCHEME_SHA256 = "32c03f44e789e9bf4c42a12d5508a5d4e27ef84c8806b46629f07c83032e512c"
SUBTYPES = ["POLE", "MMRd", "NSMP", "p53abn"]
CONTRASTS = {
    "C1": np.array([-1 / 3, -1 / 3, -1 / 3, 1.0]),
    "C2": np.array([0.5, 0.5, -1.0, 0.0]),
    "C3": np.array([1.0, -1.0, 0.0, 0.0]),
}
TARGETS = [
    ("GATA2", "C2", -1),
    ("SOX9", "C2", -1),
    ("HOXA9", "C2", -1),
    ("WT1", "C2", -1),
    ("PAX8", "C1", 1),
    ("LHX1", "C1", 1),
]
TARGET_ORDER = {x[0]: i for i, x in enumerate(TARGETS)}
CONFIGS = {
    "PRIMARY": {"cov_names": ["M4", "purity_CPE", "composition"], "mask": "mask506"},
    "SENS_nopurity": {"cov_names": ["M4", "composition"], "mask": "mask507"},
}


def utcnow() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def step_seed(step_id: str) -> int:
    return int(hashlib.sha256(f"{MASTER_SEED}:{step_id}".encode("ascii")).hexdigest()[:8], 16)


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    df.to_csv(path, sep="\t", index=False, float_format="%.17g")


def bh_fdr(pvals: dict[str, float]) -> dict[str, float]:
    items = [(k, float(v)) for k, v in pvals.items() if v is not None and np.isfinite(v)]
    items.sort(key=lambda kv: kv[1])
    out: dict[str, float] = {}
    prev = 1.0
    for i in range(len(items) - 1, -1, -1):
        key, p = items[i]
        prev = min(prev, p * len(items) / (i + 1), 1.0)
        out[key] = prev
    return out


def design(subtype: np.ndarray, covariates: list[np.ndarray]) -> np.ndarray:
    dummies = np.column_stack([(subtype == s).astype(float) for s in SUBTYPES[1:]])
    return np.column_stack([np.ones(len(subtype)), dummies, *covariates])


def fit_core(y: np.ndarray, subtype: np.ndarray, covariates: list[np.ndarray]) -> dict:
    X = design(subtype, covariates)
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    resid = y - X @ beta
    dof = len(y) - X.shape[1]
    sigma2 = float(resid @ resid) / dof
    sigma = math.sqrt(sigma2)
    means = np.array([beta[0], beta[0] + beta[1], beta[0] + beta[2], beta[0] + beta[3]])
    estimates = {c: float(L @ means) for c, L in CONTRASTS.items()}
    ds = {c: estimates[c] / sigma for c in CONTRASTS}
    R = np.zeros((3, X.shape[1]))
    R[:, 1:4] = np.eye(3)
    xtx_inv = np.linalg.inv(X.T @ X)
    rb = R @ beta
    middle = R @ (sigma2 * xtx_inv) @ R.T
    wald = float(rb @ np.linalg.inv(middle) @ rb)
    fstat = wald / 3.0
    return {
        "n": len(y), "beta": beta, "resid": resid, "sigma": sigma,
        "sigma2": sigma2, "dof": dof, "estimates": estimates, "ds": ds,
        "F": fstat, "p_F": float(stats.f.sf(fstat, 3, dof)),
    }


def bootstrap(y: np.ndarray, subtype: np.ndarray, covariates: list[np.ndarray], seed: int) -> dict:
    rng = np.random.default_rng(seed)
    n = len(y)
    est = {c: np.full(N_BOOT, np.nan) for c in CONTRASTS}
    ds = {c: np.full(N_BOOT, np.nan) for c in CONTRASTS}
    failed = 0
    for i in range(N_BOOT):
        take = rng.integers(0, n, size=n)
        ss = subtype[take]
        if len(set(ss.tolist())) < 4:
            failed += 1
            continue
        try:
            fit = fit_core(y[take], ss, [x[take] for x in covariates])
        except np.linalg.LinAlgError:
            failed += 1
            continue
        for c in CONTRASTS:
            est[c][i] = fit["estimates"][c]
            ds[c][i] = fit["ds"][c]
    out = {"failed": failed}
    for c in CONTRASTS:
        e = est[c][np.isfinite(est[c])]
        d = ds[c][np.isfinite(ds[c])]
        out[c] = {
            "est_ci_lo": float(np.percentile(e, 2.5)),
            "est_ci_hi": float(np.percentile(e, 97.5)),
            "d_ci_lo": float(np.percentile(d, 2.5)),
            "d_ci_hi": float(np.percentile(d, 97.5)),
            "n_used": len(e),
        }
    return out


def permutation(y: np.ndarray, subtype: np.ndarray, covariates: list[np.ndarray], seed: int) -> dict:
    observed = fit_core(y, subtype, covariates)
    rng = np.random.default_rng(seed)
    ge_f = 1
    ge_d = {c: 1 for c in CONTRASTS}
    for _ in range(N_PERM):
        fit = fit_core(y, rng.permutation(subtype), covariates)
        if fit["F"] >= observed["F"] - 1e-12:
            ge_f += 1
        for c in CONTRASTS:
            if abs(fit["ds"][c]) >= abs(observed["ds"][c]) - 1e-12:
                ge_d[c] += 1
    denom = N_PERM + 1
    return {
        "p_omnibus": ge_f / denom,
        "p_contrast": {c: ge_d[c] / denom for c in CONTRASTS},
    }


def parse_standard_seal(seal: Path, base: Path) -> list[dict]:
    rows = []
    for raw in seal.read_text(encoding="ascii").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        expected, rel = line.split(None, 1)
        path = (base / rel.strip()).resolve()
        actual = sha256_file(path) if path.is_file() else None
        rows.append({"seal": str(seal), "path": str(path), "expected": expected,
                     "actual": actual, "match": actual == expected})
    return rows


def parse_task028_spec_seal(seal: Path, base: Path) -> list[dict]:
    rows = []
    for raw in seal.read_text(encoding="ascii").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        first, filename = line.split(None, 1)
        if "/" in first:
            prefix, expected = first.rsplit("/", 1)
            path = base / prefix / filename.strip()
        else:
            expected = first
            path = base / filename.strip()
        actual = sha256_file(path) if path.is_file() else None
        rows.append({"seal": str(seal), "path": str(path.resolve()), "expected": expected,
                     "actual": actual, "match": actual == expected})
    return rows


def area_from_transforms(t1: np.ndarray, t2: np.ndarray, edge_rows: np.ndarray,
                         weights: np.ndarray, mor: np.ndarray, keep: np.ndarray) -> np.ndarray:
    ii = edge_rows[keep]
    w = weights[keep]
    m = mor[keep]
    wts = w / w.max()
    wtss = wts / wts.sum()
    s2 = np.sum((wtss * m)[:, None] * t2[ii], axis=0)
    s1 = np.sum(wtss[:, None] * t1[ii], axis=0)
    es = (np.abs(s2) + np.where(s1 > 0, s1, 0.0)) * np.sign(s2)
    return es * math.sqrt(float(wts @ wts))


def module_family(sealed: dict, config: str) -> dict:
    return sealed["primary"] if config == "PRIMARY" else sealed["sensitivities"][config]


def substitute_families(sealed_family: dict, changed_module: str, changed_fit: dict,
                        changed_perm: dict) -> dict:
    f1_p = {m: float(v["omnibus"]["p_F"]) for m, v in sealed_family.items()}
    before_f1 = dict(f1_p)
    f1_p[changed_module] = changed_fit["p_F"]
    f1_q = bh_fdr(f1_p)
    gated = sorted(m for m in f1_p if f1_q[m] <= Q_CUTOFF)

    f2_p: dict[str, float] = {}
    unchanged_f2_preserved = True
    for module in gated:
        for contrast in CONTRASTS:
            key = f"{module}::{contrast}"
            if module == changed_module:
                f2_p[key] = float(changed_perm["p_contrast"][contrast])
            else:
                f2_p[key] = float(sealed_family[module]["permutation"]["perm_p_contrast"][contrast])
                unchanged_f2_preserved &= (
                    f2_p[key] == float(sealed_family[module]["permutation"]["perm_p_contrast"][contrast])
                )
    f2_q = bh_fdr(f2_p)
    unchanged_f1_preserved = all(
        f1_p[m] == before_f1[m] for m in f1_p if m != changed_module
    )
    changed_keys = [f"{changed_module}::{c}" for c in CONTRASTS]
    return {
        "f1_q": f1_q, "gated": gated, "f2_p": f2_p, "f2_q": f2_q,
        "changed_f1_q": f1_q[changed_module],
        "changed_f1_gate": changed_module in gated,
        "target_f2_q": {c: f2_q.get(f"{changed_module}::{c}") for c in CONTRASTS},
        "f1_family_size": len(f1_p), "f1_substitution_count": 1,
        "f2_family_size": len(f2_p),
        "f2_changed_entries": sum(k in f2_p for k in changed_keys),
        "unchanged_f1_preserved": unchanged_f1_preserved,
        "unchanged_f2_preserved": unchanged_f2_preserved,
    }


def md_table(df: pd.DataFrame) -> str:
    header = "| " + " | ".join(df.columns) + " |"
    rule = "| " + " | ".join(["---"] * len(df.columns)) + " |"
    rows = ["| " + " | ".join(str(v).replace("|", "\\|") for v in row) + " |"
            for row in df.itertuples(index=False, name=None)]
    return "\n".join([header, rule, *rows])


def main() -> None:
    started = utcnow()
    RESULTS.mkdir(parents=True, exist_ok=True)
    MANUSCRIPT.mkdir(parents=True, exist_ok=True)
    assert sha256_file(SEED_SCHEME) == SEED_SCHEME_SHA256, "frozen seed scheme changed"

    seal_rows = []
    seal_rows += parse_task028_spec_seal(T28 / "sealed_v3/SEAL_MANIFEST.sha256", T28 / "sealed_v3")
    seal_rows += parse_standard_seal(R28 / "RESULT_SEAL.sha256", R28)
    seal_rows += parse_standard_seal(R29 / "RESULT_SEAL.sha256", R29)
    seal_df = pd.DataFrame(seal_rows)
    write_tsv(seal_df, RESULTS / "upstream_checksum_verification.tsv")
    if not bool(seal_df["match"].all()):
        raise RuntimeError("upstream checksum mismatch; analysis aborted")

    sealed = json.loads((R28 / "phase2b_models.json").read_text(encoding="ascii"))
    cptac = json.loads((R29 / "primary_results.json").read_text(encoding="ascii"))

    # No-change family reconstruction: if no module is altered, the local BH and
    # gated-family logic must reproduce every sealed F1/F2 q-value.
    family_audit_rows = []
    for config in CONFIGS:
        sf = module_family(sealed, config)
        f1_p = {m: float(v["omnibus"]["p_F"]) for m, v in sf.items()}
        f1_q = bh_fdr(f1_p)
        gated = sorted(m for m in sf if f1_q[m] <= Q_CUTOFF)
        f1_delta = max(abs(f1_q[m] - float(sf[m]["omnibus"]["BH_q_F1"])) for m in sf)
        f1_gate_match = all((m in gated) == bool(sf[m]["omnibus"]["F1_gate_pass"]) for m in sf)
        f2_p = {f"{m}::{c}": float(sf[m]["permutation"]["perm_p_contrast"][c])
                for m in gated for c in CONTRASTS}
        f2_q = bh_fdr(f2_p)
        f2_delta = max(abs(f2_q[k] - float(sf[k.split("::")[0]]["contrasts"][k.split("::")[1]]["BH_q_F2"]))
                       for k in f2_q)
        family_audit_rows.append({
            "config": config, "F1_family_size": len(f1_p), "F1_gated_module_count": len(gated),
            "F2_family_size": len(f2_p), "F1_q_max_abs_delta_vs_sealed": f1_delta,
            "F1_gate_membership_exact_match": f1_gate_match,
            "F2_q_max_abs_delta_vs_sealed": f2_delta,
            "no_change_reconstruction_pass": bool(f1_delta <= 1e-15 and f1_gate_match and f2_delta <= 1e-15),
        })
    family_audit_df = pd.DataFrame(family_audit_rows)
    write_tsv(family_audit_df, RESULTS / "family_reconstruction_audit.tsv")
    if not bool(family_audit_df["no_change_reconstruction_pass"].all()):
        raise RuntimeError("sealed no-change family reconstruction failed")

    edges = pd.read_csv(LEDGER, sep="\t", comment="#", dtype=str)
    edges = edges[edges["primary_inclusion"] == "Y"].copy()
    names = json.loads((INTER / "gene_names_v3.json").read_text(encoding="ascii"))
    expr = np.load(INTER / "log2tpm_v3.npy", mmap_mode="r")
    saved_scores = np.load(INTER / "scores_v3.npz", allow_pickle=True)
    cov = pd.read_csv(INTER / "covariates_v3.tsv", sep="\t")
    subtype_all = cov["subtype"].to_numpy()
    cov_all = {
        "M4": cov["M4_prolif"].to_numpy(float),
        "purity_CPE": cov["purity_CPE"].to_numpy(float),
        "composition": cov["composition"].to_numpy(float),
    }
    masks = {
        "mask506": cov["cpe_complete_case"].to_numpy(int) == 1,
        "mask507": np.ones(len(cov), dtype=bool),
    }

    print("[score] build exact aREA rank transforms", utcnow(), flush=True)
    means = np.asarray(expr.mean(axis=1), float)
    sds = np.asarray(expr.std(axis=1, ddof=1), float)
    finite = np.isfinite(means) & np.isfinite(sds) & (sds > 0)
    finite_rows = np.where(finite)[0]
    by_symbol: dict[str, list[int]] = {}
    for i, symbol in enumerate(names):
        by_symbol.setdefault(symbol, []).append(i)
    best: dict[str, int] = {}
    for symbol in edges["target"].unique():
        candidates = by_symbol.get(symbol, [])
        if candidates:
            row = max(candidates, key=lambda x: means[x])
            if finite[row]:
                best[symbol] = row
    union_rows = sorted(set(best.values()))
    union_index = {row: i for i, row in enumerate(union_rows)}
    reduced_pos = -np.ones(len(names), dtype=int)
    reduced_pos[finite_rows] = np.arange(len(finite_rows))
    selected_reduced = np.array([reduced_pos[r] for r in union_rows], int)
    t1 = np.empty((len(union_rows), expr.shape[1]), float)
    t2 = np.empty_like(t1)
    for sample in range(expr.shape[1]):
        z = (np.asarray(expr[finite_rows, sample], float) - means[finite_rows]) / sds[finite_rows]
        signed_rank = rankdata(z, method="average") / (len(finite_rows) + 1.0)
        abs_rank = rankdata(np.abs(z), method="average") / (len(finite_rows) + 1.0)
        t2[:, sample] = norm.ppf(signed_rank[selected_reduced])
        t1[:, sample] = norm.ppf((abs_rank[selected_reduced] + 1.0) / 2.0)

    score_audit = []
    regulons: dict[str, dict] = {}
    for tf in sorted(edges["TF"].unique()):
        listed = edges[edges["TF"] == tf].copy()
        mapped = listed[listed["target"].isin(best)].reset_index(drop=True)
        edge_rows = np.array([union_index[best[x]] for x in mapped["target"]], int)
        weights = mapped["weight"].astype(float).to_numpy()
        mor = mapped["primary_sign"].astype(float).to_numpy()
        full = area_from_transforms(t1, t2, edge_rows, weights, mor,
                                    np.ones(len(mapped), dtype=bool))
        delta = float(np.max(np.abs(full - np.asarray(saved_scores[f"M3primary__{tf}"], float))))
        score_audit.append({
            "target": tf, "listed_primary_edges": len(listed), "mapped_edges": len(mapped),
            "unique_mapped_genes": int(mapped["target"].nunique()),
            "full_score_max_abs_delta_vs_sealed": delta,
            "match_tolerance_1e-12": delta <= 1e-12,
        })
        regulons[tf] = {"mapped": mapped, "edge_rows": edge_rows, "weights": weights,
                        "mor": mor, "full": full}
    score_audit_df = pd.DataFrame(score_audit)
    write_tsv(score_audit_df, RESULTS / "full_score_reconstruction_audit.tsv")
    if not bool(score_audit_df["match_tolerance_1e-12"].all()):
        raise RuntimeError("full score reconstruction mismatch")

    point_rows: list[dict] = []
    score_by_deletion: dict[tuple[str, str], np.ndarray] = {}
    print("[point] all deletions", utcnow(), flush=True)
    for target, contrast, direction in TARGETS:
        reg = regulons[target]
        mapped = reg["mapped"]
        for gene in mapped["target"].drop_duplicates().tolist():
            dropped = mapped["target"].to_numpy() == gene
            score = area_from_transforms(reg["t1"] if "t1" in reg else t1,
                                         reg["t2"] if "t2" in reg else t2,
                                         reg["edge_rows"], reg["weights"], reg["mor"], ~dropped)
            score_by_deletion[(target, gene)] = score
            for config, spec in CONFIGS.items():
                mask = masks[spec["mask"]]
                fit = fit_core(score[mask], subtype_all[mask],
                               [cov_all[x][mask] for x in spec["cov_names"]])
                d = fit["ds"][contrast]
                boot_id = f"loo_boot__{config}__M3_{target}__drop_{gene}"
                perm_id = f"loo_perm__{config}__M3_{target}__drop_{gene}"
                point_rows.append({
                    "target": target, "config": config, "contrast": contrast,
                    "locked_direction": "negative" if direction < 0 else "positive",
                    "dropped_gene": gene, "dropped_edge_count": int(dropped.sum()),
                    "dropped_edge_weights": ";".join(mapped.loc[dropped, "weight"].tolist()),
                    "dropped_edge_signs": ";".join(mapped.loc[dropped, "primary_sign"].tolist()),
                    "remaining_edges": int((~dropped).sum()), "n": fit["n"],
                    "b": fit["estimates"][contrast], "residual_sd": fit["sigma"], "d": d,
                    "direction_pass": bool(np.sign(d) == direction),
                    "floor_pass": bool(abs(d) >= FLOOR),
                    "deterministic_gate_pass": bool(np.sign(d) == direction and abs(d) >= FLOOR),
                    "bootstrap_step_id": boot_id, "bootstrap_seed": step_seed(boot_id),
                    "permutation_step_id": perm_id, "permutation_seed": step_seed(perm_id),
                })
    point_df = pd.DataFrame(point_rows)
    all_point_pass = point_df.groupby(["target", "config"])["deterministic_gate_pass"].all().to_dict()
    point_df["target_config_all_deletions_point_pass"] = [
        bool(all_point_pass[(r.target, r.config)]) for r in point_df.itertuples(index=False)
    ]
    point_df = point_df.sort_values(["target", "config", "dropped_gene"],
                                    key=lambda s: s.map(TARGET_ORDER) if s.name == "target" else s)
    write_tsv(point_df, RESULTS / "gene_loo_all_point_rows.tsv")

    stochastic_rows: list[dict] = []
    print("[stochastic] gated deletions", utcnow(), flush=True)
    for row in point_df.to_dict("records"):
        target, config, gene = row["target"], row["config"], row["dropped_gene"]
        contrast = row["contrast"]
        module = f"M3_{target}"
        base = {
            "target": target, "config": config, "contrast": contrast,
            "dropped_gene": gene, "bootstrap_step_id": row["bootstrap_step_id"],
            "bootstrap_seed": row["bootstrap_seed"],
            "permutation_step_id": row["permutation_step_id"],
            "permutation_seed": row["permutation_seed"],
            "b": row["b"], "residual_sd": row["residual_sd"], "d": row["d"],
            "direction_pass": row["direction_pass"], "floor_pass": row["floor_pass"],
        }
        if not all_point_pass[(target, config)]:
            stochastic_rows.append({
                **base, "stochastic_status": "not_run_by_deterministic_any_deletion_gate",
                "not_run_reason": "target-config has at least one deletion failing locked direction or |d|>=0.50",
                "n_boot": N_BOOT, "n_boot_used": None, "n_boot_failed": None,
                "d_ci_lo": None, "d_ci_hi": None, "ci_excludes_0": None,
                "n_perm": N_PERM, "omnibus_F": None, "omnibus_p_analytic": None,
                "omnibus_perm_p": None, "omnibus_F1_q": None, "F1_gate_pass": None,
                "contrast_perm_p": None, "contrast_F2_q": None,
                "f1_family_size": None, "f1_substitution_count": None,
                "f2_family_size": None, "f2_changed_entries": None,
                "unchanged_f1_preserved": None, "unchanged_f2_preserved": None,
                "family_substitution_check": "not_run_by_gate", "complete_deletion_gate": False,
                "deletion_failure_reason": "deterministic target-config short-circuit",
            })
            continue
        spec = CONFIGS[config]
        mask = masks[spec["mask"]]
        y = score_by_deletion[(target, gene)][mask]
        ss = subtype_all[mask]
        cc = [cov_all[x][mask] for x in spec["cov_names"]]
        print(f"  [run] {config} {target} drop {gene}", flush=True)
        fit = fit_core(y, ss, cc)
        boot = bootstrap(y, ss, cc, int(row["bootstrap_seed"]))
        perm = permutation(y, ss, cc, int(row["permutation_seed"]))
        fam = substitute_families(module_family(sealed, config), module, fit, perm)
        ci = boot[contrast]
        ci_excl = bool(ci["d_ci_lo"] > 0 or ci["d_ci_hi"] < 0)
        f2q = fam["target_f2_q"][contrast]
        complete = bool(row["direction_pass"] and row["floor_pass"] and ci_excl and
                        fam["changed_f1_gate"] and f2q is not None and f2q <= Q_CUTOFF)
        failures = []
        if not row["direction_pass"]:
            failures.append("wrong_direction")
        if not row["floor_pass"]:
            failures.append("abs_d_below_0.50")
        if not ci_excl:
            failures.append("patient_bootstrap_CI_includes_0")
        if not fam["changed_f1_gate"]:
            failures.append("module_omnibus_F1_q_above_0.05")
        if f2q is None or f2q > Q_CUTOFF:
            failures.append("planned_contrast_full_gated_F2_q_above_0.05")
        substitution_ok = bool(fam["f1_substitution_count"] == 1 and
                               fam["unchanged_f1_preserved"] and fam["unchanged_f2_preserved"] and
                               fam["f1_family_size"] == len(module_family(sealed, config)))
        if not substitution_ok:
            raise RuntimeError("family substitution invariant failed")
        stochastic_rows.append({
            **base, "stochastic_status": "run", "not_run_reason": "",
            "n_boot": N_BOOT, "n_boot_used": ci["n_used"], "n_boot_failed": boot["failed"],
            "d_ci_lo": ci["d_ci_lo"], "d_ci_hi": ci["d_ci_hi"],
            "ci_excludes_0": ci_excl, "n_perm": N_PERM,
            "omnibus_F": fit["F"], "omnibus_p_analytic": fit["p_F"],
            "omnibus_perm_p": perm["p_omnibus"],
            "omnibus_F1_q": fam["changed_f1_q"], "F1_gate_pass": fam["changed_f1_gate"],
            "contrast_perm_p": perm["p_contrast"][contrast], "contrast_F2_q": f2q,
            "f1_family_size": fam["f1_family_size"],
            "f1_substitution_count": fam["f1_substitution_count"],
            "f2_family_size": fam["f2_family_size"],
            "f2_changed_entries": fam["f2_changed_entries"],
            "unchanged_f1_preserved": fam["unchanged_f1_preserved"],
            "unchanged_f2_preserved": fam["unchanged_f2_preserved"],
            "family_substitution_check": "PASS_one_F1_and_three_F2_values_substituted_if_gated",
            "complete_deletion_gate": complete,
            "deletion_failure_reason": ";".join(failures),
        })
    stochastic_df = pd.DataFrame(stochastic_rows)
    stochastic_df = stochastic_df.sort_values(["target", "config", "dropped_gene"],
                                               key=lambda s: s.map(TARGET_ORDER) if s.name == "target" else s)
    write_tsv(stochastic_df, RESULTS / "gene_loo_stochastic_full_gate_rows.tsv")

    summary_rows = []
    for target, contrast, direction in TARGETS:
        module = f"M3_{target}"
        for config in CONFIGS:
            p = point_df[(point_df.target == target) & (point_df.config == config)]
            s = stochastic_df[(stochastic_df.target == target) & (stochastic_df.config == config)]
            legacy = module_family(sealed, config)[module]["b3_category"][contrast]["category"]
            point_ok = bool(p["deterministic_gate_pass"].all())
            complete = bool(point_ok and s["complete_deletion_gate"].all())
            worst = p.loc[p["d"].abs().idxmin()]
            floor_fail = int((~p["floor_pass"]).sum())
            direction_fail = int((~p["direction_pass"]).sum())
            stochastic_run = int((s["stochastic_status"] == "run").sum())
            full_fail = int((s["stochastic_status"].eq("run") & ~s["complete_deletion_gate"]).sum())
            if complete:
                credit = legacy + "_CREDITED_B2.12"
                reason = "all mapped single-gene deletions pass direction, |d|>=0.50, bootstrap CI, F1, and full gated F2"
            else:
                credit = "UNCREDITED_DESCRIPTIVE_POSITIVE_GENE_LOO_FAILURE_NOT_CAT5"
                pieces = []
                if direction_fail:
                    pieces.append(f"wrong_direction_in_{direction_fail}_of_{len(p)}_deletions")
                if floor_fail:
                    pieces.append(f"abs_d_below_0.50_in_{floor_fail}_of_{len(p)}_deletions")
                if full_fail:
                    pieces.append(f"stochastic_or_family_gate_failed_in_{full_fail}_of_{stochastic_run}_run_deletions")
                reason = ";".join(pieces) or "complete_any-deletion_gate_not_satisfied"
            run_s = s[s["stochastic_status"] == "run"]
            summary_rows.append({
                "target": target, "config": config, "contrast": contrast,
                "locked_direction": "negative" if direction < 0 else "positive",
                "unique_mapped_genes_tested": len(p), "total_mapped_edges": int(regulons[target]["mapped"].shape[0]),
                "direction_fail_count": direction_fail, "floor_fail_count": floor_fail,
                "minimum_abs_d": float(p["d"].abs().min()), "maximum_abs_d": float(p["d"].abs().max()),
                "worst_deletion": worst["dropped_gene"], "worst_deletion_d": worst["d"],
                "all_deletion_point_gate_pass": point_ok,
                "stochastic_deletions_run": stochastic_run,
                "stochastic_deletions_not_run_by_gate": int((s["stochastic_status"] != "run").sum()),
                "minimum_bootstrap_CI_signed_margin_from_zero": (
                    float(np.minimum(run_s["d_ci_lo"].abs(), run_s["d_ci_hi"].abs()).min())
                    if len(run_s) else None
                ),
                "maximum_F1_q_across_run_deletions": float(run_s["omnibus_F1_q"].max()) if len(run_s) else None,
                "maximum_planned_contrast_F2_q_across_run_deletions": float(run_s["contrast_F2_q"].max()) if len(run_s) else None,
                "complete_gate_fail_count_among_run_deletions": full_fail,
                "survives_gene_LOO_any_deletion_rule": complete,
                "legacy_incomplete_B2_12_label": legacy,
                "B2_12_compliant_credit_status": credit,
                "exact_failure_or_pass_reason": reason,
            })
    summary_df = pd.DataFrame(summary_rows).sort_values(
        ["target", "config"], key=lambda s: s.map(TARGET_ORDER) if s.name == "target" else s)
    write_tsv(summary_df, RESULTS / "gene_loo_per_target_summary.tsv")

    # Scenario C follows the supplied definitions because both headline no-purity
    # rows lose B2.12 credit. The numeric purity/CPTAC magnitude comparison remains B-like.
    headline = summary_df[(summary_df.target.isin(["GATA2", "SOX9"])) &
                          (summary_df.config == "SENS_nopurity")]
    scenario = "C" if not bool(headline["survives_gene_LOO_any_deletion_rule"].all()) else "B"

    # Corrected manuscript-candidate table. Legacy labels are explicitly quarantined.
    table_rows = []
    for target, contrast, _ in TARGETS:
        primary = summary_df[(summary_df.target == target) & (summary_df.config == "PRIMARY")].iloc[0]
        nopurity = summary_df[(summary_df.target == target) & (summary_df.config == "SENS_nopurity")].iloc[0]
        pfit = sealed["primary"][f"M3_{target}"]
        nfit = sealed["sensitivities"]["SENS_nopurity"][f"M3_{target}"]
        cm = cptac["meta"][target]
        cv = cptac["per_target_verdict"][target]
        table_rows.append({
            "Target": target, "Contrast": contrast,
            "TCGA primary d": f"{pfit['contrasts'][contrast]['d']:.3f}",
            "Primary legacy label (incomplete B2.12)": primary["legacy_incomplete_B2_12_label"],
            "Primary B2.12 credit": primary["B2_12_compliant_credit_status"],
            "Primary worst LOO d": f"{primary['worst_deletion_d']:.3f}",
            "TCGA no-purity d": f"{nfit['contrasts'][contrast]['d']:.3f}",
            "No-purity legacy label (incomplete B2.12)": nopurity["legacy_incomplete_B2_12_label"],
            "No-purity B2.12 credit": nopurity["B2_12_compliant_credit_status"],
            "No-purity worst LOO d": f"{nopurity['worst_deletion_d']:.3f}",
            "CPTAC d_meta": f"{cm['d_meta']:.3f}",
            "Frozen CPTAC verdict (bytes unchanged)": cv["status"],
        })
    manuscript_df = pd.DataFrame(table_rows)
    write_tsv(manuscript_df, MANUSCRIPT / "task030_cycle2_B2_12_candidate_table.tsv")
    md = """# TASK-030 cycle-2 B2.12 candidate table

Runtime analytical candidate only. No DOCX or manuscript byte was edited. Legacy
TCGA labels reproduced the incomplete TASK-028 implementation and are not current
B2.12 credit claims. A failed positive is not CAT5; CAT5 remains available only
through the frozen H3 equivalence rule.

""" + md_table(manuscript_df) + "\n"
    (MANUSCRIPT / "task030_cycle2_B2_12_candidate_table.md").write_text(md, encoding="ascii")

    summary_md = summary_df[[
        "target", "config", "unique_mapped_genes_tested", "floor_fail_count",
        "worst_deletion", "worst_deletion_d", "survives_gene_LOO_any_deletion_rule",
        "legacy_incomplete_B2_12_label", "B2_12_compliant_credit_status",
    ]].copy()
    report = f"""# TASK-030 cycle-2 frozen gene-LOO audit

Status: DONE (experimenter computation; critic verdict not assigned here).

## Operational definition and result

"Survives gene LOO" means that the original positive CAT1/CAT2 call remains
creditable for EVERY deletion of one unique mapped target gene, with every admitted
edge for that symbol removed: the locked direction, |d|>=0.50, patient-bootstrap
95% CI excluding 0, module omnibus F1 BH q<=0.05, planned-contrast full gated F2 BH
q<=0.05, and the same M4/purity/composition adjustment configuration must all hold.
BH-of-6 was not used. A deterministic direction or floor failure is definitive;
stochastic work was short-circuited for that target-configuration while point
estimates were retained for every deletion.

Scenario {scenario} under the supplied scenario definitions: GATA2 and SOX9 lose
B2.12 credit in the official n=507 no-purity configuration because at least one
single-gene deletion fails the fixed |d| floor. The same defect occurs in PRIMARY,
so this is not a no-purity-specific loss. Holding the old incomplete labels fixed,
the previously observed magnitude comparison remains numerically B-like rather
than D-like: removing CPE does not move GATA2 or SOX9 closer to CPTAC. That numeric
description does not override Scenario C under the complete frozen category gate.

{md_table(summary_md)}

Failed positives are labelled uncredited/descriptive positives due to gene-LOO
failure, explicitly not CAT5 and not a new frozen category. PAX8 and LHX1 underwent
all 2,000 bootstraps and 2,000 permutations for every deletion and retain credit
only where every complete deletion gate passed. GATA2, SOX9, HOXA9, and WT1 had
deterministic any-deletion floor failures; their stochastic fields are explicitly
not-run-by-gate.

## Full-family substitution

For each run deletion, the altered module's analytic omnibus p replaced exactly
one member of the corresponding 21-module sealed F1 family, after which BH and F1
membership were recomputed. The F2 family was then rebuilt from all modules gated
by that altered F1 result. The altered module's three C1-C3 permutation p-values
replaced its entries when gated; every unchanged module used its corresponding
sealed PRIMARY or SENS_nopurity permutation p-values. Invariants in the stochastic
table test one F1 substitution, unchanged-value preservation, family size, and the
three altered F2 entries. This is the frozen full gated family, not BH across the
six requested targets and not a post-hoc worst-deletion substitute.

## TASK-028/TASK-029 implications

F34 is analytically affected: the four C2 positives GATA2, SOX9, HOXA9, and WT1
are legacy/incomplete-B2.12 positives and are uncredited under the complete gate in
both PRIMARY and no-purity. PAX8 and LHX1 are reported from the completed gate.
No F34 byte was changed. F35/TASK-029 numerical results, tables, and frozen verdict
bytes were not changed; however, their selection as external replication targets
inherits the F34 credit defect and therefore requires separately authorized
reassessment. This runtime audit does not silently amend F35.

## Chronology, comparability, and mandatory bounds

The TCGA no-purity model was pre-specified for TASK-028; its TASK-030 comparison
with the already-observed CPTAC gap is a retrospective audit, not a prospectively
frozen CPTAC-confounding test.

TCGA and CPTAC use rank-derived ssGSEA scores with range normalization performed
separately in TCGA and each CPTAC stratum, and M3 aREA signatures are z-standardized
across samples separately within TCGA, CPTAC Discovery, and CPTAC Confirmatory.
All-zero filtering and duplicate-symbol best-row selection are data-set specific.
Thus CPTAC is model-adapted external replication, not identical covariate replication.

Residual purity and composition confounding remain unexcluded. This is one CPTAC
cohort represented by two strata, target-level only, with method-dependent magnitude,
different subtype-classifier implementations, and relative within-tumour subtype
contrasts. No causal, per-patient biomarker, purity-removal, or cross-cohort absolute-
score claim is made. No source DOCX or existing manuscript byte was edited.

## Reproduction

All {len(seal_df)} TASK-028/TASK-029 upstream seal entries verified with zero
mismatches. All 20 reconstructed full aREA module scores matched sealed scores at
maximum absolute delta {score_audit_df['full_score_max_abs_delta_vs_sealed'].max():.3g}.
The frozen seed-scheme SHA-256 is {SEED_SCHEME_SHA256}. Exact per-deletion seeds,
family sizes, fail counts, and checks are in the TSV outputs and manifest.
"""
    (OUT / "ANALYTICAL_REPORT.md").write_text(report, encoding="ascii")

    dirty = subprocess.run(["git", "-C", str(WORKTREE), "status", "--porcelain=v1"],
                           check=True, capture_output=True, text=True).stdout
    (OUT / "git_dirty_state.txt").write_text(dirty, encoding="ascii")
    branch = subprocess.run(["git", "-C", str(WORKTREE), "branch", "--show-current"],
                            check=True, capture_output=True, text=True).stdout.strip()
    head = subprocess.run(["git", "-C", str(WORKTREE), "rev-parse", "HEAD"],
                          check=True, capture_output=True, text=True).stdout.strip()
    pip_lines = subprocess.run([sys.executable, "-m", "pip", "freeze"], check=True,
                               capture_output=True, text=True).stdout.splitlines()
    env_lines = [
        f"executable={sys.executable}", f"python={sys.version.replace(chr(10), ' ')}",
        f"platform={platform.platform()}", f"numpy={np.__version__}",
        f"pandas={pd.__version__}", f"scipy={scipy.__version__}",
        f"OPENBLAS_NUM_THREADS={os.environ.get('OPENBLAS_NUM_THREADS')}",
        f"OMP_NUM_THREADS={os.environ.get('OMP_NUM_THREADS')}", "", "pip freeze:", *pip_lines,
    ]
    (OUT / "PIP_ENVIRONMENT.txt").write_text("\n".join(env_lines) + "\n", encoding="ascii")

    commands = [
        "cd data/external/original-workspace/revgate-tcga-no-purity-verify/experiments/task030_verify_current_main/gene_loo",
        "data/external/original-workspace/task028-freeze-b-draft/verifier_v3/venv_match/bin/python -m py_compile run_task030_cycle2_gene_loo.py",
        "OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 data/external/original-workspace/task028-freeze-b-draft/verifier_v3/venv_match/bin/python run_task030_cycle2_gene_loo.py",
    ]
    (OUT / "COMMANDS.txt").write_text("\n".join(commands) + "\n", encoding="ascii")

    input_paths = [
        SEED_SCHEME, T28 / "sealed_v3/B2_MODEL_SPEC.md", T28 / "FREEZE_B_DRAFT_PACKET.md",
        T28 / "sealed_v3/SEAL_MANIFEST.sha256", LEDGER, INTER / "log2tpm_v3.npy",
        INTER / "gene_names_v3.json", INTER / "gene_ids_versioned_v3.json",
        INTER / "covariates_v3.tsv", INTER / "scores_v3.npz", R28 / "RESULT_SEAL.sha256",
        R28 / "phase2b_models.json", R28 / "REPRODUCIBILITY_MANIFEST.json",
        R29 / "RESULT_SEAL.sha256", R29 / "primary_results.json",
        R29 / "REPRODUCIBILITY_MANIFEST.json", T30C1 / "critic/CRITIC_VERDICT_task030.md",
    ]
    script_path = OUT / "run_task030_cycle2_gene_loo.py"
    durable_before_manifest = sorted(
        p for p in OUT.rglob("*") if p.is_file() and "__pycache__" not in p.parts
        and p.name not in {"REPRODUCIBILITY_MANIFEST.json", "OUTPUT_SHA256SUMS.txt", "run.log"}
    )
    completed = utcnow()
    manifest = {
        "task": "TASK-030 cycle-2 frozen gene leave-one-out audit",
        "classification": "VERIFY", "status": "done",
        "started_utc": started, "completed_utc": completed,
        "git": {"branch": branch, "HEAD": head, "dirty": bool(dirty),
                "authorized_existing_dirty_sandbox": False,
                "dirty_state_path": str(OUT / "git_dirty_state.txt"),
                "dirty_state_sha256": sha256_file(OUT / "git_dirty_state.txt"),
                "no_fetch_checkout_switch_commit_push_add": True},
        "frozen_authority": {
            "B2_MODEL_SPEC_B2_12_lines": "237-246",
            "FREEZE_B_DRAFT_PACKET_lines": "219-242",
            "operational_definition": "positive survives iff every unique mapped single-gene deletion retains locked direction, abs(d)>=0.50, patient-bootstrap CI excludes 0, module F1 BH q<=0.05, planned-contrast full gated F2 BH q<=0.05, and locked adjustment",
            "failed_positive_label_rule": "uncredited/descriptive positive; explicitly not CAT5 and not a new frozen category",
        },
        "parameters": {
            "master_seed": MASTER_SEED, "seed_hash_rule": "int(sha256('20260713:' + step_id).hexdigest()[:8],16)",
            "seed_scheme_sha256": SEED_SCHEME_SHA256, "n_boot": N_BOOT, "n_perm": N_PERM,
            "effect_floor_abs_d": FLOOR, "BH_q_cutoff": Q_CUTOFF,
            "targets": [{"target": t, "contrast": c, "direction": d} for t, c, d in TARGETS],
            "configs": {"PRIMARY": {"n": 506, "covariates": CONFIGS["PRIMARY"]["cov_names"]},
                        "SENS_nopurity": {"n": 507, "covariates": CONFIGS["SENS_nopurity"]["cov_names"]}},
            "score": "sealed signed aREA; delete every admitted mapped edge for one target symbol",
            "family_substitution": "replace altered module analytic omnibus p in full 21-module F1; recompute BH/gating; rebuild F2 over gated modules using altered module C1-C3 permutation p and sealed unchanged-module p",
            "scenario": scenario,
        },
        "random_subseeds": [
            {"target": r.target, "config": r.config, "dropped_gene": r.dropped_gene,
             "bootstrap_step_id": r.bootstrap_step_id, "bootstrap_seed": int(r.bootstrap_seed),
             "permutation_step_id": r.permutation_step_id, "permutation_seed": int(r.permutation_seed),
             "stochastic_required": bool(r.target_config_all_deletions_point_pass)}
            for r in point_df.itertuples(index=False)
        ],
        "verification": {
            "upstream_seal_entries": len(seal_df), "upstream_mismatches": int((~seal_df["match"]).sum()),
            "full_scores_reconstructed": len(score_audit_df),
            "full_score_max_abs_delta_vs_sealed": float(score_audit_df["full_score_max_abs_delta_vs_sealed"].max()),
            "all_full_scores_match_1e-12": bool(score_audit_df["match_tolerance_1e-12"].all()),
            "point_deletion_rows": len(point_df), "stochastic_rows": len(stochastic_df),
            "stochastic_runs": int((stochastic_df["stochastic_status"] == "run").sum()),
            "not_run_by_gate": int((stochastic_df["stochastic_status"] != "run").sum()),
            "family_substitution_invariants_all_pass": bool(
                (stochastic_df.loc[stochastic_df.stochastic_status == "run", "family_substitution_check"]
                 == "PASS_one_F1_and_three_F2_values_substituted_if_gated").all()),
            "sealed_no_change_family_reconstruction_all_pass": bool(
                family_audit_df["no_change_reconstruction_pass"].all()),
            "sealed_no_change_family_reconstruction": family_audit_df.to_dict("records"),
        },
        "software": {"executable": sys.executable, "python": sys.version,
                     "platform": platform.platform(), "numpy": np.__version__,
                     "pandas": pd.__version__, "scipy": scipy.__version__,
                     "pip_freeze": pip_lines,
                     "thread_environment": {"OPENBLAS_NUM_THREADS": os.environ.get("OPENBLAS_NUM_THREADS"),
                                            "OMP_NUM_THREADS": os.environ.get("OMP_NUM_THREADS")}},
        "scripts": {str(script_path): sha256_file(script_path)},
        "commands": commands,
        "input_checksums": {str(p): sha256_file(p) for p in input_paths},
        "output_checksum_scope": "durable cycle-2 artifacts excluding manifest, checksum index, transient run.log, and bytecode",
        "output_checksums_excluding_manifest_checksum_index_runlog_bytecode": {
            str(p): sha256_file(p) for p in durable_before_manifest
        },
        "execution_notes": "Current-main isolated-worktree verification rerun using the corrected frozen-family implementation and unchanged seed scheme.",
        "deviations": "none",
        "preservation": "No sealed TASK-028/TASK-029 byte, cycle-1 TASK-030 byte, revgate/src byte, revgate/docs byte, DOCX, manuscript, branch, index, or commit was altered.",
    }
    manifest_path = OUT / "REPRODUCIBILITY_MANIFEST.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="ascii")

    final_files = sorted(
        p for p in OUT.rglob("*") if p.is_file() and "__pycache__" not in p.parts
        and p.name not in {"OUTPUT_SHA256SUMS.txt", "run.log"}
    )
    lines = ["# TASK-030 cycle-2 output SHA-256 checksums",
             "# Excludes self, transient run.log, and Python bytecode cache."]
    lines += [f"{sha256_file(p)}  {p.relative_to(OUT)}" for p in final_files]
    (OUT / "OUTPUT_SHA256SUMS.txt").write_text("\n".join(lines) + "\n", encoding="ascii")
    print(f"[done] scenario={scenario}", flush=True)
    print(f"[done] outputs={OUT}", flush=True)


if __name__ == "__main__":
    main()
