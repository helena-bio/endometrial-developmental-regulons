#!/usr/bin/env python3
"""TASK A post-hoc per-class decomposition of the frozen C2 contrast.

This script is deliberately self-contained orchestration around immutable TASK-028
and TASK-029 inputs. It never writes outside the requested output directory.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import platform
import shutil
import subprocess
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("PYTHONHASHSEED", "0")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy import stats

TASK_DIR = Path(__file__).resolve().parent
WORKTREE = TASK_DIR.parents[1]
SPEC_PATH = TASK_DIR / "FROZEN_POSTHOC_SPEC.json"
CHARTER_PATH = TASK_DIR / "HYPOTHESIS_AND_ANALYSIS_CHARTER.md"
T28 = Path("data/external/original-workspace/task028-freeze-b-draft")
T29 = Path("data/external/original-workspace/task029-external-replication-feasibility")
T30 = Path("data/external/original-workspace/revgate-tcga-no-purity-verify")
ORIGINAL = Path("data/external/original-workspace/revgate")
MASTER_SEED = 20260722
N_BOOT = 2000
N_PERM = 2000
TARGETS = ["GATA2", "SOX9", "HOXA9", "WT1", "PAX8", "LHX1"]
PRIMARY = {"GATA2", "SOX9"}
SUBTYPES = ["POLE", "MMRd", "NSMP", "p53abn"]
CONTRASTS = {
    "POLE_vs_NSMP": np.array([1.0, 0.0, -1.0, 0.0]),
    "MMRd_vs_NSMP": np.array([0.0, 1.0, -1.0, 0.0]),
    "POLE_vs_MMRd": np.array([1.0, -1.0, 0.0, 0.0]),
    "FROZEN_C2": np.array([0.5, 0.5, -1.0, 0.0]),
}
MODELS = [
    "TCGA_PRIMARY_CPE_N506",
    "TCGA_NOPURITY_N507",
    "CPTAC_DISCOVERY_N95",
    "CPTAC_CONFIRMATORY_N135",
]
META_MODEL = "CPTAC_FIXED_EFFECT_META"
FAMILY_BY_MODEL = {
    "TCGA_PRIMARY_CPE_N506": "BH18_TCGA_PRIMARY",
    "TCGA_NOPURITY_N507": "BH18_TCGA_NOPURITY",
    "CPTAC_DISCOVERY_N95": "BH18_CPTAC_DISCOVERY",
    "CPTAC_CONFIRMATORY_N135": "BH18_CPTAC_CONFIRMATORY",
    META_MODEL: "BH18_CPTAC_FIXED_EFFECT_META",
}
EXPECTED_CHARTER_SHA = "0dbd83c4dd69cca9fc70eecb9a8012bf6e91bfb9b51fcc0abd10f36c84d9851a"
EXPECTED_SPEC_SHA = "044bf4c5cdb80df5825e91cf5c6c9a29b17d2e75dac4f8d67f002f1655534340"
EXPECTED_HEAD = "83503bad47b60193598b2b9ebe819c22c83e8ac1"
EXPECTED_BRANCH = "experiment/task-a-perclass-c2-sensitivity"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def run_cmd(args: list[str], cwd: Path | None = None) -> str:
    return subprocess.run(args, cwd=cwd, text=True, capture_output=True,
                          check=True).stdout


def git_snapshot(repo: Path) -> dict:
    status = run_cmd(["git", "status", "--porcelain=v1"], repo)
    diff = subprocess.run(["git", "diff", "--binary", "HEAD"], cwd=repo,
                          capture_output=True, check=True).stdout
    cached = subprocess.run(["git", "diff", "--cached", "--binary", "HEAD"],
                            cwd=repo, capture_output=True, check=True).stdout
    return {
        "path": str(repo),
        "head": run_cmd(["git", "rev-parse", "HEAD"], repo).strip(),
        "branch": run_cmd(["git", "branch", "--show-current"], repo).strip(),
        "status_porcelain": status.splitlines(),
        "status_sha256": sha256_text(status),
        "unstaged_binary_diff_sha256": hashlib.sha256(diff).hexdigest(),
        "staged_binary_diff_sha256": hashlib.sha256(cached).hexdigest(),
    }


def subseed(step_id: str) -> int:
    return int(hashlib.sha256(f"{MASTER_SEED}:{step_id}".encode()).hexdigest()[:8], 16)


def clean_json(x):
    if isinstance(x, dict):
        return {str(k): clean_json(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [clean_json(v) for v in x]
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.floating,)):
        x = float(x)
    if isinstance(x, (np.bool_,)):
        return bool(x)
    if isinstance(x, float) and not math.isfinite(x):
        return None
    return x


def write_json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(clean_json(obj), indent=2, sort_keys=True,
                               allow_nan=False) + "\n", encoding="ascii")


def write_df(path: Path, df: pd.DataFrame, sep: str = "\t") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep=sep, index=False, na_rep="NA", lineterminator="\n",
              quoting=csv.QUOTE_MINIMAL)


def verify_standard_manifest(manifest: Path, base: Path) -> dict:
    rows, bad = [], []
    for line in manifest.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#"):
            continue
        digest, rel = line.split(None, 1)
        path = (base / rel.strip()).resolve()
        got = sha256_file(path) if path.is_file() else "MISSING"
        ok = got == digest
        rows.append({"manifest": str(manifest), "path": str(path),
                     "expected_sha256": digest, "observed_sha256": got, "match": ok})
        if not ok:
            bad.append(rows[-1])
    if bad:
        raise RuntimeError(f"BLOCKED_DATA_VERSION_MISMATCH: {bad[0]}")
    return {"manifest": str(manifest), "entries": len(rows), "bad": 0, "rows": rows}


def verify_tcga_design_manifest() -> dict:
    manifest = T28 / "sealed_v3/SEAL_MANIFEST.sha256"
    rows, bad = [], []
    for line in manifest.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#"):
            continue
        first, name = line.split(None, 1)
        name = name.strip()
        if first.startswith("B1/"):
            digest = first.split("/", 1)[1]
            path = T28 / "sealed_v3/B1" / name
        else:
            digest = first
            path = T28 / "sealed_v3" / name
        got = sha256_file(path) if path.is_file() else "MISSING"
        ok = got == digest
        rows.append({"manifest": str(manifest), "path": str(path),
                     "expected_sha256": digest, "observed_sha256": got, "match": ok})
        if not ok:
            bad.append(rows[-1])
    if bad:
        raise RuntimeError(f"BLOCKED_DATA_VERSION_MISMATCH: {bad[0]}")
    return {"manifest": str(manifest), "entries": len(rows), "bad": 0, "rows": rows}


def provenance_gate(out: Path) -> dict:
    spec = json.loads(SPEC_PATH.read_text(encoding="utf-8"))
    if sha256_file(CHARTER_PATH) != EXPECTED_CHARTER_SHA:
        raise RuntimeError("BLOCKED_DATA_VERSION_MISMATCH: charter SHA256")
    if sha256_file(SPEC_PATH) != EXPECTED_SPEC_SHA:
        raise RuntimeError("BLOCKED_DATA_VERSION_MISMATCH: spec SHA256")
    head = run_cmd(["git", "rev-parse", "HEAD"], WORKTREE).strip()
    branch = run_cmd(["git", "branch", "--show-current"], WORKTREE).strip()
    if head != EXPECTED_HEAD or branch != EXPECTED_BRANCH:
        raise RuntimeError(f"BLOCKED_GIT_GUARD: head={head} branch={branch}")
    required = []
    for raw, expected in spec["required_input_sha256"].items():
        path = Path(raw)
        got = sha256_file(path) if path.is_file() else "MISSING"
        row = {"path": raw, "expected_sha256": expected,
               "observed_sha256": got, "match": got == expected}
        required.append(row)
        if not row["match"]:
            raise RuntimeError(f"BLOCKED_DATA_VERSION_MISMATCH: {row}")
    verdict = T30 / "experiments/task030_verify_current_main/critic_cycle2_targeted/CRITIC_VERDICT_TARGETED_CYCLE2.md"
    if verdict.read_text(encoding="utf-8").splitlines()[0].strip() != "SURVIVES":
        raise RuntimeError("BLOCKED_AUTHORITY: TASK-030 first line is not SURVIVES")
    seals = [verify_tcga_design_manifest()]
    seals.append(verify_standard_manifest(T28 / "execution/results_v3/RESULT_SEAL.sha256",
                                          T28 / "execution/results_v3"))
    seals.append(verify_standard_manifest(T29 / "SEAL_MANIFEST_replication.sha256", T29))
    seals.append(verify_standard_manifest(T29 / "execution/intermediate/INTERMEDIATE_SEAL.sha256",
                                          T29 / "execution/intermediate"))
    seals.append(verify_standard_manifest(T29 / "execution/results/RESULT_SEAL.sha256",
                                          T29 / "execution/results"))
    roster = T29 / "sources_task029gate/dou2023_mmc2.xlsx"
    if roster.stat().st_size != 277385:
        raise RuntimeError(f"BLOCKED_DATA_VERSION_MISMATCH: roster size={roster.stat().st_size}")
    snapshots = {"original_dirty_worktree": git_snapshot(ORIGINAL),
                 "task030_verified_worktree": git_snapshot(T30)}
    report = {"charter_sha256": sha256_file(CHARTER_PATH),
              "spec_sha256": sha256_file(SPEC_PATH), "required_inputs": required,
              "nested_seals": seals, "start_snapshots": snapshots,
              "worktree_head": head, "worktree_branch": branch}
    write_json(out / "integrity/PROVENANCE_GATE.json", report)
    all_rows = required + [r for s in seals for r in s["rows"]]
    write_df(out / "INPUT_CHECKSUMS.tsv", pd.DataFrame(all_rows))
    return report


def build_design(subtype: np.ndarray, cov: dict[str, np.ndarray], names: list[str]) -> np.ndarray:
    n = len(subtype)
    dummies = np.column_stack([(subtype == s).astype(float)
                               for s in ["MMRd", "NSMP", "p53abn"]])
    return np.column_stack([np.ones(n), dummies] + [np.asarray(cov[nm], float) for nm in names])


def emm_matrix(cov: dict[str, np.ndarray], names: list[str]) -> np.ndarray:
    ref = [float(np.mean(cov[nm])) for nm in names]
    rows = []
    for s in SUBTYPES:
        rows.append([1.0, float(s == "MMRd"), float(s == "NSMP"),
                     float(s == "p53abn"), *ref])
    return np.asarray(rows, float)


def fit_observed(y: np.ndarray, subtype: np.ndarray, cov: dict[str, np.ndarray],
                 names: list[str]) -> dict:
    X = build_design(subtype, cov, names)
    n, p = X.shape
    beta, _, rank, _ = np.linalg.lstsq(X, y, rcond=None)
    if rank != p:
        raise np.linalg.LinAlgError(f"rank deficient observed design rank={rank} p={p}")
    resid = y - X @ beta
    df = n - p
    sigma = math.sqrt(float(resid @ resid) / df)
    if not math.isfinite(sigma) or sigma <= 0:
        raise ValueError(f"invalid residual SD {sigma}")
    cov_beta = (sigma ** 2) * np.linalg.inv(X.T @ X)
    E = emm_matrix(cov, names)
    mu = E @ beta
    tcrit = float(stats.t.ppf(0.975, df))
    out = {"n": n, "p": p, "df": df, "sigma": sigma, "mu": mu, "contrasts": {}}
    for name, w in CONTRASTS.items():
        l = w @ E
        b = float(l @ beta)
        se = math.sqrt(float(l @ cov_beta @ l))
        out["contrasts"][name] = {"b": b, "se": se, "d": b / sigma,
                                    "t_lo": b - tcrit * se, "t_hi": b + tcrit * se,
                                    "wald_p": float(2 * stats.t.sf(abs(b / se), df))}
    return out


def fit_replicate(y: np.ndarray, subtype: np.ndarray, cov: dict[str, np.ndarray],
                  names: list[str]) -> tuple[float, dict[str, tuple[float, float]]]:
    X = build_design(subtype, cov, names)
    n, p = X.shape
    beta, _, rank, _ = np.linalg.lstsq(X, y, rcond=None)
    if rank != p:
        raise np.linalg.LinAlgError("rank deficient replicate")
    resid = y - X @ beta
    sigma = math.sqrt(float(resid @ resid) / (n - p))
    if not math.isfinite(sigma) or sigma <= 0:
        raise ValueError("invalid replicate residual SD")
    mu = emm_matrix(cov, names) @ beta
    return sigma, {name: (float(w @ mu), float((w @ mu) / sigma))
                   for name, w in CONTRASTS.items()}


def analyze_one(y: np.ndarray, subtype: np.ndarray, cov: dict[str, np.ndarray],
                cov_names: list[str], model: str, target: str) -> dict:
    observed = fit_observed(y, subtype, cov, cov_names)
    n = len(y)
    rng = np.random.default_rng(subseed(f"boot__{model}__{target}"))
    store = {c: {"b": [], "d": []} for c in CONTRASTS}
    failures = 0
    max_boot_b_scaled = 0.0
    max_boot_d_scaled = 0.0
    for _ in range(N_BOOT):
        ix = rng.integers(0, n, n)
        st = subtype[ix]
        if len(set(st.tolist())) != 4:
            failures += 1
            continue
        try:
            _, vals = fit_replicate(y[ix], st, {k: v[ix] for k, v in cov.items()}, cov_names)
        except (np.linalg.LinAlgError, ValueError, FloatingPointError):
            failures += 1
            continue
        for c, (b, d) in vals.items():
            store[c]["b"].append(b); store[c]["d"].append(d)
        bpn, bmn, bc2 = vals["POLE_vs_NSMP"][0], vals["MMRd_vs_NSMP"][0], vals["FROZEN_C2"][0]
        dpn, dmn, dc2 = vals["POLE_vs_NSMP"][1], vals["MMRd_vs_NSMP"][1], vals["FROZEN_C2"][1]
        max_boot_b_scaled = max(max_boot_b_scaled, abs(bc2 - .5*bpn - .5*bmn) /
                                max(1.0, abs(bc2), .5*abs(bpn)+.5*abs(bmn)))
        max_boot_d_scaled = max(max_boot_d_scaled, abs(dc2 - .5*dpn - .5*dmn) /
                                max(1.0, abs(dc2), .5*abs(dpn)+.5*abs(dmn)))
    if any(len(store[c]["d"]) < 2 for c in CONTRASTS):
        raise RuntimeError(f"NOT_EVALUABLE bootstrap {model} {target}")
    for c in CONTRASTS:
        ba = np.asarray(store[c]["b"]); da = np.asarray(store[c]["d"])
        observed["contrasts"][c].update({
            "b_boot_lo": float(np.percentile(ba, 2.5)),
            "b_boot_hi": float(np.percentile(ba, 97.5)),
            "d_boot_lo": float(np.percentile(da, 2.5)),
            "d_boot_hi": float(np.percentile(da, 97.5)),
            "d_boot_se": float(np.std(da, ddof=1)),
        })
    rng = np.random.default_rng(subseed(f"perm__{model}__{target}"))
    ge = {c: 1 for c in ["POLE_vs_NSMP", "MMRd_vs_NSMP", "POLE_vs_MMRd"]}
    obs_abs = {c: abs(observed["contrasts"][c]["d"]) for c in ge}
    perm_fail = 0
    for _ in range(N_PERM):
        st = rng.permutation(subtype)
        try:
            _, vals = fit_replicate(y, st, cov, cov_names)
        except (np.linalg.LinAlgError, ValueError, FloatingPointError):
            perm_fail += 1
            continue
        for c in ge:
            if abs(vals[c][1]) >= obs_abs[c] - 1e-12:
                ge[c] += 1
    if perm_fail:
        raise RuntimeError(f"NOT_EVALUABLE permutation failures={perm_fail} {model} {target}")
    for c in ge:
        observed["contrasts"][c]["perm_p"] = ge[c] / (N_PERM + 1)
    observed["n_boot_failed"] = failures
    observed["n_boot_used"] = N_BOOT - failures
    observed["n_perm"] = N_PERM
    observed["boot_identity_max_b_scaled_error"] = max_boot_b_scaled
    observed["boot_identity_max_d_scaled_error"] = max_boot_d_scaled
    return observed


def read_regulons() -> dict:
    path = T28 / "sealed_v3/B1/m3_primary_edge_ledger.tsv"
    df = pd.read_csv(path, sep="\t", comment="#")
    df = df[df["primary_inclusion"] == "Y"]
    out = {}
    for tf in TARGETS:
        x = df[df["TF"] == tf]
        out[tf] = {"targets": x["target"].tolist(),
                   "weights": x["weight"].astype(float).to_numpy(),
                   "mor": x["primary_sign"].astype(float).to_numpy()}
    return out


def read_module(path: Path) -> list[str]:
    return [x.strip() for x in path.read_text().splitlines()
            if x.strip() and not x.startswith("#")]


def score_cptac_stratum(stratum: str, regulons: dict) -> dict:
    sys.path.insert(0, str(T28 / "execution/scripts_v3"))
    import scoring as sealed_scoring
    inter = T29 / "execution/intermediate"
    expr = np.load(inter / f"log2tpm_{stratum}.npy")
    gene_names = json.loads((inter / "gene_names.json").read_text())
    case_order = json.loads((inter / "case_order_by_stratum.json").read_text())[stratum]
    ledger = json.loads((inter / "acq02_join_ledger.json").read_text())
    pick = ledger["pick"]
    sym_rows = defaultdict(list)
    for i, name in enumerate(gene_names):
        sym_rows[name].append(i)
    means = expr.mean(axis=1)
    def best_row(symbol):
        rows = sym_rows.get(symbol, [])
        if not rows:
            return None
        return rows[int(np.argmax(means[rows]))]
    allzero = {symbol for symbol, rows in sym_rows.items()
               if all(bool((expr[r] == 0).all()) for r in rows)}
    b1 = T28 / "sealed_v3/B1"
    estimate = {}
    for line in (T28 / "experimenter_final/sources/ESTIMATE_SI_geneset.gmt").read_text().splitlines():
        fields = line.split("\t"); estimate[fields[0]] = fields[2:]
    modules = {"M4": read_module(b1 / "M4_covariate.txt"),
               "stromal": estimate["StromalSignature"],
               "immune": estimate["ImmuneSignature"]}
    sets = {}
    for name, genes in modules.items():
        rows, seen = [], set()
        for g in genes:
            if g in allzero:
                continue
            r = best_row(g)
            if r is not None and r not in seen:
                rows.append(r); seen.add(r)
        sets[name] = rows
    cov_scores = sealed_scoring.ssgsea_matrix(expr, sets)
    # Reproduce sealed aREA once per stratum while sharing rank transforms across targets.
    z = (expr - expr.mean(axis=1, keepdims=True)) / expr.std(axis=1, ddof=1, keepdims=True)
    finite = np.all(np.isfinite(z), axis=1)
    reduced_pos = -np.ones(expr.shape[0], dtype=np.int64)
    reduced_pos[finite] = np.arange(int(finite.sum()))
    zred = z[finite]
    prepared = {}
    for tf, reg in regulons.items():
        rows, weights, mor = [], [], []
        for target, weight, sign in zip(reg["targets"], reg["weights"], reg["mor"]):
            r = best_row(target)
            if r is not None and finite[r]:
                rows.append(reduced_pos[r]); weights.append(weight); mor.append(sign)
        weights = np.asarray(weights, float); wts = weights / weights.max()
        prepared[tf] = (np.asarray(rows, int), wts / wts.sum(), np.asarray(mor, float),
                        math.sqrt(float(np.sum(wts ** 2))))
    scores = {tf: np.zeros(expr.shape[1]) for tf in TARGETS}
    nfin = zred.shape[0]
    for j in range(expr.shape[1]):
        sig = zred[:, j]
        t2all = stats.norm.ppf(stats.rankdata(sig, method="average") / (nfin + 1.0))
        t1all = stats.norm.ppf((stats.rankdata(np.abs(sig), method="average") /
                               (nfin + 1.0) + 1.0) / 2.0)
        for tf, (rows, wtss, mor, normer) in prepared.items():
            s2 = float(np.sum(wtss * mor * t2all[rows]))
            s1 = float(np.sum(wtss * t1all[rows]))
            scores[tf][j] = (abs(s2) + (s1 if s1 > 0 else 0.0)) * np.sign(s2) * normer
    return {"subtype": np.asarray([pick[c]["subtype"] for c in case_order]),
            "cov": {"M4": np.asarray(cov_scores["M4"]),
                    "composition": np.asarray(cov_scores["stromal"] + cov_scores["immune"])},
            "scores": scores, "cases": case_order,
            "targets_used": {tf: len(prepared[tf][0]) for tf in TARGETS}}


def load_all_data() -> dict:
    inter = T28 / "execution/intermediate"
    npz = np.load(inter / "scores_v3.npz", allow_pickle=True)
    covdf = pd.read_csv(inter / "covariates_v3.tsv", sep="\t")
    if list(npz["patient_order"]) != list(covdf["patient_barcode"]):
        raise RuntimeError("BLOCKED_DATA_VERSION_MISMATCH: TCGA patient order")
    scores = {tf: np.asarray(npz[f"M3primary__{tf}"], float) for tf in TARGETS}
    subtype = covdf["subtype"].to_numpy()
    cov = {"M4": covdf["M4_prolif"].to_numpy(float),
           "Aran_CPE": covdf["purity_CPE"].to_numpy(float),
           "composition": covdf["composition"].to_numpy(float)}
    mask = covdf["cpe_complete_case"].to_numpy() == 1
    data = {
        "TCGA_PRIMARY_CPE_N506": {"subtype": subtype[mask],
            "cov": {k: v[mask] for k, v in cov.items()},
            "cov_names": ["M4", "Aran_CPE", "composition"],
            "scores": {tf: v[mask] for tf, v in scores.items()}, "cohort": "TCGA", "stratum": "PRIMARY_CPE"},
        "TCGA_NOPURITY_N507": {"subtype": subtype,
            "cov": {"M4": cov["M4"], "composition": cov["composition"]},
            "cov_names": ["M4", "composition"], "scores": scores,
            "cohort": "TCGA", "stratum": "NOPURITY"},
    }
    regs = read_regulons()
    for stratum, model in [("Discovery", "CPTAC_DISCOVERY_N95"),
                           ("Confirmatory", "CPTAC_CONFIRMATORY_N135")]:
        x = score_cptac_stratum(stratum, regs)
        x.update({"cov_names": ["M4", "composition"], "cohort": "CPTAC", "stratum": stratum.upper()})
        data[model] = x
    return data


def check_counts(data: dict, spec: dict) -> list[dict]:
    rows = []
    for model in MODELS:
        observed = {s: int(np.sum(data[model]["subtype"] == s)) for s in SUBTYPES}
        observed["total"] = len(data[model]["subtype"])
        expected = spec["models"][model]["expected_counts"]
        row = {"model": model, **{f"observed_{k}": observed[k] for k in ["total", *SUBTYPES]},
               **{f"expected_{k}": expected[k] for k in ["total", *SUBTYPES]}}
        row["match"] = all(observed[k] == expected[k] for k in expected)
        rows.append(row)
        if not row["match"]:
            raise RuntimeError(f"BLOCKED_DATA_VERSION_MISMATCH subtype counts: {row}")
    return rows


def taxonomy(vals: dict, meta: bool = False) -> str:
    pn, mn, pm = [vals[x] for x in ["POLE_vs_NSMP", "MMRd_vs_NSMP", "POLE_vs_MMRd"]]
    bpn, bmn = pn["value"], mn["value"]
    prefix = "CPTAC_META_" if meta else ""
    if bpn == 0 or bmn == 0:
        return prefix + "ZERO_BOUNDARY_UNRESOLVED"
    pm_excludes = pm["lo"] > 0 or pm["hi"] < 0
    pn_excludes = pn["lo"] > 0 or pn["hi"] < 0
    mn_excludes = mn["lo"] > 0 or mn["hi"] < 0
    if bpn * bmn < 0:
        if pm_excludes and pn_excludes and mn_excludes:
            return prefix + "OPPOSITE_DIRECTION_HETEROGENEITY_SUPPORTED"
        return prefix + "OPPOSITE_DIRECTION_POINT_HETEROGENEITY_UNRESOLVED"
    if pm_excludes and abs(pm["d_pm"]) >= 0.5:
        if abs(bpn) > abs(bmn):
            return prefix + "SAME_DIRECTION_MATERIALLY_DIFFERENT_POLE_DOMINANT"
        if abs(bmn) > abs(bpn):
            return prefix + "SAME_DIRECTION_MATERIALLY_DIFFERENT_MMRD_DOMINANT"
        return prefix + "SAME_DIRECTION_MATERIALLY_DIFFERENT_TIE_ERROR"
    if pm_excludes:
        return prefix + "SAME_DIRECTION_DISTINGUISHABLE_BELOW_MATERIALITY_FLOOR"
    if pn_excludes and mn_excludes:
        return prefix + "SAME_DIRECTION_COMPATIBLE_BOTH_RESOLVED"
    return prefix + "SAME_DIRECTION_POINT_ESTIMATES_UNRESOLVED"


def bh18(pvals: list[float | None]) -> list[float | None]:
    work = np.asarray([1.0 if p is None or not math.isfinite(p) else p for p in pvals], float)
    order = np.argsort(work, kind="mergesort")
    ranked = work[order]
    qrank = np.minimum.accumulate((ranked * 18 / np.arange(1, 19))[::-1])[::-1]
    q = np.empty(18); q[order] = np.minimum(qrank, 1.0)
    return [None if p is None or not math.isfinite(p) else float(q[i]) for i, p in enumerate(pvals)]


def model_formula(model: str) -> str:
    if model == "TCGA_PRIMARY_CPE_N506":
        return "score ~ subtype + M4 + Aran_CPE + ESTIMATE_composition"
    if model == META_MODEL:
        return "CPTAC inverse-bootstrap-variance fixed-effect meta on standardized d"
    return "score ~ subtype + M4 + ESTIMATE_composition"


def target_role(target: str) -> str:
    if target in PRIMARY:
        return "PRIMARY_EXPLANATORY"
    if target in {"PAX8", "LHX1"}:
        return "SECONDARY_COMPLETENESS_PRIOR_C1_TARGET"
    return "SECONDARY_COMPLETENESS"


def execute(out: Path) -> None:
    out.mkdir(parents=True, exist_ok=False)
    prov = provenance_gate(out)
    spec = json.loads(SPEC_PATH.read_text())
    data = load_all_data()
    count_rows = check_counts(data, spec)
    write_df(out / "results/SAMPLE_COUNTS.tsv", pd.DataFrame(count_rows))
    fits = {}
    for model in MODELS:
        md = data[model]
        for target in TARGETS:
            print(f"FIT {model} {target}", flush=True)
            fits[(model, target)] = analyze_one(md["scores"][target], md["subtype"],
                                                md["cov"], md["cov_names"], model, target)
    direct_states = {}
    recon_rows = []
    for model in MODELS:
        for target in TARGETS:
            fit = fits[(model, target)]; cc = fit["contrasts"]
            vals = {c: {"value": cc[c]["b"], "lo": cc[c]["t_lo"], "hi": cc[c]["t_hi"],
                        "d_pm": cc["POLE_vs_MMRd"]["d"]} for c in CONTRASTS if c != "FROZEN_C2"}
            direct_states[(model, target)] = taxonomy(vals)
            b_re = .5*cc["POLE_vs_NSMP"]["b"] + .5*cc["MMRd_vs_NSMP"]["b"]
            d_re = .5*cc["POLE_vs_NSMP"]["d"] + .5*cc["MMRd_vs_NSMP"]["d"]
            b = cc["FROZEN_C2"]["b"]; d = cc["FROZEN_C2"]["d"]
            eb, ed = abs(b-b_re), abs(d-d_re)
            sb = max(1, abs(b), .5*abs(cc["POLE_vs_NSMP"]["b"])+.5*abs(cc["MMRd_vs_NSMP"]["b"]))
            sd = max(1, abs(d), .5*abs(cc["POLE_vs_NSMP"]["d"])+.5*abs(cc["MMRd_vs_NSMP"]["d"]))
            row = {"target": target, "model": model, "b_C2_direct": b, "b_C2_reconstructed": b_re,
                   "c2_reconstruction_abs_error": eb, "c2_scaled_error": eb/sb,
                   "d_C2_direct": d, "d_C2_reconstructed": d_re,
                   "d_c2_reconstruction_abs_error": ed, "d_c2_scaled_error": ed/sd,
                   "bootstrap_max_b_scaled_error": fit["boot_identity_max_b_scaled_error"],
                   "bootstrap_max_d_scaled_error": fit["boot_identity_max_d_scaled_error"],
                   "tolerance": 1e-12}
            row["pass"] = all(row[k] <= 1e-12 for k in ["c2_scaled_error", "d_c2_scaled_error",
                                                        "bootstrap_max_b_scaled_error", "bootstrap_max_d_scaled_error"])
            if not row["pass"]:
                raise RuntimeError(f"BLOCKING_C2_IDENTITY_FAILURE {row}")
            recon_rows.append(row)
    # Fixed-effect meta for all four contrasts, requiring both strata.
    meta_rows, meta_lookup = [], {}
    for target in TARGETS:
        for contrast in CONTRASTS:
            a = fits[("CPTAC_DISCOVERY_N95", target)]["contrasts"][contrast]
            b = fits[("CPTAC_CONFIRMATORY_N135", target)]["contrasts"][contrast]
            ses = np.asarray([a["d_boot_se"], b["d_boot_se"]], float)
            if np.any(~np.isfinite(ses)) or np.any(ses <= 0):
                raise RuntimeError(f"NOT_EVALUABLE meta SE {target} {contrast} {ses}")
            w = 1 / ses**2; ds = np.asarray([a["d"], b["d"]])
            dfe = float(np.sum(w*ds)/np.sum(w)); sefe = float(math.sqrt(1/np.sum(w)))
            lo, hi = dfe-1.96*sefe, dfe+1.96*sefe
            p = float(2*stats.norm.sf(abs(dfe/sefe)))
            row = {"target": target, "contrast": contrast, "d_FE": dfe, "SE_FE": sefe,
                   "CI_lo": lo, "CI_hi": hi, "p_raw": p,
                   "discovery_coefficient": a["b"], "discovery_coefficient_SE": a["se"],
                   "discovery_residual_SD": fits[("CPTAC_DISCOVERY_N95", target)]["sigma"],
                   "discovery_d": a["d"], "discovery_d_boot_CI_lo": a["d_boot_lo"],
                   "discovery_d_boot_CI_hi": a["d_boot_hi"], "discovery_boot_SE_d": a["d_boot_se"],
                   "confirmatory_coefficient": b["b"], "confirmatory_coefficient_SE": b["se"],
                   "confirmatory_residual_SD": fits[("CPTAC_CONFIRMATORY_N135", target)]["sigma"],
                   "confirmatory_d": b["d"], "confirmatory_d_boot_CI_lo": b["d_boot_lo"],
                   "confirmatory_d_boot_CI_hi": b["d_boot_hi"], "confirmatory_boot_SE_d": b["d_boot_se"],
                   "weight_discovery": float(w[0]), "weight_confirmatory": float(w[1]),
                   "weight_fraction_discovery": float(w[0]/w.sum()),
                   "weight_fraction_confirmatory": float(w[1]/w.sum()),
                   "opposite_direction_flag": bool(ds[0]*ds[1] < 0)}
            meta_rows.append(row); meta_lookup[(target, contrast)] = row
    meta_states = {}
    for target in TARGETS:
        vals = {c: {"value": meta_lookup[(target,c)]["d_FE"],
                    "lo": meta_lookup[(target,c)]["CI_lo"], "hi": meta_lookup[(target,c)]["CI_hi"],
                    "d_pm": meta_lookup[(target,"POLE_vs_MMRd")]["d_FE"]}
                for c in ["POLE_vs_NSMP", "MMRd_vs_NSMP", "POLE_vs_MMRd"]}
        meta_states[target] = taxonomy(vals, meta=True)
    for target in TARGETS:
        c2 = meta_lookup[(target,"FROZEN_C2")]["d_FE"]
        re = .5*meta_lookup[(target,"POLE_vs_NSMP")]["d_FE"] + .5*meta_lookup[(target,"MMRd_vs_NSMP")]["d_FE"]
        for c in CONTRASTS:
            meta_lookup[(target,c)]["meta_C2_direct_minus_reconstruction"] = c2-re
            meta_lookup[(target,c)]["interpretation_state"] = meta_states[target]
    # Build long direct + three meta rows per target.
    long_rows = []
    for model in MODELS:
        md = data[model]
        counts = {s: int(np.sum(md["subtype"] == s)) for s in SUBTYPES}
        for target in TARGETS:
            fit = fits[(model,target)]; c2 = fit["contrasts"]["FROZEN_C2"]
            for contrast in ["POLE_vs_NSMP","MMRd_vs_NSMP","POLE_vs_MMRd"]:
                x = fit["contrasts"][contrast]
                long_rows.append({"task":"TASK_A_PERCLASS_C2", "target":target,
                    "target_role":target_role(target), "model":model, "cohort":md["cohort"],
                    "stratum":md["stratum"], "contrast":contrast, "model_formula":model_formula(model),
                    "covariates":"|".join(md["cov_names"]), "analytic_n":fit["n"], "residual_df":fit["df"],
                    **{f"n_{s}":counts[s] for s in SUBTYPES},
                    "contrast_weights":json.dumps(CONTRASTS[contrast].tolist(),separators=(",",":")),
                    "coefficient":x["b"], "coefficient_SE":x["se"], "residual_SD":fit["sigma"], "d":x["d"],
                    "coefficient_t_CI_lo":x["t_lo"], "coefficient_t_CI_hi":x["t_hi"],
                    "coefficient_boot_CI_lo":x["b_boot_lo"], "coefficient_boot_CI_hi":x["b_boot_hi"],
                    "d_boot_CI_lo":x["d_boot_lo"], "d_boot_CI_hi":x["d_boot_hi"], "d_boot_SE":x["d_boot_se"],
                    "perm_p_raw":x["perm_p"], "BH18_q_descriptive":np.nan,
                    "wald_t_p_diagnostic":x["wald_p"], "point_direction":"POSITIVE" if x["b"]>0 else "NEGATIVE" if x["b"]<0 else "ZERO",
                    "n_boot_requested":N_BOOT, "n_boot_used":fit["n_boot_used"], "n_boot_failed":fit["n_boot_failed"],
                    "n_perm":N_PERM, "deterministic_seed":subseed(f"boot__{model}__{target}"),
                    "permutation_seed":subseed(f"perm__{model}__{target}"), "evaluability_status":"EVALUABLE", "evaluability_reason":"",
                    "c2_reconstruction_value":c2["b"],
                    "c2_reconstruction_abs_error":next(r["c2_reconstruction_abs_error"] for r in recon_rows if r["model"]==model and r["target"]==target),
                    "d_c2_reconstruction_abs_error":next(r["d_c2_reconstruction_abs_error"] for r in recon_rows if r["model"]==model and r["target"]==target),
                    "interpretation_state":direct_states[(model,target)],
                    "interpretation_basis":"same-fit model t CIs; direct PM; inherited |d_PM|>=0.50"})
    for target in TARGETS:
        for contrast in ["POLE_vs_NSMP","MMRd_vs_NSMP","POLE_vs_MMRd"]:
            m=meta_lookup[(target,contrast)]
            long_rows.append({"task":"TASK_A_PERCLASS_C2","target":target,"target_role":target_role(target),
                "model":META_MODEL,"cohort":"CPTAC","stratum":"DISCOVERY_PLUS_CONFIRMATORY_FE", "contrast":contrast,
                "model_formula":model_formula(META_MODEL),"covariates":"stratum-specific M4|composition", "analytic_n":230,
                "residual_df":np.nan,"n_POLE":13,"n_MMRd":72,"n_NSMP":109,"n_p53abn":36,
                "contrast_weights":json.dumps(CONTRASTS[contrast].tolist(),separators=(",",":")),
                "coefficient":np.nan,"coefficient_SE":np.nan,"residual_SD":np.nan,"d":m["d_FE"],
                "coefficient_t_CI_lo":np.nan,"coefficient_t_CI_hi":np.nan,"coefficient_boot_CI_lo":np.nan,
                "coefficient_boot_CI_hi":np.nan,"d_boot_CI_lo":m["CI_lo"],"d_boot_CI_hi":m["CI_hi"],
                "d_boot_SE":m["SE_FE"],"perm_p_raw":m["p_raw"],"BH18_q_descriptive":np.nan,
                "wald_t_p_diagnostic":np.nan,"point_direction":"POSITIVE" if m["d_FE"]>0 else "NEGATIVE" if m["d_FE"]<0 else "ZERO",
                "n_boot_requested":N_BOOT,"n_boot_used":min(fits[("CPTAC_DISCOVERY_N95",target)]["n_boot_used"],fits[("CPTAC_CONFIRMATORY_N135",target)]["n_boot_used"]),
                "n_boot_failed":fits[("CPTAC_DISCOVERY_N95",target)]["n_boot_failed"]+fits[("CPTAC_CONFIRMATORY_N135",target)]["n_boot_failed"],
                "n_perm":0,"deterministic_seed":np.nan,"permutation_seed":np.nan,"evaluability_status":"EVALUABLE","evaluability_reason":"",
                "c2_reconstruction_value":meta_lookup[(target,"FROZEN_C2")]["d_FE"],
                "c2_reconstruction_abs_error":abs(meta_lookup[(target,"FROZEN_C2")]["d_FE"]-(.5*meta_lookup[(target,"POLE_vs_NSMP")]["d_FE"]+.5*meta_lookup[(target,"MMRd_vs_NSMP")]["d_FE"])),
                "d_c2_reconstruction_abs_error":abs(meta_lookup[(target,"FROZEN_C2")]["d_FE"]-(.5*meta_lookup[(target,"POLE_vs_NSMP")]["d_FE"]+.5*meta_lookup[(target,"MMRd_vs_NSMP")]["d_FE"])),
                "interpretation_state":meta_states[target],"interpretation_basis":"contrast-specific CPTAC FE d and normal CIs; direct FE PM"})
            long_rows[-1].update({
                "meta_discovery_coefficient": m["discovery_coefficient"],
                "meta_discovery_coefficient_SE": m["discovery_coefficient_SE"],
                "meta_discovery_residual_SD": m["discovery_residual_SD"],
                "meta_discovery_d": m["discovery_d"],
                "meta_discovery_boot_SE_d": m["discovery_boot_SE_d"],
                "meta_discovery_weight": m["weight_discovery"],
                "meta_confirmatory_coefficient": m["confirmatory_coefficient"],
                "meta_confirmatory_coefficient_SE": m["confirmatory_coefficient_SE"],
                "meta_confirmatory_residual_SD": m["confirmatory_residual_SD"],
                "meta_confirmatory_d": m["confirmatory_d"],
                "meta_confirmatory_boot_SE_d": m["confirmatory_boot_SE_d"],
                "meta_confirmatory_weight": m["weight_confirmatory"],
                "meta_weight_fraction_discovery": m["weight_fraction_discovery"],
                "meta_weight_fraction_confirmatory": m["weight_fraction_confirmatory"],
                "meta_opposite_direction_flag": m["opposite_direction_flag"],
            })
    longdf=pd.DataFrame(long_rows)
    family_rows=[]
    for model in [*MODELS,META_MODEL]:
        idx=list(longdf.index[longdf["model"]==model]); assert len(idx)==18
        ps=[float(longdf.at[i,"perm_p_raw"]) if pd.notna(longdf.at[i,"perm_p_raw"]) else None for i in idx]
        qs=bh18(ps)
        for i,q in zip(idx,qs): longdf.at[i,"BH18_q_descriptive"]=q
        family_rows.append({"family":FAMILY_BY_MODEL[model],"model":model,"planned_size":18,
                            "evaluable_count":sum(p is not None for p in ps),"missing_placeholders_inserted":sum(p is None for p in ps),
                            "raw_p_source":"fixed-effect normal p" if model==META_MODEL else "permutation p"})
    for row in meta_rows:
        if row["contrast"] != "FROZEN_C2":
            q=longdf[(longdf.model==META_MODEL)&(longdf.target==row["target"])&(longdf.contrast==row["contrast"])]["BH18_q_descriptive"].iloc[0]
            row["BH18_q_descriptive"]=q
    write_df(out/"results/PER_CLASS_CONTRASTS_LONG.tsv",longdf)
    write_df(out/"results/PER_CLASS_CONTRASTS_LONG.csv",longdf,sep=",")
    write_json(out/"results/PER_CLASS_CONTRASTS.json",longdf.where(pd.notna(longdf),None).to_dict(orient="records"))
    write_df(out/"results/DESCRIPTIVE_BH18_FAMILIES.tsv",pd.DataFrame(family_rows))
    write_df(out/"results/C2_RECONSTRUCTION_CHECKS.tsv",pd.DataFrame(recon_rows))
    write_df(out/"results/CPTAC_FIXED_EFFECT_META.tsv",pd.DataFrame(meta_rows))
    taxrows=[{"target":t,"model":m,"interpretation_state":direct_states[(m,t)]} for m in MODELS for t in TARGETS]
    taxrows += [{"target":t,"model":META_MODEL,"interpretation_state":meta_states[t]} for t in TARGETS]
    write_df(out/"results/INTERPRETATION_TAXONOMY.tsv",pd.DataFrame(taxrows))
    het=longdf[longdf.target.isin(["GATA2","SOX9"])].copy()
    write_df(out/"results/GATA2_SOX9_HETEROGENEITY.tsv",het)
    write_df(out/"results/SIX_TARGET_COMPLETENESS.tsv",longdf)
    present_cols=["target","model","contrast","coefficient","coefficient_SE","residual_SD","d","d_boot_CI_lo","d_boot_CI_hi","perm_p_raw","BH18_q_descriptive","n_POLE","n_MMRd","n_NSMP","n_p53abn","interpretation_state"]
    pres=longdf[present_cols].copy(); write_df(out/"results/MANUSCRIPT_READY_POSTHOC_TABLE.tsv",pres)
    make_markdown_table(out/"results/MANUSCRIPT_READY_POSTHOC_TABLE.md",pres)
    make_heterogeneity_markdown(out/"results/GATA2_SOX9_HETEROGENEITY.md",het)
    for target in ["GATA2","SOX9"]: make_forest(out,target,longdf)
    make_report(out,longdf,pd.DataFrame(recon_rows),pd.DataFrame(meta_rows),family_rows,direct_states,meta_states,count_rows)
    end_snap={"original_dirty_worktree":git_snapshot(ORIGINAL),"task030_verified_worktree":git_snapshot(T30)}
    preservation={"start":prov["start_snapshots"],"end":end_snap,
                  "original_preserved":prov["start_snapshots"]["original_dirty_worktree"]==end_snap["original_dirty_worktree"],
                  "task030_preserved":prov["start_snapshots"]["task030_verified_worktree"]==end_snap["task030_verified_worktree"],
                  "pinned_inputs_reverified":all(sha256_file(Path(r["path"]))==r["expected_sha256"] for r in prov["required_inputs"])}
    write_json(out/"integrity/UPSTREAM_AND_DIRTY_WORKTREE_PRESERVATION.json",preservation)
    if not preservation["original_preserved"] or not preservation["task030_preserved"] or not preservation["pinned_inputs_reverified"]:
        raise RuntimeError(f"UPSTREAM_PRESERVATION_FAILURE {preservation}")
    env={"python_executable":sys.executable,"python_version":sys.version,"platform":platform.platform(),
         "packages":{"numpy":np.__version__,"pandas":pd.__version__,"scipy":scipy.__version__,"matplotlib":matplotlib.__version__},
         "environment":{"OPENBLAS_NUM_THREADS":os.environ.get("OPENBLAS_NUM_THREADS"),"OMP_NUM_THREADS":os.environ.get("OMP_NUM_THREADS"),"PYTHONHASHSEED":os.environ.get("PYTHONHASHSEED")}}
    write_json(out/"SESSION_ENVIRONMENT.json",env)
    (out/"DEPENDENCIES.lock").write_text("\n".join([f"matplotlib=={matplotlib.__version__}",f"numpy=={np.__version__}",f"pandas=={pd.__version__}",f"scipy=={scipy.__version__}"])+"\n",encoding="ascii")
    manifest={"task":"TASK_A_PERCLASS_C2","classification":"POSTHOC_EXPLANATORY_DESCRIPTIVE_NONCAUSAL",
      "status":"COMPLETE","git":{"head":EXPECTED_HEAD,"branch":EXPECTED_BRANCH},"master_seed":MASTER_SEED,
      "subseed_rule":"int(sha256('20260722:' + step_id).hexdigest()[:8],16)","n_boot":N_BOOT,"n_perm":N_PERM,
      "data_versions":spec["data_versions"],"sample_counts":count_rows,"all_21_pinned_inputs_match":True,
      "nested_seals_pass":True,"all_direct_c2_identities_pass":all(r["pass"] for r in recon_rows),
      "bh18_families":family_rows,"upstream_preservation":preservation,"deviations":"none",
      "firewall":{"posthoc":True,"descriptive":True,"noncausal":True,"frozen_verdicts_unchanged":True,
                  "manuscript_edited":False,"src_edited":False,"upstream_written":False,"thresholds_tuned":False}}
    write_json(out/"REPRODUCIBILITY_MANIFEST.json",manifest)
    (out/"SEED_SCHEME.md").write_text("# Seed scheme\n\nMaster seed: `20260722`.\n\nSubseed: `int(sha256('20260722:' + step_id).hexdigest()[:8], 16)`, with exact frozen step IDs `boot__<model>__<target>` and `perm__<model>__<target>`.\n",encoding="ascii")
    print("COMPLETE",out,flush=True)


def fmt(x, digits=3):
    if x is None or pd.isna(x): return "NA"
    return f"{float(x):.{digits}f}"


def make_markdown_table(path:Path,df:pd.DataFrame):
    lines=["Post-hoc explanatory sensitivity: per-class decomposition of frozen C2; descriptive, non-causal, and unable to change frozen categories or replication verdicts.","",
           "| Target | Model | Contrast | b (SE) | Residual SD | d [95% bootstrap CI] | Raw p | BH-18 q | Counts P/M/N/p53 | State |",
           "|---|---|---|---:|---:|---:|---:|---:|---:|---|"]
    for _,r in df.iterrows():
        lines.append(f"| {r.target} | {r.model} | {r.contrast} | {fmt(r.coefficient)} ({fmt(r.coefficient_SE)}) | {fmt(r.residual_SD)} | {fmt(r.d)} [{fmt(r.d_boot_CI_lo)}, {fmt(r.d_boot_CI_hi)}] | {fmt(r.perm_p_raw,4)} | {fmt(r.BH18_q_descriptive,4)} | {int(r.n_POLE)}/{int(r.n_MMRd)}/{int(r.n_NSMP)}/{int(r.n_p53abn)} | {r.interpretation_state} |")
    path.write_text("\n".join(lines)+"\n",encoding="ascii")


def make_heterogeneity_markdown(path:Path,df:pd.DataFrame):
    lines=["# GATA2 and SOX9 per-class heterogeneity", "", "Post-hoc, descriptive, non-causal. Direct same-fit POLE-minus-MMRd inference is used; non-significance is never equality.",""]
    lines += ["| Target | Model | Contrast | b | SE | 95% model CI | SD | d | 95% boot d CI | raw p | q18 | State |","|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|"]
    for _,r in df.iterrows():
        lines.append(f"| {r.target} | {r.model} | {r.contrast} | {fmt(r.coefficient)} | {fmt(r.coefficient_SE)} | [{fmt(r.coefficient_t_CI_lo)}, {fmt(r.coefficient_t_CI_hi)}] | {fmt(r.residual_SD)} | {fmt(r.d)} | [{fmt(r.d_boot_CI_lo)}, {fmt(r.d_boot_CI_hi)}] | {fmt(r.perm_p_raw,4)} | {fmt(r.BH18_q_descriptive,4)} | {r.interpretation_state} |")
    path.write_text("\n".join(lines)+"\n",encoding="ascii")


def make_forest(out:Path,target:str,df:pd.DataFrame):
    order=["TCGA_PRIMARY_CPE_N506","TCGA_NOPURITY_N507","CPTAC_DISCOVERY_N95","CPTAC_CONFIRMATORY_N135",META_MODEL]
    contrasts=["POLE_vs_NSMP","MMRd_vs_NSMP","POLE_vs_MMRd"]
    fig,ax=plt.subplots(figsize=(8.5,8.0)); y=[]; labels=[]; colors={contrasts[0]:"#2166ac",contrasts[1]:"#b2182b",contrasts[2]:"#4d4d4d"}
    pos=0
    for model in order:
        for c in contrasts:
            r=df[(df.target==target)&(df.model==model)&(df.contrast==c)].iloc[0]
            ax.errorbar(r.d,pos,xerr=np.array([[r.d-r.d_boot_CI_lo],[r.d_boot_CI_hi-r.d]]),fmt="o",color=colors[c],capsize=3,markersize=5)
            y.append(pos); labels.append(f"{model} | {c}"); pos+=1
        pos+=0.45
    ax.axvline(0,color="black",linewidth=.8); ax.set_yticks(y,labels); ax.invert_yaxis(); ax.set_xlabel("Standardized effect d (95% bootstrap/FE CI)"); ax.set_title(f"{target}: post-hoc per-class decomposition"); ax.grid(axis="x",alpha=.2); fig.tight_layout()
    fdir=out/"figures"; fdir.mkdir(parents=True,exist_ok=True)
    plt.rcParams["svg.hashsalt"]="taskA-perclass-c2"
    fig.savefig(fdir/f"{target}_PER_CLASS_FOREST.svg",metadata={"Date":None})
    fig.savefig(fdir/f"{target}_PER_CLASS_FOREST.png",dpi=180,metadata={"Software":"TASK_A_PERCLASS_C2"})
    fig.savefig(fdir/f"{target}_PER_CLASS_FOREST.pdf",metadata={"Creator":"TASK_A_PERCLASS_C2","CreationDate":datetime(2026,7,22,tzinfo=timezone.utc),"ModDate":datetime(2026,7,22,tzinfo=timezone.utc)})
    plt.close(fig)


def make_report(out:Path,longdf:pd.DataFrame,recon:pd.DataFrame,meta:pd.DataFrame,families,direct_states,meta_states,counts):
    lines=["# TASK A analytical report","","Status: COMPLETE producer execution; independent critic verdict pending.","",
      "This is a post-hoc explanatory sensitivity. It is descriptive and non-causal. It cannot change any frozen TCGA/CPTAC category or replication verdict, and no manuscript byte was edited.","","## Sample guards","",
      "All four exact analytic counts passed: TCGA primary 506, TCGA no-purity 507, CPTAC Discovery 95, and CPTAC Confirmatory 135. Exact subtype counts are in `results/SAMPLE_COUNTS.tsv`.","","## GATA2 and SOX9","",
      "The mechanical taxonomy below uses the direct same-fit POLE-minus-MMRd interval. A compatible or unresolved result is not equality or equivalence.",""]
    for t in ["GATA2","SOX9"]:
        lines.append(f"### {t}"); lines.append("")
        for m in MODELS: lines.append(f"- {m}: `{direct_states[(m,t)]}`")
        lines.append(f"- CPTAC fixed-effect meta: `{meta_states[t]}`"); lines.append("")
    lines += ["The equal-weight C2 is not sample-size weighted; POLE:MMRd sample proportions alone are not a class-mix artifact.","","## Algebra and multiplicity","",
      f"All {len(recon)} direct target/model coefficient and d reconstructions passed scaled tolerance 1e-12, including replicate-level bootstrap identity. Maximum observed replicate errors are recorded in `results/C2_RECONSTRUCTION_CHECKS.tsv`.",
      "CPTAC fixed-effect C2 identity is not required because residual scales and contrast-specific inverse-bootstrap-variance weights differ. Direct FE(C2) discrepancies and every weight are reported in `results/CPTAC_FIXED_EFFECT_META.tsv`.",
      "Exactly five separate descriptive BH-18 families were computed. These q values do not confer confirmatory credit.","","## Interpretation boundary","",
      "Cross-cohort differences can reflect biology, composition, acquisition, platform, classifier, or scoring context; this analysis cannot identify cause. Bulk subtype-level results are not individual-patient biomarkers. Frozen verdicts are preserved and no manuscript edit was made.",""]
    (out/"ANALYTICAL_REPORT.md").write_text("\n".join(lines),encoding="ascii")


def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--output-root",required=True,type=Path); args=ap.parse_args()
    execute(args.output_root.resolve())

if __name__=="__main__": main()
