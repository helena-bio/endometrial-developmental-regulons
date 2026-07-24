#!/usr/bin/env python3
"""
TASK-030 TCGA no-purity and purity-confounding audit.

This is a thin wrapper around the sealed TASK-028 model engine and frozen inputs.
It regenerates the TASK-028 PRIMARY and SENS_nopurity families, computes the
requested common-sample diagnostic, and joins the sealed TASK-029 comparator.
No target, model term, threshold, family, contrast, scoring rule, or verdict is
retuned or reselected. ASCII only.
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
from scipy import stats


ROOT = Path("data/external/original-workspace")
WORKTREE = ROOT / ".language-model tool/work/revgate-tcga-no-purity-verify"
T28 = ROOT / ".language-model tool/work/task028-freeze-b-draft"
T29 = ROOT / ".language-model tool/work/task029-external-replication-feasibility"
OUT = WORKTREE / "experiments/task030_verify_current_main/base_audit"
RESULTS = OUT / "results"
MANUSCRIPT = OUT / "manuscript"
GENE_LOO_SUMMARY = WORKTREE / "experiments/task030_verify_current_main/gene_loo/results/gene_loo_per_target_summary.tsv"
T28_RESULTS = T28 / "execution/results_v3"
T28_INTER = T28 / "execution/intermediate"
T28_SCRIPTS = T28 / "execution/scripts_v3"
T29_RESULTS = T29 / "execution/results"

sys.path.insert(0, str(T28_SCRIPTS))
import model_engine as ME  # noqa: E402 - sealed TASK-028 implementation


MASTER_SEED_UPSTREAM = 20260713
MASTER_SEED_WRAPPER = None
N_BOOT = 2000
N_PERM = 2000
SESOI = 0.30
FLOOR = 0.50
SUBTYPES = ["POLE", "MMRd", "NSMP", "p53abn"]
TARGETS = [
    ("GATA2", "C2", -1),
    ("SOX9", "C2", -1),
    ("HOXA9", "C2", -1),
    ("WT1", "C2", -1),
    ("PAX8", "C1", +1),
    ("LHX1", "C1", +1),
]
TARGET_ORDER = {t: i for i, (t, _, _) in enumerate(TARGETS)}


def utcnow() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def substep_seed(step_id: str) -> int:
    text = f"{MASTER_SEED_UPSTREAM}:{step_id}".encode("ascii")
    return int(hashlib.sha256(text).hexdigest()[:8], 16)


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    df.to_csv(path, sep="\t", index=False, float_format="%.17g")


def verify_standard_seal(seal_path: Path, base: Path) -> list[dict]:
    rows = []
    for line in seal_path.read_text(encoding="ascii").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        digest, rel = line.split(None, 1)
        rel = rel.strip()
        path = (base / rel).resolve()
        actual = sha256_file(path) if path.is_file() else None
        rows.append({
            "seal": str(seal_path), "path": str(path), "expected": digest,
            "actual": actual, "match": actual == digest,
        })
    return rows


def verify_task028_spec_seal(seal_path: Path, base: Path) -> list[dict]:
    """TASK-028 spec seal stores subdirectory in the first token (B1/<sha>)."""
    rows = []
    for line in seal_path.read_text(encoding="ascii").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        first, filename = line.split(None, 1)
        filename = filename.strip()
        if "/" in first:
            prefix, digest = first.rsplit("/", 1)
            path = base / prefix / filename
        else:
            digest = first
            path = base / filename
        actual = sha256_file(path) if path.is_file() else None
        rows.append({
            "seal": str(seal_path), "path": str(path.resolve()),
            "expected": digest, "actual": actual, "match": actual == digest,
        })
    return rows


def categorize(fit: dict) -> dict:
    out = {}
    gate = fit["omnibus"]["F1_gate_pass"]
    for contrast in ["C1", "C2", "C3"]:
        cc = fit["contrasts"][contrast]
        boot = fit["bootstrap"][contrast]
        d = cc["d"]
        q = cc.get("BH_q_F2")
        lo, hi = boot["d_ci_lo"], boot["d_ci_hi"]
        ci_excl0 = bool(np.isfinite(lo) and np.isfinite(hi) and (lo > 0 or hi < 0))
        cat = "DESCRIPTIVE"
        reasons = []
        if gate and q is not None and q <= 0.05 and abs(d) >= FLOOR and ci_excl0:
            cat = "CAT1_ENRICHMENT" if d > 0 else "CAT2_DEPLETION"
        else:
            if contrast == "C1" and gate and lo > -SESOI and hi < SESOI:
                cat = "CAT5_CLEAN_NULL"
                reasons.append("H3 equivalence (C1) CI-in-SESOI")
            if cat == "DESCRIPTIVE":
                if not gate:
                    reasons.append("omnibus F1 gate FAIL -> descriptive/exploratory")
                else:
                    if q is None or q > 0.05:
                        reasons.append("F2 q>0.05")
                    if abs(d) < FLOOR:
                        reasons.append("|d|<0.5")
                    if not ci_excl0:
                        reasons.append("bootstrap CI includes 0")
        out[contrast] = {
            "category": cat, "reason": "; ".join(reasons),
            "ci_excludes_0": ci_excl0,
        }
    return out


def run_family(
    config_name: str,
    modules: dict[str, np.ndarray],
    subtype: np.ndarray,
    cov_dict: dict[str, np.ndarray],
    cov_names: list[str],
    mask: np.ndarray,
) -> dict:
    """Re-run one complete sealed TASK-028 family, including bootstrap/permutation."""
    result = {}
    for module, y in modules.items():
        print(f"[regenerate] {config_name} {module}", flush=True)
        fit = ME.fit_module(y, subtype, cov_dict, cov_names, mask)
        fit["bootstrap"] = ME.bootstrap_contrasts(
            y, subtype, cov_dict, cov_names, mask,
            substep_seed(f"boot__{config_name}__{module}"), n_boot=N_BOOT,
        )
        fit["permutation"] = ME.permutation_omnibus_and_contrasts(
            y, subtype, cov_dict, cov_names, mask,
            substep_seed(f"perm__{config_name}__{module}"), n_perm=N_PERM,
        )
        result[module] = fit

    q1 = ME.bh_fdr({k: v["omnibus"]["p_F"] for k, v in result.items()})
    for module, fit in result.items():
        fit["omnibus"]["BH_q_F1"] = q1[module]
        fit["omnibus"]["F1_gate_pass"] = bool(q1[module] <= 0.05)

    p2 = {}
    for module, fit in result.items():
        if fit["omnibus"]["F1_gate_pass"]:
            for contrast in ["C1", "C2", "C3"]:
                p2[f"{module}::{contrast}"] = fit["permutation"]["perm_p_contrast"][contrast]
    q2 = ME.bh_fdr(p2)
    for key, q in q2.items():
        module, contrast = key.split("::")
        result[module]["contrasts"][contrast]["perm_p"] = p2[key]
        result[module]["contrasts"][contrast]["BH_q_F2"] = q
    for fit in result.values():
        fit["b3_category"] = categorize(fit)
    return result


def flatten_family(config_name: str, family: dict) -> pd.DataFrame:
    rows = []
    for module, fit in family.items():
        for contrast in ["C1", "C2", "C3"]:
            cc = fit["contrasts"][contrast]
            boot = fit["bootstrap"][contrast]
            om = fit["omnibus"]
            cat = fit["b3_category"][contrast]
            rows.append({
                "config": config_name, "module": module, "contrast": contrast,
                "n": fit["n"], "sigma_resid": fit["sigma_resid"],
                "omnibus_F": om["F"], "omnibus_p_F": om["p_F"],
                "omnibus_BH_q_F1": om["BH_q_F1"],
                "F1_gate_pass": om["F1_gate_pass"],
                "estimate": cc["estimate"], "d": cc["d"],
                "est_ci_lo": boot["est_ci_lo"], "est_ci_hi": boot["est_ci_hi"],
                "d_ci_lo": boot["d_ci_lo"], "d_ci_hi": boot["d_ci_hi"],
                "perm_p": cc.get("perm_p"), "BH_q_F2": cc.get("BH_q_F2"),
                "b3_category": cat["category"], "category_reason": cat["reason"],
            })
    return pd.DataFrame(rows)


def compare_regenerated(regen: pd.DataFrame, sealed_json: dict) -> pd.DataFrame:
    rows = []
    config_map = {"PRIMARY": sealed_json["primary"],
                  "SENS_nopurity": sealed_json["sensitivities"]["SENS_nopurity"]}
    numeric = [
        "n", "sigma_resid", "omnibus_F", "omnibus_p_F", "omnibus_BH_q_F1",
        "estimate", "d", "est_ci_lo", "est_ci_hi", "d_ci_lo", "d_ci_hi",
        "perm_p", "BH_q_F2",
    ]
    for r in regen.itertuples(index=False):
        sfit = config_map[r.config][r.module]
        scc = sfit["contrasts"][r.contrast]
        sb = sfit["bootstrap"][r.contrast]
        som = sfit["omnibus"]
        sealed = {
            "n": sfit["n"], "sigma_resid": sfit["sigma_resid"],
            "omnibus_F": som["F"], "omnibus_p_F": som["p_F"],
            "omnibus_BH_q_F1": som["BH_q_F1"],
            "estimate": scc["estimate"], "d": scc["d"],
            "est_ci_lo": sb["est_ci_lo"], "est_ci_hi": sb["est_ci_hi"],
            "d_ci_lo": sb["d_ci_lo"], "d_ci_hi": sb["d_ci_hi"],
            "perm_p": scc.get("perm_p"), "BH_q_F2": scc.get("BH_q_F2"),
        }
        deltas = []
        mismatches = []
        for field in numeric:
            a = getattr(r, field)
            b = sealed[field]
            if pd.isna(a) and (b is None or pd.isna(b)):
                continue
            if b is None or pd.isna(a):
                mismatches.append(field)
                continue
            delta = abs(float(a) - float(b))
            deltas.append(delta)
            if not math.isclose(float(a), float(b), rel_tol=1e-12, abs_tol=1e-12):
                mismatches.append(field)
        sealed_cat = sfit["b3_category"][r.contrast]["category"]
        if r.b3_category != sealed_cat:
            mismatches.append("b3_category")
        if bool(r.F1_gate_pass) != bool(som["F1_gate_pass"]):
            mismatches.append("F1_gate_pass")
        rows.append({
            "config": r.config, "module": r.module, "contrast": r.contrast,
            "numeric_max_abs_delta": max(deltas) if deltas else 0.0,
            "category_match": r.b3_category == sealed_cat,
            "gate_match": bool(r.F1_gate_pass) == bool(som["F1_gate_pass"]),
            "all_fields_match_tolerance_1e-12": not mismatches,
            "mismatch_fields": ";".join(mismatches),
        })
    return pd.DataFrame(rows)


def shapley_decompose(b0: float, s0: float, b1: float, s1: float) -> tuple[float, float, float]:
    """Symmetric exact decomposition of d1-d0 for d=b/s.

    Numerator = 0.5 * (b1-b0) * (1/s0 + 1/s1)
    Denominator = 0.5 * (b0+b1) * (1/s1 - 1/s0)
    """
    numerator = 0.5 * (b1 - b0) * (1.0 / s0 + 1.0 / s1)
    denominator = 0.5 * (b0 + b1) * (1.0 / s1 - 1.0 / s0)
    err = (b1 / s1 - b0 / s0) - numerator - denominator
    return numerator, denominator, err


def pooled_standardized_difference(a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    diff = float(np.mean(a) - np.mean(b))
    pooled_var = ((len(a) - 1) * np.var(a, ddof=1) + (len(b) - 1) * np.var(b, ddof=1)) / (len(a) + len(b) - 2)
    return diff, diff / math.sqrt(pooled_var)


def criterion_summary(row: dict, pred_sign: int) -> str:
    checks = {
        "omnibus_gate": bool(row["F1_gate_pass"]),
        "effect_floor": abs(row["d"]) >= FLOOR,
        "ci_excludes_zero": bool(row["d_ci_lo"] > 0 or row["d_ci_hi"] < 0),
        "frozen_F2_q": row["BH_q_F2"] <= 0.05,
        "prespecified_direction": bool(np.sign(row["d"]) == pred_sign),
    }
    failed = [k for k, v in checks.items() if not v]
    if not failed:
        return "all relevant frozen criteria remain satisfied"
    return "non-confirmatory criterion status: failed " + ", ".join(failed)


def main() -> None:
    started = utcnow()
    RESULTS.mkdir(parents=True, exist_ok=True)
    MANUSCRIPT.mkdir(parents=True, exist_ok=True)

    # Verify every upstream sealed byte before analysis.
    seal_rows = []
    seal_rows += verify_task028_spec_seal(T28 / "sealed_v3/SEAL_MANIFEST.sha256", T28 / "sealed_v3")
    seal_rows += verify_standard_seal(T28_RESULTS / "RESULT_SEAL.sha256", T28_RESULTS)
    seal_rows += verify_standard_seal(T29_RESULTS / "RESULT_SEAL.sha256", T29_RESULTS)
    seal_df = pd.DataFrame(seal_rows)
    write_tsv(seal_df, RESULTS / "upstream_checksum_verification.tsv")
    if not bool(seal_df["match"].all()):
        bad = seal_df.loc[~seal_df["match"], ["path", "expected", "actual"]]
        raise RuntimeError("upstream checksum mismatch:\n" + bad.to_string(index=False))

    scores = np.load(T28_INTER / "scores_v3.npz", allow_pickle=True)
    cov = pd.read_csv(T28_INTER / "covariates_v3.tsv", sep="\t")
    assert list(cov["patient_barcode"]) == list(scores["patient_order"])
    subtype = cov["subtype"].to_numpy()
    cov_dict = {
        "M4": cov["M4_prolif"].to_numpy(float),
        "purity_CPE": cov["purity_CPE"].to_numpy(float),
        "composition": cov["composition"].to_numpy(float),
    }
    mask506 = cov["cpe_complete_case"].to_numpy(int) == 1
    mask507 = np.ones(len(cov), dtype=bool)
    core = [str(x) for x in list(scores["CORE_TFS"])]
    modules = {"M1": np.asarray(scores["M1"], float)}
    modules.update({f"M3_{tf}": np.asarray(scores[f"M3primary__{tf}"], float) for tf in core})

    primary = run_family(
        "PRIMARY", modules, subtype, cov_dict,
        ["M4", "purity_CPE", "composition"], mask506,
    )
    nopurity = run_family(
        "SENS_nopurity", modules, subtype, cov_dict,
        ["M4", "composition"], mask507,
    )
    regen = pd.concat([
        flatten_family("PRIMARY", primary),
        flatten_family("SENS_nopurity", nopurity),
    ], ignore_index=True)
    write_tsv(regen, RESULTS / "regenerated_full_families.tsv")

    sealed_models = json.loads((T28_RESULTS / "phase2b_models.json").read_text())
    regen_audit = compare_regenerated(regen, sealed_models)
    write_tsv(regen_audit, RESULTS / "regeneration_audit.tsv")
    if not bool(regen_audit["all_fields_match_tolerance_1e-12"].all()):
        bad = regen_audit.loc[~regen_audit["all_fields_match_tolerance_1e-12"]]
        raise RuntimeError("sealed result regeneration failed:\n" + bad.to_string(index=False))

    # Common-sample no-purity is a diagnostic point-estimate model only.
    common = {}
    for target, contrast, _ in TARGETS:
        module = f"M3_{target}"
        common[module] = ME.fit_module(
            modules[module], subtype, cov_dict, ["M4", "composition"], mask506,
        )

    cptac = json.loads((T29_RESULTS / "primary_results.json").read_text())
    primary_flat = flatten_family("PRIMARY", primary)
    no_flat = flatten_family("SENS_nopurity", nopurity)

    target_rows = []
    sens_rows = []
    corr_rows = []
    decomp_rows = []
    gap_rows = []
    for target, contrast, pred_sign in TARGETS:
        module = f"M3_{target}"
        pr = primary_flat[(primary_flat.module == module) & (primary_flat.contrast == contrast)].iloc[0].to_dict()
        nr = no_flat[(no_flat.module == module) & (no_flat.contrast == contrast)].iloc[0].to_dict()
        cm = common[module]
        cm_b = cm["contrasts"][contrast]["estimate"]
        cm_s = cm["sigma_resid"]
        cm_d = cm["contrasts"][contrast]["d"]
        cmeta = cptac["meta"][target]
        cver = cptac["per_target_verdict"][target]
        direction_preserved = bool(np.sign(pr["d"]) == np.sign(nr["d"]))
        category_preserved = pr["b3_category"] == nr["b3_category"]
        reason = criterion_summary(nr, pred_sign)
        if not category_preserved:
            reason = "category changed; " + reason

        target_rows.append({
            "target": target, "contrast": contrast,
            "contrast_definition_full_subtype_model": "C1=[-1,-1,-1,+3]" if contrast == "C1" else "C2=[+1,+1,-2,0]",
            "predicted_direction": "d>0" if pred_sign > 0 else "d<0",
            "analytic_n_primary": int(pr["n"]), "analytic_n_nopurity": int(nr["n"]),
            "primary_contrast_coefficient": pr["estimate"],
            "primary_coefficient_ci_lo": pr["est_ci_lo"],
            "primary_coefficient_ci_hi": pr["est_ci_hi"],
            "primary_residual_sd": pr["sigma_resid"], "primary_d": pr["d"],
            "primary_d_ci_lo": pr["d_ci_lo"], "primary_d_ci_hi": pr["d_ci_hi"],
            "primary_raw_p": pr["perm_p"],
            "primary_raw_p_definition": "two-sided subtype-label permutation p on |d|, add-one, n_perm=2000",
            "primary_frozen_BH_q": pr["BH_q_F2"],
            "primary_frozen_family": "TASK-028 F2: all contrasts of F1-omnibus-gated modules (not BH-of-6)",
            "primary_frozen_category": pr["b3_category"],
            "nopurity_contrast_coefficient": nr["estimate"],
            "nopurity_coefficient_ci_lo": nr["est_ci_lo"],
            "nopurity_coefficient_ci_hi": nr["est_ci_hi"],
            "nopurity_residual_sd": nr["sigma_resid"], "nopurity_d": nr["d"],
            "nopurity_d_ci_lo": nr["d_ci_lo"], "nopurity_d_ci_hi": nr["d_ci_hi"],
            "nopurity_raw_p": nr["perm_p"],
            "nopurity_raw_p_definition": "two-sided subtype-label permutation p on |d|, add-one, n_perm=2000",
            "nopurity_frozen_BH_q": nr["BH_q_F2"],
            "nopurity_frozen_family": "TASK-028 sensitivity F2: all contrasts of sensitivity F1-gated modules (not BH-of-6)",
            "nopurity_frozen_category": nr["b3_category"],
            "direction_preserved": direction_preserved,
            "category_preserved": category_preserved,
            "category_change_or_criterion_reason": reason,
            "common_sample_n": int(cm["n"]),
            "common_sample_nopurity_coefficient": cm_b,
            "common_sample_nopurity_residual_sd": cm_s,
            "common_sample_nopurity_d": cm_d,
            "common_sample_label": "diagnostic only; not a frozen sensitivity or verdict-changing model",
        })

        score = modules[module][mask506]
        purity = cov_dict["purity_CPE"][mask506]
        pear = stats.pearsonr(score, purity)
        spear = stats.spearmanr(score, purity)
        corr_rows.append({
            "target": target, "n": int(mask506.sum()),
            "pearson_r_score_vs_Aran_CPE": float(pear.statistic),
            "pearson_p_two_sided": float(pear.pvalue),
            "spearman_rho_score_vs_Aran_CPE": float(spear.statistic),
            "spearman_p_two_sided": float(spear.pvalue),
        })

        sens_rows.append({
            "target": target, "contrast": contrast,
            "pearson_r_score_vs_Aran_CPE": float(pear.statistic),
            "spearman_rho_score_vs_Aran_CPE": float(spear.statistic),
            "primary_n": int(pr["n"]), "primary_coefficient": pr["estimate"],
            "primary_residual_sd": pr["sigma_resid"], "primary_d": pr["d"],
            "official_nopurity_n": int(nr["n"]),
            "official_nopurity_coefficient": nr["estimate"],
            "official_nopurity_residual_sd": nr["sigma_resid"],
            "official_nopurity_d": nr["d"],
            "common_nopurity_n": int(cm["n"]),
            "common_nopurity_coefficient": cm_b,
            "common_nopurity_residual_sd": cm_s,
            "common_nopurity_d": cm_d,
            "direction_preserved_official": direction_preserved,
            "category_primary": pr["b3_category"],
            "category_official_nopurity": nr["b3_category"],
            "category_preserved_official": category_preserved,
            "exact_outcome": reason,
        })

        stages = [
            ("purity_covariate_removal_common_n506", int(pr["n"]), int(cm["n"]), pr["estimate"], pr["sigma_resid"], cm_b, cm_s),
            ("sample_set_change_common_n506_to_official_n507", int(cm["n"]), int(nr["n"]), cm_b, cm_s, nr["estimate"], nr["sigma_resid"]),
            ("total_primary_n506_to_official_n507", int(pr["n"]), int(nr["n"]), pr["estimate"], pr["sigma_resid"], nr["estimate"], nr["sigma_resid"]),
        ]
        for stage, n0, n1, b0, s0, b1, s1 in stages:
            num, den, err = shapley_decompose(b0, s0, b1, s1)
            decomp_rows.append({
                "target": target, "contrast": contrast, "stage": stage,
                "from_n": n0, "to_n": n1, "from_coefficient": b0,
                "to_coefficient": b1, "from_residual_sd": s0,
                "to_residual_sd": s1, "from_d": b0 / s0, "to_d": b1 / s1,
                "delta_d": b1 / s1 - b0 / s0,
                "coefficient_numerator_component": num,
                "residual_sd_denominator_component": den,
                "decomposition_identity_error": err,
                "decomposition": "symmetric Shapley: num=0.5*(b1-b0)*(1/s0+1/s1); den=0.5*(b0+b1)*(1/s1-1/s0)",
            })

        gap_primary = abs(pr["d"] - cmeta["d_meta"])
        gap_no = abs(nr["d"] - cmeta["d_meta"])
        gap_rows.append({
            "target": target, "contrast": contrast,
            "tcga_primary_d": pr["d"], "tcga_official_nopurity_d": nr["d"],
            "cptac_meta_d": cmeta["d_meta"],
            "absolute_gap_primary_vs_cptac": gap_primary,
            "absolute_gap_nopurity_vs_cptac": gap_no,
            "nopurity_moves_closer_to_cptac": gap_no < gap_primary,
            "percent_magnitude_change_primary_to_nopurity": ((abs(nr["d"]) - abs(pr["d"])) / abs(pr["d"]) * 100.0) if pr["d"] != 0 else np.nan,
            "percent_absolute_gap_change": ((gap_no - gap_primary) / gap_primary * 100.0) if gap_primary != 0 else np.nan,
            "comparison_scope": "descriptive only; not a replication criterion",
        })

    target_df = pd.DataFrame(target_rows).sort_values("target", key=lambda x: x.map(TARGET_ORDER))
    sens_df = pd.DataFrame(sens_rows).sort_values("target", key=lambda x: x.map(TARGET_ORDER))
    corr_df = pd.DataFrame(corr_rows).sort_values("target", key=lambda x: x.map(TARGET_ORDER))
    decomp_df = pd.DataFrame(decomp_rows).sort_values(["target", "stage"], key=lambda x: x.map(TARGET_ORDER) if x.name == "target" else x)
    gap_df = pd.DataFrame(gap_rows).sort_values("target", key=lambda x: x.map(TARGET_ORDER))
    write_tsv(target_df, RESULTS / "tcga_primary_vs_nopurity_per_target.tsv")
    write_tsv(sens_df, RESULTS / "per_target_sensitivity_audit.tsv")
    write_tsv(corr_df, RESULTS / "purity_score_correlations.tsv")
    write_tsv(decomp_df, RESULTS / "delta_d_decomposition.tsv")
    write_tsv(gap_df, RESULTS / "descriptive_tcga_cptac_gaps.tsv")

    # Purity side-by-side group summaries on the common CPE-complete sample.
    purity = cov_dict["purity_CPE"]
    group_specs = [
        ("C2", "POLE+MMRd", np.isin(subtype, ["POLE", "MMRd"]) & mask506,
         "NSMP", (subtype == "NSMP") & mask506),
        ("C1", "p53abn", (subtype == "p53abn") & mask506,
         "pooled_non-p53abn", (subtype != "p53abn") & mask506),
    ]
    group_rows = []
    subtype_means = {s: float(np.mean(purity[(subtype == s) & mask506])) for s in SUBTYPES}
    for contrast, name_a, mask_a, name_b, mask_b in group_specs:
        a, b = purity[mask_a], purity[mask_b]
        diff, std_diff = pooled_standardized_difference(a, b)
        equal_subtype_diff = (
            0.5 * (subtype_means["POLE"] + subtype_means["MMRd"]) - subtype_means["NSMP"]
            if contrast == "C2" else
            subtype_means["p53abn"] - (subtype_means["POLE"] + subtype_means["MMRd"] + subtype_means["NSMP"]) / 3.0
        )
        group_rows.append({
            "contrast": contrast, "group_a": name_a, "group_a_n": len(a),
            "group_a_mean_Aran_CPE": float(np.mean(a)), "group_a_sd_Aran_CPE": float(np.std(a, ddof=1)),
            "group_b": name_b, "group_b_n": len(b),
            "group_b_mean_Aran_CPE": float(np.mean(b)), "group_b_sd_Aran_CPE": float(np.std(b, ddof=1)),
            "pooled_side_mean_difference_a_minus_b": diff,
            "pooled_sd_standardized_difference_Cohen_d": std_diff,
            "frozen_equal_subtype_weight_contrast_difference": equal_subtype_diff,
            "note": "Purity summary is descriptive; pooled side means differ from the equal-subtype-weight frozen model contrast where side A has multiple subtypes.",
        })
    write_tsv(pd.DataFrame(group_rows), RESULTS / "purity_group_differences.tsv")

    # Join complete B2.12 before rendering any category/replication narrative.
    # The top-level workflow runs gene-LOO first, so a clean delivery cannot
    # silently fall back to the incomplete pointwise layer.
    if not GENE_LOO_SUMMARY.is_file():
        raise FileNotFoundError(
            f"complete B2.12 summary is required before base report generation: {GENE_LOO_SUMMARY}"
        )
    loo_summary = pd.read_csv(GENE_LOO_SUMMARY, sep="\t")
    expected_loo_keys = {
        (target, config)
        for target, _, _ in TARGETS
        for config in ["PRIMARY", "SENS_nopurity"]
    }
    observed_loo_keys = set(zip(loo_summary["target"], loo_summary["config"]))
    if expected_loo_keys != observed_loo_keys:
        raise RuntimeError("complete B2.12 summary keys do not match the frozen six-target family")
    loo_by_key = loo_summary.set_index(["target", "config"])

    # Pointwise labels remain only as explicitly fenced legacy fields.
    manuscript_rows = []
    for row in target_df.to_dict("records"):
        target = row["target"]
        cm = cptac["meta"][target]
        cv = cptac["per_target_verdict"][target]
        primary_loo = loo_by_key.loc[(target, "PRIMARY")]
        nopurity_loo = loo_by_key.loc[(target, "SENS_nopurity")]
        manuscript_rows.append({
            "Target": target, "Contrast": row["contrast"],
            "TCGA primary d (95% CI)": f"{row['primary_d']:.2f} ({row['primary_d_ci_lo']:.2f}, {row['primary_d_ci_hi']:.2f})",
            "TCGA primary legacy pointwise label (incomplete B2.12)": row["primary_frozen_category"],
            "TCGA primary complete B2.12 credit": primary_loo["B2_12_compliant_credit_status"],
            "TCGA no-purity d (95% CI)": f"{row['nopurity_d']:.2f} ({row['nopurity_d_ci_lo']:.2f}, {row['nopurity_d_ci_hi']:.2f})",
            "TCGA no-purity legacy pointwise label (incomplete B2.12)": row["nopurity_frozen_category"],
            "TCGA no-purity complete B2.12 credit": nopurity_loo["B2_12_compliant_credit_status"],
            "CPTAC d_meta (95% CI)": f"{cm['d_meta']:.2f} ({cm['ci_lo']:.2f}, {cm['ci_hi']:.2f})",
            "CPTAC BH q of 6": f"{cm['BH_q_of_6']:.6g}",
            "Frozen CPTAC byte status (not TCGA B2.12 credit)": cv["status"],
        })
    manuscript_df = pd.DataFrame(manuscript_rows)
    write_tsv(manuscript_df, MANUSCRIPT / "tcga_nopurity_cptac_table.tsv")
    md_header = "| " + " | ".join(manuscript_df.columns) + " |"
    md_rule = "| " + " | ".join(["---"] * len(manuscript_df.columns)) + " |"
    md_rows = [
        "| " + " | ".join(str(v).replace("|", "\\|") for v in row) + " |"
        for row in manuscript_df.itertuples(index=False, name=None)
    ]
    md_lines = [
        "# TCGA no-purity sensitivity and CPTAC comparator",
        "",
        "Pointwise TCGA labels below are legacy/incomplete-B2.12 fields only. Complete B2.12 universal gene-LOO credit is binding. TCGA q-values use the frozen TASK-028 full gated F2 family. CPTAC q-values use the distinct frozen BH-of-6 family. CPTAC is a descriptive comparator here; no new replication criterion is introduced.",
        "",
        md_header,
        md_rule,
        *md_rows,
        "",
    ]
    (MANUSCRIPT / "tcga_nopurity_cptac_table.md").write_text("\n".join(md_lines), encoding="ascii")

    # Retain the pointwise sensitivity scenario as a legacy diagnostic, then
    # apply the binding complete B2.12 universal gene-LOO gate.
    gs = target_df[target_df.target.isin(["GATA2", "SOX9"])]
    gs_gap = gap_df[gap_df.target.isin(["GATA2", "SOX9"])]
    both_category = bool(gs["category_preserved"].all())
    both_direction = bool(gs["direction_preserved"].all())
    both_closer = bool(gs_gap["nopurity_moves_closer_to_cptac"].all())
    if both_category and both_direction and both_closer:
        pointwise_scenario = "A"
    elif both_category and both_direction:
        pointwise_scenario = "B"
    elif not both_category:
        pointwise_scenario = "C"
    else:
        pointwise_scenario = "D"

    gs_complete = loo_summary[
        loo_summary.target.isin(["GATA2", "SOX9"])
        & (loo_summary.config == "SENS_nopurity")
    ]
    gs_complete_credit = bool(gs_complete["survives_gene_LOO_any_deletion_rule"].all())
    scenario = pointwise_scenario if gs_complete_credit else "C"
    if scenario != "C":
        raise RuntimeError("binding complete B2.12 scenario must be C for current-main bytes")

    complete_rows = []
    for target, _, _ in TARGETS:
        primary_loo = loo_by_key.loc[(target, "PRIMARY")]
        nopurity_loo = loo_by_key.loc[(target, "SENS_nopurity")]
        complete_rows.append({
            "target": target,
            "primary_B2_12_credit": primary_loo["B2_12_compliant_credit_status"],
            "nopurity_B2_12_credit": nopurity_loo["B2_12_compliant_credit_status"],
            "primary_floor_fail_count": int(primary_loo["floor_fail_count"]),
            "nopurity_floor_fail_count": int(nopurity_loo["floor_fail_count"]),
        })
    complete_view = pd.DataFrame(complete_rows)
    complete_table_header = "| " + " | ".join(complete_view.columns) + " |"
    complete_table_rule = "| " + " | ".join(["---"] * len(complete_view.columns)) + " |"
    complete_table_rows = [
        "| " + " | ".join(str(v).replace("|", "\\|") for v in row) + " |"
        for row in complete_view.itertuples(index=False, name=None)
    ]
    complete_table = "\n".join(
        [complete_table_header, complete_table_rule, *complete_table_rows]
    )

    report = f"""# TASK-030 TCGA no-purity and purity-confounding audit

Status: DONE (experimenter computation; no critic verdict assigned here).

## Exact no-purity outcome

Binding integrated conclusion: Scenario {scenario} under the complete B2.12 universal gene-LOO gate. The older pointwise-only sensitivity layer is a legacy incomplete diagnostic (its mechanically calculated scenario is {pointwise_scenario}) and cannot assign B2.12 credit. GATA2, SOX9, HOXA9, and WT1 are uncredited/descriptive directional positives in both PRIMARY and SENS_nopurity; they are not CAT5 and are not credited CAT2. PAX8 and LHX1 are credited CAT1 in both configurations.

This statement applies only to the prespecified TASK-028 no-purity sensitivity. It does not summarize ABSOLUTE-purity, signed-ssGSEA, or any other scoring/purity sensitivity.

The legacy pointwise labels and numerical criteria remain unchanged, but complete credit additionally requires every mapped single-gene deletion to pass the universal gate. The complete result is generated from `gene_loo/results/gene_loo_per_target_summary.tsv`:

{complete_table}

The frozen TASK-029 GATA2/SOX9 byte statuses remain unchanged, but their evidence is target-level directional concordance only, not confirmatory external replication of a B2.12-credited TCGA discovery. HOXA9 and WT1 are evaluable-not-replicated. PAX8 and LHX1 were frozen sensitivity-only in CPTAC.

## Purity-confounding scope

The audit reports per-target Pearson and Spearman score-CPE correlations, subtype-side CPE summaries, and all six frozen target coefficients/residual SDs/d values for: primary n=506 CPE adjustment, official no-purity n=507, and a common-sample n=506 no-purity diagnostic. The common-sample model is diagnostic only and is not a frozen sensitivity or verdict-changing model. Delta-d is separated into purity-covariate removal and sample-set change, with an exact symmetric Shapley decomposition into coefficient/numerator and residual-SD/denominator components. These diagnostics do not establish that purity confounding was removed; residual purity and composition confounding cannot be excluded.

## Model equivalence assessment

At the named-variable level, the TCGA no-purity model and CPTAC primary model use the same covariate set and formula: subtype + M4 proliferation + ESTIMATE composition, with no purity term, the same sealed model engine, contrasts, score orientation, and 761 signed-edge definition. The analyses are nevertheless not literally identical covariate replication. TCGA fits one n=507 cohort and applies its own full F1/F2 frozen families. CPTAC scores and fits Discovery (n=95) and Confirmatory (n=135) separately, then inverse-variance fixed-effect meta-analyzes standardized effects and applies a distinct BH-of-6 target family plus stratum-consistency/opposite-direction rules. Covariate values are cohort/stratum-specific ranks, composition implementations are separately scored in each expression matrix, all-zero gene filtering is cohort/stratum-specific, subtype-classifier implementations differ, and the sample sources differ (TCGA vs independent CPTAC/GDC samples).

Accordingly, the TCGA sensitivity is a CPTAC-compatible specification: model-adapted external replication, not identical covariate replication.

## Statistical firewall and attestation

No target, threshold, contrast, family, scoring definition, edge set, covariate, category rule, or TASK-029 verdict was changed. No retuning and no new target selection occurred. TCGA no-purity used the frozen TASK-028 gated full F2 family, not BH across six targets. CPTAC BH-of-6 remains separate and descriptive in this audit. Sensitivity results were not used to rescue a primary result or alter F35. The C2 result remains a relative within-tumour contrast; no purity-removal or causal claim is made.

## Reproduction status

All TASK-028 and TASK-029 upstream seal entries verified. The full TASK-028 PRIMARY and SENS_nopurity families regenerated numerically under the matched TASK-028 environment; all compared fields, categories, and gates matched sealed phase2b_models.json at tolerance 1e-12. This is numeric regeneration, not a claim that newly written tables are byte-identical to upstream tables.
"""
    (OUT / "ANALYTICAL_REPORT.md").write_text(report, encoding="ascii")

    dirty = subprocess.run(
        ["git", "-C", str(WORKTREE), "status", "--porcelain=v1"],
        check=True, capture_output=True, text=True,
    ).stdout
    (OUT / "git_dirty_state.txt").write_text(dirty, encoding="ascii")
    branch = subprocess.run(
        ["git", "-C", str(WORKTREE), "branch", "--show-current"],
        check=True, capture_output=True, text=True,
    ).stdout.strip()
    head = subprocess.run(
        ["git", "-C", str(WORKTREE), "rev-parse", "HEAD"],
        check=True, capture_output=True, text=True,
    ).stdout.strip()
    pip_freeze = subprocess.run(
        [sys.executable, "-m", "pip", "freeze"], check=True,
        capture_output=True, text=True,
    ).stdout.splitlines()

    commands = [
        "cd data/external/original-workspace/revgate-tcga-no-purity-verify/experiments/task030_verify_current_main/base_audit",
        "data/external/original-workspace/task028-freeze-b-draft/verifier_v3/venv_match/bin/python -m py_compile run_task030_audit.py",
        "OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 data/external/original-workspace/task028-freeze-b-draft/verifier_v3/venv_match/bin/python run_task030_audit.py",
    ]
    (OUT / "COMMANDS.txt").write_text("\n".join(commands) + "\n", encoding="ascii")

    input_paths = [
        T28 / "sealed_v3/SEAL_MANIFEST.sha256",
        T28_RESULTS / "RESULT_SEAL.sha256",
        T28_RESULTS / "REPRODUCIBILITY_MANIFEST.json",
        T28_RESULTS / "phase2b_models.json",
        T28_RESULTS / "primary_results.tsv",
        T28_RESULTS / "sensitivity_SENS_nopurity.tsv",
        T28_INTER / "scores_v3.npz",
        T28_INTER / "covariates_v3.tsv",
        T28_SCRIPTS / "model_engine.py",
        T29 / "FROZEN_REPLICATION_DESIGN.md",
        T29_RESULTS / "RESULT_SEAL.sha256",
        T29_RESULTS / "REPRODUCIBILITY_MANIFEST.json",
        T29_RESULTS / "primary_results.json",
        T29_RESULTS / "table_primary_per_target.tsv",
        T29_RESULTS / "table_meta_fixedeffect.tsv",
        T29_RESULTS / "table_per_stratum_effects.tsv",
        GENE_LOO_SUMMARY,
    ]
    output_paths = sorted([
        p for p in OUT.rglob("*")
        if p.is_file()
        and "__pycache__" not in p.parts
        and p.name not in {"REPRODUCIBILITY_MANIFEST.json", "OUTPUT_SHA256SUMS.txt", "run.log"}
    ])
    completed = utcnow()
    manifest = {
        "task": "TASK-030 TCGA no-purity and purity-confounding audit",
        "status": "done",
        "classification": "VERIFY",
        "started_utc": started, "completed_utc": completed,
        "git": {
            "branch": branch, "HEAD": head, "authorized_existing_dirty_branch": False,
            "dirty": bool(dirty), "dirty_state_path": str(OUT / "git_dirty_state.txt"),
            "dirty_state_sha256": sha256_file(OUT / "git_dirty_state.txt"),
            "no_checkout_fetch_commit_push_add": True,
        },
        "upstream": {
            "task028_version": "Freeze B v3 PRIMARY=C",
            "task028_master_seed": MASTER_SEED_UPSTREAM,
            "task028_n_boot": N_BOOT, "task028_n_perm": N_PERM,
            "task028_SESOI": SESOI, "task028_H1_floor": FLOOR,
            "task028_seal_manifest_sha256": sha256_file(T28 / "sealed_v3/SEAL_MANIFEST.sha256"),
            "task028_result_seal_sha256": sha256_file(T28_RESULTS / "RESULT_SEAL.sha256"),
            "task029_version": "CPTAC-UCEC frozen external-replication design",
            "task029_master_seed": 20260714,
            "task029_result_seal_sha256": sha256_file(T29_RESULTS / "RESULT_SEAL.sha256"),
            "seal_entries_verified": int(len(seal_df)),
            "seal_mismatches": int((~seal_df["match"]).sum()),
        },
        "wrapper_randomness": {
            "master_seed": MASTER_SEED_WRAPPER,
            "statement": "No new randomness. Regeneration reused frozen TASK-028 master seed and deterministic sub-seed rule.",
            "sub_seed_rule": "sha256('20260713:<step_id>')[:8] -> int",
        },
        "parameters": {
            "targets": [t for t, _, _ in TARGETS],
            "contrasts": {"C1": "[-1,-1,-1,+3] normalized to [-1/3,-1/3,-1/3,+1]", "C2": "[+1,+1,-2,0] normalized to [+1/2,+1/2,-1,0]"},
            "primary_covariates": ["M4", "purity_CPE", "composition"],
            "official_nopurity_covariates": ["M4", "composition"],
            "common_sample_nopurity_covariates": ["M4", "composition"],
            "common_sample_nopurity_label": "diagnostic only; not a frozen sensitivity or verdict-changing model",
            "category_logic": "sealed TASK-028 F1 omnibus gate + |d|>=0.50 + bootstrap CI excludes 0 + frozen gated full-family F2 BH q<=0.05; C1 equivalence per sealed rule",
            "decomposition": "exact symmetric Shapley for d=b/s",
            "legacy_pointwise_scenario": pointwise_scenario,
            "binding_complete_B2_12_scenario": scenario,
            "complete_B2_12_summary_path": str(GENE_LOO_SUMMARY),
        },
        "reproduction": {
            "full_primary_and_nopurity_rows_compared": int(len(regen_audit)),
            "all_fields_match_tolerance_1e-12": bool(regen_audit["all_fields_match_tolerance_1e-12"].all()),
            "maximum_numeric_abs_delta": float(regen_audit["numeric_max_abs_delta"].max()),
            "claim": "numeric regeneration; no byte-reproduction claim",
        },
        "software": {
            "executable": sys.executable, "python": sys.version,
            "platform": platform.platform(), "numpy": np.__version__,
            "pandas": pd.__version__, "scipy": stats.__version__ if hasattr(stats, "__version__") else __import__("scipy").__version__,
            "pip_freeze": pip_freeze,
            "thread_environment": {k: os.environ.get(k) for k in ["OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS"]},
        },
        "scripts": {str(OUT / "run_task030_audit.py"): sha256_file(OUT / "run_task030_audit.py")},
        "commands": commands,
        "execution_notes": "Current-main isolated-worktree verification rerun from sealed TASK-028/TASK-029 inputs.",
        "input_checksums": {str(p): sha256_file(p) for p in input_paths},
        "output_checksum_scope": "All durable analysis artifacts. Excludes self-referential manifest/checksum file, transient tee run.log, and Python bytecode cache.",
        "output_checksums_excluding_self_referential_and_transient_files": {str(p): sha256_file(p) for p in output_paths},
        "deviations": "Cycle-2 narrative correction integrates the complete B2.12 universal gene-LOO gate; scientific parameters and result bytes are unchanged.",
        "attestation": "No retuning, no new target selection, no changes to frozen families/categories/contrasts/score orientation/edges/covariates/verdicts.",
    }
    manifest_path = OUT / "REPRODUCIBILITY_MANIFEST.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="ascii")

    final_files = sorted([
        p for p in OUT.rglob("*")
        if p.is_file()
        and "__pycache__" not in p.parts
        and p.name not in {"OUTPUT_SHA256SUMS.txt", "run.log"}
    ])
    seal_lines = [
        "# TASK-030 output checksums",
        "# Excludes OUTPUT_SHA256SUMS.txt (self), transient tee run.log, and Python bytecode cache.",
    ]
    for p in final_files:
        seal_lines.append(f"{sha256_file(p)}  {p.relative_to(OUT)}")
    (OUT / "OUTPUT_SHA256SUMS.txt").write_text("\n".join(seal_lines) + "\n", encoding="ascii")

    print(f"[done] scenario={scenario}", flush=True)
    print(f"[done] outputs={OUT}", flush=True)


if __name__ == "__main__":
    main()
