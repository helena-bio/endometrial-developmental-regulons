#!/usr/bin/env python3
"""Independent TASK B Phase-2 cycle-6 reconstruction.

This implementation reads only the frozen inputs and specification.  It does
not import or execute producer run_phase2.py or any producer result-generating
function.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path("data/external/original-workspace/revgate-task-b-grade-histology")
OUT = ROOT / "experiments/taskB_grade_histology/phase2_critic_cycle6"
PRODUCER = ROOT / "experiments/taskB_grade_histology/phase2_execution_cycle6/run1"
SCORES = Path("data/external/original-workspace/task028-freeze-b-draft/execution/intermediate/scores_v3.npz")
COVARIATES = Path("data/external/original-workspace/task028-freeze-b-draft/execution/intermediate/covariates_v3.tsv")
ORDER = Path("data/external/original-workspace/task028-freeze-b-draft/execution/intermediate/patient_order_v3.json")
LEDGER = ROOT / "experiments/taskB_grade_histology/phase1_execution/run1/LINKAGE_LEDGER.tsv"

TARGETS = ("GATA2", "SOX9", "HOXA9", "WT1", "PAX8", "LHX1")
SUBTYPES = ("POLE", "MMRd", "NSMP", "p53abn")
MASTER_SEED = 20260713
BOOTSTRAPS = 2000
PERMUTATIONS = 2000


def derived_seed(step_id: str) -> int:
    payload = f"{MASTER_SEED}|TASKB_PHASE2|{step_id}".encode("ascii")
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "big", signed=False)


def matrix(subtype: np.ndarray, continuous: dict[str, np.ndarray],
           clinical: np.ndarray | None) -> tuple[np.ndarray, list[str]]:
    columns = [np.ones(subtype.size)]
    names = ["intercept"]
    for level in SUBTYPES[1:]:
        columns.append((subtype == level).astype(float))
        names.append("subtype_" + level)
    for name, values in continuous.items():
        columns.append(np.asarray(values, dtype=float))
        names.append(name)
    if clinical is not None:
        columns.append(np.asarray(clinical, dtype=float))
        names.append("clinical_binary")
    return np.column_stack(columns), names


def ols(y: np.ndarray, x: np.ndarray) -> dict:
    n, p = x.shape
    beta_vector, _, _, singular = np.linalg.lstsq(x, y, rcond=None)
    tolerance = np.finfo(float).eps * max(n, p) * singular[0]
    rank = int(np.sum(singular > tolerance))
    if rank != p:
        raise np.linalg.LinAlgError("rank deficient")
    residual = y - x @ beta_vector
    residual_df = n - p
    residual_sd = math.sqrt(float(residual @ residual) / residual_df)
    xtx_inverse = np.linalg.inv(x.T @ x)
    contrast = np.zeros(p)
    contrast[1] = 0.5
    contrast[2] = -1.0
    coefficient = float(contrast @ beta_vector)
    se = math.sqrt(float(residual_sd**2 * (contrast @ xtx_inverse @ contrast)))
    t_value = coefficient / se
    raw_p = float(2.0 * stats.t.sf(abs(t_value), residual_df))
    t_critical = float(stats.t.ppf(0.975, residual_df))
    return {
        "vector": beta_vector,
        "residual": residual,
        "inverse": xtx_inverse,
        "singular": singular,
        "rank": rank,
        "rank_tolerance": tolerance,
        "residual_df": residual_df,
        "residual_sd": residual_sd,
        "coefficient": coefficient,
        "se": se,
        "t": t_value,
        "raw_p": raw_p,
        "ci_lo": coefficient - t_critical * se,
        "ci_hi": coefficient + t_critical * se,
        "d": coefficient / residual_sd,
    }


def resample_single(y: np.ndarray, subtype: np.ndarray,
                    continuous: dict[str, np.ndarray],
                    clinical: np.ndarray | None, seed: int) -> dict:
    rng = np.random.default_rng(seed)
    coefficients: list[float] = []
    ds: list[float] = []
    n = y.size
    for _ in range(BOOTSTRAPS):
        take = rng.integers(0, n, size=n)
        if np.unique(subtype[take]).size != 4:
            continue
        try:
            xb, _ = matrix(
                subtype[take],
                {name: values[take] for name, values in continuous.items()},
                None if clinical is None else clinical[take],
            )
            fit = ols(y[take], xb)
        except (np.linalg.LinAlgError, FloatingPointError):
            continue
        coefficients.append(fit["coefficient"])
        ds.append(fit["d"])
    return {
        "valid": len(ds),
        "coefficient_ci": np.percentile(coefficients, [2.5, 97.5]),
        "d_ci": np.percentile(ds, [2.5, 97.5]),
    }


def permutation_probability(y: np.ndarray, subtype: np.ndarray,
                            continuous: dict[str, np.ndarray],
                            clinical: np.ndarray | None, observed_d: float,
                            seed: int) -> float:
    rng = np.random.default_rng(seed)
    exceedances = 0
    for _ in range(PERMUTATIONS):
        permuted = rng.permutation(subtype)
        xp, _ = matrix(permuted, continuous, clinical)
        try:
            permuted_fit = ols(y, xp)
        except np.linalg.LinAlgError:
            continue
        exceedances += abs(permuted_fit["d"]) >= abs(observed_d) - 1e-12
    return (1 + exceedances) / (PERMUTATIONS + 1)


def paired_resample(y: np.ndarray, subtype: np.ndarray,
                    continuous: dict[str, np.ndarray], clinical: np.ndarray,
                    seed: int) -> dict:
    rng = np.random.default_rng(seed)
    draws: list[list[float]] = []
    n = y.size
    for _ in range(BOOTSTRAPS):
        take = rng.integers(0, n, size=n)
        if np.unique(subtype[take]).size != 4:
            continue
        cov = {name: values[take] for name, values in continuous.items()}
        try:
            base = ols(y[take], matrix(subtype[take], cov, None)[0])
            adjusted = ols(y[take], matrix(subtype[take], cov, clinical[take])[0])
        except (np.linalg.LinAlgError, FloatingPointError):
            continue
        signed_delta = adjusted["d"] - base["d"]
        magnitude_attenuation = abs(base["d"]) - abs(adjusted["d"])
        percentage = (
            100 * magnitude_attenuation / abs(base["d"])
            if abs(base["d"]) >= 0.05 and np.sign(base["d"]) == np.sign(adjusted["d"])
            else np.nan
        )
        scale = base["coefficient"] * (
            1 / adjusted["residual_sd"] - 1 / base["residual_sd"]
        )
        coefficient = (
            adjusted["coefficient"] - base["coefficient"]
        ) / adjusted["residual_sd"]
        draws.append([
            signed_delta,
            abs(signed_delta),
            magnitude_attenuation,
            percentage,
            adjusted["coefficient"] - base["coefficient"],
            scale,
            coefficient,
        ])
    array = np.asarray(draws)
    return {
        "valid": array.shape[0],
        "percentiles": np.nanpercentile(array, [2.5, 97.5], axis=0),
    }


def classify(base_d: float, adjusted_d: float, unstable: bool) -> str:
    if unstable:
        return "unstable"
    if abs(base_d) < 0.05:
        return "base_near_null_unclassifiable"
    if np.sign(base_d) != np.sign(adjusted_d):
        return "sign_reversal"
    if abs(adjusted_d) > abs(base_d) + 1e-12:
        return "amplified"
    if abs(adjusted_d) < 0.05:
        return "near_null_or_disappears_descriptively"
    attenuation = abs(base_d) - abs(adjusted_d)
    percent = 100 * attenuation / abs(base_d)
    if attenuation >= 0.10 and percent >= 20:
        return "same_direction_materially_attenuated"
    return "largely_retained"


def diagnostics(y: np.ndarray, x: np.ndarray, subtype: np.ndarray,
                clinical: np.ndarray | None, node: str, fit: dict) -> tuple[dict, list[dict]]:
    n, p = x.shape
    predictors = x[:, 1:]
    standardized = (predictors - predictors.mean(axis=0)) / predictors.std(axis=0)
    singular = np.linalg.svd(standardized, compute_uv=False)
    condition = float(singular[0] / singular[-1])
    correlation = np.corrcoef(predictors, rowvar=False)
    correlation_inverse = np.linalg.pinv(correlation)
    max_vif = float(np.diag(correlation_inverse).max())
    subtype_indices = [0, 1, 2]
    other_indices = list(range(3, correlation.shape[0]))
    determinant = float(np.linalg.det(correlation))
    if not other_indices:
        subtype_gvif = 1.0
    elif determinant <= 0:
        subtype_gvif = float("inf")
    else:
        det_subtype = float(np.linalg.det(correlation[np.ix_(subtype_indices, subtype_indices)]))
        det_other = float(np.linalg.det(correlation[np.ix_(other_indices, other_indices)]))
        subtype_gvif = float((det_subtype * det_other / determinant) ** (1 / 6))

    inverse = fit["inverse"]
    hat = np.sum(x * (x @ inverse), axis=1)
    clipped_hat = np.clip(hat, 0, 1 - 1e-14)
    residual = fit["residual"]
    residual_sd = fit["residual_sd"]
    residual_df = fit["residual_df"]
    studentized_internal = residual / (residual_sd * np.sqrt(1 - clipped_hat))
    studentized_external = studentized_internal * np.sqrt(
        (residual_df - 1)
        / np.maximum(1e-15, residual_df - studentized_internal**2)
    )
    cooks = (
        (residual**2 / (p * residual_sd**2))
        * clipped_hat
        / (1 - clipped_hat) ** 2
    )
    dfbeta = (
        (inverse @ x.T).T
        * residual[:, None]
        / (
            residual_sd
            * np.sqrt(np.diag(inverse))[None, :]
            * (1 - clipped_hat)[:, None]
        )
    )
    max_abs_dfbeta = np.max(np.abs(dfbeta), axis=1)
    flagged = np.where(
        (cooks > 4 / residual_df)
        | (max_abs_dfbeta > 2 / math.sqrt(n))
        | (np.abs(studentized_external) > 3)
    )[0]

    maximum_leave_one_out_change = 0.0
    direction_change = False
    for index in flagged:
        keep = np.arange(n) != index
        try:
            leave_one_out = ols(y[keep], x[keep])
        except np.linalg.LinAlgError:
            maximum_leave_one_out_change = float("inf")
            direction_change = True
            continue
        maximum_leave_one_out_change = max(
            maximum_leave_one_out_change,
            abs(leave_one_out["d"] - fit["d"]),
        )
        direction_change |= np.sign(leave_one_out["d"]) != np.sign(fit["d"])

    reasons: list[str] = []
    if n < max(80, 10 * p):
        reasons.append("support_n")
    if residual_df < 30:
        reasons.append("residual_df")
    if set(subtype) != set(SUBTYPES):
        reasons.append("missing_subtype")
    subtype_counts = {level: int(np.sum(subtype == level)) for level in SUBTYPES}
    required_subtype_n = 5 if node == "histology_matched" else 10 if node == "endometrioid_grade_matched" else 0
    if required_subtype_n and min(subtype_counts.values()) < required_subtype_n:
        reasons.append("subtype_support")

    cramers_v = None
    max_clinical_correlation = None
    clinical_counts: dict[str, int] = {}
    if clinical is not None:
        clinical_counts = {str(level): int(np.sum(clinical == level)) for level in (0.0, 1.0)}
        if min(clinical_counts.values()) < 20:
            reasons.append("clinical_overall_support")
        table = np.array([
            [np.sum((subtype == level) & (clinical == value)) for value in (0.0, 1.0)]
            for level in SUBTYPES
        ])
        if np.sum(np.sum(table >= 5, axis=0) >= 2) < 2:
            reasons.append("clinical_cross_subtype_support")
        for column in range(2):
            if table[:, column].sum() and table[:, column].max() / table[:, column].sum() >= 0.9:
                reasons.append("clinical_level_concentration")
        chi_square = stats.chi2_contingency(table, correction=False)[0]
        cramers_v = float(math.sqrt(chi_square / (n * min(table.shape[0] - 1, table.shape[1] - 1))))
        max_clinical_correlation = float(np.nanmax(np.abs(correlation[:3, -1])))
        frequency_ratio = max(clinical_counts.values()) / max(1, min(clinical_counts.values()))
        percent_unique = 200 / n
        if frequency_ratio > 19 and percent_unique <= 10:
            reasons.append("near_zero_variance")
        if cramers_v > 0.9:
            reasons.append("cramers_v")
        if max_clinical_correlation >= 0.95:
            reasons.append("dummy_correlation")
    if fit["rank"] < p:
        reasons.append("rank")
    if condition > 100:
        reasons.append("condition_number_not_estimable")
    elif condition > 30:
        reasons.append("condition_number_unstable")
    if max_vif > 10 or subtype_gvif > 10:
        reasons.append("vif_fail")
    elif max_vif > 5 or subtype_gvif > 5:
        reasons.append("vif_unstable")
    fraction_high_hat = float(np.mean(hat > 3 * p / n))
    if float(hat.max()) >= 0.5 or fraction_high_hat > 0.05:
        reasons.append("leverage")
    if direction_change:
        reasons.append("single_case_direction_change")
    if maximum_leave_one_out_change >= 0.2:
        reasons.append("single_case_delta_d")
    if not reasons:
        status = "PASS"
    elif any(reason in reasons for reason in ("rank", "condition_number_not_estimable", "vif_fail")):
        status = "NOT_ESTIMABLE"
    else:
        status = "UNSTABLE_NOT_INTERPRETED"

    summary = {
        "n": n,
        "p": p,
        "rank": fit["rank"],
        "residual_df": residual_df,
        "condition_number": condition,
        "max_vif": max_vif,
        "subtype_gvif_adjusted": subtype_gvif,
        "max_hat": float(hat.max()),
        "fraction_hat_gt_3p_n": fraction_high_hat,
        "n_hat_gt_2p_n": int(np.sum(hat > 2 * p / n)),
        "n_cooks_flags": int(np.sum(cooks > 4 / residual_df)),
        "n_dfbeta_flags": int(np.sum(max_abs_dfbeta > 2 / math.sqrt(n))),
        "n_studentized_flags": int(np.sum(np.abs(studentized_external) > 3)),
        "n_influence_union": int(flagged.size),
        "max_single_case_abs_d_change": maximum_leave_one_out_change,
        "single_case_direction_change": bool(direction_change),
        "cramers_v": cramers_v,
        "max_clinical_dummy_abs_correlation": max_clinical_correlation,
        "status": status,
        "reasons": ";".join(sorted(set(reasons))) or "none",
        "subtype_counts": json.dumps(subtype_counts, sort_keys=True, separators=(",", ":")),
        "clinical_counts": json.dumps(clinical_counts, sort_keys=True, separators=(",", ":")),
    }
    influence = [
        {
            "row_index": int(index),
            "hat": float(hat[index]),
            "cooks_d": float(cooks[index]),
            "max_abs_dfbeta": float(max_abs_dfbeta[index]),
            "externally_studentized": float(studentized_external[index]),
        }
        for index in flagged
    ]
    return summary, influence


def read_producer_tsv(name: str) -> list[dict[str, str]]:
    with open(PRODUCER / name, newline="", encoding="ascii") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(name: str, rows: list[dict]) -> None:
    if not rows:
        return
    with open(OUT / name, "w", newline="", encoding="ascii") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def compare_rows(independent: list[dict], producer: list[dict[str, str]],
                 keys: tuple[str, ...], name: str) -> dict:
    left = {tuple(str(row[key]) for key in keys): row for row in independent}
    right = {tuple(row[key] for key in keys): row for row in producer}
    result = {
        "table": name,
        "independent_rows": len(independent),
        "producer_rows": len(producer),
        "missing_keys": sorted(set(right) - set(left)),
        "extra_keys": sorted(set(left) - set(right)),
        "numeric_fields_compared": 0,
        "text_fields_compared": 0,
        "mismatch_count": 0,
        "max_absolute_difference": 0.0,
        "max_relative_difference": 0.0,
        "max_absolute_location": None,
        "max_relative_location": None,
    }
    for key in sorted(set(left) & set(right)):
        a = left[key]
        b = right[key]
        for field, bv in b.items():
            if field in keys:
                continue
            av = a.get(field)
            if bv == "NA" or bv == "":
                result["text_fields_compared"] += 1
                expected = None if bv in ("NA", "") else bv
                if av is not None and str(av) not in ("", "NA"):
                    result["mismatch_count"] += 1
                continue
            try:
                bf = float(bv)
                af = float(av)
            except (TypeError, ValueError):
                result["text_fields_compared"] += 1
                if str(av) != bv:
                    result["mismatch_count"] += 1
                continue
            result["numeric_fields_compared"] += 1
            absolute = abs(af - bf)
            relative = absolute / max(abs(bf), 1e-300)
            if absolute > result["max_absolute_difference"]:
                result["max_absolute_difference"] = absolute
                result["max_absolute_location"] = [*key, field, af, bf]
            if relative > result["max_relative_difference"]:
                result["max_relative_difference"] = relative
                result["max_relative_location"] = [*key, field, af, bf]
            if not np.isclose(af, bf, rtol=1e-9, atol=1e-11):
                result["mismatch_count"] += 1
    return result


def main() -> None:
    patient_order = json.loads(ORDER.read_text(encoding="ascii"))
    covariates = pd.read_csv(COVARIATES, sep="\t", dtype={"patient_barcode": str})
    assert covariates["patient_barcode"].tolist() == patient_order
    archive = np.load(SCORES, allow_pickle=True)
    assert archive["patient_order"].tolist() == patient_order
    outcomes = {
        target: np.asarray(archive[f"M3primary__{target}"], dtype=float)
        for target in TARGETS
    }
    archive.close()

    linkage = pd.read_csv(LEDGER, sep="\t", dtype=str, keep_default_na=False)
    tcga = linkage[linkage["cohort"] == "TCGA"].drop_duplicates("patient_id").set_index("patient_id")
    tcga = tcga.reindex(patient_order)
    assert tcga.shape[0] == len(patient_order) and not tcga.index.hasnans
    histology = tcga["histology_proposed"].to_numpy(dtype=str)
    grade = tcga["grade_proposed"].to_numpy(dtype=str)
    subtype = covariates["subtype"].to_numpy(dtype=str)

    base_definitions = {
        "primary_cpe": (
            ("M4_prolif", "purity_CPE", "composition"),
            covariates["cpe_complete_case"].astype(str).str.lower().isin(["1", "true"]).to_numpy(),
        ),
        "no_purity": (
            ("M4_prolif", "composition"),
            np.ones(len(patient_order), dtype=bool),
        ),
    }
    nodes = {
        "full_base": (lambda roster: roster.copy(), None, ("base",)),
        "histology_matched": (
            lambda roster: roster & np.isin(histology, ["endometrioid", "non_endometrioid"]),
            (histology == "non_endometrioid").astype(float),
            ("base_refit", "base_plus_binary_histology"),
        ),
        "endometrioid_grade_matched": (
            lambda roster: roster & (histology == "endometrioid") & np.isin(grade, ["low_grade", "high_grade"]),
            (grade == "low_grade").astype(float),
            ("base_refit", "base_plus_binary_grade"),
        ),
    }

    models: list[dict] = []
    diagnostic_rows: list[dict] = []
    influence_rows: list[dict] = []
    decompositions: list[dict] = []
    seeds: list[dict] = []
    stored_fits: dict[tuple[str, str, str, str], tuple[dict, dict]] = {}
    roster_rows: list[dict] = []

    for base_name, (covariate_names, roster) in base_definitions.items():
        for node_name, (filter_function, clinical_all, variants) in nodes.items():
            selected = filter_function(roster)
            indices = np.flatnonzero(selected)
            node_subtype = subtype[indices]
            continuous = {
                name: covariates[name].to_numpy(dtype=float)[indices]
                for name in covariate_names
            }
            node_clinical = None if clinical_all is None else clinical_all[indices]
            for index in indices:
                roster_rows.append({
                    "base_specification": base_name,
                    "node": node_name,
                    "row_index": int(index),
                    "patient_barcode": patient_order[index],
                    "subtype": subtype[index],
                })
            for target in TARGETS:
                y = outcomes[target][indices]
                for variant in variants:
                    use_clinical = node_clinical if variant.startswith("base_plus") else None
                    x, _ = matrix(node_subtype, continuous, use_clinical)
                    fit = ols(y, x)
                    step_prefix = "|".join(("TCGA-UCEC", base_name, node_name, target, "C2"))
                    bootstrap_seed = derived_seed(step_prefix + "|bootstrap|" + variant)
                    permutation_seed = derived_seed(step_prefix + "|permutation|" + variant)
                    seeds.extend([
                        {"model_id": f"{base_name}|{node_name}|{variant}|{target}", "procedure": "bootstrap", "seed": bootstrap_seed},
                        {"model_id": f"{base_name}|{node_name}|{variant}|{target}", "procedure": "permutation", "seed": permutation_seed},
                    ])
                    bootstrap = resample_single(
                        y, node_subtype, continuous, use_clinical, bootstrap_seed
                    )
                    permutation_p = permutation_probability(
                        y, node_subtype, continuous, use_clinical, fit["d"], permutation_seed
                    )
                    diagnostic, influence = diagnostics(
                        y, x, node_subtype, use_clinical, node_name, fit
                    )
                    model_id = f"{base_name}|{node_name}|{variant}|{target}"
                    stored_fits[(base_name, node_name, variant, target)] = (fit, diagnostic)
                    models.append({
                        "cohort": "TCGA-UCEC",
                        "base_specification": base_name,
                        "node": node_name,
                        "model": variant,
                        "target": target,
                        "contrast": "C2",
                        "n": len(indices),
                        "p_design": x.shape[1],
                        "coefficient": fit["coefficient"],
                        "se": fit["se"],
                        "t": fit["t"],
                        "residual_df": fit["residual_df"],
                        "ci_lo": fit["ci_lo"],
                        "ci_hi": fit["ci_hi"],
                        "raw_p": fit["raw_p"],
                        "direction": "negative" if fit["coefficient"] < 0 else "positive" if fit["coefficient"] > 0 else "zero",
                        "residual_sd": fit["residual_sd"],
                        "d": fit["d"],
                        "d_ci_lo": bootstrap["d_ci"][0],
                        "d_ci_hi": bootstrap["d_ci"][1],
                        "coefficient_boot_ci_lo": bootstrap["coefficient_ci"][0],
                        "coefficient_boot_ci_hi": bootstrap["coefficient_ci"][1],
                        "bootstrap_attempts": BOOTSTRAPS,
                        "bootstrap_valid": bootstrap["valid"],
                        "diagnostic_permutations": PERMUTATIONS,
                        "permutation_p": permutation_p,
                        "interpretability": diagnostic["status"],
                    })
                    diagnostic_rows.append({"model_id": model_id, **diagnostic})
                    for influence_row in influence:
                        influence_rows.append({
                            "model_id": model_id,
                            "patient_barcode": patient_order[indices[influence_row["row_index"]]],
                            **influence_row,
                        })

            if node_name != "full_base":
                adjusted_variant = variants[1]
                for target in TARGETS:
                    base_fit, base_diagnostic = stored_fits[(base_name, node_name, "base_refit", target)]
                    adjusted_fit, adjusted_diagnostic = stored_fits[(base_name, node_name, adjusted_variant, target)]
                    delta = adjusted_fit["d"] - base_fit["d"]
                    attenuation = abs(base_fit["d"]) - abs(adjusted_fit["d"])
                    scale_contribution = base_fit["coefficient"] * (
                        1 / adjusted_fit["residual_sd"] - 1 / base_fit["residual_sd"]
                    )
                    coefficient_contribution = (
                        adjusted_fit["coefficient"] - base_fit["coefficient"]
                    ) / adjusted_fit["residual_sd"]
                    decomposition_error = abs(
                        scale_contribution + coefficient_contribution - delta
                    )
                    paired_seed = derived_seed(
                        "|".join(("TCGA-UCEC", base_name, node_name, target, "C2", "paired_bootstrap"))
                    )
                    seeds.append({
                        "model_id": f"{base_name}|{node_name}|{target}",
                        "procedure": "paired_bootstrap",
                        "seed": paired_seed,
                    })
                    paired = paired_resample(
                        outcomes[target][indices],
                        node_subtype,
                        continuous,
                        node_clinical,
                        paired_seed,
                    )
                    unstable = (
                        base_diagnostic["status"] != "PASS"
                        or adjusted_diagnostic["status"] != "PASS"
                        or paired["valid"] < 1900
                    )
                    percent = (
                        None
                        if abs(base_fit["d"]) < 0.05
                        or np.sign(base_fit["d"]) != np.sign(adjusted_fit["d"])
                        else 100 * attenuation / abs(base_fit["d"])
                    )
                    decompositions.append({
                        "cohort": "TCGA-UCEC",
                        "base_specification": base_name,
                        "node": node_name,
                        "target": target,
                        "n": len(indices),
                        "base_beta": base_fit["coefficient"],
                        "adjusted_beta": adjusted_fit["coefficient"],
                        "base_residual_sd": base_fit["residual_sd"],
                        "adjusted_residual_sd": adjusted_fit["residual_sd"],
                        "base_d": base_fit["d"],
                        "adjusted_d": adjusted_fit["d"],
                        "signed_delta_d": delta,
                        "absolute_delta_d": abs(delta),
                        "magnitude_attenuation": attenuation,
                        "percent_attenuation": percent,
                        "beta_change": adjusted_fit["coefficient"] - base_fit["coefficient"],
                        "residual_scale_contribution": scale_contribution,
                        "coefficient_contribution": coefficient_contribution,
                        "decomposition_error": decomposition_error,
                        "paired_bootstrap_attempts": BOOTSTRAPS,
                        "paired_bootstrap_valid": paired["valid"],
                        "signed_delta_d_ci_lo": paired["percentiles"][0, 0],
                        "signed_delta_d_ci_hi": paired["percentiles"][1, 0],
                        "taxonomy": classify(base_fit["d"], adjusted_fit["d"], unstable),
                        "interpretability": "UNSTABLE_NOT_INTERPRETED" if unstable else "PASS",
                    })

    write_tsv("INDEPENDENT_MODEL_RESULTS.tsv", models)
    write_tsv("INDEPENDENT_MODEL_DIAGNOSTICS.tsv", diagnostic_rows)
    write_tsv("INDEPENDENT_INFLUENCE_RECORDS.tsv", influence_rows)
    write_tsv("INDEPENDENT_MATCHED_DECOMPOSITIONS.tsv", decompositions)
    write_tsv("INDEPENDENT_SEEDS.tsv", seeds)
    write_tsv("INDEPENDENT_ROSTERS.tsv", roster_rows)

    comparisons = [
        compare_rows(
            models,
            read_producer_tsv("MODEL_RESULTS.tsv"),
            ("cohort", "base_specification", "node", "model", "target", "contrast"),
            "MODEL_RESULTS.tsv",
        ),
        compare_rows(
            diagnostic_rows,
            read_producer_tsv("MODEL_DIAGNOSTICS.tsv"),
            ("model_id",),
            "MODEL_DIAGNOSTICS.tsv",
        ),
        compare_rows(
            influence_rows,
            read_producer_tsv("INFLUENCE_RECORDS.tsv"),
            ("model_id", "patient_barcode", "row_index"),
            "INFLUENCE_RECORDS.tsv",
        ),
        compare_rows(
            decompositions,
            read_producer_tsv("MATCHED_DECOMPOSITIONS.tsv"),
            ("cohort", "base_specification", "node", "target"),
            "MATCHED_DECOMPOSITIONS.tsv",
        ),
    ]
    summary = {
        "implementation_independence": "No producer script imported, executed, copied, or called.",
        "model_rows": len(models),
        "diagnostic_rows": len(diagnostic_rows),
        "influence_rows": len(influence_rows),
        "decomposition_rows": len(decompositions),
        "seed_rows": len(seeds),
        "roster_rows": len(roster_rows),
        "taxonomy_counts": dict(pd.Series([row["taxonomy"] for row in decompositions]).value_counts()),
        "diagnostic_status_counts": dict(pd.Series([row["status"] for row in diagnostic_rows]).value_counts()),
        "bootstrap_valid_min": min(row["bootstrap_valid"] for row in models),
        "paired_bootstrap_valid_min": min(row["paired_bootstrap_valid"] for row in decompositions),
        "maximum_decomposition_error": max(row["decomposition_error"] for row in decompositions),
        "maximum_condition_number": max(row["condition_number"] for row in diagnostic_rows),
        "maximum_vif": max(row["max_vif"] for row in diagnostic_rows),
        "maximum_adjusted_subtype_gvif": max(row["subtype_gvif_adjusted"] for row in diagnostic_rows),
        "maximum_hat": max(row["max_hat"] for row in diagnostic_rows),
        "maximum_single_case_abs_d_change": max(row["max_single_case_abs_d_change"] for row in diagnostic_rows),
        "single_case_direction_change_count": sum(row["single_case_direction_change"] for row in diagnostic_rows),
        "comparisons": comparisons,
    }
    (OUT / "INDEPENDENT_REPRODUCTION_SUMMARY.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True, default=lambda x: int(x) if isinstance(x, np.integer) else float(x)) + "\n",
        encoding="ascii",
    )
    print(json.dumps(summary, sort_keys=True, default=lambda x: int(x) if isinstance(x, np.integer) else float(x)))


if __name__ == "__main__":
    main()
