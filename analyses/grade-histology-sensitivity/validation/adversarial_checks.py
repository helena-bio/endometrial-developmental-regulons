#!/usr/bin/env python3
"""Adversarial checks for matched roster, taxonomy, and case-deletion fragility."""

import csv
import json
from collections import Counter

import numpy as np
import pandas as pd

import independent_reconstruction as independent


def frozen_taxonomy(base_d, adjusted_d):
    return independent.classify(base_d, adjusted_d, False)


def main():
    order = json.loads(independent.ORDER.read_text(encoding="ascii"))
    covariates = pd.read_csv(independent.COVARIATES, sep="\t", dtype={"patient_barcode": str})
    archive = np.load(independent.SCORES, allow_pickle=True)
    outcomes = {
        target: np.asarray(archive[f"M3primary__{target}"], dtype=float)
        for target in independent.TARGETS
    }
    archive.close()
    linkage = pd.read_csv(independent.LEDGER, sep="\t", dtype=str, keep_default_na=False)
    tcga = linkage[linkage["cohort"] == "TCGA"].drop_duplicates("patient_id").set_index("patient_id").reindex(order)
    histology = tcga["histology_proposed"].to_numpy(dtype=str)
    grade = tcga["grade_proposed"].to_numpy(dtype=str)
    subtype = covariates["subtype"].to_numpy(dtype=str)
    rosters = {
        "primary_cpe": covariates["cpe_complete_case"].astype(str).str.lower().isin(["1", "true"]).to_numpy(),
        "no_purity": np.ones(len(order), dtype=bool),
    }
    covariate_names = {
        "primary_cpe": ("M4_prolif", "purity_CPE", "composition"),
        "no_purity": ("M4_prolif", "composition"),
    }
    node_definitions = {
        "histology_matched": (
            np.isin(histology, ["endometrioid", "non_endometrioid"]),
            (histology == "non_endometrioid").astype(float),
        ),
        "endometrioid_grade_matched": (
            (histology == "endometrioid") & np.isin(grade, ["low_grade", "high_grade"]),
            (grade == "low_grade").astype(float),
        ),
    }
    rows = []
    for base_name, base_roster in rosters.items():
        for node_name, (clinical_roster, clinical_all) in node_definitions.items():
            indices = np.flatnonzero(base_roster & clinical_roster)
            node_subtype = subtype[indices]
            continuous = {
                name: covariates[name].to_numpy(float)[indices]
                for name in covariate_names[base_name]
            }
            clinical = clinical_all[indices]
            for target in independent.TARGETS:
                y = outcomes[target][indices]
                base = independent.ols(y, independent.matrix(node_subtype, continuous, None)[0])
                adjusted = independent.ols(y, independent.matrix(node_subtype, continuous, clinical)[0])
                original_taxonomy = frozen_taxonomy(base["d"], adjusted["d"])
                loo_taxonomies = Counter()
                loo_base_d = []
                loo_adjusted_d = []
                loo_attenuation = []
                loo_percent = []
                for remove in range(len(indices)):
                    keep = np.arange(len(indices)) != remove
                    cov = {name: values[keep] for name, values in continuous.items()}
                    base_loo = independent.ols(
                        y[keep], independent.matrix(node_subtype[keep], cov, None)[0]
                    )
                    adjusted_loo = independent.ols(
                        y[keep], independent.matrix(node_subtype[keep], cov, clinical[keep])[0]
                    )
                    loo_taxonomies[frozen_taxonomy(base_loo["d"], adjusted_loo["d"])] += 1
                    loo_base_d.append(base_loo["d"])
                    loo_adjusted_d.append(adjusted_loo["d"])
                    attenuation = abs(base_loo["d"]) - abs(adjusted_loo["d"])
                    loo_attenuation.append(attenuation)
                    loo_percent.append(100 * attenuation / abs(base_loo["d"]))
                rows.append({
                    "base_specification": base_name,
                    "node": node_name,
                    "target": target,
                    "n": len(indices),
                    "original_taxonomy": original_taxonomy,
                    "leave_one_out_taxonomy_counts": dict(loo_taxonomies),
                    "leave_one_out_taxonomy_change_n": len(indices) - loo_taxonomies[original_taxonomy],
                    "leave_one_out_base_sign_change_n": sum(np.sign(value) != np.sign(base["d"]) for value in loo_base_d),
                    "leave_one_out_adjusted_sign_change_n": sum(np.sign(value) != np.sign(adjusted["d"]) for value in loo_adjusted_d),
                    "leave_one_out_magnitude_attenuation_min": min(loo_attenuation),
                    "leave_one_out_magnitude_attenuation_max": max(loo_attenuation),
                    "leave_one_out_percent_attenuation_min": min(loo_percent),
                    "leave_one_out_percent_attenuation_max": max(loo_percent),
                })

    model_rows = list(csv.DictReader(
        open(independent.PRODUCER / "MODEL_RESULTS.tsv", encoding="ascii"), delimiter="\t"
    ))
    decomposition_rows = list(csv.DictReader(
        open(independent.PRODUCER / "MATCHED_DECOMPOSITIONS.tsv", encoding="ascii"), delimiter="\t"
    ))
    report = {
        "leave_one_out": rows,
        "leave_one_out_any_c2_sign_change": any(
            row["leave_one_out_base_sign_change_n"] or row["leave_one_out_adjusted_sign_change_n"]
            for row in rows
        ),
        "leave_one_out_taxonomy_change_rows": [
            {
                "base_specification": row["base_specification"],
                "node": row["node"],
                "target": row["target"],
                "count": row["leave_one_out_taxonomy_change_n"],
                "taxonomies": row["leave_one_out_taxonomy_counts"],
            }
            for row in rows if row["leave_one_out_taxonomy_change_n"]
        ],
        "all_model_c2_estimates_negative": all(float(row["d"]) < 0 for row in model_rows),
        "producer_nodes": sorted(set(row["node"] for row in model_rows)),
        "producer_cohorts": sorted(set(row["cohort"] for row in model_rows)),
        "producer_contrasts": sorted(set(row["contrast"] for row in model_rows)),
        "producer_taxonomies": dict(Counter(row["taxonomy"] for row in decomposition_rows)),
        "q_value_columns": sorted({
            field for row in model_rows + decomposition_rows for field in row
            if field.lower().startswith("q") or "fdr" in field.lower()
        }),
        "threshold_margin_primary": [
            {
                "base_specification": row["base_specification"],
                "node": row["node"],
                "target": row["target"],
                "magnitude_attenuation": float(row["magnitude_attenuation"]),
                "percent_attenuation": float(row["percent_attenuation"]),
                "margin_to_0_10": 0.10 - float(row["magnitude_attenuation"]),
                "margin_to_20_percent": 20 - float(row["percent_attenuation"]),
                "would_pass_relaxed_0_08_and_15_percent": (
                    float(row["magnitude_attenuation"]) >= 0.08
                    and float(row["percent_attenuation"]) >= 15
                ),
            }
            for row in decomposition_rows
            if row["target"] in ("GATA2", "SOX9")
        ],
    }
    (independent.OUT / "ADVERSARIAL_CHECKS.json").write_text(
        json.dumps(
            report,
            indent=2,
            sort_keys=True,
            default=lambda value: int(value) if isinstance(value, np.integer) else float(value),
        ) + "\n",
        encoding="ascii",
    )
    print(json.dumps({
        "any_sign_change": report["leave_one_out_any_c2_sign_change"],
        "taxonomy_change_rows": report["leave_one_out_taxonomy_change_rows"],
        "all_model_c2_estimates_negative": report["all_model_c2_estimates_negative"],
        "q_value_columns": report["q_value_columns"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
