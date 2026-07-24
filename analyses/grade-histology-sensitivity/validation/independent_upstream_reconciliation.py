#!/usr/bin/env python3
"""Independently reconcile cycle-6 full-base point estimates to pinned upstream."""

import csv
import hashlib
import json
from pathlib import Path


ROOT = Path("data/external/original-workspace/revgate-task-b-grade-histology")
CRITIC = ROOT / "experiments/taskB_grade_histology/phase2_critic_cycle6"
CURRENT = CRITIC / "INDEPENDENT_MODEL_RESULTS.tsv"
PRIMARY = Path("data/external/original-workspace/task028-freeze-b-draft/execution/results_v3/primary_results.tsv")
NOPURITY = Path("data/external/original-workspace/task028-freeze-b-draft/execution/results_v3/sensitivity_SENS_nopurity.tsv")
TASK030 = Path("data/external/original-workspace/task030-tcga-nopurity-audit/results/tcga_primary_vs_nopurity_per_target.tsv")
TARGETS = ("GATA2", "SOX9", "HOXA9", "WT1", "PAX8", "LHX1")
EXPECTED = {
    PRIMARY: "50df517a55744c12cac1db62a40b123976b1e4dc7efc2b806fe9f6d2a1608f9f",
    NOPURITY: "be4a89f723cc1e8f62405ea32f4acb6dccf386cf7cb06f26627a6120478c43d7",
    TASK030: "f49685bf10cf6f8c8302ddaa81fe4f8c0d60ee36107e7d7f01c89692a52e5399",
}


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(path: Path):
    with open(path, newline="", encoding="ascii") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    for path, expected in EXPECTED.items():
        assert digest(path) == expected
    current = {
        (row["base_specification"], row["target"]): row
        for row in rows(CURRENT)
        if row["node"] == "full_base" and row["model"] == "base"
    }
    assert len(current) == 12
    comparisons = []

    def check(source, base, target, quantity, upstream):
        field = {"coefficient": "coefficient", "residual_sd": "residual_sd", "d": "d"}[quantity]
        observed = float(current[(base, target)][field])
        reference = float(upstream)
        delta = observed - reference
        comparisons.append({
            "source": source,
            "base_specification": base,
            "target": target,
            "quantity": quantity,
            "independent_cycle6": observed,
            "upstream": reference,
            "signed_delta": delta,
            "absolute_delta": abs(delta),
            "status": "MATCH" if abs(delta) <= 1e-12 else "MISMATCH",
        })

    primary = {
        row["module"].removeprefix("M3_"): row
        for row in rows(PRIMARY)
        if row["contrast"] == "C2" and row["module"] in {f"M3_{target}" for target in TARGETS}
    }
    assert set(primary) == set(TARGETS)
    for target in TARGETS:
        check("TASK028_PRIMARY", "primary_cpe", target, "coefficient", primary[target]["estimate"])
        check("TASK028_PRIMARY", "primary_cpe", target, "residual_sd", primary[target]["sigma_resid"])
        check("TASK028_PRIMARY", "primary_cpe", target, "d", primary[target]["d"])

    no_purity = {
        row["module"].removeprefix("M3_"): row
        for row in rows(NOPURITY)
        if row["contrast"] == "C2" and row["module"] in {f"M3_{target}" for target in TARGETS}
    }
    assert set(no_purity) == set(TARGETS)
    for target in TARGETS:
        check("TASK028_NOPURITY", "no_purity", target, "coefficient", no_purity[target]["estimate"])
        check("TASK028_NOPURITY", "no_purity", target, "d", no_purity[target]["d"])

    task030 = {
        row["target"]: row for row in rows(TASK030)
        if row["contrast"] == "C2" and row["target"] in TARGETS
    }
    assert set(task030) == {"GATA2", "SOX9", "HOXA9", "WT1"}
    for target, row in sorted(task030.items()):
        for quantity, primary_field, no_purity_field in (
            ("coefficient", "primary_contrast_coefficient", "nopurity_contrast_coefficient"),
            ("residual_sd", "primary_residual_sd", "nopurity_residual_sd"),
            ("d", "primary_d", "nopurity_d"),
        ):
            check("TASK030_AUDIT", "primary_cpe", target, quantity, row[primary_field])
            check("TASK030_AUDIT", "no_purity", target, quantity, row[no_purity_field])

    with open(CRITIC / "INDEPENDENT_UPSTREAM_RECONCILIATION.tsv", "w", newline="", encoding="ascii") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(comparisons[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(comparisons)
    summary = {
        "check_count": len(comparisons),
        "mismatch_count": sum(row["status"] != "MATCH" for row in comparisons),
        "max_absolute_delta": max(row["absolute_delta"] for row in comparisons),
        "task028_targets": list(TARGETS),
        "task030_targets": sorted(task030),
        "input_hashes": {str(path): digest(path) for path in EXPECTED},
    }
    (CRITIC / "INDEPENDENT_UPSTREAM_RECONCILIATION.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="ascii"
    )
    print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()
