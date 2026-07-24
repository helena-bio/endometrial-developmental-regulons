#!/usr/bin/env python3
"""Point-estimate-only upstream reconciliation after producer-output freeze."""
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent
TARGETS = ["GATA2", "SOX9", "HOXA9", "WT1", "PAX8", "LHX1"]
CURRENT = ROOT / "run1" / "MODEL_RESULTS.tsv"
TASK028_PRIMARY = Path("data/external/original-workspace/task028-freeze-b-draft/execution/results_v3/primary_results.tsv")
TASK028_NOPURITY = Path("data/external/original-workspace/task028-freeze-b-draft/execution/results_v3/sensitivity_SENS_nopurity.tsv")
TASK030 = Path("data/external/original-workspace/task030-tcga-nopurity-audit/results/tcga_primary_vs_nopurity_per_target.tsv")
EXPECTED_HASHES = {
    str(TASK028_PRIMARY): "50df517a55744c12cac1db62a40b123976b1e4dc7efc2b806fe9f6d2a1608f9f",
    str(TASK028_NOPURITY): "be4a89f723cc1e8f62405ea32f4acb6dccf386cf7cb06f26627a6120478c43d7",
    str(TASK030): "f49685bf10cf6f8c8302ddaa81fe4f8c0d60ee36107e7d7f01c89692a52e5399",
}
TOLERANCE = 1e-12


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_tsv(path):
    with open(path, encoding="ascii", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


for path_text, expected in EXPECTED_HASHES.items():
    actual = sha256(Path(path_text))
    if actual != expected:
        raise RuntimeError(f"upstream hash mismatch: {path_text}")

current_rows = [
    row
    for row in read_tsv(CURRENT)
    if row["node"] == "full_base" and row["target"] in TARGETS and row["contrast"] == "C2"
]
current = {(row["base_specification"], row["target"]): row for row in current_rows}
if len(current) != 12:
    raise RuntimeError(f"expected 12 current full-base rows, found {len(current)}")

checks = []


def add_check(source, target, base_spec, quantity, upstream_value):
    current_field = {"coefficient": "coefficient", "residual_sd": "residual_sd", "d": "d"}[quantity]
    current_value = float(current[(base_spec, target)][current_field])
    upstream_value = float(upstream_value)
    delta = current_value - upstream_value
    checks.append(
        {
            "source": source,
            "target": target,
            "contrast": "C2",
            "base_specification": base_spec,
            "quantity": quantity,
            "current_value": format(current_value, ".17g"),
            "upstream_value": format(upstream_value, ".17g"),
            "signed_delta": format(delta, ".17g"),
            "absolute_delta": format(abs(delta), ".17g"),
            "tolerance": format(TOLERANCE, ".17g"),
            "status": "MATCH" if abs(delta) <= TOLERANCE else "MISMATCH",
        }
    )


task028_primary_rows = {
    row["module"].removeprefix("M3_"): row
    for row in read_tsv(TASK028_PRIMARY)
    if row["module"] in {f"M3_{target}" for target in TARGETS} and row["contrast"] == "C2"
}
if set(task028_primary_rows) != set(TARGETS):
    raise RuntimeError("Task-028 primary C2 target set mismatch")
for target in TARGETS:
    row = task028_primary_rows[target]
    add_check("TASK028_PRIMARY", target, "primary_cpe", "coefficient", row["estimate"])
    add_check("TASK028_PRIMARY", target, "primary_cpe", "residual_sd", row["sigma_resid"])
    add_check("TASK028_PRIMARY", target, "primary_cpe", "d", row["d"])

task028_nopurity_rows = {
    row["module"].removeprefix("M3_"): row
    for row in read_tsv(TASK028_NOPURITY)
    if row["module"] in {f"M3_{target}" for target in TARGETS} and row["contrast"] == "C2"
}
if set(task028_nopurity_rows) != set(TARGETS):
    raise RuntimeError("Task-028 no-purity C2 target set mismatch")
for target in TARGETS:
    row = task028_nopurity_rows[target]
    add_check("TASK028_NOPURITY", target, "no_purity", "coefficient", row["estimate"])
    add_check("TASK028_NOPURITY", target, "no_purity", "d", row["d"])

task030_rows = {
    row["target"]: row
    for row in read_tsv(TASK030)
    if row["target"] in TARGETS and row["contrast"] == "C2"
}
if set(task030_rows) != {"GATA2", "SOX9", "HOXA9", "WT1"}:
    raise RuntimeError("Task-030 pinned table has an unexpected C2 target set")
for target, row in sorted(task030_rows.items()):
    for quantity, primary_field, nopurity_field in (
        ("coefficient", "primary_contrast_coefficient", "nopurity_contrast_coefficient"),
        ("residual_sd", "primary_residual_sd", "nopurity_residual_sd"),
        ("d", "primary_d", "nopurity_d"),
    ):
        add_check("TASK030_AUDIT", target, "primary_cpe", quantity, row[primary_field])
        add_check("TASK030_AUDIT", target, "no_purity", quantity, row[nopurity_field])

fields = [
    "source",
    "target",
    "contrast",
    "base_specification",
    "quantity",
    "current_value",
    "upstream_value",
    "signed_delta",
    "absolute_delta",
    "tolerance",
    "status",
]
with open(ROOT / "UPSTREAM_POINT_ESTIMATE_RECONCILIATION.tsv", "w", encoding="ascii", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(checks)

max_delta = max(float(row["absolute_delta"]) for row in checks)
summary = {
    "status": "MATCH" if all(row["status"] == "MATCH" for row in checks) else "MISMATCH",
    "timing": "performed only after both producer runs were byte-compared, checksummed, and made read-only",
    "scope": "point estimates only; C2 full_base rows; no category or verdict mutation",
    "current_file": str(CURRENT),
    "current_file_sha256": sha256(CURRENT),
    "upstream_files": [
        {"path": path, "sha256": expected}
        for path, expected in sorted(EXPECTED_HASHES.items())
    ],
    "check_count": len(checks),
    "mismatch_count": sum(row["status"] != "MATCH" for row in checks),
    "max_absolute_delta": max_delta,
    "tolerance": TOLERANCE,
    "task030_target_scope": ["GATA2", "SOX9", "HOXA9", "WT1"],
    "task030_scope_note": "The pinned Task-030 table contains four C2 targets; the pinned Task-028 primary and no-purity tables provide all six frozen targets.",
    "new_output_mutation": False,
    "upstream_verdict_or_category_mutation": False,
}
with open(ROOT / "UPSTREAM_POINT_ESTIMATE_RECONCILIATION.json", "w", encoding="ascii", newline="\n") as handle:
    json.dump(summary, handle, sort_keys=True, indent=2, allow_nan=False)
    handle.write("\n")

access_rows = [
    {
        "utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "action": "READ_FILTERED_C2_POINT_ROWS_AFTER_OUTPUT_FREEZE",
        "path": path,
        "sha256": expected,
        "scope": "GATA2,SOX9,HOXA9,WT1,PAX8,LHX1 where available; C2 only",
    }
    for path, expected in sorted(EXPECTED_HASHES.items())
]
with open(ROOT / "UPSTREAM_RECONCILIATION_ACCESS_LOG.tsv", "w", encoding="ascii", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=["utc", "action", "path", "sha256", "scope"],
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(access_rows)

print(json.dumps(summary, sort_keys=True))
