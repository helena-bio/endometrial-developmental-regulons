#!/usr/bin/env python3
"""Frozen deny-pattern scan for the offline clinical audit program."""

import json
import re
from pathlib import Path


HERE = Path(__file__).resolve().parent
TARGET = HERE / "audit_phase1.py"
PATTERNS = {
    "forbidden_external_worktree": r"(?i)(revgate-task-a-perclass|task030-|data/external/original-workspace/revgate/)",
    "forbidden_scientific_artifact_keyword": r"(?i)(regulon|score_array|primary_results|per_target|coefficient|effect_size|p_value|q_value|confidence_interval)",
    "forbidden_outcome_keyword": r"(?i)(overall_survival|progression|recurrence|days_to_death|vital_status)",
    "forbidden_target_scientific_pair": r"(?i)(gata2|sox9|hoxa9|wt1|pax8|lhx1).*(result|score|coefficient|effect|figure)",
    "binary_array_load": r"(?i)(numpy|np)\.load\s*\([^\n]*(npy|npz)",
    "parquet_read": r"(?i)(read_parquet|parquet\.ParquetFile)",
    "recursive_discovery": r"(?i)(rglob\s*\(|os\.walk\s*\(|glob\s*\(.*\*\*)",
    "network_module": r"(?i)(urllib|requests|httpx|socket|aiohttp)",
}


def main():
    text = TARGET.read_text(encoding="utf-8")
    findings = {name: [m.group(0) for m in re.finditer(pattern, text)] for name, pattern in PATTERNS.items()}
    findings = {name: values for name, values in findings.items() if values}
    report = {"script": str(TARGET), "sha256_checked_separately_in_manifest": True, "patterns": PATTERNS, "findings": findings, "status": "PASS" if not findings else "REJECT"}
    (HERE / "STATIC_SCAN_REPORT.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if findings:
        raise SystemExit("STATIC_SCAN_REJECT")
    print("STATIC_SCAN_PASS")


if __name__ == "__main__":
    main()
