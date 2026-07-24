#!/usr/bin/env python3
"""Independent audit of producer and critic OS-level traces."""

from __future__ import annotations

import csv
import json
import pathlib
import re


ROOT = pathlib.Path("data/external/original-workspace/revgate-task-b-grade-histology")
EXEC = ROOT / "experiments/taskB_grade_histology/phase1_execution"
OUT = ROOT / "experiments/taskB_grade_histology/phase1_critic"
PRODUCER_TRACES = [
    EXEC / "acquisition_round1_open_trace.log",
    EXEC / "acquisition_round1_retry_open_trace.log",
    EXEC / "acquisition_round2_open_trace.log",
    EXEC / "run1_complete_open_trace.log",
    EXEC / "run2_complete_open_trace.log",
]
CRITIC_TRACES = [OUT / "acquire_round1_open_trace.log", OUT / "acquire_round2_open_trace.log", OUT / "independent_analysis_open_trace.log"]
FORBIDDEN = re.compile(r"revgate-task-a-perclass|/task030-|data/external/original-workspace/revgate/|/results/|/manuscript/|/final-manuscript/|B3_CLAIM_MATRIX|scores?[^/]*\.(?:npy|npz)|primary_results|per_target|model[^/]*_results|ANALYTICAL_REPORT|CRITIC_VERDICT|MANUSCRIPT_READY|\.docx|\.pdf|/Figure|/figure", re.I)


def audit_trace(path: pathlib.Path, analysis: bool) -> dict:
    lines = path.read_text(errors="replace").splitlines()
    quoted = []
    external_network = []
    for line in lines:
        if "openat(" in line or "open(" in line or "creat(" in line:
            quoted.extend(re.findall(r'"([^"\\]*(?:\\.[^"\\]*)*)"', line))
        if analysis and ("connect(" in line or "sendto(" in line or "recvfrom(" in line):
            if ("AF_INET" in line or "AF_INET6" in line) and "127.0.0.1" not in line and "::1" not in line and "0.0.0.0" not in line:
                external_network.append(line)
    forbidden = sorted({p for p in quoted if FORBIDDEN.search(p)})
    return {"path": str(path), "line_count": len(lines), "quoted_open_path_count": len(quoted), "forbidden_open_paths": forbidden, "external_network_events": external_network}


def main() -> None:
    producer = [audit_trace(p, i >= 3) for i, p in enumerate(PRODUCER_TRACES)]
    critic = [audit_trace(p, i == 2) for i, p in enumerate(CRITIC_TRACES)]
    with (EXEC / "OPEN_FILE_AUDIT.tsv").open(newline="", encoding="utf-8") as handle:
        parsed = list(csv.DictReader(handle, delimiter="\t"))
    summary = json.loads((EXEC / "OPEN_FILE_AUDIT_SUMMARY.json").read_text())
    static = json.loads((EXEC / "STATIC_SCAN_REPORT.json").read_text())
    lock = json.loads((EXEC / "INPUTS_LOCK.json").read_text())
    commands = (EXEC / "COMMANDS.txt").read_text()
    environment = json.loads((EXEC / "ENVIRONMENT.json").read_text())
    result = {
        "producer_raw_trace_stages": producer,
        "producer_stage_count": len(producer),
        "producer_forbidden_open_count": sum(len(x["forbidden_open_paths"]) for x in producer),
        "producer_external_analysis_network_count": sum(len(x["external_network_events"]) for x in producer),
        "producer_parsed_open_rows": len(parsed),
        "producer_parsed_expected_1325": len(parsed) == 1325,
        "producer_summary_status": summary["status"],
        "producer_summary_denied_count": len(summary["denied_open_events"]),
        "producer_summary_external_inet_count": len(summary["external_inet_events"]),
        "producer_static_scan": static,
        "producer_env_i_commands": commands.count("env -i"),
        "producer_bwrap_unshare_net_commands": commands.count("bwrap --unshare-net"),
        "producer_numpy": environment.get("numpy"),
        "producer_failed_attempts": lock.get("failed_attempts"),
        "critic_raw_traces": critic,
        "critic_forbidden_open_count": sum(len(x["forbidden_open_paths"]) for x in critic),
        "critic_external_analysis_network_count": sum(len(x["external_network_events"]) for x in critic),
        "early_preflight_assessment": "Current PREFLIGHT_SCHEMA_REPORT is sanitized (dynamic identifiers redacted). The disclosed superseded version contained linkage identifiers only; no retained initial bytes exist, so exact prior content cannot be independently proven. Identifiers are permitted linkage metadata, not grade/histology or molecular outcomes.",
        "gdc_case_mapping_assessment": "SHA matches allowlist but top-level schema is a broad mapping definition, not a case crosswalk; producer denied and did not use it.",
        "status": "PASS" if all([sum(len(x["forbidden_open_paths"]) for x in producer) == 0, sum(len(x["external_network_events"]) for x in producer) == 0, len(parsed) == 1325, sum(len(x["forbidden_open_paths"]) for x in critic) == 0, sum(len(x["external_network_events"]) for x in critic) == 0]) else "FAIL",
    }
    (OUT / "FIREWALL_OPEN_TRACE_AUDIT.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="ascii")
    print(json.dumps({k: result[k] for k in ["status", "producer_forbidden_open_count", "producer_external_analysis_network_count", "producer_parsed_open_rows", "critic_forbidden_open_count", "critic_external_analysis_network_count"]}, sort_keys=True))


if __name__ == "__main__":
    main()
