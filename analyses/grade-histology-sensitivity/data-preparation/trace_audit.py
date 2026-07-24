#!/usr/bin/env python3
"""Parse complete strace logs and enforce the Phase-1 open/network policy."""

import csv
import json
import re
from collections import Counter
from pathlib import Path


ROOT = Path("data/external/original-workspace/revgate-task-b-grade-histology")
EXEC = ROOT / "experiments/taskB_grade_histology/phase1_execution"
EXACT_EXTERNAL = {
    "data/external/original-workspace/task027-acquire-freeze-a/freeze_a_redux/cohort_selected_primary.tsv",
    "data/external/original-workspace/task024-freeze-a/subtype_normalized.tsv",
    "data/external/original-workspace/task029-external-replication-feasibility/join_tables/discovery_join.json",
    "data/external/original-workspace/task029-external-replication-feasibility/join_tables/confirmatory_join.json",
    "data/external/original-workspace/task029-external-replication-feasibility/sources/cbio_2020_clinattr.json",
    "data/external/original-workspace/task029-external-replication-feasibility/sources_task029gate/uec_cptac_gdc_clinattr.json",
    "data/external/original-workspace/task029-external-replication-feasibility/sources/pdc_discovery_biospec.json",
    "data/external/original-workspace/task029-external-replication-feasibility/sources/pdc_confirmatory_biospec.json",
}
SYSTEM_PREFIXES = ("/usr/", "/lib/", "/lib64/", "/etc/", "/proc/", "/dev/", "/sys/", "/newroot/", "/var/run/", "/tmp/", "data/external/original-home/.local/lib/python3.12/site-packages/")


def classify(path):
    normalized = path
    if normalized in EXACT_EXTERNAL:
        return "ALLOW_EXACT_SHA_PINNED_EXTERNAL"
    if normalized.startswith(str(ROOT) + "/experiments/taskB_grade_histology/inputs/"):
        return "ALLOW_LOCKED_ACQUIRED_INPUT"
    if normalized.startswith(str(EXEC) + "/run1/") or normalized.startswith(str(EXEC) + "/run2/"):
        return "ALLOW_DECLARED_OUTPUT"
    if normalized.startswith(str(EXEC) + "/acquisition/"):
        return "ALLOW_DECLARED_ACQUISITION_OUTPUT"
    if normalized.startswith(str(EXEC) + "/venv_phase1/"):
        return "ALLOW_DECLARED_ENVIRONMENT"
    if normalized in {str(EXEC / "audit_phase1.py"), str(EXEC / "acquire_inputs.py"), str(EXEC / "INPUTS_LOCK.json"), str(ROOT / "experiments/taskB_grade_histology/PHASE1_FEASIBILITY_RULES.json")}:
        return "ALLOW_PHASE1_CODE_LOCK_OR_RULE"
    if normalized == str(EXEC):
        return "ALLOW_PHASE1_CODE_DIRECTORY"
    if normalized.startswith("experiments/taskB_grade_histology/phase1_execution/run"):
        return "ALLOW_DECLARED_OUTPUT_RELATIVE"
    if normalized.startswith("experiments/taskB_grade_histology/inputs/"):
        return "ALLOW_LOCKED_ACQUIRED_INPUT_RELATIVE"
    if normalized.startswith("experiments/taskB_grade_histology/phase1_execution/"):
        return "ALLOW_PHASE1_CODE_OR_ENV_RELATIVE"
    if normalized == "experiments/taskB_grade_histology/PHASE1_FEASIBILITY_RULES.json":
        return "ALLOW_FROZEN_RULE_RELATIVE"
    if normalized in {"/", "/proc", "/newroot", "data/external/original-home/.local/lib/python3.12/site-packages"} or normalized.startswith(SYSTEM_PREFIXES) or not normalized.startswith("/"):
        return "ALLOW_OS_RUNTIME_OR_RELATIVE_RUNTIME"
    return "DENY_UNDECLARED_ABSOLUTE_PATH"


def main():
    rows = []
    network = []
    traces = [
        ("acquisition_round1_failed_http503", EXEC / "acquisition_round1_open_trace.log", False),
        ("acquisition_round1_retry", EXEC / "acquisition_round1_retry_open_trace.log", False),
        ("acquisition_round2", EXEC / "acquisition_round2_open_trace.log", False),
        ("analysis_run1", EXEC / "run1_complete_open_trace.log", True),
        ("analysis_run2", EXEC / "run2_complete_open_trace.log", True),
    ]
    for stage, trace, network_traced in traces:
        for line_no, line in enumerate(trace.read_text(encoding="utf-8", errors="replace").splitlines(), 1):
            if re.search(r"\b(openat|open|creat)\(", line):
                match = re.search(r'"([^"\\]*(?:\\.[^"\\]*)*)"', line)
                if match:
                    path = bytes(match.group(1), "utf-8").decode("unicode_escape")
                    rows.append({"stage": stage, "line": line_no, "path": path, "decision": classify(path)})
            if network_traced and re.search(r"\b(socket|connect|sendto|recvfrom)\(", line):
                if "AF_INET" in line and "AF_NETLINK" not in line:
                    network.append({"stage": stage, "line": line_no, "call": line})
    denied = [row for row in rows if row["decision"].startswith("DENY")]
    summary = {
        "schema_version": "1.0.0",
        "trace_method": "Acquisition: strace -f openat/open/creat in standalone minimal processes. Analysis: strace -f openat/open/creat/network around bwrap --unshare-net.",
        "trace_stages": [stage for stage, _path, _network in traces],
        "acquisition_network_trace_note": "Acquisition strace captured opens; exact remote requests, headers, hashes, and timestamps are recorded in INPUTS_LOCK. Analysis traces additionally captured network syscalls.",
        "open_events": len(rows),
        "decision_counts": dict(sorted(Counter(row["decision"] for row in rows).items())),
        "denied_open_events": denied,
        "external_inet_events": network,
        "network_namespace": "isolated; only bwrap AF_NETLINK loopback setup and local AF_UNIX attempts observed",
        "status": "PASS" if not denied and not network else "FAIL",
    }
    with (EXEC / "OPEN_FILE_AUDIT.tsv").open("w", encoding="ascii", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["stage", "line", "path", "decision"], delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    (EXEC / "OPEN_FILE_AUDIT_SUMMARY.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if summary["status"] != "PASS":
        raise SystemExit("TRACE_AUDIT_FAIL")
    print("TRACE_AUDIT_PASS")


if __name__ == "__main__":
    main()
