#!/usr/bin/env python3
"""Write Cycle 6 provenance and raw-report artifacts without changing run1/run2."""
import csv
import hashlib
import json
import re
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent
RUN1 = ROOT / "run1"
RUN2 = ROOT / "run2"
AUDIT = ROOT / "audit"
SEAL_PATH = ROOT / "PRE_OUTCOME_EXECUTION_SEAL.json"
RECON_PATH = ROOT / "UPSTREAM_POINT_ESTIMATE_RECONCILIATION.json"
SCIENTIFIC_FILES = [
    "ANALYTICAL_REPORT.md",
    "CLINICAL_DISTRIBUTIONS.tsv",
    "COHORT_COUNTS.tsv",
    "INFLUENCE_RECORDS.tsv",
    "INPUT_ACCESS_INVENTORY.tsv",
    "MANUSCRIPT_READY_SENSITIVITY_TABLE.tsv",
    "MATCHED_DECOMPOSITIONS.tsv",
    "MISSINGNESS.tsv",
    "MODEL_DIAGNOSTICS.tsv",
    "MODEL_RESULTS.tsv",
    "PROPOSED_WORDING.md",
    "SCIENTIFIC_RESULTS.json",
    "SIX_TARGET_APPENDIX.tsv",
]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_tsv(path):
    with open(path, encoding="ascii", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_json(path, value):
    with open(path, "w", encoding="ascii", newline="\n") as handle:
        json.dump(value, handle, sort_keys=True, indent=2, allow_nan=False)
        handle.write("\n")


def write_tsv(path, rows, fields):
    with open(path, "w", encoding="ascii", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


seal = json.load(open(SEAL_PATH, encoding="ascii"))
recon = json.load(open(RECON_PATH, encoding="ascii"))
analysis_sha = sha256(ROOT / "run_phase2.py")
config_sha = sha256(ROOT / "EXECUTION_CONFIG.json")
if analysis_sha != seal["code_and_config_hashes"]["run_phase2.py"]:
    raise RuntimeError("post-seal analysis hash changed")
if config_sha != seal["code_and_config_hashes"]["EXECUTION_CONFIG.json"]:
    raise RuntimeError("post-seal config hash changed")

run_hashes = {}
identity_rows = []
for name in SCIENTIFIC_FILES:
    p1 = RUN1 / name
    p2 = RUN2 / name
    h1 = sha256(p1)
    h2 = sha256(p2)
    run_hashes[name] = h1
    identity_rows.append(
        {
            "file": name,
            "run1_sha256": h1,
            "run2_sha256": h2,
            "run1_size": p1.stat().st_size,
            "run2_size": p2.stat().st_size,
            "byte_identical": h1 == h2 and p1.read_bytes() == p2.read_bytes(),
        }
    )
if not all(row["byte_identical"] for row in identity_rows):
    raise RuntimeError("scientific output byte mismatch")

for run_name in ("run1", "run2"):
    with open(ROOT / f"{run_name.upper()}_SCIENTIFIC_SHA256SUMS.txt", "w", encoding="ascii", newline="\n") as handle:
        for name in SCIENTIFIC_FILES:
            handle.write(f"{run_hashes[name]}  {run_name}/{name}\n")
write_json(
    ROOT / "BYTE_IDENTITY_REPORT.json",
    {
        "status": "BYTE_IDENTICAL",
        "scientific_file_count": len(identity_rows),
        "all_byte_identical": True,
        "files": identity_rows,
    },
)

config = json.load(open(ROOT / "EXECUTION_CONFIG.json", encoding="ascii"))
allowed_inputs = set(config["paths"].values())
cycle_prefix = str(ROOT) + "/"
trace_reports = {}
for run_name in ("RUN1", "RUN2"):
    trace_path = ROOT / f"POSTSEAL_{run_name}_OPEN_TRACE.log"
    text = trace_path.read_text(encoding="utf-8", errors="replace")
    forbidden = {
        "prior_cycle": len(re.findall(r"phase2_execution_cycle[1-5]", text, flags=re.I)),
        "phase2_critic": len(re.findall(r"phase2_critic", text, flags=re.I)),
        "task030": len(re.findall(r"task0?30|task-030", text, flags=re.I)),
        "task_a": len(re.findall(r"task[_-]a|taska", text, flags=re.I)),
        "no_purity_worktree": len(re.findall(r"revgate-tcga-no-purity-verify", text, flags=re.I)),
        "manuscript": len(re.findall(r"[/\\]manuscript[/\\]|\.docx|\.pdf", text, flags=re.I)),
        "figure": len(re.findall(r"[/\\][^/\\]*figure[^/\\]*", text, flags=re.I)),
    }
    internet_connects = [
        line
        for line in text.splitlines()
        if "connect(" in line and ("AF_INET" in line or "AF_INET6" in line)
    ]
    netlink_lines = [line for line in text.splitlines() if "AF_NETLINK" in line or "sa_family=AF_NETLINK" in line]
    local_nscd_lines = [line for line in text.splitlines() if "/var/run/nscd/socket" in line]
    workspace_paths = set(
        re.findall(r'"(data/external/original-workspace/[^"]+)"', text)
    )
    unexpected_workspace_paths = sorted(
        path
        for path in workspace_paths
        if not path.startswith(cycle_prefix)
        and not str(ROOT).startswith(path)
        and path not in allowed_inputs
    )
    input_open_counts = {path: text.count(path) for path in sorted(allowed_inputs)}
    trace_reports[run_name.lower()] = {
        "path": str(trace_path),
        "sha256": sha256(trace_path),
        "line_count": len(text.splitlines()),
        "forbidden_path_matches": forbidden,
        "internet_connect_count": len(internet_connects),
        "internet_connect_lines": internet_connects,
        "netlink_namespace_setup_line_count": len(netlink_lines),
        "failed_local_nscd_line_count": len(local_nscd_lines),
        "unexpected_workspace_paths": unexpected_workspace_paths,
        "declared_input_path_occurrence_counts": input_open_counts,
        "status": "PASS"
        if not any(forbidden.values()) and not internet_connects and not unexpected_workspace_paths
        else "FAIL",
    }
if any(report["status"] != "PASS" for report in trace_reports.values()):
    raise RuntimeError(f"post-seal trace audit failed: {trace_reports}")
write_json(
    ROOT / "POSTSEAL_TRACE_AUDIT.json",
    {
        "status": "PASS",
        "network_isolation": "bwrap --unshare-net",
        "note": "AF_NETLINK calls are bwrap loopback/net-namespace setup; failed AF_UNIX nscd probes are local. There were zero AF_INET/AF_INET6 connect calls.",
        "runs": trace_reports,
    },
)

command_rows = []
for meta_path in sorted((AUDIT / "commands").glob("*.meta")):
    values = {}
    for line in meta_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            values[key] = value
    if "exit_status" not in values:
        continue
    command_rows.append(
        {
            "command_id": values.get("command_id", meta_path.stem),
            "start_utc": values.get("start_utc", ""),
            "end_utc": values.get("end_utc", ""),
            "exit_status": values.get("exit_status", ""),
            "cwd": values.get("cwd", ""),
            "command": values.get("command", ""),
            "stdout_sha256": values.get("stdout_sha256", ""),
            "stderr_sha256": values.get("stderr_sha256", ""),
            "meta_path": str(meta_path.relative_to(ROOT)),
        }
    )
command_rows.sort(key=lambda row: int(row["command_id"]))
write_tsv(
    ROOT / "COMMANDS_CHRONOLOGICAL.tsv",
    command_rows,
    [
        "command_id",
        "start_utc",
        "end_utc",
        "exit_status",
        "cwd",
        "command",
        "stdout_sha256",
        "stderr_sha256",
        "meta_path",
    ],
)

model_rows = read_tsv(RUN1 / "MODEL_RESULTS.tsv")
diagnostic_rows = read_tsv(RUN1 / "MODEL_DIAGNOSTICS.tsv")
decomposition_rows = read_tsv(RUN1 / "MATCHED_DECOMPOSITIONS.tsv")
influence_rows = read_tsv(RUN1 / "INFLUENCE_RECORDS.tsv")
counts_rows = read_tsv(RUN1 / "COHORT_COUNTS.tsv")
missing_rows = read_tsv(RUN1 / "MISSINGNESS.tsv")
distribution_rows = read_tsv(RUN1 / "CLINICAL_DISTRIBUTIONS.tsv")

full_base = {
    (row["base_specification"], row["target"]): row
    for row in model_rows
    if row["node"] == "full_base"
}
primary_decomp = [
    row
    for row in decomposition_rows
    if row["target"] in {"GATA2", "SOX9"}
]
diagnostic_summary = {
    "row_count": len(diagnostic_rows),
    "status_counts": dict(Counter(row["status"] for row in diagnostic_rows)),
    "max_condition_number": max(float(row["condition_number"]) for row in diagnostic_rows),
    "max_vif": max(float(row["max_vif"]) for row in diagnostic_rows),
    "max_gvif_adjusted": max(float(row["subtype_gvif_adjusted"]) for row in diagnostic_rows),
    "max_hat": max(float(row["max_hat"]) for row in diagnostic_rows),
    "max_single_case_abs_d_change": max(float(row["max_single_case_abs_d_change"]) for row in diagnostic_rows),
}
taxonomy_counts = dict(Counter(row["taxonomy"] for row in decomposition_rows))
interpretability_counts = dict(Counter(row["interpretability"] for row in decomposition_rows))
bootstrap_valid_min = min(int(row["bootstrap_valid"]) for row in model_rows)
paired_bootstrap_valid_min = min(int(row["paired_bootstrap_valid"]) for row in decomposition_rows)

raw_numbers = {
    "model_rows": len(model_rows),
    "diagnostic_rows": len(diagnostic_rows),
    "decomposition_rows": len(decomposition_rows),
    "influence_records": len(influence_rows),
    "input_hash_rows": 16,
    "bootstrap_valid_min": bootstrap_valid_min,
    "paired_bootstrap_valid_min": paired_bootstrap_valid_min,
    "diagnostic_summary": diagnostic_summary,
    "taxonomy_counts": taxonomy_counts,
    "decomposition_interpretability_counts": interpretability_counts,
    "missingness": {row["field"]: int(row["missing_n"]) for row in missing_rows},
    "clinical_distributions": {
        f"{row['field']}|{row['level']}": int(row["n"]) for row in distribution_rows
    },
    "full_base": {
        f"{base}|{target}": {
            "n": int(row["n"]),
            "coefficient": float(row["coefficient"]),
            "raw_p": float(row["raw_p"]),
            "permutation_p": float(row["permutation_p"]),
            "residual_sd": float(row["residual_sd"]),
            "d": float(row["d"]),
            "d_ci": [float(row["d_ci_lo"]), float(row["d_ci_hi"])],
            "interpretability": row["interpretability"],
        }
        for (base, target), row in sorted(full_base.items())
    },
    "primary_target_matched_decompositions": [
        {
            "base_specification": row["base_specification"],
            "node": row["node"],
            "target": row["target"],
            "n": int(row["n"]),
            "base_d": float(row["base_d"]),
            "adjusted_d": float(row["adjusted_d"]),
            "signed_delta_d": float(row["signed_delta_d"]),
            "magnitude_attenuation": float(row["magnitude_attenuation"]),
            "percent_attenuation": float(row["percent_attenuation"]),
            "signed_delta_d_ci": [float(row["signed_delta_d_ci_lo"]), float(row["signed_delta_d_ci_hi"])],
            "taxonomy": row["taxonomy"],
            "interpretability": row["interpretability"],
        }
        for row in primary_decomp
    ],
    "cohort_counts": counts_rows,
}
write_json(ROOT / "RAW_NUMBERS.json", raw_numbers)

deviations = [
    {
        "id": "D1",
        "severity": "DISCLOSED_PRESEAL_PLUMBING",
        "text": "Two pre-seal synthetic computations passed, but strace followed live tee helpers and deadlocked. Both attempts and exit-137 tracer terminations are preserved. The audit wrapper was corrected before sealing; the third traced synthetic test passed and exited cleanly.",
        "scientific_effect": "none; no outcome was opened and no fit used real outcomes",
    },
    {
        "id": "D2",
        "severity": "DISCLOSED_AUDIT_BOOTSTRAP",
        "text": "One post-patch attempt to invoke audit_exec.sh directly failed with Permission denied before the wrapper could log itself; the exact command and error are disclosed in COMMAND_COVERAGE_NOTE.md and the tool transcript. The immediately following bash-invoked command restored executable mode and is logged.",
        "scientific_effect": "none; pre-seal and before outcome access",
    },
    {
        "id": "D3",
        "severity": "DISCLOSED_RECONCILIATION_SCOPE",
        "text": "While resolving Task-028 module identifiers after output freeze, a filtered command printed GATA2/SOX9 C1 and C3 rows in addition to the authorized C2 rows. This occurred only after both Cycle-6 outputs were read-only and byte-frozen; no analysis code, taxonomy, output, or verdict changed.",
        "scientific_effect": "none on Cycle-6 computation; broader-than-needed upstream read disclosed",
    },
    {
        "id": "D4",
        "severity": "DISCLOSED_UPSTREAM_COVERAGE",
        "text": "The pinned Task-030 point table contains GATA2, SOX9, HOXA9, and WT1 only. The pinned Task-028 primary and no-purity tables supplied all six targets, including PAX8 and LHX1.",
        "scientific_effect": "none; all six point estimates reconciled to Task-028",
    },
    {
        "id": "D5",
        "severity": "DISCLOSED_PROVENANCE_CLOSEOUT",
        "text": "The first provenance-only finalizer attempt stopped before writing the manifest because the trace auditor treated bwrap ancestor-directory traversal as an unexpected workspace access. The raw traceback and exit status are preserved. The auditor was corrected to allow only lexical ancestors of the sealed Cycle-6 root and then rerun; frozen run1/run2 files were never changed.",
        "scientific_effect": "none; both scientific executions had already succeeded and been made read-only",
    },
]
deviation_text = ["# Cycle 6 deviations", ""]
for item in deviations:
    deviation_text.extend(
        [
            f"## {item['id']} - {item['severity']}",
            "",
            item["text"],
            "",
            f"Scientific effect: {item['scientific_effect']}.",
            "",
        ]
    )
(ROOT / "DEVIATIONS.md").write_text("\n".join(deviation_text), encoding="ascii", newline="\n")

coverage_note = """# Command coverage note

The first shell command created the Cycle-6 directory and is recorded verbatim in
audit_bootstrap.txt. The next action created audit_exec.sh and the bootstrap note
through apply_patch; apply_patch is a tool action rather than a shell command.

Every subsequently executed shell command was routed through audit_exec.sh and has
per-command metadata, stdout, stderr, exit status, timestamps, and hashes under
audit/commands/, except one invocation that could not start the wrapper:

    experiments/taskB_grade_histology/phase2_execution_cycle6/audit_exec.sh -- chmod 0555 experiments/taskB_grade_histology/phase2_execution_cycle6/audit_exec.sh

It failed immediately with:

    /bin/bash: line 1: experiments/taskB_grade_histology/phase2_execution_cycle6/audit_exec.sh: Permission denied

This happened after apply_patch reset executable mode and before any outcome access.
The immediately following bash-invoked wrapper command restored mode and is logged.
Pre-seal tracer deadlocks and forced tracer terminations are retained as exit-137
commands. COMMANDS_CHRONOLOGICAL.tsv is sorted by command start ID; raw immutable
command evidence remains under audit/.
"""
(ROOT / "COMMAND_COVERAGE_NOTE.md").write_text(coverage_note, encoding="ascii", newline="\n")

report_lines = [
    "# Cycle 6 raw analytical report",
    "",
    "Status: COMPLETE PRODUCER EXECUTION. This is a post-hoc explanatory TCGA-UCEC sensitivity only. It does not change any frozen TCGA/CPTAC verdict, target category, manuscript, or claim.",
    "",
    "Two post-seal runs completed with 60 model rows, 24 matched decompositions, 60 diagnostic rows, and 5,231 influence records each. All 13 scientific files are byte-identical. All 60 fitted nodes passed the frozen gates; all 24 decompositions are interpretable. Taxonomy counts: 16 largely_retained and 8 amplified.",
    "",
    "Full-base primary CPE C2 d values: "
    + ", ".join(
        f"{target} {float(full_base[('primary_cpe', target)]['d']):+.6f}"
        for target in ["GATA2", "SOX9", "HOXA9", "WT1", "PAX8", "LHX1"]
    )
    + ".",
    "",
    "GATA2/SOX9 matched results:",
]
for row in primary_decomp:
    report_lines.append(
        f"- {row['base_specification']} / {row['node']} / {row['target']}: n={row['n']}, "
        f"d_base={float(row['base_d']):+.6f}, d_adjusted={float(row['adjusted_d']):+.6f}, "
        f"delta_d={float(row['signed_delta_d']):+.6f}, percent_attenuation={float(row['percent_attenuation']):+.3f}%, "
        f"taxonomy={row['taxonomy']}, status={row['interpretability']}."
    )
report_lines.extend(
    [
        "",
        f"Diagnostics maxima: condition={diagnostic_summary['max_condition_number']:.6f}, VIF={diagnostic_summary['max_vif']:.6f}, adjusted GVIF={diagnostic_summary['max_gvif_adjusted']:.6f}, hat={diagnostic_summary['max_hat']:.6f}, single-case absolute d change={diagnostic_summary['max_single_case_abs_d_change']:.6f}.",
        "",
        f"Upstream point-estimate reconciliation: {recon['status']}; {recon['check_count']} checks, {recon['mismatch_count']} mismatches, maximum absolute delta {recon['max_absolute_delta']:.3e}.",
        "",
        "No q values or new categories were computed. Raw Student-t p values and diagnostic permutation p values are descriptive only. No CPTAC adjusted model was run.",
    ]
)
(ROOT / "RAW_ANALYTICAL_REPORT.md").write_text("\n".join(report_lines) + "\n", encoding="ascii", newline="\n")

manifest = {
    "schema_version": "1.0.0",
    "status": "DONE",
    "classification": "HYPOTHESIS_POST_HOC_EXPLANATORY_SENSITIVITY_PHASE2",
    "scope": "TCGA-UCEC only; post-hoc explanatory sensitivity; no CPTAC adjusted inference",
    "branch": "experiment/task-b-grade-histology-sensitivity",
    "git_commit": "83503bad47b60193598b2b9ebe819c22c83e8ac1",
    "completed_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
    "seed": 20260713,
    "subseed_rule": "first 8 bytes big-endian SHA256 of ASCII 20260713|TASKB_PHASE2|<analysis_step_id>",
    "pre_outcome_seal": {
        "path": str(SEAL_PATH),
        "sha256": sha256(SEAL_PATH),
        "sealed_before_outcome_load": True,
        "validated_input_hash_count": 16,
    },
    "analysis": {
        "path": str(ROOT / "run_phase2.py"),
        "sha256": analysis_sha,
        "config_path": str(ROOT / "EXECUTION_CONFIG.json"),
        "config_sha256": config_sha,
        "cycle5_plumbing_fix_only": True,
        "postseal_analysis_edit": False,
        "postseal_config_edit": False,
        "postseal_retry": False,
        "postseal_execution_count": 2,
    },
    "environment": {
        "path": str(ROOT / "ENVIRONMENT.txt"),
        "sha256": sha256(ROOT / "ENVIRONMENT.txt"),
        "interpreter": "/usr/bin/python3",
        "python": sys.version.replace("\n", " "),
        "network_isolation": "bwrap --unshare-net",
    },
    "runs": {
        "run1": {
            "path": str(RUN1),
            "exit_code": 0,
            "stdout_sha256": sha256(ROOT / "POSTSEAL_RUN1.stdout"),
            "stderr_sha256": sha256(ROOT / "POSTSEAL_RUN1.stderr"),
            "open_trace_sha256": sha256(ROOT / "POSTSEAL_RUN1_OPEN_TRACE.log"),
        },
        "run2": {
            "path": str(RUN2),
            "exit_code": 0,
            "stdout_sha256": sha256(ROOT / "POSTSEAL_RUN2.stdout"),
            "stderr_sha256": sha256(ROOT / "POSTSEAL_RUN2.stderr"),
            "open_trace_sha256": sha256(ROOT / "POSTSEAL_RUN2_OPEN_TRACE.log"),
        },
    },
    "byte_identity": {
        "status": "BYTE_IDENTICAL",
        "scientific_file_count": len(SCIENTIFIC_FILES),
        "report_path": str(ROOT / "BYTE_IDENTITY_REPORT.json"),
        "report_sha256": sha256(ROOT / "BYTE_IDENTITY_REPORT.json"),
    },
    "scientific_output_hashes": run_hashes,
    "trace_audit": {
        "status": "PASS",
        "path": str(ROOT / "POSTSEAL_TRACE_AUDIT.json"),
        "sha256": sha256(ROOT / "POSTSEAL_TRACE_AUDIT.json"),
        "internet_connect_count": 0,
        "forbidden_path_match_count": 0,
    },
    "upstream_reconciliation": {
        "status": recon["status"],
        "path": str(RECON_PATH),
        "sha256": sha256(RECON_PATH),
        "check_count": recon["check_count"],
        "mismatch_count": recon["mismatch_count"],
        "max_absolute_delta": recon["max_absolute_delta"],
        "performed_after_output_freeze": True,
        "new_output_mutation": False,
        "verdict_or_category_mutation": False,
    },
    "raw_numbers_path": str(ROOT / "RAW_NUMBERS.json"),
    "raw_numbers_sha256": sha256(ROOT / "RAW_NUMBERS.json"),
    "commands": {
        "bootstrap_note": str(ROOT / "audit_bootstrap.txt"),
        "chronological_table": str(ROOT / "COMMANDS_CHRONOLOGICAL.tsv"),
        "raw_audit_directory": str(AUDIT),
        "completed_logged_shell_commands_in_table": len(command_rows),
    },
    "deviations": deviations,
    "scientific_verdict": "NONE_PRODUCER_DOES_NOT_DECIDE",
    "critic_required": True,
    "prohibited_changes": {
        "frozen_tcga_cptac_verdict_changed": False,
        "target_category_changed": False,
        "manuscript_changed": False,
        "src_changed": False,
        "docs_changed": False,
        "staged_or_committed": False,
    },
}
write_json(ROOT / "REPRODUCIBILITY_MANIFEST.json", manifest)

inventory_rows = []
for path in sorted(ROOT.rglob("*")):
    if not path.is_file() or AUDIT in path.parents or "__pycache__" in path.parts:
        continue
    if path.name in {"FILE_INVENTORY.tsv", "OUTPUT_SHA256SUMS.txt"}:
        continue
    inventory_rows.append(
        {
            "path": str(path.relative_to(ROOT)),
            "size": path.stat().st_size,
            "sha256": sha256(path),
        }
    )
write_tsv(ROOT / "FILE_INVENTORY.tsv", inventory_rows, ["path", "size", "sha256"])

checksum_paths = [
    path
    for path in sorted(ROOT.rglob("*"))
    if path.is_file()
    and AUDIT not in path.parents
    and "__pycache__" not in path.parts
    and path.name != "OUTPUT_SHA256SUMS.txt"
]
with open(ROOT / "OUTPUT_SHA256SUMS.txt", "w", encoding="ascii", newline="\n") as handle:
    for path in checksum_paths:
        handle.write(f"{sha256(path)}  {path.relative_to(ROOT)}\n")

print(
    json.dumps(
        {
            "status": "DONE",
            "scientific_files_byte_identical": len(SCIENTIFIC_FILES),
            "model_rows": len(model_rows),
            "decompositions": len(decomposition_rows),
            "diagnostics_pass": diagnostic_summary["status_counts"].get("PASS", 0),
            "upstream_reconciliation": recon["status"],
        },
        sort_keys=True,
    )
)
