#!/usr/bin/env python3
"""Create the Cycle 6 pre-outcome execution seal without loading outcomes."""
import csv
import difflib
import hashlib
import importlib.util
import json
import platform
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
import statsmodels

ROOT = Path(__file__).resolve().parent
CONFIG = ROOT / "EXECUTION_CONFIG.json"
ANALYSIS = ROOT / "run_phase2.py"
PRIOR_ANALYSIS = ROOT.parent / "phase2_execution_cycle5" / "run_phase2.py"
TRACE = ROOT / "PRESEAL_SYNTHETIC_OPEN_TRACE_SUCCESS.log"
SYNTHETIC_STDOUT = ROOT / "PRESEAL_SYNTHETIC_SUCCESS.stdout"
SYNTHETIC_STDERR = ROOT / "PRESEAL_SYNTHETIC_SUCCESS.stderr"
FREEZE_ROOT = ROOT.parent / "phase2_freeze"
SEAL = ROOT / "PRE_OUTCOME_EXECUTION_SEAL.json"


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def write_json(path, value):
    with open(path, "w", encoding="ascii", newline="\n") as handle:
        json.dump(value, handle, sort_keys=True, indent=2, allow_nan=False)
        handle.write("\n")


if SEAL.exists():
    raise RuntimeError("pre-outcome seal already exists")
for output_dir in (ROOT / "run1", ROOT / "run2"):
    if output_dir.exists():
        raise RuntimeError(f"post-seal output directory already exists: {output_dir}")

module_spec = importlib.util.spec_from_file_location("cycle6_analysis_preseal", ANALYSIS)
analysis_module = importlib.util.module_from_spec(module_spec)
module_spec.loader.exec_module(analysis_module)
expected = dict(analysis_module.EXPECTED)
config = json.load(open(CONFIG, encoding="ascii"))
paths = config["paths"]
if set(paths) != set(expected):
    raise RuntimeError("config key set differs from frozen expected input key set")

roles = {
    "scores": "OUTCOME_HASH_ONLY_PRESEAL_LOAD_POSTSEAL",
    "covariates": "COVARIATE_HASH_ONLY_PRESEAL_READ_POSTSEAL",
    "order": "ORDER_HASH_ONLY_PRESEAL_READ_POSTSEAL",
}
inventory = []
for key in sorted(paths):
    path = Path(paths[key])
    actual = sha256(path)
    if actual != expected[key]:
        raise RuntimeError(f"hash mismatch: {key}")
    inventory.append(
        {
            "key": key,
            "path": str(path),
            "access": roles.get(key, "READ_PRESEAL_AND_HASH_POSTSEAL"),
            "expected_sha256": expected[key],
            "actual_sha256": actual,
            "status": "MATCH",
        }
    )

with open(ROOT / "INPUT_ACCESS_ALLOWLIST.tsv", "w", encoding="ascii", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=["key", "path", "access", "expected_sha256", "actual_sha256", "status"],
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(inventory)

freeze_expected = {
    "PHASE2_ANALYSIS_CHARTER.md": "195aa495733a5a3ace0e3bda8ce3a11db8eff588477f16d10fe9eb9c868d1950",
    "PHASE2_FROZEN_SPEC.json": "14343c81ed4a9211e26117d5489ef24d3d158a028121b4c78e633e956a24ef35",
    "PHASE2_INPUT_HASHES.tsv": "20a6b6335854e0f6a220b808668241370db6d496f0de38033569e079dd5c5c3f",
    "PHASE2_RESEARCH_ACCESS_LOG.tsv": "163af9e1bc2c472f9a9835e546b44a60fdfffd71c44897231549ad4bfe5bb07a",
}
freeze_inventory = {}
for name, expected_hash in sorted(freeze_expected.items()):
    actual = sha256(FREEZE_ROOT / name)
    if actual != expected_hash:
        raise RuntimeError(f"freeze hash mismatch: {name}")
    freeze_inventory[name] = actual

analysis_text = ANALYSIS.read_text(encoding="ascii")
config_text = CONFIG.read_text(encoding="ascii")
forbidden_patterns = {
    "task030": r"task0?30|task-030",
    "task_a": r"task[_-]a|taska",
    "manuscript_path": r"[/\\]manuscript[/\\]",
    "docx": r"\.docx\b",
    "pdf": r"\.pdf\b",
    "figure_path": r"[/\\][^/\\]*figure[^/\\]*",
    "task028_results": r"task028[^'\"]*[/\\]results",
    "no_purity_worktree": r"revgate-tcga-no-purity-verify",
}
static_findings = {}
for label, pattern in forbidden_patterns.items():
    matches = []
    for source_name, text in (("run_phase2.py", analysis_text), ("EXECUTION_CONFIG.json", config_text)):
        if re.search(pattern, text, flags=re.IGNORECASE):
            matches.append(source_name)
    static_findings[label] = matches
if any(static_findings.values()):
    raise RuntimeError(f"forbidden static path/reference found: {static_findings}")

trace_text = TRACE.read_text(encoding="utf-8", errors="replace")
trace_forbidden_patterns = {
    "scores": r"scores_v3\.npz",
    "covariates": r"covariates_v3\.tsv",
    "order": r"patient_order_v3\.json",
    "task030": r"task0?30|task-030",
    "task_a": r"task[_-]a|taska",
    "results": r"[/\\]results(?:[/\\]|[\"'])",
    "manuscript": r"manuscript|\.docx|\.pdf",
    "figure": r"[/\\][^/\\]*figure[^/\\]*",
}
trace_findings = {
    label: len(re.findall(pattern, trace_text, flags=re.IGNORECASE))
    for label, pattern in trace_forbidden_patterns.items()
}
if any(trace_findings.values()):
    raise RuntimeError(f"preseal synthetic trace touched forbidden paths: {trace_findings}")
network_syscall_lines = [
    line
    for line in trace_text.splitlines()
    if re.search(r"\b(connect|sendto|recvfrom|socket)\(", line)
]

synthetic_stdout = SYNTHETIC_STDOUT.read_text(encoding="ascii")
synthetic_stderr = SYNTHETIC_STDERR.read_text(encoding="ascii")
if synthetic_stdout != "SYNTHETIC_TEST_PASS\n" or synthetic_stderr:
    raise RuntimeError("successful synthetic test artifacts do not have exact expected content")

prior_lines = PRIOR_ANALYSIS.read_text(encoding="ascii").splitlines(keepends=True)
current_lines = analysis_text.splitlines(keepends=True)
code_diff = "".join(
    difflib.unified_diff(
        prior_lines,
        current_lines,
        fromfile="cycle5/run_phase2.py",
        tofile="cycle6/run_phase2.py",
    )
)
(ROOT / "CYCLE5_TO_CYCLE6_CODE_DIFF.patch").write_text(code_diff, encoding="ascii", newline="\n")

static_audit = {
    "status": "PASS",
    "analysis_path": str(ANALYSIS),
    "config_path": str(CONFIG),
    "static_forbidden_findings": static_findings,
    "synthetic_trace_path": str(TRACE),
    "synthetic_trace_line_count": len(trace_text.splitlines()),
    "synthetic_trace_forbidden_findings": trace_findings,
    "synthetic_trace_network_syscall_count": len(network_syscall_lines),
    "synthetic_trace_network_syscall_lines": network_syscall_lines,
    "network_namespace": "bwrap --unshare-net",
}
write_json(ROOT / "STATIC_ACCESS_SCAN.json", static_audit)

environment_lines = [
    f"Interpreter: {sys.executable}",
    f"Interpreter SHA-256: {sha256(sys.executable)}",
    f"Python: {sys.version.replace(chr(10), ' ')}",
    f"NumPy: {np.__version__}",
    f"SciPy: {scipy.__version__}",
    f"pandas: {pd.__version__}",
    f"statsmodels: {statsmodels.__version__}",
    f"Platform: {platform.system()} {platform.machine()}",
    "Locale discipline: ASCII output",
    "Network discipline: bwrap --unshare-net",
]
(ROOT / "ENVIRONMENT.txt").write_text("\n".join(environment_lines) + "\n", encoding="ascii", newline="\n")

failed_attempts = []
for trace_name, state in (
    ("PRESEAL_SYNTHETIC_OPEN_TRACE.log", "synthetic_pass_then_tracer_live_tee_deadlock"),
    ("PRESEAL_SYNTHETIC_OPEN_TRACE_CLEAN.log", "synthetic_pass_then_tracer_live_tee_deadlock"),
):
    path = ROOT / trace_name
    failed_attempts.append(
        {
            "trace": trace_name,
            "state": state,
            "preserved": path.exists(),
            "sha256": sha256(path) if path.exists() else None,
        }
    )

code_hashes = {
    "run_phase2.py": sha256(ANALYSIS),
    "EXECUTION_CONFIG.json": sha256(CONFIG),
    "prepare_seal.py": sha256(Path(__file__)),
    "audit_exec.sh": sha256(ROOT / "audit_exec.sh"),
    "trace_exec.sh": sha256(ROOT / "trace_exec.sh"),
}
seal = {
    "schema_version": "1.0.0",
    "sealed_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
    "sealed_before_outcome_load": True,
    "classification": "HYPOTHESIS_POST_HOC_EXPLANATORY_SENSITIVITY_PHASE2",
    "branch": "experiment/task-b-grade-histology-sensitivity",
    "git_commit": "83503bad47b60193598b2b9ebe819c22c83e8ac1",
    "outcome_firewall": {
        "prior_numerical_outputs_opened": False,
        "outcome_file_preseal_access": "SHA256_ONLY",
        "covariate_file_preseal_access": "SHA256_ONLY",
        "patient_order_preseal_access": "SHA256_ONLY",
        "synthetic_test": "PASS",
        "synthetic_serialization_fixture": "Python and NumPy bool/int/uint/float/string scalars, 0d/1d/2d arrays, tuples, nested mappings, and actual synthetic diagnostics",
        "synthetic_trace_audit": "PASS",
        "static_access_scan": "PASS",
        "network_namespace": "bwrap --unshare-net",
    },
    "plumbing_change_only": {
        "cycle5_script_sha256": sha256(PRIOR_ANALYSIS),
        "cycle6_script_sha256": code_hashes["run_phase2.py"],
        "diff_path": str(ROOT / "CYCLE5_TO_CYCLE6_CODE_DIFF.patch"),
        "diff_sha256": sha256(ROOT / "CYCLE5_TO_CYCLE6_CODE_DIFF.patch"),
        "change": "recursive JSON container/scalar normalization plus exhaustive synthetic serialization validation; no numeric computation changed",
    },
    "preseal_attempt_history": failed_attempts
    + [
        {
            "trace": TRACE.name,
            "state": "PASS",
            "preserved": True,
            "sha256": sha256(TRACE),
        }
    ],
    "code_and_config_hashes": code_hashes,
    "input_allowlist_path": str(ROOT / "INPUT_ACCESS_ALLOWLIST.tsv"),
    "input_allowlist_sha256": sha256(ROOT / "INPUT_ACCESS_ALLOWLIST.tsv"),
    "validated_input_hash_count": len(inventory),
    "input_hashes": {row["key"]: row["actual_sha256"] for row in inventory},
    "freeze_hashes": freeze_inventory,
    "environment_path": str(ROOT / "ENVIRONMENT.txt"),
    "environment_sha256": sha256(ROOT / "ENVIRONMENT.txt"),
    "static_access_scan_path": str(ROOT / "STATIC_ACCESS_SCAN.json"),
    "static_access_scan_sha256": sha256(ROOT / "STATIC_ACCESS_SCAN.json"),
    "synthetic_stdout_sha256": sha256(SYNTHETIC_STDOUT),
    "synthetic_stderr_sha256": sha256(SYNTHETIC_STDERR),
    "seed": 20260713,
    "bootstrap_attempts_per_fit": 2000,
    "diagnostic_permutations_per_fit": 2000,
    "postseal_rule": "exactly two executions from absent clean run1/run2 directories; stop on either failure; no code/config/input change or retry",
    "scientific_files": [
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
    ],
}
write_json(SEAL, seal)
print(json.dumps({"status": "SEALED", "inputs": len(inventory), "analysis_sha256": code_hashes["run_phase2.py"]}, sort_keys=True))
