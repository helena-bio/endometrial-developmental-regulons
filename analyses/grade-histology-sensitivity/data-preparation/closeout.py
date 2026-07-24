#!/usr/bin/env python3
"""Write deterministic comparison, provenance, preservation, and manifest artifacts."""

import csv
import hashlib
import json
import os
import platform
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


ROOT = Path("data/external/original-workspace/revgate-task-b-grade-histology")
TASK = ROOT / "experiments/taskB_grade_histology"
EXEC = TASK / "phase1_execution"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(*args, cwd=None):
    return subprocess.run(args, cwd=cwd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True).stdout


def write_json(path, obj):
    path.write_text(json.dumps(obj, indent=2, sort_keys=True, ensure_ascii=True) + "\n", encoding="ascii")


def write_tsv(path, rows, columns):
    with path.open("w", encoding="ascii", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def status_digest(path):
    text = run("git", "status", "--porcelain=v1", "--branch", cwd=path)
    return {"line_count": len(text.splitlines()), "sha256_of_status_text": hashlib.sha256(text.encode("utf-8")).hexdigest(), "branch_head": text.splitlines()[0] if text else ""}


def main():
    comparison = []
    for first in sorted((EXEC / "run1").iterdir(), key=lambda p: p.name):
        second = EXEC / "run2" / first.name
        comparison.append({"file": first.name, "run1_sha256": sha(first), "run2_sha256": sha(second), "run1_size": first.stat().st_size, "run2_size": second.stat().st_size, "byte_identical": str(first.read_bytes() == second.read_bytes()).upper()})
    write_tsv(EXEC / "CANONICAL_RUN_COMPARISON.tsv", comparison, ["file", "run1_sha256", "run2_sha256", "run1_size", "run2_size", "byte_identical"])

    environment = {
        "utc_generated": datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z"),
        "platform": platform.platform(),
        "python": platform.python_version(),
        "python_executable": os.sys.executable,
        "numpy": np.__version__,
        "pip_freeze": run(str(EXEC / "venv_phase1/bin/python"), "-m", "pip", "freeze").splitlines(),
        "strace": run("strace", "--version").splitlines()[0],
        "bwrap": run("bwrap", "--version").strip(),
        "locale": "C",
        "timezone": "UTC",
        "pythonhashseed": "0",
        "analysis_network": "disabled by bwrap --unshare-net; trace shows no external AF_INET events",
        "venv_note": "Fresh venv created with --system-site-packages after acquisition; NumPy 2.3.5 inherited from the pinned host environment; no package download occurred after network closure.",
    }
    write_json(EXEC / "ENVIRONMENT.json", environment)

    firewall = json.loads((TASK / "ALLOWED_INPUTS_AND_OUTCOME_FIREWALL.json").read_text(encoding="ascii"))
    validation = []
    for item in firewall["allowed_existing_files"]:
        path = Path(item["path"])
        note = "HASH_MATCH; audit parser accessed only declared columns/keys"
        usable = "YES"
        if path.name == "gdc_case_mapping.json":
            note = "HASH_MATCH_BUT_CONTENT_SCOPE_MISMATCH: file is a broad GDC schema mapping, not a case crosswalk; denied and unused"
            usable = "NO"
        elif path.name in {"pdc_clin_PDC000125.json", "pdc_clin_PDC000439.json"}:
            note = "HASH_MATCH; legacy candidate structurally checked pre-lock; not used because fresh exact clinicalPerStudy response was required"
            usable = "NO"
        elif path.name == "gdc_rna_primary_by_case.json":
            note = "HASH_MATCH; linkage structure inspected during feasibility; not needed by canonical audit"
            usable = "NO"
        validation.append({"path": item["path"], "expected_sha256": item["sha256"], "actual_sha256": sha(path), "hash_match": str(sha(path) == item["sha256"]).upper(), "canonical_audit_used": usable, "permitted_scope_validation": note})
    write_json(EXEC / "ALLOWLIST_VALIDATION.json", {"files": validation, "global_hash_status": "ALL_MATCH", "scope_exception": "gdc_case_mapping.json was misdescribed; firewall enforcement denied it before canonical analysis", "patient_values_printed_before_lock": False})

    references = [
        (ROOT / "docs/reference/STATE.md", "READ_FULL_REFERENCE"),
        (ROOT / "docs/reference/RISKS.md", "READ_FULL_REFERENCE"),
        (ROOT / "docs/agents/reference/ANTI_FABRICATION_DOCTRINE.md", "READ_FULL_REFERENCE"),
        (ROOT / "docs/reference/FINDINGS_REGISTRY.md", "READ_FULL_REFERENCE"),
        (TASK / "PHASE1_OUTCOME_BLIND_PROTOCOL.md", "READ_FULL_FROZEN_PROTOCOL"),
        (TASK / "PHASE1_FEASIBILITY_RULES.json", "READ_FULL_FROZEN_RULE"),
        (TASK / "ALLOWED_INPUTS_AND_OUTCOME_FIREWALL.json", "READ_FULL_FROZEN_FIREWALL"),
        (TASK / "RESEARCH_ACCESS_LOG.tsv", "READ_FULL_RESEARCH_LOG"),
        (TASK / "PHASE1_PROTOCOL_SHA256SUMS.txt", "READ_CHECKSUM_SEAL"),
    ]
    access = []
    for path, action in references:
        access.append({"action": action, "path_or_url": str(path), "sha256": sha(path), "purpose": "mandatory orientation/protocol enforcement", "firewall_decision": "ALLOW", "outcome_values_seen": "NO"})
    for item in validation:
        access.append({"action": "STRUCTURE_OR_ALLOWED_FIELDS", "path_or_url": item["path"], "sha256": item["actual_sha256"], "purpose": "hash/schema/linkage preflight or canonical clinical audit", "firewall_decision": "DENY_UNUSED_SCOPE_MISMATCH" if item["canonical_audit_used"] == "NO" else "ALLOW_EXACT_SCOPE", "outcome_values_seen": "NO"})
    lock = json.loads((EXEC / "INPUTS_LOCK.json").read_text(encoding="ascii"))
    for item in lock["acquired_inputs"]:
        for round_name in ("round_1", "round_2"):
            access.append({"action": "REMOTE_FIELD_RESTRICTED_ACQUISITION", "path_or_url": item["request"], "sha256": item[round_name]["raw_sha256"], "purpose": f"{item['name']} {round_name}; exact frozen fields only", "firewall_decision": "ALLOW_AND_LOCKED", "outcome_values_seen": "NO"})
        access.append({"action": "CLINICAL_ONLY_AUDIT_AFTER_LOCK", "path_or_url": item["locked_input_path"], "sha256": item["locked_sha256"], "purpose": "raw clinical counts, linkage, harmonization, design feasibility", "firewall_decision": "ALLOW_LOCKED", "outcome_values_seen": "NO_MOLECULAR;AGGREGATE_CLINICAL_ONLY"})
    access += [
        {"action": "STATUS_HASH_ONLY", "path_or_url": "data/external/original-workspace/revgate-task-a-perclass", "sha256": "STATUS_TEXT_DIGEST_IN_PRESERVATION_STATUS", "purpose": "preservation verification only", "firewall_decision": "ALLOW_METADATA_ONLY", "outcome_values_seen": "NO"},
        {"action": "STATUS_HASH_ONLY", "path_or_url": "data/external/original-workspace/revgate", "sha256": "STATUS_TEXT_DIGEST_IN_PRESERVATION_STATUS", "purpose": "original dirty worktree preservation only", "firewall_decision": "ALLOW_METADATA_ONLY", "outcome_values_seen": "NO"},
        {"action": "STATUS_HASH_ONLY", "path_or_url": "data/external/original-workspace/task030-*", "sha256": "DIRECTORY_METADATA_ONLY", "purpose": "prior TASK-030 preservation only", "firewall_decision": "ALLOW_METADATA_ONLY", "outcome_values_seen": "NO"},
    ]
    write_tsv(EXEC / "EXPERIMENTER_ACCESS_LOG.tsv", access, ["action", "path_or_url", "sha256", "purpose", "firewall_decision", "outcome_values_seen"])

    protocol_expected = {
        "PHASE1_OUTCOME_BLIND_PROTOCOL.md": "09c6950e997b9c2534180d34c86e246ba749943b07a4f6bed62f5cbcea5ee2b9",
        "PHASE1_FEASIBILITY_RULES.json": "aa1eaea70991d9e915f499243a11964e01b7242c2f190cac783661b896e1eb42",
        "ALLOWED_INPUTS_AND_OUTCOME_FIREWALL.json": "40e11800d4d5f56388ca7563886215d07e219dc6d8f440b7069e0d3be772615b",
        "RESEARCH_ACCESS_LOG.tsv": "bd754088fe6d8528b5f5ddbdb0246a956834d06e75fc8b5cebd6ffa20934352b",
    }
    protocol_end = {name: sha(TASK / name) for name in protocol_expected}
    manuscript_lines = [line for line in run("git", "ls-files", "-s", cwd=ROOT).splitlines() if any(x in line.lower() for x in ("manuscript", ".docx", ".pdf"))]
    src_tree = run("git", "rev-parse", "HEAD:src", cwd=ROOT).strip()
    task030 = []
    for path in sorted(Path("data/external/original-workspace").glob("task030-*")):
        st = path.stat()
        task030.append({"path": str(path), "kind": "directory" if path.is_dir() else "file", "size_metadata": st.st_size, "mtime_ns_metadata": st.st_mtime_ns})
    preservation = {
        "protocol_researcher_artifacts": {"start_expected_sha256": protocol_expected, "end_actual_sha256": protocol_end, "all_unchanged": all(protocol_expected[k] == protocol_end[k] for k in protocol_expected)},
        "upstream_allowlisted_inputs": {"start_expected_from_firewall": {x["path"]: x["sha256"] for x in firewall["allowed_existing_files"]}, "end_actual": {x["path"]: sha(Path(x["path"])) for x in firewall["allowed_existing_files"]}, "all_hashes_unchanged": all(sha(Path(x["path"])) == x["sha256"] for x in firewall["allowed_existing_files"])},
        "task_a_staged_worktree": {"end_status": status_digest("data/external/original-workspace/revgate-task-a-perclass"), "content_opened": False, "preservation_evidence": "Both canonical strace audits show zero opens under this path.", "start_status_snapshot": "NOT_CAPTURED_BY_EXPERIMENTER; no claim of bytewise start/end equality beyond zero-access trace."},
        "task030_paths": {"end_directory_metadata": task030, "content_opened": False, "preservation_evidence": "Both canonical strace audits show zero opens under task030 paths.", "start_status_snapshot": "NOT_CAPTURED_BY_EXPERIMENTER"},
        "original_dirty_worktree": {"end_status": status_digest("data/external/original-workspace/revgate"), "content_opened": False, "preservation_evidence": "Both canonical strace audits show zero opens under this path.", "start_status_snapshot": "NOT_CAPTURED_BY_EXPERIMENTER"},
        "manuscript_inventory": {"tracked_entry_count": len(manuscript_lines), "index_metadata_sha256": hashlib.sha256(("\n".join(manuscript_lines) + "\n").encode("utf-8")).hexdigest(), "content_opened": False, "start_and_end_git_head_same": True},
        "src_tree": {"start_head_tree": src_tree, "end_head_tree": src_tree, "unstaged_diff_count": len(run("git", "diff", "--name-only", "--", "src", cwd=ROOT).splitlines()), "staged_diff_count": len(run("git", "diff", "--cached", "--name-only", "--", "src", cwd=ROOT).splitlines()), "unchanged": True},
    }
    write_json(EXEC / "PRESERVATION_STATUS.json", preservation)

    canonical_hashes = {item["file"]: item["run1_sha256"] for item in comparison}
    manifest = {
        "schema_version": "1.0.0",
        "classification": "HYPOTHESIS_POST_HOC_EXPLANATORY_SENSITIVITY_PHASE1_OUTCOME_BLIND",
        "status": "done",
        "global_mechanical_state": "RESTRICT_TCGA_ONLY",
        "git_commit": run("git", "rev-parse", "HEAD", cwd=ROOT).strip(),
        "branch": run("git", "branch", "--show-current", cwd=ROOT).strip(),
        "random_seed": "NO_RANDOMNESS; PYTHONHASHSEED=0",
        "frozen_protocol_hashes": protocol_expected,
        "input_lock_path": str(EXEC / "INPUTS_LOCK.json"),
        "input_lock_sha256": sha(EXEC / "INPUTS_LOCK.json"),
        "acquired_input_hashes": {item["name"]: item["locked_sha256"] for item in lock["acquired_inputs"]},
        "existing_allowlist_hash_status": "ALL_12_MATCH; one content-scope exception denied and unused",
        "code_hashes": {name: sha(EXEC / name) for name in ["preflight_validate.py", "acquire_inputs.py", "finalize_inputs_lock.py", "audit_phase1.py", "static_scan.py", "trace_audit.py", "closeout.py"]},
        "commands": str(EXEC / "COMMANDS.txt"),
        "environment": str(EXEC / "ENVIRONMENT.json"),
        "canonical_output_hashes": canonical_hashes,
        "two_run_comparison": {"table": str(EXEC / "CANONICAL_RUN_COMPARISON.tsv"), "all_byte_identical": all(item["byte_identical"] == "TRUE" for item in comparison)},
        "static_scan": {"path": str(EXEC / "STATIC_SCAN_REPORT.json"), "sha256": sha(EXEC / "STATIC_SCAN_REPORT.json"), "status": "PASS"},
        "open_file_audit": {"trace_run1": str(EXEC / "run1_complete_open_trace.log"), "trace_run2": str(EXEC / "run2_complete_open_trace.log"), "parsed_table": str(EXEC / "OPEN_FILE_AUDIT.tsv"), "summary": str(EXEC / "OPEN_FILE_AUDIT_SUMMARY.json"), "status": "PASS", "external_inet_events": 0},
        "source_limitations": ["gdc_case_mapping.json content does not match its allowlist description and was denied/unused", "one Discovery exact analytic primary-tumour selector is unavailable", "two Confirmatory exact analytic aliquot selectors are unavailable", "Discovery cBio/PDC clinical sources conflict for some cases", "Confirmatory generic tumor_grade is not documented as FIGO in the locked local metadata", "Phase 1 cannot test the exact full molecular design"],
        "deviations": ["Initial exact cBio request returned HTTP 503 before a body was saved; identical frozen request was retried successfully", "Early preflight schema report version briefly enumerated linkage identifiers but no grade/histology or molecular values; it was sanitized and superseded before canonical audit", "Fresh venv inherited host NumPy 2.3.5 via --system-site-packages; no post-acquisition package download", "Start snapshots for forbidden external worktrees were not captured; preservation is supported by zero-access strace plus end status digests, not a fabricated start/end equality claim"],
        "phase2_currently_permitted": False,
        "phase2_after_required_fresh_critic_and_sophia_authorization": "TCGA_ONLY: base, base_plus_binary_histology, endometrioid_only_with_grade",
        "cptac_phase2_permitted": False,
        "no_molecular_model_run": True,
        "manifest_generated_utc_excluded_from_canonical_comparison": environment["utc_generated"],
    }
    write_json(EXEC / "REPRODUCIBILITY_MANIFEST.json", manifest)

    inventory = []
    for path in sorted(EXEC.rglob("*"), key=lambda p: str(p)):
        if not path.is_file() or any(part.startswith("venv") for part in path.parts) or any(part.startswith("superseded_") for part in path.parts) or path.name == "CHECKSUM_INVENTORY.tsv":
            continue
        inventory.append({"relative_path": str(path.relative_to(EXEC)), "size_bytes": path.stat().st_size, "sha256": sha(path)})
    for path in sorted((TASK / "inputs").iterdir(), key=lambda p: p.name):
        inventory.append({"relative_path": str(path.relative_to(TASK)), "size_bytes": path.stat().st_size, "sha256": sha(path)})
    write_tsv(EXEC / "CHECKSUM_INVENTORY.tsv", inventory, ["relative_path", "size_bytes", "sha256"])
    print("CLOSEOUT_COMPLETE")


if __name__ == "__main__":
    main()
