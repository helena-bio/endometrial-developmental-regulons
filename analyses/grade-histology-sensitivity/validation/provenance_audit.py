#!/usr/bin/env python3
"""Independent checksum, chronology, trace, and preservation audit for cycle 6."""

from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import stat
import subprocess
from pathlib import Path


WORKTREE = Path("data/external/original-workspace/revgate-task-b-grade-histology")
CYCLE = WORKTREE / "experiments/taskB_grade_histology/phase2_execution_cycle6"
CRITIC = WORKTREE / "experiments/taskB_grade_histology/phase2_critic_cycle6"
FREEZE = WORKTREE / "experiments/taskB_grade_histology/phase2_freeze"
EXPECTED_EXTERNAL = {
    Path("data/external/original-workspace/revgate"): "17fcf6c711c9b376b1fef426440e7583526547ef456f6bada618ef0d2970bc22",
    Path("data/external/original-workspace/revgate-task-a-perclass"): "84d57032f00ad4846d9e8c8223e9d1b1fd1d6ad9dcbeb85443fff58672ef59c3",
    Path("data/external/original-workspace/revgate-tcga-no-purity-verify"): "16a3e9e9cf8bc9b3ae225c8e25c964939eb1a8d7a8212f844f652f45923e5210",
}


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def hash_lines(path: Path, base: Path) -> dict:
    checked = []
    with open(path, encoding="ascii") as handle:
        for line_number, line in enumerate(handle, 1):
            line = line.rstrip("\n")
            expected, relative = line.split("  ", 1)
            target = base / relative
            actual = sha(target) if target.is_file() else None
            checked.append({
                "line": line_number,
                "path": relative,
                "expected": expected,
                "actual": actual,
                "match": expected == actual,
            })
    return {
        "line_count": len(checked),
        "match_count": sum(row["match"] for row in checked),
        "mismatch_count": sum(not row["match"] for row in checked),
        "rows": checked,
    }


def inventory(path: Path, base: Path) -> dict:
    checked = []
    with open(path, newline="", encoding="ascii") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            target = base / row["path"]
            exists = target.is_file()
            actual_size = target.stat().st_size if exists else None
            actual_sha = sha(target) if exists else None
            checked.append({
                "path": row["path"],
                "exists": exists,
                "expected_size": int(row["size"]),
                "actual_size": actual_size,
                "expected_sha256": row["sha256"],
                "actual_sha256": actual_sha,
                "match": exists and actual_size == int(row["size"]) and actual_sha == row["sha256"],
            })
    return {
        "row_count": len(checked),
        "match_count": sum(row["match"] for row in checked),
        "mismatch_count": sum(not row["match"] for row in checked),
        "rows": checked,
    }


def phase2_input_hashes() -> dict:
    checked = []
    with open(FREEZE / "PHASE2_INPUT_HASHES.tsv", newline="", encoding="ascii") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            target = Path(row["path"])
            actual = sha(target) if target.is_file() else None
            checked.append({
                "class": row["class"],
                "path": row["path"],
                "expected": row["sha256"],
                "actual": actual,
                "match": row["sha256"] == actual,
            })
    return {
        "row_count": len(checked),
        "match_count": sum(row["match"] for row in checked),
        "mismatch_count": sum(not row["match"] for row in checked),
        "rows": checked,
    }


def input_allowlist() -> dict:
    checked = []
    with open(CYCLE / "INPUT_ACCESS_ALLOWLIST.tsv", newline="", encoding="ascii") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            target = Path(row["path"])
            actual = sha(target) if target.is_file() else None
            checked.append({
                "key": row["key"],
                "path": row["path"],
                "expected": row["expected_sha256"],
                "recorded_actual": row["actual_sha256"],
                "independent_actual": actual,
                "match": actual == row["expected_sha256"] == row["actual_sha256"] and row["status"] == "MATCH",
            })
    return {
        "row_count": len(checked),
        "match_count": sum(row["match"] for row in checked),
        "mismatch_count": sum(not row["match"] for row in checked),
        "rows": checked,
    }


def command_audit() -> dict:
    with open(CYCLE / "audit/command_log.tsv", encoding="ascii") as handle:
        log_rows = [line.rstrip("\n").split("\t", 4) for line in handle if line.strip()]
    ids = [row[0] for row in log_rows]
    expected_ids = [f"{index:04d}" for index in range(1, len(log_rows) + 1)]
    raw_checks = []
    for command_id, start, end, exit_status, command in log_rows:
        prefix = CYCLE / "audit/commands" / command_id
        meta_path = prefix.with_suffix(".meta")
        stdout_path = prefix.with_suffix(".stdout")
        stderr_path = prefix.with_suffix(".stderr")
        meta = {}
        for line in meta_path.read_text(encoding="ascii").splitlines():
            key, value = line.split("=", 1)
            meta[key] = value
        raw_checks.append({
            "command_id": command_id,
            "meta_exists": meta_path.is_file(),
            "stdout_exists": stdout_path.is_file(),
            "stderr_exists": stderr_path.is_file(),
            "log_meta_identity": (
                meta.get("command_id") == command_id
                and meta.get("start_utc") == start
                and meta.get("end_utc") == end
                and meta.get("exit_status") == exit_status
                and meta.get("command") == command
            ),
            "stdout_hash_match": meta.get("stdout_sha256") == sha(stdout_path),
            "stderr_hash_match": meta.get("stderr_sha256") == sha(stderr_path),
        })
    with open(CYCLE / "COMMANDS_CHRONOLOGICAL.tsv", newline="", encoding="ascii") as handle:
        table = list(csv.DictReader(handle, delimiter="\t"))
    table_ids = [row["command_id"] for row in table]
    return {
        "raw_command_count": len(log_rows),
        "raw_ids_contiguous": ids == expected_ids,
        "raw_per_command_all_match": all(
            all(value for key, value in row.items() if key != "command_id")
            for row in raw_checks
        ),
        "raw_checks": raw_checks,
        "chronological_table_count": len(table),
        "chronological_table_ids": table_ids,
        "chronological_table_missing_raw_ids": sorted(set(ids) - set(table_ids)),
        "chronological_table_extra_ids": sorted(set(table_ids) - set(ids)),
        "permission_error_disclosed_only_outside_wrapper": True,
        "permission_error_note": "The direct wrapper invocation could not start and therefore has no command ID; COMMAND_COVERAGE_NOTE.md preserves the exact command/error.",
    }


def quoted_absolute_paths(trace: Path) -> set[Path]:
    found = set()
    text = trace.read_text(encoding="utf-8", errors="replace")
    for value in re.findall(r'"(/srv/[^"]+)"', text):
        found.add(Path(value))
    return found


def lexical_ancestors(path: Path) -> set[Path]:
    ancestors = set()
    current = path
    while str(current).startswith("data/external/original-workspace"):
        ancestors.add(current)
        if current == current.parent:
            break
        current = current.parent
    return ancestors


def trace_audit() -> dict:
    seal = json.loads((CYCLE / "PRE_OUTCOME_EXECUTION_SEAL.json").read_text(encoding="ascii"))
    config = json.loads((CYCLE / "EXECUTION_CONFIG.json").read_text(encoding="ascii"))
    exact_inputs = {Path(path) for path in config["paths"].values()}
    allowed_cycle = lexical_ancestors(CYCLE)

    preseal = {}
    forbidden_preseal = (
        "scores_v3.npz", "covariates_v3.tsv", "patient_order_v3.json",
        "phase2_execution_cycle1", "phase2_execution_cycle2", "phase2_execution_cycle3",
        "phase2_execution_cycle4", "phase2_execution_cycle5", "phase2_critic",
        "task030", "task-030", "task_a", "task-a", "/manuscript/", ".docx", ".pdf",
    )
    for name in (
        "PRESEAL_SYNTHETIC_OPEN_TRACE.log",
        "PRESEAL_SYNTHETIC_OPEN_TRACE_CLEAN.log",
        "PRESEAL_SYNTHETIC_OPEN_TRACE_SUCCESS.log",
    ):
        text = (CYCLE / name).read_text(encoding="utf-8", errors="replace").lower()
        preseal[name] = {
            "line_count": len(text.splitlines()),
            "forbidden_occurrence_count": sum(text.count(token.lower()) for token in forbidden_preseal),
            "internet_connect_count": len(re.findall(r"connect\([^\n]*(?:AF_INET|AF_INET6)", text)),
        }

    postseal = {}
    forbidden_postseal = (
        "phase2_execution_cycle1", "phase2_execution_cycle2", "phase2_execution_cycle3",
        "phase2_execution_cycle4", "phase2_execution_cycle5", "phase2_critic",
        "task030", "task-030", "task_a", "task-a", "revgate-tcga-no-purity-verify",
        "/manuscript/", ".docx", ".pdf",
    )
    for name in ("POSTSEAL_RUN1_OPEN_TRACE.log", "POSTSEAL_RUN2_OPEN_TRACE.log"):
        trace = CYCLE / name
        text = trace.read_text(encoding="utf-8", errors="replace")
        lower = text.lower()
        workspace_paths = {
            path for path in quoted_absolute_paths(trace)
            if str(path).startswith("data/external/original-workspace")
        }
        unexpected = sorted(
            str(path) for path in workspace_paths
            if path not in exact_inputs
            and path not in allowed_cycle
            and not str(path).startswith(str(CYCLE) + "/")
            and str(path) != "data/external/original-workspace/.language-model tool/wo"
        )
        postseal[name] = {
            "line_count": len(text.splitlines()),
            "forbidden_occurrence_count": sum(lower.count(token.lower()) for token in forbidden_postseal),
            "internet_connect_count": len(re.findall(r"connect\([^\n]*(?:AF_INET|AF_INET6)", text)),
            "unexpected_workspace_paths": unexpected,
            "declared_inputs_seen": {
                str(path): str(path) in text for path in sorted(exact_inputs, key=str)
            },
        }

    upstream_trace = CYCLE / "UPSTREAM_RECONCILIATION_OPEN_TRACE.log"
    upstream_text = upstream_trace.read_text(encoding="utf-8", errors="replace")
    upstream_allowed = {
        Path("data/external/original-workspace/task028-freeze-b-draft/execution/results_v3/primary_results.tsv"),
        Path("data/external/original-workspace/task028-freeze-b-draft/execution/results_v3/sensitivity_SENS_nopurity.tsv"),
        Path("data/external/original-workspace/task030-tcga-nopurity-audit/results/tcga_primary_vs_nopurity_per_target.tsv"),
    }
    upstream_workspace_paths = {
        path for path in quoted_absolute_paths(upstream_trace)
        if str(path).startswith("data/external/original-workspace")
    }
    upstream_unexpected = sorted(
        str(path) for path in upstream_workspace_paths
        if path not in upstream_allowed
        and path not in allowed_cycle
        and not str(path).startswith(str(CYCLE) + "/")
        and str(path) != "data/external/original-workspace/.language-model tool/wo"
    )
    return {
        "preseal": preseal,
        "preseal_all_clean": all(
            row["forbidden_occurrence_count"] == 0 and row["internet_connect_count"] == 0
            for row in preseal.values()
        ),
        "postseal": postseal,
        "postseal_all_clean": all(
            row["forbidden_occurrence_count"] == 0
            and row["internet_connect_count"] == 0
            and not row["unexpected_workspace_paths"]
            for row in postseal.values()
        ),
        "upstream": {
            "line_count": len(upstream_text.splitlines()),
            "internet_connect_count": len(re.findall(r"connect\([^\n]*(?:AF_INET|AF_INET6)", upstream_text)),
            "authorized_inputs_seen": {str(path): str(path) in upstream_text for path in upstream_allowed},
            "unexpected_workspace_paths": upstream_unexpected,
        },
        "seal_claims_preoutcome": seal["outcome_firewall"],
    }


def preservation() -> dict:
    external = []
    for path, expected in EXPECTED_EXTERNAL.items():
        status_bytes = subprocess.check_output(
            ["git", "-C", str(path), "status", "--porcelain", "--untracked-files=no"]
        )
        observed = hashlib.sha256(status_bytes).hexdigest()
        external.append({
            "path": str(path),
            "expected": expected,
            "observed": observed,
            "match": expected == observed,
            "head": subprocess.check_output(["git", "-C", str(path), "rev-parse", "HEAD"], text=True).strip(),
            "branch": subprocess.check_output(["git", "-C", str(path), "branch", "--show-current"], text=True).strip(),
        })
    unstaged_src_docs = subprocess.check_output(
        ["git", "diff", "--name-only", "--", "src", "docs"], cwd=WORKTREE, text=True
    ).splitlines()
    staged_src_docs = subprocess.check_output(
        ["git", "diff", "--cached", "--name-only", "--", "src", "docs"], cwd=WORKTREE, text=True
    ).splitlines()
    changed = subprocess.check_output(
        ["git", "status", "--porcelain"], cwd=WORKTREE, text=True
    ).splitlines()
    manuscript_like = [
        row for row in changed
        if re.search(r"manuscript|\\.docx$|\\.pdf$", row, re.IGNORECASE)
    ]
    return {
        "external": external,
        "external_all_match": all(row["match"] for row in external),
        "unstaged_src_docs_paths": unstaged_src_docs,
        "staged_src_docs_paths": staged_src_docs,
        "src_paths_changed": [path for path in unstaged_src_docs + staged_src_docs if path.startswith("src/")],
        "manuscript_docx_pdf_status_paths": manuscript_like,
        "note": "The listed staged docs changes predate cycle 6; their Phase-1/freeze bytes are separately hash-pinned and reverified.",
    }


def run_identity() -> dict:
    names1 = sorted(path.name for path in (CYCLE / "run1").iterdir() if path.is_file())
    names2 = sorted(path.name for path in (CYCLE / "run2").iterdir() if path.is_file())
    rows = []
    for name in sorted(set(names1) | set(names2)):
        one = CYCLE / "run1" / name
        two = CYCLE / "run2" / name
        rows.append({
            "file": name,
            "run1_sha256": sha(one) if one.is_file() else None,
            "run2_sha256": sha(two) if two.is_file() else None,
            "byte_identical": one.is_file() and two.is_file() and one.read_bytes() == two.read_bytes(),
            "run1_mode": stat.S_IMODE(one.stat().st_mode) if one.is_file() else None,
            "run2_mode": stat.S_IMODE(two.stat().st_mode) if two.is_file() else None,
        })
    return {
        "run1_count": len(names1),
        "run2_count": len(names2),
        "identity_count": sum(row["byte_identical"] for row in rows),
        "all_read_only": all(row["run1_mode"] == 0o444 and row["run2_mode"] == 0o444 for row in rows),
        "run1_directory_mode": stat.S_IMODE((CYCLE / "run1").stat().st_mode),
        "run2_directory_mode": stat.S_IMODE((CYCLE / "run2").stat().st_mode),
        "rows": rows,
    }


def main():
    result = {
        "branch": subprocess.check_output(["git", "branch", "--show-current"], cwd=WORKTREE, text=True).strip(),
        "head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=WORKTREE, text=True).strip(),
        "output_sha256s": hash_lines(CYCLE / "OUTPUT_SHA256SUMS.txt", CYCLE),
        "run1_sha256s": hash_lines(CYCLE / "RUN1_SCIENTIFIC_SHA256SUMS.txt", CYCLE),
        "run2_sha256s": hash_lines(CYCLE / "RUN2_SCIENTIFIC_SHA256SUMS.txt", CYCLE),
        "file_inventory": inventory(CYCLE / "FILE_INVENTORY.tsv", CYCLE),
        "phase2_input_hashes": phase2_input_hashes(),
        "input_allowlist": input_allowlist(),
        "run_identity": run_identity(),
        "commands": command_audit(),
        "traces": trace_audit(),
        "preservation": preservation(),
    }
    (CRITIC / "PROVENANCE_AUDIT.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="ascii"
    )
    concise = {
        "branch": result["branch"],
        "head": result["head"],
        "output_checks": [
            result["output_sha256s"]["match_count"],
            result["output_sha256s"]["line_count"],
        ],
        "inventory_checks": [
            result["file_inventory"]["match_count"],
            result["file_inventory"]["row_count"],
        ],
        "phase2_inputs": [
            result["phase2_input_hashes"]["match_count"],
            result["phase2_input_hashes"]["row_count"],
        ],
        "sealed_allowlist": [
            result["input_allowlist"]["match_count"],
            result["input_allowlist"]["row_count"],
        ],
        "run_identity": [
            result["run_identity"]["identity_count"],
            result["run_identity"]["run1_count"],
        ],
        "raw_commands": result["commands"]["raw_command_count"],
        "command_table_missing": result["commands"]["chronological_table_missing_raw_ids"],
        "traces_clean": [
            result["traces"]["preseal_all_clean"],
            result["traces"]["postseal_all_clean"],
            not result["traces"]["upstream"]["unexpected_workspace_paths"],
        ],
        "external_preservation": result["preservation"]["external_all_match"],
    }
    print(json.dumps(concise, sort_keys=True))


if __name__ == "__main__":
    main()
