#!/usr/bin/env python3
"""Seal double-fetched clinical-only inputs before any value-level audit."""

import hashlib
import json
import shutil
from pathlib import Path


ROOT = Path("data/external/original-workspace/revgate-task-b-grade-histology")
TASK = ROOT / "experiments/taskB_grade_histology"
EXEC = TASK / "phase1_execution"
ACQ = EXEC / "acquisition"


def file_sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    first = json.loads((ACQ / "round_1_metadata.json").read_text(encoding="ascii"))
    second = json.loads((ACQ / "round_2_metadata.json").read_text(encoding="ascii"))
    if [x["name"] for x in first] != [x["name"] for x in second]:
        raise SystemExit("ROUND_NAME_MISMATCH")
    locked = []
    for a, b in zip(first, second):
        if a["raw_sha256"] != b["raw_sha256"]:
            if a["canonical_sha256"] != b["canonical_sha256"]:
                raise SystemExit(f"DOUBLE_FETCH_MISMATCH:{a['name']}")
            raise SystemExit(f"CANONICALIZATION_REQUIRED_NOT_USED:{a['name']}")
        source = Path(a["raw_path"])
        target = Path(a["final_path"])
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
        if file_sha(target) != a["raw_sha256"]:
            raise SystemExit(f"COPY_HASH_MISMATCH:{a['name']}")
        locked.append({
            "name": a["name"],
            "authority": a["authority"],
            "method": a["method"],
            "request": a["request"],
            "graphql_query": a.get("graphql_query"),
            "requested_attribute_id": a.get("requested_attribute_id"),
            "pdc_study_id": a.get("pdc_study_id"),
            "study_uuid_version_pin": a.get("study_uuid_version_pin"),
            "study_version_pin": a.get("study_version_pin"),
            "round_1": {k: a[k] for k in ("retrieval_utc", "status", "headers", "raw_path", "raw_size", "raw_sha256", "canonical_sha256", "record_count")},
            "round_2": {k: b[k] for k in ("retrieval_utc", "status", "headers", "raw_path", "raw_size", "raw_sha256", "canonical_sha256", "record_count")},
            "comparison": "RAW_BYTE_IDENTICAL",
            "locked_input_path": str(target),
            "locked_sha256": file_sha(target),
        })
    preflight = json.loads((EXEC / "PREFLIGHT_SCHEMA_REPORT.json").read_text(encoding="ascii"))
    lock = {
        "schema_version": "1.0.0",
        "phase": "TASK_B_PHASE1_OUTCOME_BLIND",
        "base_commit": "83503bad47b60193598b2b9ebe819c22c83e8ac1",
        "protocol_files": {
            "PHASE1_OUTCOME_BLIND_PROTOCOL.md": "09c6950e997b9c2534180d34c86e246ba749943b07a4f6bed62f5cbcea5ee2b9",
            "PHASE1_FEASIBILITY_RULES.json": "aa1eaea70991d9e915f499243a11964e01b7242c2f190cac783661b896e1eb42",
            "ALLOWED_INPUTS_AND_OUTCOME_FIREWALL.json": "40e11800d4d5f56388ca7563886215d07e219dc6d8f440b7069e0d3be772615b",
            "RESEARCH_ACCESS_LOG.tsv": "bd754088fe6d8528b5f5ddbdb0246a956834d06e75fc8b5cebd6ffa20934352b",
        },
        "checksum_file_format_warning": "PHASE1_PROTOCOL_SHA256SUMS.txt has two leading metadata lines before four checksum-format lines; four hashes verified explicitly; frozen file not edited.",
        "existing_allowlisted_inputs": preflight,
        "acquired_inputs": locked,
        "failed_attempts": [{"round": 1, "attempt": 1, "endpoint": first[0]["request"], "http_status": 503, "response_body_saved": False, "resolution": "Exact frozen request retried without changing fields or endpoint."}],
        "lock_rule": "Written before any patient value, raw category frequency, or aggregate clinical distribution was printed or computed.",
    }
    out = EXEC / "INPUTS_LOCK.json"
    out.write_text(json.dumps(lock, indent=2, sort_keys=True) + "\n", encoding="ascii")
    print("INPUTS_LOCK_WRITTEN")


if __name__ == "__main__":
    main()
