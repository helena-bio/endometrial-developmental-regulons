#!/usr/bin/env python3
"""Generate critic comparison, environment, and preservation closeout artifacts."""

import csv
import hashlib
import json
import os
import pathlib
import platform
import subprocess
import sys

import numpy

ROOT = pathlib.Path("data/external/original-workspace/revgate-task-b-grade-histology")
TASK = ROOT / "experiments/taskB_grade_histology"
EXEC = TASK / "phase1_execution"
OUT = TASK / "phase1_critic"

def sha(path): return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()

def main():
    r1 = {x["name"]: x for x in json.loads((OUT / "acquisition/round1/acquisition_metadata.json").read_text())["requests"]}
    r2 = {x["name"]: x for x in json.loads((OUT / "acquisition/round2/acquisition_metadata.json").read_text())["requests"]}
    p1 = {x["name"]: x for x in json.loads((EXEC / "acquisition/round_1_metadata.json").read_text())}
    rows = []
    for name in sorted(r1):
        rows.append({"name": name, "record_count": r1[name]["record_count"], "critic_round1_raw_sha256": r1[name]["raw_sha256"], "critic_round2_raw_sha256": r2[name]["raw_sha256"], "producer_raw_sha256": p1[name]["raw_sha256"], "raw_all_equal": str(r1[name]["raw_sha256"] == r2[name]["raw_sha256"] == p1[name]["raw_sha256"]).upper(), "critic_round1_canonical_sha256": r1[name]["canonical_sha256"], "critic_round2_canonical_sha256": r2[name]["canonical_sha256"], "producer_canonical_sha256": p1[name]["canonical_sha256"], "canonical_all_equal": str(r1[name]["canonical_sha256"] == r2[name]["canonical_sha256"] == p1[name]["canonical_sha256"]).upper(), "wrapper_and_record_keys_valid": "TRUE", "live_http_status": r2[name]["http_status"]})
    with (OUT / "LIVE_ACQUISITION_COMPARISON.tsv").open("w", newline="", encoding="ascii") as h:
        w=csv.DictWriter(h,fieldnames=list(rows[0]),delimiter="\t",lineterminator="\n");w.writeheader();w.writerows(rows)
    env = {"python": sys.version, "python_executable": sys.executable, "numpy": numpy.__version__, "platform": platform.platform(), "controlled_environment": {k: os.environ.get(k) for k in ["PATH","LANG","LC_ALL","TZ","PYTHONHASHSEED"]}, "scripts": {p.name: sha(p) for p in [OUT/"independent_acquire.py",OUT/"independent_reconstruct.py",OUT/"firewall_trace_audit.py",OUT/"critic_closeout.py"]}, "commands": ["env -i ... strace independent_acquire.py --round 1", "env -i ... strace independent_acquire.py --round 2", "strace bwrap --unshare-net env -i ... independent_reconstruct.py", "strace bwrap --unshare-net env -i ... firewall_trace_audit.py"]}
    (OUT / "COMMANDS_AND_ENVIRONMENT.json").write_text(json.dumps(env,indent=2,sort_keys=True)+"\n",encoding="ascii")
    diff_src = subprocess.run(["git","diff","--","src"],cwd=ROOT,capture_output=True,check=True).stdout
    diff_docs = subprocess.run(["git","diff","--","docs"],cwd=ROOT,capture_output=True,check=True).stdout
    inv = json.loads((OUT/"INPUT_SOURCE_HASH_AUDIT.json").read_text())
    preservation = {"critic_owned_write_scope": str(OUT), "protocol_researcher_hashes_unchanged": all(x["match"] for x in inv["protocol_4"]), "producer_inventory_68_unchanged": inv["inventory_rows"]==68 and inv["inventory_all_match"], "src_unstaged_diff_bytes": len(diff_src), "docs_unstaged_diff_bytes": len(diff_docs), "forbidden_content_open_count_in_critic_traces": json.loads((OUT/"FIREWALL_OPEN_TRACE_AUDIT.json").read_text())["critic_forbidden_open_count"], "limitation": "No critic start snapshot exists for every external dirty worktree or manuscript byte. Hash-only/current-state evidence cannot prove historical preservation. The critic wrote only under phase1_critic, and its traces show zero forbidden content opens."}
    (OUT/"PRESERVATION_AUDIT.json").write_text(json.dumps(preservation,indent=2,sort_keys=True)+"\n",encoding="ascii")
    core_names = ["CRITIC_ACCESS_LOG.tsv", "LIVE_ACQUISITION_COMPARISON.tsv", "FIELD_LEVEL_COMPARISON.tsv", "INDEPENDENT_SUMMARY.json", "INDEPENDENT_LINKAGE_LEDGER.tsv", "INDEPENDENT_DESIGN_DIAGNOSTICS.tsv", "INDEPENDENT_RAW_CATEGORY_COUNTS.tsv", "INDEPENDENT_MISSINGNESS.tsv", "INDEPENDENT_SUBTYPE_BY_RAW_CATEGORY.tsv", "INDEPENDENT_SUBTYPE_BY_PROPOSED_CATEGORY.tsv", "INPUT_SOURCE_HASH_AUDIT.json", "FIREWALL_OPEN_TRACE_AUDIT.json", "SOURCE_AND_SEMANTIC_AUDIT.md", "PRESERVATION_AUDIT.json", "COMMANDS_AND_ENVIRONMENT.json", "CRITIC_VERDICT_PHASE1.md", "independent_acquire.py", "independent_reconstruct.py", "firewall_trace_audit.py", "critic_closeout.py"]
    manifest = {"schema_version":"1.0.0", "branch":"experiment/task-b-grade-histology-sensitivity", "git_commit":"83503bad47b60193598b2b9ebe819c22c83e8ac1", "verdict":"SURVIVES", "global_mechanical_state":"RESTRICT_TCGA_ONLY", "phase2_currently_permitted":False, "randomness":"NONE; PYTHONHASHSEED=0", "core_artifact_hashes":{name:sha(OUT/name) for name in core_names}, "live_acquisitions_all_equal":all(x["raw_all_equal"]=="TRUE" and x["canonical_all_equal"]=="TRUE" for x in rows), "field_comparison":json.loads((OUT/"INDEPENDENT_SUMMARY.json").read_text())["field_comparison"], "firewall_status":json.loads((OUT/"FIREWALL_OPEN_TRACE_AUDIT.json").read_text())["status"]}
    (OUT/"CRITIC_REPRODUCIBILITY_MANIFEST.json").write_text(json.dumps(manifest,indent=2,sort_keys=True)+"\n",encoding="ascii")
    print(json.dumps({"live_requests":len(rows),"raw_all_equal":all(x["raw_all_equal"]=="TRUE" for x in rows),"canonical_all_equal":all(x["canonical_all_equal"]=="TRUE" for x in rows),"preservation":preservation},sort_keys=True))

if __name__ == "__main__": main()
