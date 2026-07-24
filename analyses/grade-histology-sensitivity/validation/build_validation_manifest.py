#!/usr/bin/env python3
"""Build the cycle-6 critic reproducibility manifest and checksum inventory."""

import hashlib
import json
import platform
import subprocess
import sys
from pathlib import Path

import numpy
import pandas
import scipy


ROOT = Path("data/external/original-workspace/revgate-task-b-grade-histology")
HERE = ROOT / "experiments/taskB_grade_histology/phase2_critic_cycle6"
CYCLE = ROOT / "experiments/taskB_grade_histology/phase2_execution_cycle6"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    reproduction = json.loads((HERE / "INDEPENDENT_REPRODUCTION_SUMMARY.json").read_text())
    provenance = json.loads((HERE / "PROVENANCE_AUDIT.json").read_text())
    upstream = json.loads((HERE / "INDEPENDENT_UPSTREAM_RECONCILIATION.json").read_text())
    adversarial = json.loads((HERE / "ADVERSARIAL_CHECKS.json").read_text())
    excluded = {
        "CRITIC_REPRODUCIBILITY_MANIFEST.json",
        "CRITIC_SHA256SUMS.txt",
        "CRITIC_VERDICT_CYCLE6.md",
        "__pycache__",
    }
    artifacts = []
    for path in sorted(HERE.iterdir()):
        if path.name in excluded or not path.is_file():
            continue
        artifacts.append({
            "path": path.name,
            "size": path.stat().st_size,
            "sha256": sha(path),
        })
    manifest = {
        "classification": "HYPOTHESIS_POST_HOC_EXPLANATORY_SENSITIVITY_PHASE2",
        "role": "fresh_independent_binding_critic_cycle6",
        "verdict": "SURVIVES",
        "branch": subprocess.check_output(["git", "branch", "--show-current"], cwd=ROOT, text=True).strip(),
        "git_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "master_seed": 20260713,
        "subseed_rule": "first 8 SHA256 bytes, unsigned big-endian, of ASCII 20260713|TASKB_PHASE2|<analysis_step_id>",
        "commands": [
            "/usr/bin/python3 experiments/taskB_grade_histology/phase2_critic_cycle6/independent_reconstruction.py",
            "/usr/bin/python3 experiments/taskB_grade_histology/phase2_critic_cycle6/independent_upstream_reconciliation.py",
            "/usr/bin/python3 experiments/taskB_grade_histology/phase2_critic_cycle6/provenance_audit.py",
            "(cd experiments/taskB_grade_histology/phase2_critic_cycle6 && /usr/bin/python3 adversarial_checks.py)",
            "/usr/bin/python3 experiments/taskB_grade_histology/phase2_critic_cycle6/build_critic_manifest.py",
        ],
        "environment": {
            "interpreter": sys.executable,
            "python": sys.version.replace("\n", " "),
            "numpy": numpy.__version__,
            "pandas": pandas.__version__,
            "scipy": scipy.__version__,
            "platform": platform.platform(),
        },
        "sealed_inputs": {
            "scores_v3.npz": sha("data/external/original-workspace/task028-freeze-b-draft/execution/intermediate/scores_v3.npz"),
            "covariates_v3.tsv": sha("data/external/original-workspace/task028-freeze-b-draft/execution/intermediate/covariates_v3.tsv"),
            "patient_order_v3.json": sha("data/external/original-workspace/task028-freeze-b-draft/execution/intermediate/patient_order_v3.json"),
            "LINKAGE_LEDGER.tsv": sha(ROOT / "experiments/taskB_grade_histology/phase1_execution/run1/LINKAGE_LEDGER.tsv"),
            "PHASE2_FROZEN_SPEC.json": sha(ROOT / "experiments/taskB_grade_histology/phase2_freeze/PHASE2_FROZEN_SPEC.json"),
            "cycle6_run1_MODEL_RESULTS.tsv": sha(CYCLE / "run1/MODEL_RESULTS.tsv"),
            "cycle6_run1_MATCHED_DECOMPOSITIONS.tsv": sha(CYCLE / "run1/MATCHED_DECOMPOSITIONS.tsv"),
            "cycle6_run1_MODEL_DIAGNOSTICS.tsv": sha(CYCLE / "run1/MODEL_DIAGNOSTICS.tsv"),
            "cycle6_run1_INFLUENCE_RECORDS.tsv": sha(CYCLE / "run1/INFLUENCE_RECORDS.tsv"),
        },
        "results": {
            "model_rows": reproduction["model_rows"],
            "decomposition_rows": reproduction["decomposition_rows"],
            "diagnostic_rows": reproduction["diagnostic_rows"],
            "influence_rows": reproduction["influence_rows"],
            "field_mismatch_count": sum(row["mismatch_count"] for row in reproduction["comparisons"]),
            "maximum_independent_producer_absolute_delta": max(row["max_absolute_difference"] for row in reproduction["comparisons"]),
            "upstream_checks": upstream["check_count"],
            "upstream_mismatches": upstream["mismatch_count"],
            "upstream_maximum_absolute_delta": upstream["max_absolute_delta"],
            "all_c2_estimates_negative": adversarial["all_model_c2_estimates_negative"],
            "leave_one_out_any_c2_sign_change": adversarial["leave_one_out_any_c2_sign_change"],
            "output_sha256_checks": provenance["output_sha256s"]["line_count"],
            "output_sha256_mismatches": provenance["output_sha256s"]["mismatch_count"],
            "phase2_input_hash_checks": provenance["phase2_input_hashes"]["row_count"],
            "phase2_input_hash_mismatches": provenance["phase2_input_hashes"]["mismatch_count"],
        },
        "artifacts": artifacts,
        "scope": "TCGA-UCEC post-hoc explanatory compatibility only; no CPTAC adjusted inference; no verdict/category/manuscript change.",
    }
    manifest_path = HERE / "CRITIC_REPRODUCIBILITY_MANIFEST.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="ascii")
    print(json.dumps({"manifest_sha256": sha(manifest_path), "artifact_count": len(artifacts)}, sort_keys=True))


if __name__ == "__main__":
    main()
