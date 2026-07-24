#!/usr/bin/env python3
"""Outcome-blind byte and schema preflight; never emits patient values."""

import csv
import hashlib
import json
from pathlib import Path


FILES = {
    "data/external/original-workspace/task027-acquire-freeze-a/freeze_a_redux/cohort_selected_primary.tsv": ("e6f9d390730bff4248379b2280b12e2914ec800e2f56516b466f2151b169f3be", "tsv"),
    "data/external/original-workspace/task024-freeze-a/subtype_normalized.tsv": ("17a899a4869fccf806859342454bde440027786a48bdd50f1d1eace810b1ef9e", "tsv"),
    "data/external/original-workspace/task029-external-replication-feasibility/join_tables/discovery_join.json": ("315a9d29f1add19fd90db85188762937ac2e49262d64a51813aca5fb5ccb48ad", "json"),
    "data/external/original-workspace/task029-external-replication-feasibility/join_tables/confirmatory_join.json": ("95b9f754970e7f0691faddf6e350aca6cfe465c5c0c98600a4ed34ebbfe43743", "json"),
    "data/external/original-workspace/task029-external-replication-feasibility/sources/cbio_2020_clinattr.json": ("25161910befe62fea99ba5d754e0e43f8d2d80413538b8bfcac0c52aa5e58e6e", "json"),
    "data/external/original-workspace/task029-external-replication-feasibility/sources_task029gate/uec_cptac_gdc_clinattr.json": ("96794b79b83a663a9caf5dc51dce3d71a29983040a2b4125f972a100c7627628", "json"),
    "data/external/original-workspace/task029-external-replication-feasibility/sources/pdc_clin_PDC000125.json": ("61d05e1298312fbdf98e23b027b93e9b0310612baf1fc9fb6fe088913bc08291", "json"),
    "data/external/original-workspace/task029-external-replication-feasibility/sources/pdc_clin_PDC000439.json": ("6ff6418b4f201d7f2db3b1c69617d3435a01b4877cf37a06f4c56a4fb7ba4012", "json"),
    "data/external/original-workspace/task029-external-replication-feasibility/sources/pdc_discovery_biospec.json": ("c98f51f1ec9916252a210c40fae93adfe2ae42e998decdac5adbf82c8bd3b74e", "json"),
    "data/external/original-workspace/task029-external-replication-feasibility/sources/pdc_confirmatory_biospec.json": ("bc159c1d0500cce635b7631a069c814a9f3c802f947dc3c9bf828573e0b7aae9", "json"),
    "data/external/original-workspace/task029-external-replication-feasibility/sources_task029gate/gdc_case_mapping.json": ("4e5fd1736d2b5a17f6cd66929953a6893f4858a0a480480faa26ccafd2e798c8", "json"),
    "data/external/original-workspace/task029-external-replication-feasibility/sources_task029gate/gdc_rna_primary_by_case.json": ("87514a428a9cf4810621ff7c46f5c6696f78bda2e3fecf66e7fe9920cbe9fbda", "json"),
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def keysets(obj, depth=0, dynamic_map=False):
    """Return container shapes without emitting dynamic identifier keys."""
    out = []
    if isinstance(obj, dict):
        keys = sorted(obj)
        looks_dynamic = dynamic_map or (len(keys) > 20)
        out.append({
            "depth": depth,
            "kind": "dict",
            "key_count": len(keys),
            "keys": "REDACTED_DYNAMIC_IDENTIFIERS" if looks_dynamic else keys,
        })
        for value in obj.values():
            if isinstance(value, (dict, list)):
                out.extend(keysets(value, depth + 1, looks_dynamic))
    elif isinstance(obj, list):
        out.append({"depth": depth, "kind": "list", "length": len(obj)})
        seen = set()
        for value in obj:
            if isinstance(value, dict):
                keys = tuple(sorted(value))
                if keys not in seen:
                    seen.add(keys)
                    out.append({
                        "depth": depth + 1,
                        "kind": "record",
                        "keys": "REDACTED_DYNAMIC_IDENTIFIERS" if len(keys) > 20 else list(keys),
                    })
                for nested in value.values():
                    if isinstance(nested, (dict, list)):
                        out.extend(keysets(nested, depth + 2))
            elif isinstance(value, list):
                out.extend(keysets(value, depth + 1))
    return out


def main():
    report = []
    for name, (expected, kind) in FILES.items():
        path = Path(name)
        actual = sha(path)
        if actual != expected:
            raise SystemExit(f"HASH_MISMATCH:{name}")
        item = {"path": name, "sha256": actual, "hash_ok": True, "kind": kind}
        if kind == "tsv":
            with path.open("r", encoding="utf-8", newline="") as handle:
                item["columns"] = next(csv.reader(handle, delimiter="\t"))
        else:
            with path.open("r", encoding="utf-8") as handle:
                item["container_schema"] = keysets(json.load(handle))
        report.append(item)
    out = Path(__file__).resolve().parent / "PREFLIGHT_SCHEMA_REPORT.json"
    out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="ascii")
    print("PREFLIGHT_OK")


if __name__ == "__main__":
    main()
