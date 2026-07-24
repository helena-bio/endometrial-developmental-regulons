#!/usr/bin/env python3
"""Independent exact field-restricted Phase-1 acquisitions.

This program intentionally contains the six frozen requests explicitly and does
not import producer code.  It validates wrapper and record keys before examining
or reporting values.  It prints metadata and hashes only.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import pathlib
import urllib.request


CBIO = [
    (
        "tcga_cbio_grade",
        "https://www.cbioportal.org/api/studies/ucec_tcga_pan_can_atlas_2018/clinical-data?clinicalDataType=SAMPLE&attributeId=GRADE&pageSize=100000&pageNumber=0&projection=SUMMARY",
        "GRADE",
        {"clinicalAttributeId", "patientId", "sampleId", "studyId", "uniquePatientKey", "uniqueSampleKey", "value"},
    ),
    (
        "tcga_cbio_icdo_histology",
        "https://www.cbioportal.org/api/studies/ucec_tcga_pan_can_atlas_2018/clinical-data?clinicalDataType=PATIENT&attributeId=ICD_O_3_HISTOLOGY&pageSize=100000&pageNumber=0&projection=SUMMARY",
        "ICD_O_3_HISTOLOGY",
        {"clinicalAttributeId", "patientId", "studyId", "uniquePatientKey", "value"},
    ),
    (
        "cptac_discovery_cbio_grade",
        "https://www.cbioportal.org/api/studies/ucec_cptac_2020/clinical-data?clinicalDataType=PATIENT&attributeId=HISTOLOGIC_GRADE_FIGO&pageSize=100000&pageNumber=0&projection=SUMMARY",
        "HISTOLOGIC_GRADE_FIGO",
        {"clinicalAttributeId", "patientId", "studyId", "uniquePatientKey", "value"},
    ),
    (
        "cptac_discovery_cbio_histology",
        "https://www.cbioportal.org/api/studies/ucec_cptac_2020/clinical-data?clinicalDataType=PATIENT&attributeId=HISTOLOGIC_TYPE&pageSize=100000&pageNumber=0&projection=SUMMARY",
        "HISTOLOGIC_TYPE",
        {"clinicalAttributeId", "patientId", "studyId", "uniquePatientKey", "value"},
    ),
]

PDC = [
    (
        "pdc_PDC000125",
        '{ clinicalPerStudy(study_id: "c935c587-0cd1-11e9-a064-0a9c39d33490") { case_id case_submitter_id morphology primary_diagnosis tumor_grade } }',
    ),
    (
        "pdc_PDC000439",
        '{ clinicalPerStudy(study_id: "401b6a4e-e36e-4bd2-be2a-2926eaa44d88") { case_id case_submitter_id morphology primary_diagnosis tumor_grade } }',
    ),
]
PDC_KEYS = {"case_id", "case_submitter_id", "morphology", "primary_diagnosis", "tumor_grade"}


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def canonical(obj: object) -> bytes:
    return (json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=True) + "\n").encode("ascii")


def get(url: str) -> tuple[bytes, dict[str, str], int]:
    req = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": "RevGate-Phase1-Independent-Critic/1.0"})
    with urllib.request.urlopen(req, timeout=120) as res:
        return res.read(), dict(res.headers.items()), res.status


def post(query: str) -> tuple[bytes, dict[str, str], int]:
    body = json.dumps({"query": query}, separators=(",", ":")).encode("ascii")
    req = urllib.request.Request(
        "https://pdc.cancer.gov/graphql",
        data=body,
        method="POST",
        headers={"Accept": "application/json", "Content-Type": "application/json", "User-Agent": "RevGate-Phase1-Independent-Critic/1.0"},
    )
    with urllib.request.urlopen(req, timeout=120) as res:
        return res.read(), dict(res.headers.items()), res.status


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--round", required=True, type=int)
    ap.add_argument("--out", required=True, type=pathlib.Path)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    summary: list[dict[str, object]] = []

    for name, url, attribute, allowed_keys in CBIO:
        raw, headers, status = get(url)
        obj = json.loads(raw)
        if not isinstance(obj, list):
            raise RuntimeError(f"{name}: wrapper is not a list")
        for i, record in enumerate(obj):
            if not isinstance(record, dict) or set(record) != allowed_keys:
                raise RuntimeError(f"{name}: record {i} key mismatch")
        if any(record["clinicalAttributeId"] != attribute for record in obj):
            raise RuntimeError(f"{name}: unexpected clinicalAttributeId")
        can = canonical(obj)
        (args.out / f"{name}.raw.json").write_bytes(raw)
        (args.out / f"{name}.canonical.json").write_bytes(can)
        summary.append({"name": name, "http_status": status, "record_count": len(obj), "raw_bytes": len(raw), "raw_sha256": sha(raw), "canonical_sha256": sha(can), "response_headers": headers})

    for name, query in PDC:
        raw, headers, status = post(query)
        obj = json.loads(raw)
        if not isinstance(obj, dict) or set(obj) != {"data"}:
            raise RuntimeError(f"{name}: top-level wrapper key mismatch")
        data = obj["data"]
        if not isinstance(data, dict) or set(data) != {"clinicalPerStudy"}:
            raise RuntimeError(f"{name}: data wrapper key mismatch")
        records = data["clinicalPerStudy"]
        if not isinstance(records, list):
            raise RuntimeError(f"{name}: records are not a list")
        for i, record in enumerate(records):
            if not isinstance(record, dict) or set(record) != PDC_KEYS:
                raise RuntimeError(f"{name}: record {i} key mismatch")
        can = canonical(obj)
        (args.out / f"{name}.raw.json").write_bytes(raw)
        (args.out / f"{name}.canonical.json").write_bytes(can)
        summary.append({"name": name, "http_status": status, "record_count": len(records), "raw_bytes": len(raw), "raw_sha256": sha(raw), "canonical_sha256": sha(can), "response_headers": headers})

    out = {"round": args.round, "requests": summary}
    (args.out / "acquisition_metadata.json").write_bytes(canonical(out))
    print(json.dumps({"round": args.round, "all_key_valid": True, "request_count": len(summary), "metadata_sha256": sha(canonical(out))}, sort_keys=True))


if __name__ == "__main__":
    main()
