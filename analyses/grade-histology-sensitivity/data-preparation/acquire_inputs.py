#!/usr/bin/env python3
"""Acquire the six exact frozen clinical-only queries twice and lock bytes."""

import argparse
import hashlib
import json
import os
import shutil
import urllib.request
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path("data/external/original-workspace/revgate-task-b-grade-histology")
TASK = ROOT / "experiments/taskB_grade_histology"
EXEC = TASK / "phase1_execution"
INPUT = TASK / "inputs"

CBIO = [
    ("tcga_cbio_grade", "https://www.cbioportal.org/api/studies/ucec_tcga_pan_can_atlas_2018/clinical-data?clinicalDataType=SAMPLE&attributeId=GRADE&pageSize=100000&pageNumber=0&projection=SUMMARY", "GRADE", {"clinicalAttributeId", "patientId", "sampleId", "studyId", "uniquePatientKey", "uniqueSampleKey", "value"}, INPUT / "tcga_cbio_grade.json"),
    ("tcga_cbio_icdo_histology", "https://www.cbioportal.org/api/studies/ucec_tcga_pan_can_atlas_2018/clinical-data?clinicalDataType=PATIENT&attributeId=ICD_O_3_HISTOLOGY&pageSize=100000&pageNumber=0&projection=SUMMARY", "ICD_O_3_HISTOLOGY", {"clinicalAttributeId", "patientId", "studyId", "uniquePatientKey", "value"}, INPUT / "tcga_cbio_icdo_histology.json"),
    ("cptac_discovery_cbio_grade", "https://www.cbioportal.org/api/studies/ucec_cptac_2020/clinical-data?clinicalDataType=PATIENT&attributeId=HISTOLOGIC_GRADE_FIGO&pageSize=100000&pageNumber=0&projection=SUMMARY", "HISTOLOGIC_GRADE_FIGO", {"clinicalAttributeId", "patientId", "studyId", "uniquePatientKey", "value"}, INPUT / "cptac_discovery_grade.json"),
    ("cptac_discovery_cbio_histology", "https://www.cbioportal.org/api/studies/ucec_cptac_2020/clinical-data?clinicalDataType=PATIENT&attributeId=HISTOLOGIC_TYPE&pageSize=100000&pageNumber=0&projection=SUMMARY", "HISTOLOGIC_TYPE", {"clinicalAttributeId", "patientId", "studyId", "uniquePatientKey", "value"}, INPUT / "cptac_discovery_histology.json"),
]

PDC = [
    ("pdc_PDC000125", "PDC000125", "c935c587-0cd1-11e9-a064-0a9c39d33490", '{ clinicalPerStudy(study_id: "c935c587-0cd1-11e9-a064-0a9c39d33490") { case_id case_submitter_id morphology primary_diagnosis tumor_grade } }', INPUT / "pdc_PDC000125_grade_histology.json"),
    ("pdc_PDC000439", "PDC000439", "401b6a4e-e36e-4bd2-be2a-2926eaa44d88", '{ clinicalPerStudy(study_id: "401b6a4e-e36e-4bd2-be2a-2926eaa44d88") { case_id case_submitter_id morphology primary_diagnosis tumor_grade } }', INPUT / "pdc_PDC000439_grade_histology.json"),
]

PDC_KEYS = {"case_id", "case_submitter_id", "morphology", "primary_diagnosis", "tumor_grade"}


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def utcnow():
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def fetch(req):
    with urllib.request.urlopen(req, timeout=120) as response:
        return response.read(), {k.lower(): v for k, v in response.headers.items()}, response.status


def validate_cbio(raw, attribute):
    payload = json.loads(raw)
    if not isinstance(payload, list):
        raise RuntimeError("CBIO_WRAPPER_NOT_LIST")
    for row in payload:
        if set(row) != next(x[3] for x in CBIO if x[2] == attribute):
            raise RuntimeError(f"CBIO_EXTRA_OR_MISSING_KEYS:{attribute}:{sorted(row)}")
        if row["clinicalAttributeId"] != attribute:
            raise RuntimeError(f"CBIO_MIXED_ATTRIBUTE:{attribute}")
    return len(payload)


def validate_pdc(raw):
    payload = json.loads(raw)
    if set(payload) != {"data"} or set(payload["data"]) != {"clinicalPerStudy"}:
        raise RuntimeError("PDC_WRAPPER_KEYS_MIXED")
    rows = payload["data"]["clinicalPerStudy"]
    if not isinstance(rows, list):
        raise RuntimeError("PDC_RECORDS_NOT_LIST")
    for row in rows:
        if set(row) != PDC_KEYS:
            raise RuntimeError(f"PDC_EXTRA_OR_MISSING_KEYS:{sorted(row)}")
    return len(rows)


def canonical(raw):
    return json.dumps(json.loads(raw), ensure_ascii=True, sort_keys=True, separators=(",", ":")).encode("ascii") + b"\n"


def one_round(round_no):
    raw_dir = EXEC / "acquisition" / f"raw_round_{round_no}"
    raw_dir.mkdir(parents=True, exist_ok=False)
    records = []
    for name, url, attribute, _keys, final_path in CBIO:
        req = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": "RevGate-TaskB-Phase1/1.0"})
        raw, headers, status = fetch(req)
        count = validate_cbio(raw, attribute)
        path = raw_dir / f"{name}.json"
        path.write_bytes(raw)
        records.append({"name": name, "authority": "cBioPortal public API", "request": url, "method": "GET", "status": status, "headers": headers, "retrieval_utc": utcnow(), "raw_path": str(path), "raw_size": len(raw), "raw_sha256": sha(raw), "canonical_sha256": sha(canonical(raw)), "record_count": count, "requested_attribute_id": attribute, "study_version_pin": name.split("_cbio")[0], "final_path": str(final_path)})
    for name, pdc_id, study_uuid, query, final_path in PDC:
        body = json.dumps({"query": query}, separators=(",", ":")).encode("ascii")
        req = urllib.request.Request("https://pdc.cancer.gov/graphql", data=body, headers={"Accept": "application/json", "Content-Type": "application/json", "User-Agent": "RevGate-TaskB-Phase1/1.0"}, method="POST")
        raw, headers, status = fetch(req)
        count = validate_pdc(raw)
        path = raw_dir / f"{name}.json"
        path.write_bytes(raw)
        records.append({"name": name, "authority": "Proteomic Data Commons public API", "request": "https://pdc.cancer.gov/graphql", "graphql_query": query, "method": "POST", "status": status, "headers": headers, "retrieval_utc": utcnow(), "raw_path": str(path), "raw_size": len(raw), "raw_sha256": sha(raw), "canonical_sha256": sha(canonical(raw)), "record_count": count, "pdc_study_id": pdc_id, "study_uuid_version_pin": study_uuid, "final_path": str(final_path)})
    return records


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--round", type=int, choices=(1, 2), required=True)
    args = parser.parse_args()
    INPUT.mkdir(parents=True, exist_ok=True)
    records = one_round(args.round)
    meta = EXEC / "acquisition" / f"round_{args.round}_metadata.json"
    meta.write_text(json.dumps(records, indent=2, sort_keys=True) + "\n", encoding="ascii")
    print(f"ACQUISITION_ROUND_{args.round}_STRUCTURE_OK")


if __name__ == "__main__":
    main()
