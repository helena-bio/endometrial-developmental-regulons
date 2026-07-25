#!/usr/bin/env python3
"""
TASK-029 ACQUISITION step 01 -- authoritative GDC CPTAC-3 UCEC STAR-Counts file query.

Queries GDC (live, Data Release pinned by status probe) for CPTAC-3 uterine
Transcriptome-Profiling / Gene-Expression-Quantification / STAR-Counts files.
Records file_id, file_name, md5sum, file_size, cases.submitter_id, and full
sample metadata (sample_submitter_id, sample_type, tissue_type, aliquot).

NO expression value is read here -- metadata only. Deterministic (no randomness).
Writes the raw GDC response and query manifest to the configured acquisition directory.

Git HEAD unchanged: 83503bad47b60193598b2b9ebe819c22c83e8ac1.
"""
import json
import os
import sys
import hashlib
import urllib.request
from datetime import datetime, timezone

BASE = os.path.dirname(
    os.path.dirname(os.path.abspath(__file__))
)
WORK = os.environ.get(
    "CPTAC_WORK_DIR",
    os.path.join(BASE, "work"),
)
OUT = os.environ.get(
    "CPTAC_ACQUISITION_OUTPUT",
    os.path.join(WORK, "acquisition"),
)
os.makedirs(OUT, exist_ok=True)
GDC_FILES = "https://api.gdc.cancer.gov/files"
GDC_STATUS = "https://api.gdc.cancer.gov/status"


def utcnow():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def post_json(url, payload):
    data = json.dumps(payload).encode()
    req = urllib.request.Request(url, data=data,
                                 headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.loads(r.read().decode())


def get_json(url):
    with urllib.request.urlopen(url, timeout=60) as r:
        return json.loads(r.read().decode())


status = get_json(GDC_STATUS)
release = status["data_release"]
print("GDC data_release:", release, utcnow(), flush=True)

# Full expand of sample + aliquot metadata for exact join resolution.
payload = {
    "filters": {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "cases.project.project_id",
                                     "value": ["CPTAC-3"]}},
            {"op": "in", "content": {"field": "cases.primary_site",
                                     "value": ["Uterus, NOS", "Corpus uteri"]}},
            {"op": "in", "content": {"field": "data_category",
                                     "value": ["Transcriptome Profiling"]}},
            {"op": "in", "content": {"field": "data_type",
                                     "value": ["Gene Expression Quantification"]}},
            {"op": "in", "content": {"field": "analysis.workflow_type",
                                     "value": ["STAR - Counts"]}},
        ],
    },
    "format": "JSON",
    "size": "5000",
    "fields": ",".join([
        "file_id", "file_name", "md5sum", "file_size", "data_format",
        "state", "access", "created_datetime", "updated_datetime",
        "analysis.workflow_type",
        "cases.submitter_id", "cases.case_id", "cases.project.project_id",
        "cases.primary_site",
        "cases.samples.submitter_id", "cases.samples.sample_id",
        "cases.samples.sample_type", "cases.samples.sample_type_id",
        "cases.samples.tissue_type",
        "cases.samples.portions.analytes.aliquots.submitter_id",
        "cases.samples.portions.analytes.aliquots.aliquot_id",
    ]),
}

resp = post_json(GDC_FILES, payload)
hits = resp["data"]["hits"]
pag = resp["data"]["pagination"]
print("returned hits:", len(hits), "total:", pag["total"], utcnow(), flush=True)
assert len(hits) == pag["total"], "pagination incomplete; increase size"

manifest = {
    "queried_at": utcnow(),
    "gdc_data_release": release,
    "gdc_status": status,
    "query_payload": payload,
    "n_files": len(hits),
    "endpoint": GDC_FILES,
}
with open(os.path.join(OUT, "acq01_gdc_star_counts_raw.json"), "w") as f:
    json.dump(resp, f, indent=2)
with open(os.path.join(OUT, "acq01_gdc_query_manifest.json"), "w") as f:
    json.dump(manifest, f, indent=2)

# checksum of the raw response for provenance
h = hashlib.sha256(json.dumps(resp, sort_keys=True).encode()).hexdigest()
print("raw response sha256 (sorted-keys):", h, flush=True)
print("wrote acq01 outputs to", OUT, flush=True)
