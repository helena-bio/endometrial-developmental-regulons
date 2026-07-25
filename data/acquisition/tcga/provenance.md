TASK-027 -- PROVENANCE CHAIN
============================
Every byte in this substrate is traceable to a pinned GDC identifier.

SOURCE
  GDC (api.gdc.cancer.gov), project TCGA-UCEC, workflow STAR - Counts,
  augmented_star_gene_counts.tsv, reference GENCODE v36.
  Data Release 45.0 - December 04, 2025 (major 45 / minor 0 / release_date
  2025-12-04), status commit 8f7c2a51ab0084b216ad1b62a3fae8b945439c53, tag 8.5.0.
  Verified live 2026-07-12 via GET https://api.gdc.cancer.gov/status
  (release_status.json, sha256 e0217ad3931ecad57c34fd662dc76b55dbb72056945a1a6da4e34a4a461cabb5).
  No release drift vs the pinned checkpoint.

PIN -> MANIFEST -> DOWNLOAD -> VERIFY chain
  1. PIN. The exact 542 files (507 primary + 35 normal), each with UUID + md5 +
     size + barcode, were pinned in the Vladimir-approved task025 checkpoint:
       primary: would_download_files.tsv (sha256 0801754019bf... -> 0801754019fe0b69d88980c650c4f4d604558ec05afd334b714315356acca8af)
       normal : gdc_manifest_full.tsv rows sample_type=='Solid Tissue Normal'
                (sha256 85899f606b64c9f53c0115670eb8997889c47ad833e09c6041a66cab8c5234cb)
     Both SHA-256s were re-verified before use (match).

  2. MANIFEST. scripts/build_download_manifest.py merged the two pins into ONE
     download_manifest.tsv (sha256 60a713aee692a2f3e7e0f1027c61961f60de43ddf2455e1fb9bab1f5397c35f6)
     and asserted counts/bytes/uniqueness BEFORE any download:
       507 primary + 35 normal = 542 ; sum 2,287,429,373 bytes ;
       all uuid/md5 non-empty and unique. PASS.

  3. DOWNLOAD. scripts/download_and_verify.py fetched each UUID from
     https://api.gdc.cancer.gov/data/<uuid> into originals/ with retry/backoff,
     confirming HTTP 200 and byte length == pinned size for every file.

  4. VERIFY. For every file the computed md5 was asserted == the pinned md5
     (an md5 mismatch would have been a HARD STOP with no substitution) and a
     sha256 was recorded. checksums.tsv (sha256
     efa25ed34aa2bae46c9f3f48eb663ef4d606132c1dff851adf3a12b4877465ed):
       542/542 md5_match true ; 542/542 size_match true ;
       sum(computed_size) 2,287,429,373 (exact) ; 542 distinct sha256.

CONTENT-ANCHOR (barcode binding is not merely trusted)
  STAR-Counts file bodies carry NO TCGA barcode (intrinsic to GDC). The
  uuid->aliquot_barcode binding is GDC metadata. It is content-anchored two ways:
    (a) each file's computed md5 == the pinned md5 GDC bound to that UUID/barcode;
    (b) each served Content-Disposition filename == the GDC-recorded file_name for
        that file_id in gdc_files_raw.json (sha256
        6254d90e56dec110098a293b00ba69780fe443e52b4e88654340bc91824bc0cc).
  So the bytes held ARE the file GDC associated with that UUID (hence that
  barcode). Cohort membership was re-derived from these ANCHORED materialized
  files, not from the metadata prediction.

FROZEN RULES / TABLES (applied, not re-derived)
  SELECTION_RULE.txt (sha256 b92d7478dbc2c5059160624d15ad0e492af045e23cbe2ed2a4e8813283b4ad87)
    applied verbatim: '01' primary restriction + lexicographically-lowest-barcode
    patient dedup.
  subtype_normalized.tsv (sha256 17a899a4869fccf806859342454bde440027786a48bdd50f1d1eace810b1ef9e)
    read-only JOIN; subtypes NOT re-derived/re-mapped/re-called.

DETERMINISM
  No randomness enters any step. seed = N/A. Re-running the three scripts on the
  same originals/ reproduces every count and checksum bit-for-bit.

GIT
  HEAD start == end == 83503bad47b60193598b2b9ebe819c22c83e8ac1 (branch
  experiment/embryonic-atlas). No git write op performed. Artifacts live outside
  the revgate git tree.

ENVIRONMENT
  /usr/bin/python3 3.12.3, standard library only (no third-party packages, no
  venv). curl 8.5.0 for the GDC /status GET and a one-file connectivity/md5 test.
