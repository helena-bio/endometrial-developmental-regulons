TASK-027 -- ACQUIRE PINNED GDC DR45.0 COHORT + FREEZE A FROM SCRATCH (HUMAN SUMMARY)
====================================================================================
Directed by Vladimir 2026-07-12: Toil/option-A is OUT; acquire the PINNED GDC DR45
cohort (PRIMARY=C, single GDC STAR pipeline), re-run the INPUT feasibility freeze
(Freeze A) from scratch on the materialized files, then STOP before Freeze B.
NO biological scoring. Git must remain unchanged.

VERDICT: COMPLETE. Acquisition clean (542/542 md5+size exact). Freeze A redux =
GREEN (smallest mandatory arm = POLE = 49 >= 40). Reported AS DERIVED; no number
forced. Freeze B NOT entered. This document does NOT interpret whether the result
is "good".

------------------------------------------------------------------
0. INTEGRITY (the hard constraints)
------------------------------------------------------------------
git HEAD start == end == 83503bad47b60193598b2b9ebe819c22c83e8ac1
  branch experiment/embryonic-atlas. NO git write op of any kind (no
  add/commit/checkout/branch/reset). The acquisition artifacts were generated in an isolated workspace outside
  the analysis repository.
NO expression VALUE (tpm/count/fpkm) read into any score/contrast/module/
  ranking/fold-change/statistic. Structural checks only (columns, 4 STAR summary
  rows, GENCODE-v36 gene-row count). Grep-auditable: no script sums, ranks, or
  contrasts a numeric expression column.
The acquisition step did not modify the analysis source tree.

------------------------------------------------------------------
1. FROZEN INPUTS -- SHA-256 VERIFIED BEFORE USE (all match)
------------------------------------------------------------------
SELECTION_RULE.txt        b92d7478...ad87  OK  (applied verbatim)
subtype_normalized.tsv    17a899a4...1ef9e OK  (read-only JOIN; not re-derived)
would_download_files.tsv  08017540...caf   OK  (507 primary pins)
gdc_manifest_full.tsv     85899f60...234b  OK  (35 normal pins)
gdc_files_raw.json        6254d90e...c0cc  OK  (content-anchor cross-check)

------------------------------------------------------------------
2. RELEASE IDENTITY (STEP 2) -- NO DRIFT
------------------------------------------------------------------
Live GET https://api.gdc.cancer.gov/status returned EXACTLY:
  data_release "Data Release 45.0 - December 04, 2025"
  major 45 / minor 0 / release_date 2025-12-04
  commit 8f7c2a51ab0084b216ad1b62a3fae8b945439c53  tag 8.5.0
Identical to the pinned checkpoint. The pinned UUIDs all resolved and served
byte-exact. Raw JSON: release_status.json.

------------------------------------------------------------------
3. DOWNLOAD MANIFEST (STEP 1) -- ASSERTIONS PASS
------------------------------------------------------------------
download_manifest.tsv: 507 primary + 35 normal = 542 rows.
  n_primary==507, n_normal==35, n_total==542                 PASS
  sum(file_size_bytes)==2,287,429,373                        PASS
    primary 2,139,373,015 ; normal 148,056,358               PASS
  all file_uuid/md5sum non-empty and unique                  PASS

------------------------------------------------------------------
4. MATERIALIZATION + CHECKSUMS (STEP 3/4) -- 542/542 CLEAN
------------------------------------------------------------------
originals/ holds 542 files (<uuid>__augmented_star_gene_counts.tsv).
checksums.tsv:
  files verified          542 / 542
  md5_match true          542 / 542   (each computed md5 == pinned md5)
  size_match true         542 / 542   (each byte length == pinned size)
  sha256 recorded         542 / 542   (all distinct)
  sum(computed_size)      2,287,429,373  (EXACT)
  HARD STOPS              NONE (0 md5 mismatch, 0 size mismatch, 0 non-200)
An md5 mismatch would have been a HARD STOP with no substitution; none occurred.
Originals are PRESERVED as the permanent, provenance-complete substrate.

------------------------------------------------------------------
5. FREEZE A FROM SCRATCH (STEP 5) -- ON THE MATERIALIZED PRIMARY FILES (gate = 507 primary only)
------------------------------------------------------------------
Cohort membership was re-derived FROM THE ACTUAL DOWNLOADED FILES: every on-disk
UUID is confirmed present, md5_match=true, and CONTENT-ANCHORED (the served GDC
filename equals the GDC-recorded file_name for that file_id, so the bytes held are
the file GDC bound to that UUID/barcode). STAR-Counts bodies carry no barcode --
the uuid->barcode binding is metadata-intrinsic to GDC; it is anchored via md5 +
GDC filename, not merely trusted from the prediction.

FROZEN SELECTION RULE applied verbatim to the materialized primary barcodes:
  STEP1 keep sample-type code '01' : 507 kept, 0 excluded non-01
  STEP2 one per 12-char patient, keep lexicographically-lowest full aliquot barcode:
        507 independent patients; 0 patients with >1 primary; 0 files dropped.
  DEDUP PROVENANCE (disclosed, not a deviation): the tie-break fired 0 times HERE
  because the pinned 507 set is already one-per-patient. task025 applied the SAME
  frozen rule to the FULL 553-primary cohort; its 4 tie patients (TCGA-BK-A0CA,
  -A139, -A26L, -A0CC; 8 files dropped) each appear in the 507 as the KEPT
  lowest-barcode winner (verified). All 4 are NSMP/p53abn, none POLE -> the
  tie-break lifts no marginal arm (consistent with a blind rule).

READ-ONLY JOIN to the frozen subtype table:
  matched               507
  primary-no-subtype      0
  subtype-no-primary      0    (every pinned subtype patient has a primary file)

PER-SUBTYPE (matched, independent patients) -- DERIVED, NOT FORCED:
  POLE     49
  MMRd    148
  NSMP    147
  p53abn  163
  smallest mandatory arm = 49 (POLE)

STRUCTURAL INTEGRITY (STEP 5d) -- 507/507 valid:
  every primary file: '# gene-model: GENCODE v36' header; the 9 expected data
  columns (gene_id,gene_name,gene_type,unstranded,stranded_first,stranded_second,
  tpm_unstranded,fpkm_unstranded,fpkm_uq_unstranded); the 4 STAR summary rows
  (N_unmapped,N_multimapping,N_noFeature,N_ambiguous); gene-row count == 60660.
  Distinct gene-row counts across all 507 files: {60660}. STRUCTURAL ONLY.

GATE BAND (STEP 5e, pre-registration bands): GREEN
  GREEN >=40 in every mandatory arm ; YELLOW smallest 20-39 ; RED any <20.
  Smallest mandatory arm = 49 (POLE) -> GREEN.
  POLE = 49 >= 20: YES. All 49 counted POLE patients have a materialized,
  md5-verified, GENCODE-v36-valid primary file that survives '01'+dedup+join
  (pole_all_files_md5_and_struct_valid = True; per-patient detail in
  freeze_a_redux/pole_patient_detail.json).
  Contrast (for the record only, not an interpretation): the original on-host
  Xena/TOIL Freeze A (task024) was RED with POLE=16. The pinned GDC cohort
  materializes POLE=49. No number was tuned to reach any target.

------------------------------------------------------------------
6. THE 35 NORMALS (STEP 5f) -- VERIFY + PRESERVE ONLY
------------------------------------------------------------------
35 Solid-Tissue-Normal STAR-Counts files verified+preserved; all md5 true.
Case barcodes + any same-patient subtype label recorded in
freeze_a_redux/normals_preserved.tsv.
CONSTRAINT (flagged explicitly): adjacent/peri-tumoral Solid-Tissue-Normal,
potentially affected by field effect + surgical ascertainment + small N.
RESERVED for a DEFERRED limited-B sensitivity. NOT in the gate, selection,
thresholds, or labels. NOT used now.

------------------------------------------------------------------
7. HALT (STEP 7)
------------------------------------------------------------------
STOPPED before Freeze B. No module/score/contrast/fold-change computed. Nothing
staged, nothing committed. git HEAD end == start == 83503bad...8ac1.

PROHIBITIONS ENFORCED: no subtype merges; no histology fallback; no
recurrent/metastatic in the cohort (0 excluded non-01, the 1 recurrent file is
not in the pinned 507); no multiple aliquots as independent observations; no
post-hoc redesign after inspection; no threshold tuning; POLE not forced.

ARTIFACTS: MANIFEST.json (machine), this file, PROVENANCE.md, download_manifest.tsv,
release_status.json, checksums.tsv, originals/ (542 files), freeze_a_redux/
(summary.json, join_matched.tsv, cohort_selected_primary.tsv, dedup_ties.tsv,
dropped_by_collapse.tsv, excluded_not_primary01.tsv, primary_no_subtype.tsv,
subtype_no_primary.tsv, structural_integrity_primary.tsv, pole_patient_detail.json,
normals_preserved.tsv), scripts/.
