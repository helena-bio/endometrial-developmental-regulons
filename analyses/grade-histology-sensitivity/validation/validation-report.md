SURVIVES

I cannot break the cycle-6 scientific or provenance package. No FATAL or
relevant unmitigated MAJOR objection remains. The prior chronology MAJOR is
resolved by a new pre-outcome seal and a complete raw command record, and the
cycle-5 NumPy-bool serialization failure is closed. Two MINOR reporting defects
remain and are stated below; neither can change the numerical result or its
strict explanatory boundary.

## Binding decision

- Classification is `HYPOTHESIS_POST_HOC_EXPLANATORY_SENSITIVITY_PHASE2`.
- Branch and HEAD are exact:
  `experiment/task-b-grade-histology-sensitivity` at
  `83503bad47b60193598b2b9ebe819c22c83e8ac1`.
- This pass establishes only a reproducible TCGA-UCEC post-hoc explanatory
  sensitivity package. It does not alter any frozen TCGA/CPTAC verdict, target
  category, or claim.
- Producer output freeze is valid. All outcome-bearing code, configuration,
  inputs, thresholds, taxonomy, and model nodes were fixed before the first real
  outcome load. I found no post-seal scientific edit or retry.

## Correction-by-correction resolution

1. **Prior completed-package chronology MAJOR: RESOLVED.** Cycle 6 is a new
   clean-run directory. `run_phase2.py` SHA-256
   `4f5d4959a5b4bf6964ff7858eefb8cf181ced2de93fbc98d00c6c5058265c2aa`,
   config SHA-256
   `45d1f85dd7e23f9cfc3447625cbb97ee94f41b122ffba09cd8a06a7c447e01f4`,
   and all 16 authorized inputs were sealed at 2026-07-23T10:47:05Z.
   The seal SHA-256 is
   `a829e3664ab6c73e8e95c805d055d0e35b90852b77fc8224cc6fb33ff5962cc0`.
   Real outcome access first occurs in the traced run 1 after that seal.
   There is no unpreserved real-outcome failure, no post-outcome code repair,
   and no hidden retry suggested by the command IDs, file birth/modify times,
   process records, or open traces.

2. **Cycle-5 NumPy-bool failure: RESOLVED.** The exact cycle-5 script/config/
   seal hashes regenerate as
   `15b957cd206f4c5c3a3ce041f0968b97cd862f0d65d9c0f6af663bb9375f70f3`,
   `45d1f85dd7e23f9cfc3447625cbb97ee94f41b122ffba09cd8a06a7c447e01f4`,
   and `4aad1333990d7b006c8bca8c05ede35b1a4086af4d7d5c66a507a4d0d1689fb9`.
   The actual cycle-5 to cycle-6 diff is byte-identical to
   `CYCLE5_TO_CYCLE6_CODE_DIFF.patch`. The production change is only recursive
   JSON container/scalar normalization; the additional changes are synthetic
   serialization coverage and temporary-file hygiene. All ten complete
   pre-JSON cycle-5 scientific TSV artifacts are byte-identical to cycle 6,
   showing no numerical change. Cycle-6 `SCIENTIFIC_RESULTS.json` parses with
   60 models, 24 decompositions, and 60 diagnostics, and its formerly failing
   `single_case_direction_change` values deserialize as Python booleans.

3. **Prior wording MINOR: NOT RESOLVED.** This remains presentation-only and is
   not a scientific gate failure; exact deficiencies are below.

## Independent numerical reproduction

I reconstructed the analysis from the sealed score NPZ, covariate TSV, patient
order, Phase-1 linkage ledger, and frozen Phase-2 specification. I did not
import, execute, copy, or call producer `run_phase2.py` or any producer
result-generating function.

- Coverage: 60/60 OLS rows, 24/24 matched decompositions, 60/60 diagnostic
  rows, 5,231/5,231 flagged influence rows, 144 deterministic seed derivations,
  and all six exact node rosters.
- Field comparison: 23,468 numeric and 540 textual fields. There are zero
  mismatches. Model, diagnostic, and decomposition fields are exactly equal;
  the maximum absolute influence-field difference is
  `1.3877787807814457e-17`, with maximum relative difference
  `6.362744937852376e-16`.
- All 2,000 bootstraps are valid for every model row. All 2,000 paired
  bootstraps are valid for every decomposition. The SHA-derived seeds reproduce
  the producer confidence intervals and diagnostic permutation p values.
- Maximum decomposition algebra error is
  `1.214306433183765e-16`, well below the frozen `1e-12` tolerance.
- All 60 designs have rank equal to p and status PASS. Diagnostic maxima
  regenerate exactly: condition number `4.400066109474329`, VIF
  `3.7163612061832887`, adjusted subtype GVIF `1.1544599442039043`, hat
  `0.2026518224135735`, and single-case absolute d change
  `0.04213103957486519`. There are zero leave-one-case C2 direction changes.
- Taxonomy reproduces as 16 `largely_retained` and 8 `amplified`. No q value or
  confirmatory category is present.

Exact rosters and subtype counts regenerate:

- primary CPE: full 506 (49/148/146/163), histology matched 505
  (49/148/146/162), endometrioid-grade matched 380 (44/142/138/56);
- no purity: full 507 (49/148/147/163), histology matched 506
  (49/148/147/162), endometrioid-grade matched 381 (44/142/139/56).

The sole primary/no-purity roster difference is the frozen CPE-incomplete NSMP
patient `TCGA-BS-A0TG`. Histology is 382 endometrioid, 124
non-endometrioid, and 1 missing; grade is 207 low, 291 high, and 9 missing.
There is no patient deletion beyond the frozen complete-case/node filters, no
imputation, no subtype repair, no grade/histology collapse beyond the Phase-1
map, and no unsupported model.

## Chronology, leakage, and network audit

The immutable raw audit contains 106 contiguous command IDs. Every command has
matching metadata, start/end time, exit status, stdout/stderr bytes, and
stdout/stderr SHA-256.

- 10:34:12Z: command logging starts after the disclosed directory creation and
  apply-patch bootstrap.
- 10:36:50Z to 10:38:52Z: cycle-5 code/config are copied and the JSON-only
  change is applied. No real outcome is opened.
- 10:39:48Z and 10:42:06Z: two synthetic computations print PASS, but tracer
  tee helpers deadlock; the retained commands terminate exit 137. Their traces
  contain no score, covariate-row, patient-order, prior-result, or network
  access.
- Between commands 0045 and 0046, direct invocation of the wrapper fails with
  Permission denied after apply-patch reset its execute bit. The exact command
  and error are disclosed. The wrapper cannot start and cannot load an outcome.
  The next bash-invoked wrapper command restores executable mode.
- 10:44:42Z: the third traced synthetic test succeeds cleanly. All three
  pre-seal traces independently have zero forbidden occurrences and zero
  internet connects.
- 10:47:04Z to 10:47:05Z: pre-outcome seal creation validates 16/16 input
  hashes. Real outcome, covariate rows, and patient-order contents are still not
  loaded; their bytes are SHA-read only.
- 10:47:53Z to 10:49:00Z: post-seal run 1 succeeds, stderr empty.
- 10:49:28Z to 10:50:37Z: post-seal run 2 succeeds, stderr empty.
- 10:50:51Z to 10:51:40Z: both runs are compared, all 13 canonical files per
  run are checksummed and byte-identical, files are made mode 0444, and run
  directories mode 0555.
- 10:51:48Z onward: upstream discovery, hash/header checks, and filtered reads
  occur only after output freeze. The disclosed C1/C3 over-read occurs at
  command 0072, after the scientific bytes are read-only.
- 10:55:25Z: traced upstream reconciliation reads only the three pinned result
  files and the frozen cycle-6 output. There are zero internet connects.
- 11:01:44Z: provenance-only finalizer attempt fails on an over-strict lexical
  ancestor check. Its traceback is preserved. `finalize_cycle6.py` is repaired
  after scientific freeze and rerun successfully at 11:02:32Z. Run1/run2 bytes
  and code/config remain unchanged.

Both post-seal analysis traces have zero forbidden prior-cycle, Task A, Task
030, manuscript/DOCX/PDF, or unpinned-result opens. Their only network-related
calls are local AF_UNIX lookup and bwrap AF_NETLINK namespace setup; internet
connect count is zero. The `bwrap --unshare-net` boundary is effective.

No outcome-bearing or prior numerical result was opened before the seal. No
code, config, mapping, threshold, contrast, target set, taxonomy, or model node
changed after outcome access.

## Deviations D1-D5

- **D1: HARMLESS DISCLOSURE.** Both exit-137 events follow successful synthetic
  computation and are tracer/tee termination failures. No real outcome access
  or scientific fit occurred.
- **D2: HARMLESS DISCLOSURE.** The permission error occurs before the seal and
  before outcome access; the wrapper cannot execute. It has no scientific
  effect.
- **D3: MINOR, fully bounded.** GATA2/SOX9 C1 and C3 upstream rows were printed
  while resolving Task-028 identifiers. This is broader than required. It
  happens only after 13/13 scientific files are hashed and read-only; neither
  the scientific report nor proposed wording can change, and no C1/C3 value is
  incorporated. It cannot affect cycle-6 C2 estimates or taxonomy.
- **D4: HARMLESS DISCLOSURE.** Task 030 contains four C2 targets; pinned
  Task 028 supplies all six. Independent reconciliation still covers all
  required full-base point estimates.
- **D5: MINOR provenance packaging defect, no scientific effect.** The current
  finalizer is provenance-only and post-freeze. Its first traceback and retry
  are present in raw commands 0101-0103, but the pre-fix finalizer bytes were
  not separately preserved. Also, `COMMANDS_CHRONOLOGICAL.tsv` ends at 0102;
  the raw immutable log continues through 0106. Thus the summary table is not a
  complete command index, although the authoritative raw evidence is complete.
  This should be corrected in future packaging, but it cannot reopen the
  already read-only scientific runs.

## Upstream, checksum, and preservation results

- Independent upstream reconciliation: 54/54 checks match, zero mismatch,
  maximum absolute delta `4.440892098500626e-15`. Task-028 primary and
  no-purity source hashes are
  `50df517a55744c12cac1db62a40b123976b1e4dc7efc2b806fe9f6d2a1608f9f`
  and
  `be4a89f723cc1e8f62405ea32f4acb6dccf386cf7cb06f26627a6120478c43d7`;
  Task-030 source hash is
  `f49685bf10cf6f8c8302ddaa81fe4f8c0d60ee36107e7d7f01c89692a52e5399`.
- `OUTPUT_SHA256SUMS.txt`: 69/69 lines independently verified.
- `FILE_INVENTORY.tsv`: 68/68 path, size, and hash rows verified.
- Per-run scientific inventories: 13/13 each verified; 13/13 run pairs are
  byte-identical and read-only.
- Phase-2 frozen input table: 26/26 hashes verified. Cycle-6 sealed allowlist:
  16/16 verified. Phase-1 reviewer, decision, mapping, linkage, missingness, and
  distribution bytes match their frozen hashes.
- External tracked-status fingerprints independently remain exact:
  original dirty clone
  `17fcf6c711c9b376b1fef426440e7583526547ef456f6bada618ef0d2970bc22`,
  Task A
  `84d57032f00ad4846d9e8c8223e9d1b1fd1d6ad9dcbeb85443fff58672ef59c3`,
  and no-purity/Task-030 source worktree
  `16a3e9e9cf8bc9b3ae225c8e25c964939eb1a8d7a8212f844f652f45923e5210`.
- There is no `src/**` diff attributable to this cycle. Existing staged docs
  predate cycle 6; all Phase-1/freeze bytes used here retain their pins. No
  manuscript, DOCX, PDF, figure, staged finding, commit, push, publication, or
  freeze action occurred in cycle 6.

## Statistical, taxonomy, and manufacture-a-positive attack

All 60 C2 estimates are negative. Rank, residual-df, support, condition, VIF,
GVIF, Cramer's V, dummy-correlation, near-zero variance, leverage, Cook,
DFBETA, externally studentized residual, and leave-one-case rules regenerate.
Flags trigger audit only; no case was deleted.

I recomputed every one-patient deletion for all 24 decompositions (10,632 paired
refits). There is no base or adjusted C2 sign change. All GATA2/SOX9 grade
deletions remain below the joint frozen material-attenuation rule. Four of 505
primary-CPE SOX9 histology deletions change the descriptive label from
`amplified` to `largely_retained` because the observed delta is nearly zero.
This demonstrates that `amplified` must not be interpreted as an established
increase; its paired interval spans zero.

GATA2 endometrioid-grade attenuation is close to the frozen materiality
threshold: magnitude attenuation is `0.0893`/`0.0900` and percent attenuation
is `17.08%`/`17.63%`. Relaxing the frozen cutoffs from 0.10 and 20% to 0.08 and
15% would manufacture a materially-attenuated label. Such a threshold change
is forbidden and was not used. Under the pre-frozen rule, both rows remain
`largely_retained`; this is a taxonomy result, not equivalence or absence of
attenuation.

No non-significant p value is used to claim equality, disappearance, or a clean
null. Raw Student-t and diagnostic permutation p values receive no
multiplicity or confirmatory credit.

## Wording assessment and residual MINOR defects

`PROPOSED_WORDING.md` preserves post-hoc, TCGA-only, matched-row, non-causal,
no-purity-independence, no-target-promotion, and no-verdict-change boundaries.
It does **not** fully resolve the prior MINOR wording defect:

1. **MINOR M1:** Results does not state the observed result: histology produces
   only tiny same-direction magnitude increases for GATA2/SOX9, with paired
   intervals spanning zero, while endometrioid-only grade attenuates GATA2 by
   about 17-18% and SOX9 by about 9-10%.
2. **MINOR M2:** Methods/Limitations does not explicitly say that CPTAC-adjusted
   inference, all-histology grade, grade-plus-histology, and grade-stratified
   models are forbidden/not supported.

The charter, frozen spec, machine tables, raw analytical report, and reviewer
boundary do state these limits. Therefore these are presentation defects, not
an undisclosed model or scientific overclaim. The proposed text must not be
called manuscript-ready until they are repaired.

## Adversarial checklist C1-C8

- **C1 circularity: holds.** This is explicitly post-hoc. Model nodes, mappings,
  C2 contrast, 0.05/0.10/20% taxonomy rules, seeds, and diagnostics were frozen
  before cycle-6 outcome loading. No observed result selected a threshold.
- **C2 leakage: holds.** Pre-seal traces contain no outcome/prior-result open.
  Real outcomes are first loaded from the sealed code in traced run 1.
  Upstream results are opened only after read-only output freeze.
- **C3 ecological fallacy: holds within scope.** Outcomes and clinical factors
  are patient-linked, but inference is only a TCGA subtype-level descriptive
  compatibility sensitivity. There is no individual-patient, causal, or
  treatment claim.
- **C4 reproducibility: holds.** Two producer runs are 13/13 byte-identical;
  independent reconstruction and upstream reconciliation regenerate as
  quantified above.
- **C5 source/manuscript drift: holds.** Current C2 point fits match pinned
  Task-028/030 sources; no manuscript or frozen verdict is changed.
- **C6 hidden terms: holds.** Full four-subtype models are fitted. The p53abn
  contrast weight of zero is explicit, not a dropped covariate. Primary CPE and
  no-purity specifications are explicit; no term is silently zeroed.
- **C7 statistics: holds.** Exact matched N, Student-t OLS inference, paired
  bootstrap, diagnostic permutation, design gates, influence diagnostics,
  decomposition, and frozen taxonomy all reproduce. No q value or
  non-significance-as-equivalence claim is present.
- **C8 overreach: holds.** The package is post-hoc, descriptive, non-causal,
  TCGA-only, and no-verdict-change. No CPTAC adjusted model, new category,
  target promotion, purity-independence claim, journal selection, or causal
  claim appears.

## Exact result boundary

The survived statement is:

> In this post-hoc TCGA-UCEC sensitivity, binary histology adjustment does not
> materially attenuate the frozen GATA2 or SOX9 C2 contrast under the
> pre-specified taxonomy; the tiny magnitude increases have intervals spanning
> zero. Within endometrioid carcinoma, binary grade adjustment attenuates GATA2
> magnitude by about 17-18% and SOX9 by about 9-10%, retains negative direction,
> and remains below the joint frozen 0.10 and 20% material-attenuation rule.

This is descriptive compatibility, not proof that histology or grade does or
does not causally explain the association, not statistical equivalence, not
purity independence, not CPTAC adjusted inference, and not an individual
biomarker or treatment result. Frozen TCGA/CPTAC verdicts and target categories
remain unchanged.

reviewer reproducibility manifest:
`VALIDATION_MANIFEST.json`, SHA-256
`2c910e2c987c74674c074aa5f080af43bffa8c93a40ba7e37dd7184546b2f56b`.
`CRITIC_SHA256SUMS.txt` inventories the reviewer report, manifest, code, and
independent outputs. I wrote only inside `phase2_critic_cycle6/` and did not
repair producer artifacts.
