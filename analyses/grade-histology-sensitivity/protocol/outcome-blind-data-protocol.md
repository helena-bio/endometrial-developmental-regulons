# TASK B Phase 1: outcome-blind grade and histology feasibility protocol

Freeze timestamp: 2026-07-22T17:25:15Z

Classification: HYPOTHESIS / post-hoc explanatory sensitivity / Phase 1 only

Base commit: 83503bad47b60193598b2b9ebe819c22c83e8ac1

Branch: experiment/task-b-grade-histology-sensitivity

## Decision and scope

This protocol freezes an outcome-blind clinical-annotation audit. Phase 1 may identify, link,
count, and assess the estimability of grade and histology. It must not load, compute, quote, or
infer a regulon/module score, target coefficient, standardized effect, interval, significance
value, direction, result category, figure, manuscript statement, or Task A scientific result.

This is not confirmatory. It cannot modify the frozen TCGA or CPTAC verdicts, cannot create a
new target claim, and cannot earn multiplicity credit. Its only decision is whether later
post-hoc explanatory sensitivity models are mechanically supportable.

Structured researcher hypothesis:

- Claim: a source-pinned, outcome-blind audit can determine whether grade and histology are
  semantically harmonizable and estimable as explanatory covariates in the frozen TCGA and
  CPTAC designs.
- Prediction: clean clinical-only sources, exact patient/tumour linkage, adequate cross-class
  support, and stable design diagnostics will yield a mechanical feasibility state; otherwise
  the audit will stop at a narrower or blocked state without accessing molecular outcomes.
- Test method: deterministic extraction of only approved clinical fields; raw-label and
  linkage ledgers; frozen mapping; missingness and subtype-by-clinical tables; design-matrix
  diagnostics; two byte-identical runs; independent reviewer reproduction.
- Success criterion: the mechanical state rules in `PHASE1_FEASIBILITY_RULES.json` are applied
  exactly. No favorable biological result is part of success.

## Authoritative files

- `PHASE1_FEASIBILITY_RULES.json` contains the machine-readable mapping, support,
  diagnostics, state, Phase-2 hierarchy, and output definitions.
- `ALLOWED_INPUTS_AND_OUTCOME_FIREWALL.json` contains the exact pre-existing file allowlist,
  hashes, field-restricted acquisition specifications, deny classes, process controls, and
  attestations.
- `RESEARCH_ACCESS_LOG.tsv` records every non-reference file touched during protocol design.

On conflict, the JSON rules are binding. Any ambiguity is resolved conservatively: no read,
no inference, and `BLOCKED`.

## Phase-1 objective

For TCGA-UCEC and each frozen CPTAC-UCEC stratum, the experimenter will determine, without
opening molecular outcomes:

1. the exact official grade and histology source, release/version, field name, coding system,
   row identifier, and patient/sample/aliquot key;
2. whether the annotation is intrinsic clinical/pathology metadata rather than an
   expression-inferred or molecularly inferred label;
3. whether the record links to the same frozen analytic patient and same primary tumour;
4. raw category frequencies, duplicate and conflict counts, and missingness overall and by
   cohort, stratum, and frozen subtype;
5. subtype by raw category and subtype by proposed category tables;
6. whether source meanings support the frozen non-outcome-driven mappings; and
7. whether the resulting clinical/subtype design passes the pre-fixed support and stability
   rules.

No clinical distribution was inspected when these rules were written.

## Required audit artifacts

The Phase-1 experimenter must produce, at minimum:

- `INPUTS_LOCK.json`: exact URL/query, study release/version, retrieval UTC, byte size,
  response headers where available, local path, and SHA-256 for every field-restricted input;
- `SOURCE_FIELD_DICTIONARY.tsv`: cohort/stratum, authority, release, source field, definition,
  coding vocabulary, identifier columns, clinical-versus-inferred status, and evidence URL;
- `RAW_CATEGORY_COUNTS.tsv`: raw token counts including explicit missing tokens, never only a
  cleaned table;
- `LINKAGE_LEDGER.tsv`: one row per frozen analytic patient, source row and sample/aliquot,
  duplicate/conflict status, linkage status, and exclusion reason;
- `MISSINGNESS.tsv`: overall and by cohort, stratum, subtype, and proposed clinical category;
- `SUBTYPE_BY_RAW_CATEGORY.tsv` and `SUBTYPE_BY_PROPOSED_CATEGORY.tsv`;
- `HARMONIZATION_MAP.tsv`: every observed raw token, its proposed value, authority for the
  mapping, and reason; original tokens remain immutable;
- `DESIGN_DIAGNOSTICS.tsv`: N, p, rank, residual df, singular values, condition number,
  VIF/adjusted GVIF, Cramer's V, near-zero variance, leverage, and support-rule results for
  every planned node;
- `PHASE1_DECISION.json`: per variable, cohort, stratum, and model node status plus the global
  mechanical state;
- `REPRODUCIBILITY_MANIFEST.json`: code SHA-256, input hashes, commands, environment, seeds,
  timestamps, output hashes, access-log hash, and both-run comparison; and
- `OUTCOME_FIREWALL_ATTESTATION.txt`: all required attestations verbatim.

Raw counts are clinical-only. They are required even when a mapping fails, because a failure
must be auditable and must not be hidden by collapsing categories.

## Clean-source acquisition gate

Existing broad clinical exports and supplements are not assumed safe merely because they
contain grade or histology. The saved broad cBioPortal clinical-value file, broad GDC case
export, and the Dou confirmatory supplement contain fields beyond this audit's scope. They are
not allowlisted for patient-value extraction.

The experimenter must use the exact field-restricted acquisitions declared in
`ALLOWED_INPUTS_AND_OUTCOME_FIREWALL.json`:

- TCGA PanCanAtlas cBioPortal `GRADE` and `ICD_O_3_HISTOLOGY` only;
- CPTAC Discovery cBioPortal `HISTOLOGIC_GRADE_FIGO` and `HISTOLOGIC_TYPE` only; and
- PDC `clinicalPerStudy` for the two pinned CPTAC study UUIDs, requesting only case linkage,
  morphology, primary diagnosis, and tumour grade, then joining the separately pinned
  biospecimen roster.

Each response must contain only declared keys, be written to the exact input path, and be
SHA-locked before any patient value or frequency is printed. A response with an extra survival,
recurrence, treatment-response, expression, molecular-score, or target-result field is mixed
and must not be read. If a clean field-restricted response cannot supply both annotation and
exact linkage, the variable/stratum is `BLOCKED`; a broad supplement is not a fallback.

## Linkage, duplicate, and exclusion rules

The analytic unit is one frozen patient and one frozen primary analytic tumour. TCGA must match
the frozen cohort roster and frozen four-class subtype table by exact participant/case ID, then
verify the analytic sample. CPTAC must retain Discovery and Confirmatory as separate
donor-independent strata and match case, sample, aliquot, specimen type, and tumour role. A
technical reference is never a patient.

Only exact matches are allowed. There is no nearest sample, first available sample, cross-case
substitution, or post-hoc choice among aliquots. Byte-identical duplicates are collapsed and
ledgered. A source-designated final central-pathology record may outrank a preliminary record
only when that status is defined before values are inspected. Same-precedence conflicts become
missing. No majority vote, value preference, or subtype-based repair is allowed.

Exclusions are limited to:

- not in the frozen analytic roster;
- not the frozen primary analytic tumour;
- technical reference or normal material;
- exact linkage unavailable;
- duplicate/conflicting same-precedence annotation;
- source-defined unclassified/unknown/not reported; or
- failure of a frozen support or design rule.

Every exclusion remains in the ledger with one primary reason and all secondary flags. No
patient is excluded because of a molecular result or because exclusion improves a later model.

## Outcome-independent harmonization hierarchy

### Grade

Keep the raw value first. Normalize only punctuation, case, whitespace, and an explicit grade
prefix. Use an official source definition or coding dictionary; never infer grade from subtype,
histology, stage, expression, or molecular phenotype.

The candidate binary grade mapping is FIGO grades 1-2 versus FIGO grade 3. It is valid for
endometrioid carcinoma when the source documents FIGO grading. It is not valid to manufacture
an across-histology grade by assigning serous, clear cell, mixed, undifferentiated carcinoma,
or carcinosarcoma to grade 3. If the cohort source does not document a common grading construct
for all retained histologies, grade adjustment is restricted to the endometrioid subset.
Unknown, not reported, indeterminate, and conflicting grade remain missing. The raw three-level
grade table is descriptive only and cannot replace the binary freeze after outcomes.

### Histology

The primary mapping is endometrioid versus non-endometrioid. The source must explicitly name
the histology or provide a documented ICD-O morphology code. Serous, clear cell, mixed,
undifferentiated carcinoma, carcinosarcoma, and other explicit non-endometrioid types map to
non-endometrioid. Mixed histology is never reassigned to a pure component.

A separate three-level endometrioid/serous/other-non-endometrioid mapping is allowed only if all
pre-fixed support rules pass: each level has at least 20 observations; serous and
other-non-endometrioid each have at least 5 observations in at least two subtypes. Otherwise the
binary mapping is used. Categories are never collapsed to obtain a favorable molecular effect.

### Cross-cohort meaning

Equivalent labels do not establish equivalent definitions. The audit must compare coding
authority, specimen context, pathology-review status, grading system, and release. A compatible
but non-identical definition is a labelled cohort-specific sensitivity. An undocumented or
incompatible definition is `NOT_HARMONIZABLE`. TCGA and CPTAC clinical coefficients are never
pooled merely because their cleaned labels have the same spelling.

## Missing data and model support

Every model uses complete cases for its own required variables. There is no clinical or
molecular imputation, missingness indicator, best-available composite, or value borrowed from
another cohort. In later attenuation comparisons, the base model must be refit on the exact
complete-case rows of the adjusted model.

The fixed support rules are in `PHASE1_FEASIBILITY_RULES.json`. In summary:

- additive adjustment: complete-case N at least `max(80, 10*p)`, at least 30 residual degrees
  of freedom, all four subtype levels, at least 5 per subtype, each retained clinical level at
  least 20 overall and represented by at least 5 observations in at least two subtypes;
- endometrioid-only: complete-case N at least `max(80, 10*p)`, all four subtype levels, and at
  least 10 per subtype;
- optional grade-stratified: at least two grade strata must each meet N at least
  `max(80, 10*p)`, all four subtype levels, and at least 10 per subtype; otherwise it is not run;
- a retained clinical level cannot derive 90% or more of its observations from one subtype; and
- a subtype missing fraction above one half, or missingness that destroys the cell rules,
  fails that variable/cohort model.

These rules are for continuous molecular outcomes and design stability; they are not
event-per-variable rules and are not tuned to observed effects.

## Rank, collinearity, and instability rules

For each planned clinical/subtype design, construct the exact model matrix with a declared
reference level. Full column rank is mandatory using the frozen SVD tolerance. Residual degrees
of freedom must be at least 30. Compute condition number on standardized non-intercept columns:
at most 30 passes, greater than 30 through 100 is unstable, and greater than 100 is not
estimable.

One-degree VIF and multi-degree `GVIF^(1/(2*df))` at most 5 pass. Values greater than 5 through
10 are unstable and sensitivity-only; greater than 10 fail. Cramer's V above 0.80 warns and
above 0.90 fails. Complete nesting, structural separation, any clinical/subtype dummy absolute
correlation at least 0.95, or the frozen near-zero-variance rule fails the parameterization.

Phase 1 computes leverage from the clinical/subtype design only. Hat values above `2*p/n` are
flagged; at least 0.50 for any row, or more than 5% above `3*p/n`, fails. No row is automatically
deleted.

The expression-derived frozen covariate arrays are behind the outcome firewall, so Phase 1
cannot honestly certify the exact full molecular model. Phase 2 must repeat rank, condition,
VIF, leverage, and outcome-dependent influence checks on the exact full design. Cook's D,
DFBETA, studentized-residual, sign-change, and single-case delta-d rules are frozen now but are
not computed in Phase 1.

## Mechanical decision states

Apply the exact definitions and precedence in `PHASE1_FEASIBILITY_RULES.json`:

- `PASS`: clean and linked sources plus mandatory hierarchy nodes pass in TCGA and required
  CPTAC strata; optional grade-stratified failure is recorded but does not prevent PASS.
- `RESTRICT_TCGA_ONLY`: TCGA supports at least one requested adjustment, but CPTAC is blocked or
  not harmonizable for it.
- `PARTIAL_MODEL_HIERARCHY`: at least one frozen node is estimable, but one or more of +both,
  endometrioid-only, or a clinical variable fails. Only supported nodes may proceed.
- `NOT_HARMONIZABLE`: clean annotations exist but their meanings cannot be mapped without
  inference, or every admissible mapping fails support/overlap.
- `BLOCKED`: clean source, hash, release, linkage, or firewall enforcement is absent. Nothing
  downstream may start.

The audit reports per-variable/per-stratum states before the global state. No state is a
scientific finding.

## Frozen Phase-2 hierarchy, if separately authorized

Phase 2 may run only after a fresh reviewer verifies Phase 1. No coefficient is computed now.

For every target and frozen contrast, use this order without adding a rescue model:

1. exact frozen base model;
2. base plus grade, if supported;
3. base plus histology, if supported;
4. base plus grade plus histology, only if both are supported and the exact design passes;
5. endometrioid-only, only if supported; and
6. optional grade-stratified models, only if the separate stringent rule passes.

TCGA must run both the frozen primary Aran-CPE purity model and the otherwise identical frozen
no-purity model. CPTAC must use its frozen primary no-purity specification, with Discovery and
Confirmatory fitted separately and combined only by the frozen meta-analysis structure. A
clinical addition is additive; no grade-by-subtype or histology-by-subtype interaction is
introduced in this task.

Grade and histology can be confounders, mediators, consequences of subtype biology, or proxies
for pathology and sampling. The adjusted model is not a uniquely causal model. A change after
adjustment cannot establish why the change occurred.

## Frozen target scope

The eventual primary explanatory scope is GATA2 and SOX9. HOXA9, WT1, PAX8, and LHX1 are
completeness-only. There are no new targets, categories, claims, multiplicity credit, or frozen
verdict changes. A supported adjustment can only describe sensitivity of an existing estimate.

## Required Phase-2 outputs and attenuation definitions

For every cohort, stratum, model, target, and frozen contrast, report analytic and complete-case
N, missingness, clinical distributions, coefficient, standard error, confidence interval,
p value, direction, residual SD, standardized d and interval, VIF/adjusted GVIF, condition
number, rank, residual degrees of freedom, and leverage/influence summary. Report raw and
harmonized clinical tables alongside molecular output so sample changes are visible.

All comparisons use a matched-N base refit. With `beta` the frozen subtype contrast and `sigma`
the residual SD:

- `d = beta / sigma`;
- signed delta d is `d_adjusted - d_base`;
- absolute delta d is `abs(d_adjusted - d_base)`;
- magnitude attenuation is `abs(d_base) - abs(d_adjusted)`;
- percentage attenuation is
  `100 * (abs(d_base) - abs(d_adjusted)) / abs(d_base)`;
- coefficient change is `beta_adjusted - beta_base`;
- residual-scale contribution to signed delta d is
  `beta_base * (1/sigma_adjusted - 1/sigma_base)`; and
- coefficient contribution is `(beta_adjusted - beta_base) / sigma_adjusted`.

The last two terms exactly sum to signed delta d. If `abs(d_base) < 0.05`, percentage attenuation
is undefined. If `abs(beta_base) < 1e-8`, coefficient percentage attenuation is undefined. A
direction flip is labelled sign reversal and percentage attenuation is not reported. Increased
magnitude is retained as negative attenuation and labelled amplification. Unstable models are
reported but not interpreted.

Attenuation is descriptive. It is compatible with confounding, mediation, overlapping subtype
definitions, altered residual scale, selection by complete cases, or instability. It does not
prove a causal explanation.

## Outcome firewall and refusal rule

Phase-1 scripts must be statically scanned, run in a fresh minimal environment, and monitored by
an OS-level file-open audit. Network is allowed only for the exact clinical acquisition calls;
the analysis run is offline. Every opened file must match the allowlist/hash or be a declared
Phase-1 output. Recursive discovery is forbidden.

The experimenter must attest that no score array, expression matrix, target result, prior
coefficient/effect table, result taxonomy, figure, manuscript, Task A artifact, or TASK-030
result was accessed. The open-file trace, access log, input lock, code hash, output hashes, and
attestation are part of the manifest.

If clean separation cannot be enforced, stop. Do not use a mixed file, do not print a partial
distribution, and do not improvise a parser after seeing values. Return `BLOCKED`.

## Primary and official context

The clinical constructs are not automatically interchangeable:

- The TCGA endometrial carcinoma study integrated endometrioid, serous, and mixed morphology
  and showed that morphology and genomic subtype overlap imperfectly, including molecularly
  serous-like high-grade endometrioid tumours. This supports considering grade/histology but
  argues against treating either as a substitute for subtype. Source: TCGA Research Network,
  Nature 2013, DOI 10.1038/nature12113,
  https://www.nature.com/articles/nature12113.
- The CPTAC Discovery and Confirmatory studies are separately collected proteogenomic cohorts
  and document their own pathology and genomic-subtype implementations. A shared disease name
  does not establish identical pathology fields or review. Sources: Dou et al., Cell 2020,
  DOI 10.1016/j.cell.2020.01.026,
  https://pmc.ncbi.nlm.nih.gov/articles/PMC7233456/; Dou et al., Cancer Cell 2023,
  DOI 10.1016/j.ccell.2023.07.007,
  https://pmc.ncbi.nlm.nih.gov/articles/PMC10631452/.
- ISGYP recommends FIGO grading for endometrioid carcinoma and notes that serous, clear cell,
  undifferentiated carcinoma, and carcinosarcoma are intrinsically high-grade tumour types,
  not FIGO grade-3 endometrioid tumours. This is the reason not to force non-endometrioid types
  into the grade-3 bin. Source: Soslow et al., International Journal of Gynecological Pathology
  2019, https://pmc.ncbi.nlm.nih.gov/articles/PMC6295928/.
- The official FIGO 2023 staging paper distinguishes low-grade endometrioid (grades 1-2),
  high-grade endometrioid (grade 3), and aggressive non-endometrioid histological types. Source:
  Berek et al., International Journal of Gynecology and Obstetrics 2023,
  https://pmc.ncbi.nlm.nih.gov/articles/PMC10482588/.
- Expert re-review of TCGA high-grade endometrioid carcinoma found nontrivial diagnostic
  variability, especially near the serous/high-grade endometrioid boundary. This is strong
  counter-evidence to naive cross-cohort harmonization. Source: McConechy et al., Modern
  Pathology 2016, https://pmc.ncbi.nlm.nih.gov/articles/PMC4934379/.
- GDC defines `primary_diagnosis`, `morphology`, and `tumor_grade` as distinct diagnosis
  properties and updates permissible values across data-dictionary releases. Exact release and
  field provenance therefore matter. Sources: NCI GDC Data Dictionary,
  https://gdc.cancer.gov/about-data/data-dictionary; release notes,
  https://docs.gdc.cancer.gov/Data_Dictionary/Release_Notes/Data_Dictionary_Release_Notes/.
- PDC exposes versioned clinical data and a field-restricted `clinicalPerStudy` query including
  case linkage, morphology, primary diagnosis, and tumour grade. Source: PDC API documentation,
  https://pdc.cancer.gov/pdc/publicapi-documentation and PDC usage documentation,
  https://pdc.cancer.gov/pdc-docs/usage.

Strongest case against the hypothesis: grade is not one universal construct across endometrial
histologies; pathology labels can be release-, specimen-, and reviewer-dependent; molecular
subtypes and morphology are associated; small or empty subtype-by-clinical cells may make
adjustment collinear; and confirmatory CPTAC grade may not be available through a clean,
linkable field-restricted source. Any of these can make the honest result `BLOCKED`,
`NOT_HARMONIZABLE`, or a partial/TCGA-only hierarchy.

No source predicts that adjustment will attenuate, preserve, or strengthen a molecular effect.

## Exact experimenter handoff

The Phase-1 experimenter owns only `experiments/taskB_grade_histology/phase1_execution/` and
may create `inputs/` under the task directory. The experimenter must:

1. verify branch and HEAD, verify this protocol and both JSON hashes, and stop on mismatch;
2. run the exact field-restricted acquisition in a clean process, lock inputs before values,
   and refuse broad/mixed fallbacks;
3. verify all pre-existing allowlist hashes before use;
4. write one deterministic ASCII-only audit script using only clinical fields, linkage, roster,
   and subtype labels;
5. emit every required clinical-only artifact, including raw categories and conflict ledgers;
6. apply mappings and decision states mechanically, never by narrative judgment after counts;
7. run from a fresh output directory twice with identical inputs and require byte-identical
   canonical outputs; timestamps belong only in a manifest field excluded from canonical
   comparison;
8. record commands, package versions, locale, timezone, Python hash seed, input/output hashes,
   open-file trace, and access log; and
9. stop after Phase 1. Do not load any frozen molecular score or run a molecular model.

Phase 2 cannot begin until a fresh reviewer independently verifies the firewall, official source
and release, exact linkage, clinical counts and missingness, mappings, diagnostics, two-run
reproduction, and mechanical state, and explicitly states that no outcome was accessed.

## Freeze attestation

At this timestamp, no grade/histology distribution, subtype-by-grade/histology table, or
molecular outcome was inspected to choose these rules. The protocol is a design package, not a
clinical audit, scientific finding, or verdict.
