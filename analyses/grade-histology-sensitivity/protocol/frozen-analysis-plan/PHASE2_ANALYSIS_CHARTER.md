# TASK B Phase 2 outcome-blind analysis charter

Freeze UTC: 2026-07-22T18:37:50Z

Classification: HYPOTHESIS / post-hoc explanatory sensitivity

Branch: `experiment/task-b-grade-histology-sensitivity`

Base commit: `83503bad47b60193598b2b9ebe819c22c83e8ac1`

Binding Phase-1 state: `RESTRICT_TCGA_ONLY`

Binding Phase-1 reviewer verdict: `SURVIVES`; `I cannot break it`.

## Claim, prediction, test, and success criterion

Claim: in TCGA-UCEC, the frozen C2 regulon-score contrast for GATA2 and SOX9 can be
re-estimated under the two frozen base specifications and compared on exactly matched
patients before and after the only Phase-1-supported histology or grade additions.

Prediction: if the existing subtype contrast is largely independent of the supported
clinical annotation in this descriptive model, matched-N adjustment will leave its
direction and magnitude largely retained; if the annotation accounts for overlapping
variation, the contrast may attenuate, approach zero, reverse, amplify, or become
unstable. No direction or category is predicted in advance.

Test method: fit frozen OLS models to the frozen equal-weight C2 contrast, with exact
matched-N base refits; calculate coefficient, residual-standardized d, paired-bootstrap
uncertainty, parametric raw p, diagnostic permutation p, and the pre-frozen attenuation
decomposition; enforce the full-design support, rank, collinearity, leverage, residual,
and single-case influence gates before interpretation.

Success criterion: two deterministic runs produce byte-identical scientific outputs,
all input hashes match, all required nodes are reported without repair, and every node
is labelled interpretable or unstable solely by the frozen gates. Success is execution
and honest reporting, not attenuation, retention, a p value, a direction, a target
claim, or a changed upstream verdict.

## Scope

- Cohort: TCGA-UCEC only. There is no CPTAC adjusted analysis and no CPTAC fallback.
- Primary explanatory targets: GATA2 and SOX9.
- Completeness-only targets: HOXA9, WT1, PAX8, and LHX1.
- Outcome: each target's already frozen TASK-028 signed regulon score. No new score,
  target, edge set, category, per-class inference, or outcome transformation is allowed.
- Contrast: C2 only, with weights on `[POLE, MMRd, NSMP, p53abn]` equal to
  `[0.5, 0.5, -1.0, 0.0]`. Positive means the equal-weight mean of POLE and MMRd is
  higher than NSMP. p53abn remains in the full model with contrast weight zero.
- This package creates no confirmatory family, q value, multiplicity credit, target
  category, or verdict. All p values are raw and descriptive.

## Frozen base specifications

The primary CPE base is fit on the frozen CPE-complete roster and contains intercept,
four-level subtype, `M4_prolif`, `purity_CPE`, and `composition`. The no-purity base is
fit on the frozen full roster and contains intercept, subtype, `M4_prolif`, and
`composition`. The two rosters may differ only through the frozen CPE-completeness
indicator. Continuous covariates are used as stored; they are not re-scaled for fitting.
EMMs use fitted-sample means, although C2 cancels the common covariate contribution.

Both base specifications run exactly these nodes:

1. `full_base`: the frozen base reference on its full eligible roster.
2. `histology_matched`: histology-complete rows only, with the base refit on exactly
   those rows and compared with base plus binary histology.
3. `endometrioid_grade_matched`: frozen endometrioid histology and nonmissing grade
   rows only, with the base refit on exactly those rows and compared with base plus
   binary FIGO grade.

Forbidden are all-histology grade, grade plus histology, grade-stratified models,
interactions, three-level histology or grade rescue, alternative category collapse,
patient deletion to improve fit, covariate deletion, imputation, and all CPTAC models.

Clinical bytes, missingness, conflicts, linkage, and categories come only from the
sealed Phase-1 run1 ledger and harmonization specification. Histology uses
endometrioid as reference and non-endometrioid as the treatment-coded level. Within
the endometrioid-only node, grade uses high_grade as reference and low_grade as the
treatment-coded level, matching the Phase-1 declared diagnostics.

The Phase-1 reviewer MINOR M1 is binding: CPTAC diagnostic rows that reused base-design
placeholders are ignored. They are not literal fitted-node diagnostics and must not be
copied into any Phase-2 report.

## OLS and uncertainty

Subtype treatment coding uses POLE as reference, with dummy order MMRd, NSMP, p53abn.
For a design with coefficient order `[intercept, MMRd, NSMP, p53abn, ...]`, C2 is
`0.5 * beta_MMRd - beta_NSMP`. OLS uses `numpy.linalg.lstsq` and requires full rank.
Residual variance is `SSE/(n-p)`. The contrast SE is the square root of
`L sigma2 (X'X)^-1 L'`; the primary raw p is the two-sided Student t p value with
`n-p` df; the coefficient CI is the two-sided 95% t CI. Residual SD is the square root
of residual variance, and `d = beta_C2 / residual_SD`.

Patient-level bootstrap uses 2,000 attempted resamples with replacement. Each attempt
resamples complete patients and all their fields jointly, refits the exact model, and
is invalid if a subtype is absent, rank is deficient, or the fit is non-finite. No
retry replaces an invalid attempt. The 95% d and coefficient CIs are percentile CIs
over valid attempts. Fewer than 1,900 valid attempts makes bootstrap uncertainty
unavailable and the node unstable. Matched comparisons use the same resample indices
for base and adjusted fits, yielding paired percentile CIs for all attenuation and
decomposition quantities.

Subtype-label permutation is a diagnostic only: 2,000 permutations, with outcome,
clinical factor, and covariates fixed, and C2 evaluated by absolute d. The add-one
two-sided diagnostic p is `(1 + count(abs(d_perm) >= abs(d_obs))) / 2001`. No
permutation p replaces the primary t p and neither p creates a category.

Master seed is 20260713. Each stochastic call uses a deterministic unsigned 64-bit
subseed equal to the first eight bytes, big-endian, of SHA-256 over the ASCII string
`20260713|TASKB_PHASE2|<analysis_step_id>`. The step id is the canonical pipe-joined
tuple cohort, base specification, node, target, contrast, and procedure. Run number,
clock time, process id, and output path never enter the seed.

## Full-design gates

Every exact fitted design, not merely the clinical-only Phase-1 design, repeats:

- support: `n >= max(80, 10*p)`, residual df at least 30, all four subtypes;
  additive histology requires at least 5 per subtype and both clinical levels at least
  20 overall and at least 5 in at least two subtypes; endometrioid grade requires at
  least 10 per subtype and both grade levels at least 20 overall and represented by at
  least 5 in at least two subtypes; no level may derive 90% or more from one subtype;
- rank: SVD tolerance `machine_epsilon * max(n,p) * largest_singular_value`; rank below
  p fails;
- condition number on standardized non-intercept columns: at most 30 passes, above 30
  through 100 is unstable, and above 100 is not estimable;
- one-df VIF and adjusted multi-df `GVIF^(1/(2*df))`: at most 5 passes, above 5 through
  10 is unstable, and above 10 fails;
- clinical-by-subtype Cramer's V: above 0.80 warns and above 0.90 fails; complete
  nesting, structural separation, or absolute dummy correlation at least 0.95 fails;
- clinical near-zero variance: frequency ratio above 19 together with percent unique
  at most 10 fails;
- leverage: flag hat above `2*p/n`; any hat at least 0.50 or more than 5% of rows above
  `3*p/n` fails;
- influence flags: Cook's D above `4/(n-p)`, absolute DFBETA above `2/sqrt(n)`, and
  absolute externally studentized residual above 3 are reported and trigger the
  prespecified leave-one-case audit; they do not authorize deletion;
- single-case gate: any leave-one-case fit that changes C2 direction or changes
  absolute d by at least 0.20 makes the node unstable and not interpretable.

No failing node is simplified or repaired. Its estimates and diagnostics remain in
the audit outputs with `UNSTABLE_NOT_INTERPRETED` or `NOT_ESTIMABLE` status.

## Matched-N attenuation and exact decomposition

For base and adjusted fits on identical rows, let `beta_b`, `sigma_b`, and `d_b` be
the matched-N base quantities and `beta_a`, `sigma_a`, and `d_a` the adjusted ones:

- `d_b = beta_b / sigma_b`; `d_a = beta_a / sigma_a`;
- signed delta d = `d_a - d_b`;
- absolute delta d = `abs(d_a - d_b)`;
- magnitude attenuation = `abs(d_b) - abs(d_a)`;
- percent attenuation = `100 * (abs(d_b) - abs(d_a)) / abs(d_b)`;
- beta change = `beta_a - beta_b`;
- residual-scale contribution = `beta_b * (1/sigma_a - 1/sigma_b)`;
- coefficient contribution = `(beta_a - beta_b) / sigma_a`.

The two contributions must sum to signed delta d within absolute tolerance `1e-12`;
otherwise execution stops as an implementation error. If `abs(d_b) < 0.05`, percent
attenuation is NA. If `abs(beta_b) < 1e-8`, coefficient-percent attenuation is NA.
A sign reversal makes percent attenuation NA. Increased magnitude is retained as
negative attenuation and labelled amplification. A comparison is forbidden until the
base is refit on exactly the adjusted-model rows.

## Neutral descriptive taxonomy

Taxonomy is assigned only after all gates, has no statistical or biological claim,
and uses this precedence:

1. `unstable`: any design, uncertainty, or single-case gate fails.
2. `base_near_null_unclassifiable`: `abs(d_b) < 0.05`; report absolute quantities only.
3. `sign_reversal`: nonzero base and adjusted d have opposite signs.
4. `amplified`: same direction and `abs(d_a) > abs(d_b)` beyond `1e-12`.
5. `near_null_or_disappears_descriptively`: same direction or adjusted d is exactly
   zero, `abs(d_b) >= 0.05`, and `abs(d_a) < 0.05`.
6. `same_direction_materially_attenuated`: same direction, `abs(d_a) >= 0.05`,
   magnitude attenuation at least 0.10, and percent attenuation at least 20%.
7. `largely_retained`: every other stable same-direction comparison.

The 0.05 threshold is inherited from the frozen percentage-attenuation edge case. The
0.10 and 20% joint materiality rule is descriptive: 0.10 is one fifth of TASK-028's
pre-outcome half-SD biological-effect floor, while the relative condition prevents a
fixed small change from dominating a small base. These thresholds cannot create or
revise a claim. Non-significance alone is never disappearance, equivalence, or null.

## Outputs and reproducibility

Required outputs include cohort/subtype counts, missingness, clinical distributions,
all six target estimates, coefficient/SE/t CI/raw p/direction/residual SD/d and CI,
exact diagnostics and influence records, GATA2/SOX9 matched-N attenuation and paired
decomposition, a six-target appendix, machine-readable TSV/CSV/JSON, a sensitivity
table suitable for later manuscript drafting, missingness/distribution tables,
scripts, environment, commands, manifest, checksums, and proposed Results, Methods,
and Limitations wording. The wording is a proposal only; no manuscript is edited.

The producer runs twice in fresh deterministic environments. All scientific output
files must be byte-identical across runs. Run-specific timestamps and paths are kept
only in separate provenance manifests and excluded from the scientific comparison.
An independent reviewer must reproduce the analysis and explicitly say whether it
survives. Nothing is documented as established in this task.

The `full_base` C2 values are compared with sealed upstream values only after both new
runs and their checksums are complete. No upstream result file is currently authorized
for researcher or producer discovery. Sophia must provide an exact path and SHA before
a separate reconciliation process; otherwise reconciliation is reported
`BLOCKED_PENDING_EXPLICIT_UPSTREAM_PIN`. A mismatch is reported without changing any
new result, taxonomy, interpretation, target status, category, or verdict.

## Interpretation firewall and literature frame

Clinical adjustment is not a unique causal estimand. Grade and histology can be
confounders, mediators, subtype consequences, pathology proxies, or markers of sample
composition. ESTIMATE addresses stromal/immune admixture but not complete cell-type
composition; CPE is an estimated purity construct. The analysis cannot establish
purity independence, causality, a target/category/verdict, an individual-patient
biomarker, treatment response, journal selection, or a manuscript change.

Supporting context: TCGA's integrated endometrial analysis showed overlap but not
identity between morphology and molecular classes (Nature 2013,
https://doi.org/10.1038/nature12113). ISGYP supports binary FIGO grading for
endometrioid carcinoma and cautions against conflating other histotypes with grade-3
endometrioid disease (Int J Gynecol Pathol 2019,
https://pmc.ncbi.nlm.nih.gov/articles/PMC6295928/). Aran et al. describe CPE as a
consensus purity estimate across methods (Nat Commun 2015,
https://pmc.ncbi.nlm.nih.gov/articles/PMC4671203/). Yoshihara et al. define ESTIMATE's
expression-derived stromal and immune scores (Nat Commun 2013,
https://www.nature.com/articles/ncomms3612). The same sources are counter-evidence to
strong causal interpretation: pathology and molecular subtype overlap, and purity or
composition scores are proxies rather than direct causal measurements.

## Experimenter firewall

The exact permitted and forbidden paths are machine-readable in
`PHASE2_FROZEN_SPEC.json`. The producer may read only the pinned TASK-028 design/code
and arrays, pinned Phase-1 run1 clinical artifacts, this freeze package, and files it
creates under the declared Phase-2 execution root. Recursive discovery is forbidden.
TASK-030, Task A, the no-purity verification worktree, manuscripts, figures, all other
TASK-028 results, and unpinned upstream results are forbidden.

## Outcome-access attestation

During this research freeze I did not load or inspect `scores_v3.npz`; I SHA-hashed it
only. I did not inspect any row of `covariates_v3.tsv` or any content of
`patient_order_v3.json`; the covariate file received one purpose-built header-only,
no-row schema inspection. I did not access any regulon score value, target coefficient,
d, CI, p, q, direction, result table, result-bearing script, manuscript, figure, Task A
artifact, TASK-030 artifact, or no-purity verification artifact. The clean TASK-028
model engine contained no target-specific or expected-outcome literals. No molecular
result or verdict appears in this package.
