# TASK A: post-hoc per-class decomposition of the frozen C2 contrast

Status: FROZEN ANALYSIS CHARTER BEFORE ANY TASK-A OUTCOME COMPUTATION

Classification: HYPOTHESIS / POST-HOC EXPLANATORY SENSITIVITY

This charter is explanatory, not confirmatory. It does not replace, amend, or
re-credit the frozen TASK-028 C2 contrast, the frozen TASK-029 CPTAC verdicts,
or the TASK-030 verification. Nothing in TASK A may change a target, category,
threshold, regulon, edge, score, covariate, patient rule, multiplicity family,
replication verdict, or manuscript byte.

## 1. Claim, scope, and success

### Claim

For each of the six frozen regulon targets, the equal-weight frozen C2 contrast

    C2 = 0.5 * (POLE - NSMP) + 0.5 * (MMRd - NSMP)

can be decomposed, without refitting a different biological model, into POLE
versus NSMP, MMRd versus NSMP, and POLE versus MMRd components. The decomposition
may show compatible same-direction class effects, a same-direction magnitude
difference, one class dominating the common direction, opposite-direction
heterogeneity, or unresolved estimates. Any of those outcomes is informative.

The primary explanatory targets are GATA2 and SOX9. HOXA9, WT1, PAX8, and LHX1
are completeness targets and must be analyzed and reported by the identical
rules. PAX8 and LHX1 were frozen C1 targets; their TASK-A C2 decompositions are
new descriptive sensitivities only and cannot alter their C1 status.

### Null and counter-outcomes

The charter explicitly permits all of the following:

1. Both POLE-versus-NSMP and MMRd-versus-NSMP estimates point in the same
   direction and their difference is compatible with zero.
2. Both point in the same direction but their magnitudes differ directly under
   the POLE-minus-MMRd contrast.
3. Only one class estimate is resolved while the other remains imprecise. This
   is not proof that the unresolved class has no effect.
4. POLE and MMRd point in opposite directions, so their equal-weight average
   masks class heterogeneity.
5. All estimates are imprecise or near zero. This is an honest unresolved/null-
   compatible explanatory result, not a failure of execution.
6. Patterns differ across TCGA primary, TCGA no-purity, CPTAC Discovery, and
   CPTAC Confirmatory. Such differences may be biological, compositional, or
   technical; TASK A cannot identify a causal source.

### Exact success criterion

TASK A succeeds as an explanatory sensitivity if and only if all of the
following hold, regardless of the signs, magnitudes, p values, or q values:

1. The frozen upstream bytes and data versions in Section 12 verify exactly;
   otherwise execution refuses to run.
2. One unchanged frozen design is fit per target within each of TCGA primary,
   TCGA no-purity, CPTAC Discovery, and CPTAC Confirmatory, and all three new
   contrasts are derived from that one fit.
3. Every expected target-by-contrast row is present, or is marked NOT_EVALUABLE
   with a mechanical reason fixed independently of the outcome.
4. The four within-cohort/stratum unstandardized and standardized C2 identities
   pass the numerical tolerance in Section 8 for every target.
5. Each of the five descriptive 18-test families in Section 6 contains exactly
   18 tests and uses no outcome-dependent additions, deletions, or regrouping.
6. All required tables, forests, manifests, commands, environment records,
   checksums, and two deterministic run comparisons are delivered.
7. The report states that TASK A is post-hoc, descriptive, non-causal, and
   unable to change any frozen verdict or manuscript byte.

No positive result is required. Explanatory resolution, including a clear null,
heterogeneity, or honest non-evaluability, is the success condition.

## 2. Frozen targets, cohorts, and models

Target order is fixed:

    GATA2, SOX9, HOXA9, WT1, PAX8, LHX1

Model/cohort order is fixed:

1. `TCGA_PRIMARY_CPE_N506`: frozen TCGA primary purity-adjusted model; score as
   outcome; four-level subtype factor; covariates M4 proliferation, Aran CPE
   purity, and frozen ESTIMATE composition; CPE complete cases only.
2. `TCGA_NOPURITY_N507`: frozen TCGA no-purity sensitivity; score as outcome;
   same four-level subtype factor; covariates M4 and frozen ESTIMATE composition;
   the frozen full cohort.
3. `CPTAC_DISCOVERY_N95`: frozen CPTAC V2 no-purity model, fit in Discovery only;
   score as outcome; same four labels; covariates M4 and frozen ESTIMATE
   composition.
4. `CPTAC_CONFIRMATORY_N135`: the same frozen CPTAC V2 model, fit in Confirmatory
   only.
5. `CPTAC_FIXED_EFFECT_META`: inverse-variance fixed-effect synthesis of
   standardized effects from the two separately fit CPTAC strata. Raw expression
   is never pooled.

The exact analytic sample and subtype counts must be regenerated and printed.
Expected upstream counts, which are guards rather than values to be optimized,
are:

| model | total | POLE | MMRd | NSMP | p53abn |
| --- | ---: | ---: | ---: | ---: | ---: |
| TCGA primary | 506 | 49 | 148 | 146 | 163 |
| TCGA no-purity | 507 | 49 | 148 | 147 | 163 |
| CPTAC Discovery | 95 | 7 | 25 | 43 | 20 |
| CPTAC Confirmatory | 135 | 6 | 47 | 66 | 16 |

A mismatch is a data-version or cohort-definition failure. It is not repaired
by dropping, adding, relabeling, imputing, or substituting patients.

## 3. Exact estimands and sign conventions

Let the adjusted subtype marginal means from the single fitted target model be

    mu = (mu_POLE, mu_MMRd, mu_NSMP, mu_p53abn).

At the frozen model's reference covariate values, derive all three TASK-A
contrasts simultaneously from this same `mu` and coefficient covariance matrix:

| name | weights on [POLE, MMRd, NSMP, p53abn] | estimand |
| --- | --- | --- |
| `POLE_vs_NSMP` | `[+1, 0, -1, 0]` | `b_PN = mu_POLE - mu_NSMP` |
| `MMRd_vs_NSMP` | `[0, +1, -1, 0]` | `b_MN = mu_MMRd - mu_NSMP` |
| `POLE_vs_MMRd` | `[+1, -1, 0, 0]` | `b_PM = mu_POLE - mu_MMRd` |

The frozen C2 contrast remains:

    b_C2 = 0.5 * mu_POLE + 0.5 * mu_MMRd - mu_NSMP
         = 0.5 * b_PN + 0.5 * b_MN.

The POLE:MMRd sample proportions never enter these weights. The 0.5/0.5 weights
are equal biological subtype weights inherited from the frozen design.

Sign convention is always first-named class minus second-named class:

- Positive `b_PN` means higher regulon activity in POLE than NSMP.
- Positive `b_MN` means higher regulon activity in MMRd than NSMP.
- Positive `b_PM` means higher regulon activity in POLE than MMRd.
- Negative values mean the reverse. No target-specific sign flipping is allowed.

For each observed fit, let `s` be that fit's OLS residual SD. Report

    d_PN = b_PN / s
    d_MN = b_MN / s
    d_PM = b_PM / s
    d_C2 = b_C2 / s.

Because all contrasts within a cohort/stratum use the same fit and the same
residual SD,

    d_C2 = 0.5 * d_PN + 0.5 * d_MN

must hold numerically within that fit. A separately refit POLE/NSMP or MMRd/NSMP
subset is forbidden: it would change the residual scale, covariate reference,
and estimand and would destroy the exact identity.

## 4. Estimation and uncertainty

For every target in every directly fit cohort/stratum:

1. Fit the frozen OLS design exactly once.
2. Derive the three coefficient contrasts and their covariance-based SEs from
   that fit.
3. Report the two-sided 95% model-based t interval for each unstandardized
   coefficient using the fit's residual degrees of freedom.
4. Report `d = coefficient / residual_SD`.
5. Use a patient/case bootstrap with replacement, 2,000 replicates, refitting
   the same full model and deriving all three contrasts in every replicate.
   Report percentile 95% CIs for coefficient and d, bootstrap SE(d), bootstrap
   attempts, usable replicates, and mechanical failures. The analytic unit is
   the patient/case, never an aliquot.
6. Use the frozen-style two-sided subtype-label permutation test with 2,000
   permutations: permute the four-class label over the fitted cases, keep score
   and covariates fixed, refit the same model, and count
   `abs(d_perm) >= abs(d_observed) - 1e-12`. Use the add-one numerator and
   denominator. This permutation p is the single raw p entering BH-18.
7. Model-based t p values may be emitted as diagnostics but must be explicitly
   named `wald_t_p_diagnostic`, never substituted selectively for the frozen
   `perm_p_raw` family input.

For CPTAC fixed-effect meta-analysis, synthesize each target/contrast separately
on the standardized d scale:

    weight_h = 1 / bootstrap_SE(d_h)^2
    d_FE = sum_h(weight_h * d_h) / sum_h(weight_h)
    SE_FE = sqrt(1 / sum_h(weight_h)), h in {Discovery, Confirmatory}.

Report the two-sided normal 95% CI and raw p, both stratum effects and CIs, both
weights, the normalized weight fraction, and an opposite-direction flag. If a
bootstrap SE is zero/non-finite or a stratum is not evaluable, the corresponding
meta estimate is NOT_EVALUABLE; no fallback SE or one-stratum meta is allowed.

## 5. Required result fields

The long-form result table must contain, at minimum:

- task, target, target_role, model, cohort, stratum, contrast;
- model_formula, covariates, analytic_n, residual_df;
- n_POLE, n_MMRd, n_NSMP, n_p53abn;
- contrast_weights, coefficient, coefficient_SE, residual_SD, d;
- coefficient_t_CI_lo, coefficient_t_CI_hi;
- coefficient_boot_CI_lo, coefficient_boot_CI_hi;
- d_boot_CI_lo, d_boot_CI_hi, d_boot_SE;
- perm_p_raw, BH18_q_descriptive, point_direction;
- n_boot_requested, n_boot_used, n_boot_failed;
- n_perm, deterministic_seed;
- evaluability_status, evaluability_reason;
- c2_reconstruction_value, c2_reconstruction_abs_error,
  d_c2_reconstruction_abs_error;
- interpretation_state and interpretation_basis.

Every row must retain full floating-point precision in TSV/JSON. Rounded values
may appear only in human-readable Markdown and figures.

For `CPTAC_FIXED_EFFECT_META`, unstandardized coefficient and residual SD fields
are NA by design because the two strata do not share an outcome scale. The meta
row must instead contain `d_FE`, `SE_FE`, its CI/p/q, and both stratum-specific
coefficients, SEs, residual SDs, d values, and weights. It is forbidden to invent
a pooled residual SD or back-calculate a pooled unstandardized coefficient.

The manuscript-ready table must show, for each target and model, the POLE-NSMP,
MMRd-NSMP, and POLE-MMRd coefficient, SE, residual SD, d, 95% d CI, raw p,
descriptive BH-18 q, subtype counts, and the pre-specified interpretation state.
Its caption must begin with `Post-hoc explanatory sensitivity` and state that it
cannot change the frozen categories or replication verdicts.

GATA2 and SOX9 each require a forest plot with side-by-side POLE-NSMP and
MMRd-NSMP standardized effects and 95% bootstrap CIs for all four directly fit
cohort/stratum models, plus separate POLE-MMRd rows. CPTAC FE estimates may be
shown in a visually distinct meta row. Plot order, x-axis sign, and labels must
be identical between the two genes. No significance stars.

## 6. Fixed descriptive multiplicity

There are exactly five disjoint, descriptive 18-test families:

1. `BH18_TCGA_PRIMARY`: 6 targets x 3 contrasts.
2. `BH18_TCGA_NOPURITY`: 6 targets x 3 contrasts.
3. `BH18_CPTAC_DISCOVERY`: 6 targets x 3 contrasts.
4. `BH18_CPTAC_CONFIRMATORY`: 6 targets x 3 contrasts.
5. `BH18_CPTAC_FIXED_EFFECT_META`: 6 targets x 3 contrasts.

Within each family, Benjamini-Hochberg correction is applied to all 18 two-sided
raw p values at once. Direct-fit families use `perm_p_raw`; the CPTAC meta family
uses the two-sided fixed-effect normal p. No family is gated by the frozen F1
omnibus, the prior target status, a nominal p, a direction, or evaluability in a
different cohort. If a row is mechanically not evaluable, retain its row with
NA p/q and report both the planned family size (18) and the evaluable count. For
the finite rows, compute BH after placing a value of 1.0 in each mechanically
missing p-value slot, so the denominator remains exactly 18; restore the missing
rows to NA after the calculation. Do not shrink the scientific family or relabel
the remaining tests as a smaller pre-specified family.

The q values are descriptive only. They confer no TASK-028 category, B2.12
credit, TASK-029 replication, or manuscript claim. The frozen TCGA F1/F2 families
and frozen CPTAC BH-of-6 family are not recomputed or replaced.

## 7. Interpretation taxonomy fixed before results

Taxonomy is assigned separately for each target and directly fit model using
the signs and model-based 95% t intervals of `b_PN`, `b_MN`, and the direct
same-model `b_PM`. The standardized `d_PM` is used only for the inherited
materiality floor. Bootstrap intervals remain mandatory uncertainty outputs but
do not replace the fixed taxonomy interval after results are seen.

The pre-existing frozen biological effect floor `abs(d) >= 0.50` is reused for
`d_PM`; no new threshold is invented. This is a descriptive materiality marker,
not confirmatory credit. No equality or equivalence threshold is used. In
particular, a POLE-minus-MMRd CI that includes zero does not prove equality.

Assign exactly one state in this order:

1. `NOT_EVALUABLE`: any required class/design term is not estimable, residual SD
   is non-positive/non-finite, or required uncertainty cannot be computed.
2. `ZERO_BOUNDARY_UNRESOLVED`: either `b_PN` or `b_MN` is exactly zero at stored
   precision. Do not force a same/opposite sign label.
3. `OPPOSITE_DIRECTION_HETEROGENEITY_SUPPORTED`: `b_PN * b_MN < 0`, the 95%
   coefficient CI for `b_PM` excludes zero, and both class-versus-NSMP 95%
   coefficient CIs exclude zero. This is descriptive heterogeneity, not causal.
4. `OPPOSITE_DIRECTION_POINT_HETEROGENEITY_UNRESOLVED`: `b_PN * b_MN < 0` but
   the full support condition in state 3 is not met.
5. `SAME_DIRECTION_MATERIALLY_DIFFERENT_POLE_DOMINANT`: `b_PN * b_MN > 0`, the
   95% coefficient CI for `b_PM` excludes zero, `abs(d_PM) >= 0.50`, and POLE is
   farther from NSMP in the common direction: `abs(b_PN) > abs(b_MN)`.
6. `SAME_DIRECTION_MATERIALLY_DIFFERENT_MMRD_DOMINANT`: same as state 5, but
   `abs(b_MN) > abs(b_PN)`.
7. `SAME_DIRECTION_DISTINGUISHABLE_BELOW_MATERIALITY_FLOOR`: `b_PN * b_MN > 0`,
   the 95% coefficient CI for `b_PM` excludes zero, but `abs(d_PM) < 0.50`.
8. `SAME_DIRECTION_COMPATIBLE_BOTH_RESOLVED`: `b_PN * b_MN > 0`, the 95%
   coefficient CI for `b_PM` includes zero, and both class-versus-NSMP 95%
   coefficient CIs exclude zero in their common direction.
9. `SAME_DIRECTION_POINT_ESTIMATES_UNRESOLVED`: `b_PN * b_MN > 0` and none of
   states 5-8 applies.

For states 5 and 6, if stored magnitudes are exactly tied, use
`SAME_DIRECTION_MATERIALLY_DIFFERENT_TIE_ERROR`; this logically inconsistent
combination must be investigated, not broken by an arbitrary tie rule.

The words `compatible`, `concordant`, and `unresolved` are mandatory when the
direct POLE-minus-MMRd interval includes zero. The words `equal`, `same effect`,
`no difference`, and `equivalent` are forbidden unless a separately powered,
pre-specified equivalence analysis exists; none exists in TASK A.

For `CPTAC_FIXED_EFFECT_META`, apply the same ordered taxonomy with `d_PN`,
`d_MN`, and `d_PM` and their fixed-effect normal 95% CIs substituted for the
unstandardized coefficients and t intervals. The direct FE POLE-minus-MMRd
estimate remains the basis for a meta-level class difference. Label the result
`CPTAC_META_<state>` to distinguish it from either stratum. Because the three
meta contrasts have contrast-specific weights, this label summarizes
standardized evidence and does not imply the within-stratum C2 identity at the
meta level.

### Across-cohort descriptive synthesis

Do not force one global biological label. Report the four direct-fit states and
the CPTAC meta state side by side. An across-cohort phrase may be used only as:

- `cross-cohort pattern compatible`: all evaluable direct-fit class effects have
  the same sign pattern, without claiming equality;
- `cross-cohort magnitude variation`: sign pattern is retained but at least one
  direct PM interval supports a magnitude difference;
- `cross-cohort sign heterogeneity`: the class sign pattern differs between
  cohorts/strata;
- `cross-cohort unresolved`: remaining cases.

## 8. Reconstruction and algebra checks

For each target and each of the four directly fitted models, calculate the
frozen C2 estimate directly from `[0.5, 0.5, -1, 0]` and reconstruct it from the
two class-versus-NSMP estimates.

Define

    scale_b = max(1, abs(b_C2), 0.5*abs(b_PN) + 0.5*abs(b_MN))
    err_b = abs(b_C2 - (0.5*b_PN + 0.5*b_MN))
    scale_d = max(1, abs(d_C2), 0.5*abs(d_PN) + 0.5*abs(d_MN))
    err_d = abs(d_C2 - (0.5*d_PN + 0.5*d_MN)).

Pass requires both:

    err_b <= 1e-12 * scale_b
    err_d <= 1e-12 * scale_d.

Any failure is a blocking code/design error. Do not loosen the tolerance after
seeing results.

Bootstrap percentile endpoints need not obey the identity because quantiles of
linear combinations are not combinations of marginal quantiles. The identity
must nevertheless hold replicate by replicate before marginal quantiles are
taken; report the maximum replicate-level scaled error under the same tolerance.

CPTAC FE estimates are not required to obey the C2 identity. Each stratum uses
its own residual SD, and each contrast receives its own inverse-variance weights.
Thus direct FE(C2) generally differs from
`0.5*FE(POLE-NSMP) + 0.5*FE(MMRd-NSMP)`. Report the difference and the three sets
of weights as a meta-analysis consequence, not a failed reconstruction. Never
alter weights to force equality.

## 9. GATA2/SOX9 heterogeneity assessment

For GATA2 and SOX9, the report must lead with a side-by-side table of:

- POLE-NSMP coefficient, SE, 95% coefficient CI, residual SD, d, 95% d CI;
- MMRd-NSMP coefficient, SE, 95% coefficient CI, residual SD, d, 95% d CI;
- direct POLE-MMRd coefficient, SE, 95% coefficient CI, d, 95% d CI, raw p,
  and descriptive BH-18 q;
- subtype counts and evaluability;
- the pre-specified state from Section 7.

The POLE-minus-MMRd coefficient and CI from the same linear model are the only
formal basis for saying one class differs from the other. Comparing one nominal
p value below 0.05 with another above 0.05 is forbidden. Sample mix by itself is
not an artifact and cannot explain the equal-weight C2 mechanically.

## 10. Fixed firewall

The following are prohibited:

1. Calling C2 a sample-size-weighted contrast or calling POLE:MMRd sample mix
   alone an artifact.
2. Refitting class-pair subsets, reweighting C2, or using cohort proportions as
   contrast weights.
3. Adding, dropping, relabeling, or regrouping a target, subtype, patient,
   covariate, regulon edge, score component, test, or family in response to an
   observed result.
4. Retuning the 0.50 effect floor, any category threshold, a covariate, or a
   scoring method.
5. Treating TASK-A q values as confirmatory or using them to change TASK-028,
   TASK-029, or TASK-030 credit/verdict fields.
6. Choosing between purity/no-purity, Discovery/Confirmatory, or class-specific
   results based on which looks stronger.
7. Treating a non-significant class estimate or POLE-minus-MMRd contrast as proof
   of absence, equality, or equivalence.
8. Making a causal, mechanistic, treatment-predictive, or individual-patient
   biomarker claim from this subtype-level bulk-expression analysis.
9. Editing any manuscript, DOCX, PDF, upstream artifact, `src/**`, frozen seal,
   or verified TASK-030 file.

The permitted claim is limited to descriptive subtype-level decomposition of a
bulk regulon-activity contrast under frozen, model-adapted designs.

## 11. Primary-source literature context

The literature motivates testing heterogeneity; it does not predict the sign of
any of the six target-specific contrasts.

### Evidence supporting biological heterogeneity

1. The TCGA primary study defined POLE-ultramutated and MSI-hypermutated
   endometrial cancers as distinct genomic classes arising from different
   mutational processes, while both differ from copy-number-low/NSMP tumors.
   This supports decomposing a pooled POLE+MMRd contrast rather than assuming
   interchangeability. Source: Cancer Genome Atlas Research Network,
   `Integrated genomic characterization of endometrial carcinoma`, Nature 2013,
   DOI 10.1038/nature12113, https://doi.org/10.1038/nature12113.
2. Howitt et al. observed a higher predicted neoantigen burden in POLE tumors
   than MSI tumors, while both classes had more immune infiltration than
   microsatellite-stable tumors. This directly supports the possibility of a
   shared direction with different magnitude. Source: JAMA Oncology 2015,
   DOI 10.1001/jamaoncol.2015.2151,
   https://doi.org/10.1001/jamaoncol.2015.2151.
3. van Gool et al. experimentally and computationally linked POLE proofreading
   mutations to an antitumor immune response in endometrial cancer. Source:
   Clinical Cancer Research 2015, DOI 10.1158/1078-0432.CCR-15-0057,
   https://doi.org/10.1158/1078-0432.CCR-15-0057.
4. Stelloo et al. found marked clinicopathological and outcome differences
   among molecular subgroups in independent PORTEC trial cohorts, supporting
   the view that subtype labels can capture more than a single pooled axis.
   Source: Clinical Cancer Research 2016,
   DOI 10.1158/1078-0432.CCR-15-2878,
   https://doi.org/10.1158/1078-0432.CCR-15-2878.

### Counter-evidence and reasons not to over-interpret heterogeneity

1. Howitt et al. also found that POLE and MSI tumors shared increased T-cell
   infiltration and checkpoint expression relative to microsatellite-stable
   tumors. Eggink et al. independently profiled high-risk tumors and likewise
   identified both POLE-mutant and MSI classes as immune-active. Convergent
   downstream biology therefore makes same-direction compatible effects
   plausible; different mutational mechanisms do not imply a difference for
   every regulon. Sources: Howitt et al., DOI above; Eggink et al., OncoImmunology
   2017, DOI 10.1080/2162402X.2016.1264565,
   https://doi.org/10.1080/2162402X.2016.1264565.
2. The PORTEC-3 primary molecular analysis showed strong prognostic separation
   between POLEmut and MMRd groups, but prognosis and treatment response are not
   the regulon-activity estimands in TASK A. It is counter-evidence against
   collapsing the classes clinically, not proof of target-specific expression
   heterogeneity. Source: Leon-Castillo et al., Journal of Clinical Oncology
   2020, DOI 10.1200/JCO.20.00549,
   https://doi.org/10.1200/JCO.20.00549.
3. Dou et al. characterized independent CPTAC endometrial carcinoma cohorts
   using multi-omic platforms and distinct subtype implementations. These
   studies support external biological comparison, but not identical-covariate
   or identical-classifier replication across TCGA and CPTAC. Sources: Cell
   2020, DOI 10.1016/j.cell.2020.01.026,
   https://doi.org/10.1016/j.cell.2020.01.026; Cancer Cell 2023,
   DOI 10.1016/j.ccell.2023.07.007,
   https://doi.org/10.1016/j.ccell.2023.07.007.
4. Aran et al. demonstrated that tumor purity systematically influences bulk
   molecular measurements, and Yoshihara et al. showed that ESTIMATE derives
   stromal/immune admixture from expression itself. Residual cellular-composition
   differences can therefore change apparent subtype effects, especially when
   purity estimators differ or are unavailable. Sources: Nature Communications
   2015, DOI 10.1038/ncomms9971,
   https://doi.org/10.1038/ncomms9971; Nature Communications 2013,
   DOI 10.1038/ncomms3612, https://doi.org/10.1038/ncomms3612.
5. The SEQC consortium's multisite, cross-platform primary experiment showed
   that RNA-seq precision and reproducibility depend on platform, site, and data
   processing choices. GDC harmonization and within-stratum scoring reduce but
   do not erase cohort acquisition, library, classifier, and sample-composition
   differences. Source: Nature Biotechnology 2014,
   DOI 10.1038/nbt.2957, https://doi.org/10.1038/nbt.2957.

### Explicit inference boundary

It is a sourced fact that POLE and MMRd arise from distinct genomic mechanisms,
that both can share immune-active phenotypes, and that purity/platform/cohort
features can influence bulk expression. It is an inference, not a sourced fact,
that any of these mechanisms causes a TASK-A GATA2, SOX9, HOXA9, WT1, PAX8, or
LHX1 class-specific contrast. TASK A can describe compatibility and
heterogeneity; it cannot assign cause.

## 12. Exact experimenter handoff and byte guards

The experimenter must verify the following paths and SHA-256 values before any
TASK-A outcome computation. One mismatch causes `BLOCKED_DATA_VERSION_MISMATCH`;
there is no silent substitution, redownload, regeneration, or nearest-version
fallback.

### Verified TASK-030 authority

| absolute path | required SHA-256 |
| --- | --- |
| `data/external/original-workspace/revgate-tcga-no-purity-verify/experiments/task030_verify_current_main/REPRODUCIBILITY_MANIFEST.json` | `f3f7d520e3479dd1675e6f220581c0a44a386b7357e9022da7c410d9fbd81c4a` |
| `data/external/original-workspace/revgate-tcga-no-purity-verify/experiments/task030_verify_current_main/SHA256SUMS.txt` | `7a3146f8f631b7d1cbc70601edd7715e19aef10b105dd98f8510041fa19c71f9` |
| `data/external/original-workspace/revgate-tcga-no-purity-verify/experiments/task030_verify_current_main/results/SIX_TARGET_RESULTS.tsv` | `e7445a4df6a48ca05a2a73460687a6fd8399e08b3e7c66843702f7385b465e13` |
| `data/external/original-workspace/revgate-tcga-no-purity-verify/experiments/task030_verify_current_main/results/COVARIATE_COMPATIBILITY.tsv` | `25ab052b2ee1ba79bee7307ffb2f4cdd56bf1170f673e22e607320128c524072` |
| `data/external/original-workspace/revgate-tcga-no-purity-verify/experiments/task030_verify_current_main/critic_cycle2_targeted/VALIDATION_REPORT_TARGETED_CYCLE2.md` | `e3a25f3607805e0e8d8335498807656b5688c1234b898c8076ede6a59f3bfe4d` |

The final targeted reviewer verdict must contain the first-line verdict
`SURVIVES`. Its conclusion and all TASK-030 scientific bytes remain read-only.

### TCGA frozen inputs and model authority

| absolute path | required SHA-256 |
| --- | --- |
| `data/external/original-workspace/task028-freeze-b-draft/sealed_v3/B2_MODEL_SPEC.md` | `682099a03923778d170ff4bd46e8b475405419fe1fadb3259ce4c2e75b43f2c1` |
| `data/external/original-workspace/task028-freeze-b-draft/sealed_v3/SEAL_MANIFEST.sha256` | `0d1bfa4f808c1b9b4c63c8b7a3774c47ed43c537de5d1c2f35d24b9389e3020e` |
| `data/external/original-workspace/task028-freeze-b-draft/execution/intermediate/scores_v3.npz` | `956b7a968f35d1da57719e3bc7392decbdd6e292ac31d4dec100b6632518b0cb` |
| `data/external/original-workspace/task028-freeze-b-draft/execution/intermediate/covariates_v3.tsv` | `cc46495c0abc6877182c9bff8d4ba44f11599d8bc5bc7a8cb6f8331aa1ec3c9a` |
| `data/external/original-workspace/task028-freeze-b-draft/execution/scripts_v3/model_engine.py` | `4d7bfb997c9dc188daf9e584a791ff490dfcfdd7fa365c20bca8a2c066f42065` |
| `data/external/original-workspace/task028-freeze-b-draft/execution/results_v3/phase2b_models.json` | `376835bfcc90821ce1d0a40e3614b83300fa1667f7d02210ebe56bfe34958c60` |
| `data/external/original-workspace/task028-freeze-b-draft/execution/results_v3/RESULT_SEAL.sha256` | `49406f449d876185b247040e886e03c52db7c5c7ab67117ef039198a67cade5e` |

The TCGA result seal is standard `sha256sum -c` format and must pass from
`execution/results_v3/`. The sealed-v3 design manifest uses a historical custom
layout: a row `B1/<digest>  <name>` means hash `sealed_v3/B1/<name>` against
`<digest>`; ordinary `<digest>  <relative-path>` rows are checked normally.
Verify every non-comment entry with that exact convention and require zero
missing/mismatched entries. Do not rewrite the seal into a new format. The
experimenter may import or wrap the sealed model code but must not modify it.
New TASK-A code lives only under `experiments/taskA_perclass_c2/`.

### CPTAC frozen inputs and authority

| absolute path | required SHA-256 |
| --- | --- |
| `data/external/original-workspace/task029-external-replication-feasibility/FROZEN_REPLICATION_DESIGN.md` | `3647de824054f1a01597e8161ad04aeb3f00ea12efbd714c01b892caf8fcea43` |
| `data/external/original-workspace/task029-external-replication-feasibility/SEAL_MANIFEST_replication.sha256` | `54ea5d95c12e83203b2977d70655d5a368c00f66df17006ae02ff78bfca06743` |
| `data/external/original-workspace/task029-external-replication-feasibility/execution/intermediate/INTERMEDIATE_SEAL.sha256` | `1ac085f8fd2939fa526ea1fa4ed027a312ab322b78c75870f2bfcccd4b32309f` |
| `data/external/original-workspace/task029-external-replication-feasibility/execution/intermediate/log2tpm_Discovery.npy` | `e108c1d33ab2a214959a8c1a0a8e276c33d47853127d3c66968622fd735274c7` |
| `data/external/original-workspace/task029-external-replication-feasibility/execution/intermediate/log2tpm_Confirmatory.npy` | `bd50b7ad762ff4d5573153878d264245424ef11fa8516ebc2462f7ec454d4b02` |
| `data/external/original-workspace/task029-external-replication-feasibility/execution/intermediate/acq02_join_ledger.json` | `7bfa84133e84070086782f23c2e30a6f2c15487ee4bf039e6d8eab7b187b4713` |
| `data/external/original-workspace/task029-external-replication-feasibility/execution/scripts/exec01_score_model_meta.py` | `6d2459227196868822f778cb9abfbb84e8ea04358d6298625c7fc086dc40b8fd` |
| `data/external/original-workspace/task029-external-replication-feasibility/execution/results/REPRODUCIBILITY_MANIFEST.json` | `94bfcb5c6a12e840a0363ceb250dd94c70dab2d5be4570095837d69bb4db3f9f` |
| `data/external/original-workspace/task029-external-replication-feasibility/execution/results/RESULT_SEAL.sha256` | `8653ebbe1f884d335a9efb8dfa7654e6e2e225b71d4a3395a975d69cfa96ad4a` |

Run the standard CPTAC design, intermediate, and result seal checks from their
owning directories and require zero mismatch. Required data version remains GDC
Data Release 45.0, 2025-12-04,
commit `8f7c2a51ab0084b216ad1b62a3fae8b945439c53`, STAR Counts, GENCODE v36,
60,660-gene universe. The Dou-2023 confirmatory roster remains pinned to SHA-256
`2ea92c4279918c3a6158b24ebfa1e0e36ffc2876d50fd2ed3f07265a389e0e31`
and size 277,385 bytes.

### Git and environment guards

- Required analysis worktree:
  `data/external/original-workspace/revgate-task-a-perclass`.
- Required branch: `experiment/task-a-perclass-c2-sensitivity`.
- Required base HEAD: `83503bad47b60193598b2b9ebe819c22c83e8ac1`.
- Original dirty worktree, verified TASK-030 worktree, and TASK-028/TASK-029
  upstream directories are read-only.
- The experimenter records start/end branch, HEAD, and full porcelain hash.
- No commit, push, stage, manuscript edit, or `src/**` edit.

### Deterministic seeds

TASK-A master seed is fixed at `20260722`. For every stochastic unit, derive

    int(sha256("20260722:" + step_id).hexdigest()[:8], 16)

where `step_id` is exactly

    boot__<model>__<target>
    perm__<model>__<target>

with model and target names exactly as fixed above. A single bootstrap draw for
a target/model must calculate all three contrasts, preserving their joint
algebra. Set `OPENBLAS_NUM_THREADS=1`, `OMP_NUM_THREADS=1`, and
`PYTHONHASHSEED=0`. Record the Python executable and all imported package
versions. Run the complete analysis twice and require byte-identical scientific
TSV/JSON/SVG output.

## 13. Required deliverable names

All new files are confined to `experiments/taskA_perclass_c2/`:

- `HYPOTHESIS_AND_ANALYSIS_CHARTER.md` (this file);
- `FROZEN_POSTHOC_SPEC.json`;
- `run_taskA_perclass_c2.py`;
- `COMMANDS.txt`;
- `SEED_SCHEME.md`;
- `SESSION_ENVIRONMENT.json`;
- `DEPENDENCIES.lock`;
- `INPUT_CHECKSUMS.tsv`;
- `REPRODUCIBILITY_MANIFEST.json`;
- `ANALYTICAL_REPORT.md`;
- `results/PER_CLASS_CONTRASTS_LONG.tsv`;
- `results/PER_CLASS_CONTRASTS.json`;
- `results/DESCRIPTIVE_BH18_FAMILIES.tsv`;
- `results/C2_RECONSTRUCTION_CHECKS.tsv`;
- `results/CPTAC_FIXED_EFFECT_META.tsv`;
- `results/INTERPRETATION_TAXONOMY.tsv`;
- `results/MANUSCRIPT_READY_POSTHOC_TABLE.tsv`;
- `results/MANUSCRIPT_READY_POSTHOC_TABLE.md`;
- `figures/GATA2_PER_CLASS_FOREST.svg`;
- `figures/SOX9_PER_CLASS_FOREST.svg`;
- `integrity/DETERMINISTIC_RERUN_COMPARISON.tsv`;
- `SHA256SUMS.txt`.

The manifest records every input/data version, every script hash, the seed and
subseed rule, parameters, sample counts, software, start/end Git state, algebra
checks, family sizes, two-run determinism, deviations, and the firewall
attestation. No result is citable without this manifest.

## 14. Experimenter return shape

Return raw results and provenance, not an advocacy verdict:

- `status`: COMPLETE or BLOCKED with exact reason;
- all deliverable paths;
- analytic and subtype counts;
- raw coefficients, SEs, residual SDs, d values, CIs, p/q values;
- all taxonomy labels generated mechanically from Section 7;
- all identity/meta caveat checks;
- reproducibility manifest, command log, input and output checksums;
- explicit statements that no upstream or manuscript byte changed and no
  frozen category or verdict was reassigned.
