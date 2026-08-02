# Limitations

This repository records the completed computational analyses and their
documented reproducibility boundary. It does not remove the scientific,
technical, and data limitations of the underlying study.

## Observational and cohort-level design

The analyses are based on observational bulk-tumour transcriptomic data.

They support associations between molecular subtype and signed target-expression
patterns within the documented models. They do not establish:

- causal transcription-factor activity;
- direct transcription-factor binding;
- a mechanistic regulatory pathway;
- a developmental-state transition;
- therapeutic dependence;
- treatment response.

The unit of analysis is the tumour or patient record represented in each
cohort. Cohort-level associations must not be interpreted as direct
single-cell or patient-level mechanisms.

## Bulk RNA and cellular composition

Bulk RNA measurements combine signals from malignant, stromal, immune, and
other cells present in the tumour sample.

The primary TCGA model includes broad composition terms and tumour purity.
The no-purity and clinical-composition sensitivity analyses test selected
model dependencies, but they do not eliminate all possible composition
effects or unmeasured confounding.

ESTIMATE-derived stromal and immune scores are broad summaries rather than
cell-type-resolved measurements.

No malignant-cell-resolved RNA, chromatin, protein-activity, or direct-binding
measurement is available in the reported analysis.

## Regulon-score interpretation

The developmental regulon scores summarise signed expression patterns across
mapped target genes.

They are not direct measurements of:

- transcription-factor protein abundance;
- post-translational activation;
- nuclear localisation;
- DNA occupancy;
- enhancer or promoter activity;
- regulatory causality.

A regulon score may also be influenced by target-network coverage, gene
expression detectability, cohort-specific variance, and the realised gene
universe.

## Target-set robustness

Universal single-target-deletion robustness is a stringent property of the
specified target set and locked decision criteria.

Failure of the universal deletion gate does not erase the pointwise
association. It means the analysis cannot assign universal robustness credit
under the prespecified rule.

Passing the deletion gate does not establish external transportability,
causality, clinical validity, or independence from all possible target-set
definitions.

## External evaluation

The CPTAC-UCEC analysis is model-adapted rather than an identical replication
of the TCGA model.

The cohorts differ in:

- cohort structure;
- subtype implementation and annotation;
- platform and source-file composition;
- realised score distributions;
- residual variance;
- available covariates;
- multiplicity families;
- external fixed-effect meta-analysis.

Their standardized effects are therefore compatible in design logic but are
not interchangeable estimands.

The CPTAC POLE strata are small. This limits precision in class-specific and
direct POLE-versus-MMRd comparisons.

GATA2 and SOX9 meet the frozen CPTAC target-level conditions, but their
supported interpretation is directional concordance rather than confirmatory
external replication because the source TCGA signals lack universal deletion
credit.

## Underpowered external questions

C1 and external M1 equivalence testing were classified as underpowered before
the outcome-bearing CPTAC expression read.

Results from those questions remain sensitivity or unresolved evidence. They
must not be promoted because a point estimate appears favourable.

An underpowered result is not equivalent to evidence of no effect, no
difference, or equivalence.

## Clinical-variable sensitivity

Formal grade and histology sensitivity modelling is restricted to TCGA-UCEC.

The outcome-blind audit permitted:

- binary histology adjustment;
- endometrioid-restricted binary grade adjustment.

It excluded:

- CPTAC clinical modelling;
- all-histology grade models;
- joint grade-plus-histology models;
- grade-stratified models;
- interaction models;
- rescue models.

These exclusions were frozen before outcome read. They limit the clinical
composition claims that can be made.

The supported statements are descriptive sensitivity results, not causal
adjustment estimates.

Grade may be a confounder, a correlate of differentiation state, or partly
downstream of molecular subtype. The adjusted model is therefore not a
uniquely correct causal estimand.

Residual clinical-composition confounding in CPTAC remains possible.

## Post-hoc analyses

The per-class C2 decomposition and grade/histology sensitivity analyses are
post-hoc and descriptive.

They were separately frozen before their own outcome reads, but they cannot:

- alter a prior TCGA category;
- grant deletion-robustness credit;
- rescue a failed external result;
- replace a frozen CPTAC verdict;
- establish equality or equivalence.

A confidence interval spanning zero is not proof of equality.

## Multiple testing and thresholds

The study uses prespecified omnibus gates, effect-size floors, confidence
interval rules, and Benjamini-Hochberg families.

The resulting categories depend on those frozen criteria.

Changing an effect floor, multiplicity family, contrast weight, confidence
rule, or deletion criterion creates a different analysis and may change the
category assigned to a target.

No single p value or q value summarises the complete evidence profile.

## Fetal Mullerian module

The broad fetal Mullerian epithelial module is non-confirmatory under the
prespecified rules.

This does not establish a clean absence of every developmental signal.

Likewise, selected regulon differences do not establish uniform fetal
reactivation, reacquisition of a specific developmental trajectory, or
equivalence of the broad fetal programme across classes.

## Data availability

The repository does not redistribute every upstream data object used during
execution.

Large source datasets must be reacquired from their responsible repositories
under the applicable access and redistribution terms.

A complete rerun depends on continued availability of:

- the named source releases;
- source identifiers and mappings;
- compatible clinical annotations;
- the documented dependency environment.

Future source releases may differ from the frozen study inputs.

## Reproducibility limits

Matching repository checksums establishes byte identity of published files.

Independent reconstruction and deterministic reruns establish numerical and
provenance integrity within their documented scope.

They do not establish:

- biological correctness;
- absence of modelling bias;
- causal validity;
- clinical validity;
- generalisability to every cohort or platform;
- byte-identical floating-point output on every future environment.

A rerun using substituted inputs or altered software must document the
deviation.

## Clinical and translational boundary

The repository does not contain or validate a clinically deployable
biomarker.

The reported signals are not approved or established for:

- diagnosis;
- prognosis;
- treatment selection;
- patient stratification;
- clinical reporting;
- therapeutic targeting.

Further work would require malignant-cell-resolved measurements, independent
cohort evaluation, mechanistic experiments, predefined clinical endpoints,
and appropriate clinical validation.

## Repository boundary

The repository structure, documentation, figures, and public release records
do not expand the scientific claim beyond the frozen protocols,
machine-readable results, manuscript, and supplementary material.

Where presentation wording and a frozen result record appear to conflict, the
scientific record and explicit claim boundaries take precedence.
