# TASK-030 current-main analytical verification

Status: COMPLETE (experimenter raw computation; no reviewer verdict assigned here).

## Six-target numerical result

TCGA primary is the Aran-CPE-adjusted model on n=506: POLE 49, MMRd 148, NSMP 146, p53abn 163. TCGA no-purity is the pre-specified full-cohort sensitivity on n=507: POLE 49, MMRd 148, NSMP 147, p53abn 163. CPTAC consists of separately fitted Discovery n=95 (7/25/43/20) and Confirmatory n=135 (6/47/66/16) strata in POLE/MMRd/NSMP/p53abn order.

TCGA raw p is the add-one two-sided subtype-label permutation p on absolute d with 2,000 permutations. For each TCGA configuration, F1 is BH across 21 module omnibus tests; 17 modules gate into F2; F2 is BH across all 51 gated module-by-C1/C2/C3 tests. It is not BH across these six targets. CPTAC uses a distinct frozen BH-of-6 target family after inverse-variance fixed-effect meta-analysis of the two strata.

| target | contrast | tcga_primary_n | tcga_primary_d | tcga_primary_d_ci_lo | tcga_primary_d_ci_hi | tcga_primary_raw_permutation_p | tcga_primary_full_gated_F2_BH_q | tcga_primary_B2_12_credit_status | tcga_nopurity_n | tcga_nopurity_d | tcga_nopurity_d_ci_lo | tcga_nopurity_d_ci_hi | tcga_nopurity_raw_permutation_p | tcga_nopurity_full_gated_F2_BH_q | tcga_nopurity_B2_12_credit_status | cptac_meta_d | cptac_exact_frozen_verdict | cptac_replication_credit_status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GATA2 | C2 | 506 | -0.5530308599103358 | -0.8553990734860628 | -0.2865297115768852 | 0.0004997501249375 | 0.0019605581824472 | UNCREDITED_DESCRIPTIVE_POSITIVE_GENE_LOO_FAILURE_NOT_CAT5 | 507 | -0.5390942111232796 | -0.8424245104295747 | -0.2843641636819755 | 0.0004997501249375 | 0.0018205183122724 | UNCREDITED_DESCRIPTIVE_POSITIVE_GENE_LOO_FAILURE_NOT_CAT5 | -0.8772657732024451 | externally_replicated | TARGET_LEVEL_DIRECTIONAL_CONCORDANCE_ONLY_NOT_CONFIRMATORY_EXTERNAL_REPLICATION_OF_B2_12_CREDITED_TCGA_DISCOVERY |
| SOX9 | C2 | 506 | -0.5271109871580488 | -0.8049003736154324 | -0.2865659764584904 | 0.0004997501249375 | 0.0019605581824472 | UNCREDITED_DESCRIPTIVE_POSITIVE_GENE_LOO_FAILURE_NOT_CAT5 | 507 | -0.5089779877688975 | -0.7751319516749414 | -0.2604085787232398 | 0.0004997501249375 | 0.0018205183122724 | UNCREDITED_DESCRIPTIVE_POSITIVE_GENE_LOO_FAILURE_NOT_CAT5 | -0.7886590700390439 | externally_replicated | TARGET_LEVEL_DIRECTIONAL_CONCORDANCE_ONLY_NOT_CONFIRMATORY_EXTERNAL_REPLICATION_OF_B2_12_CREDITED_TCGA_DISCOVERY |
| HOXA9 | C2 | 506 | -0.64601938358519 | -0.9211556693192235 | -0.3817487893633889 | 0.0004997501249375 | 0.0019605581824472 | UNCREDITED_DESCRIPTIVE_POSITIVE_GENE_LOO_FAILURE_NOT_CAT5 | 507 | -0.614332050975297 | -0.8936548768618777 | -0.3542522721139036 | 0.0004997501249375 | 0.0018205183122724 | UNCREDITED_DESCRIPTIVE_POSITIVE_GENE_LOO_FAILURE_NOT_CAT5 | -0.3328812474576446 | evaluable--not_replicated | EVALUABLE_NOT_REPLICATED |
| WT1 | C2 | 506 | -0.5909520933744161 | -0.8780963850878176 | -0.3490396482363203 | 0.0004997501249375 | 0.0019605581824472 | UNCREDITED_DESCRIPTIVE_POSITIVE_GENE_LOO_FAILURE_NOT_CAT5 | 507 | -0.5982549355797853 | -0.8614831264931119 | -0.3330639388039594 | 0.0004997501249375 | 0.0018205183122724 | UNCREDITED_DESCRIPTIVE_POSITIVE_GENE_LOO_FAILURE_NOT_CAT5 | -0.1746230284630709 | evaluable--not_replicated | EVALUABLE_NOT_REPLICATED |
| PAX8 | C1 | 506 | 0.6930742631907577 | 0.4734053139677799 | 0.9241928046490812 | 0.0004997501249375 | 0.0019605581824472 | CAT1_ENRICHMENT_CREDITED_B2.12 | 507 | 0.6864853341124657 | 0.4713836272632972 | 0.9248342856375392 | 0.0004997501249375 | 0.0018205183122724 | CAT1_ENRICHMENT_CREDITED_B2.12 | 0.795493738253157 | underpowered_sensitivity | SENSITIVITY_ONLY_DIRECTIONALLY_CONSISTENT_NOT_PRIMARY_REPLICATION_VERDICT |
| LHX1 | C1 | 506 | 0.9616207324967918 | 0.7344898379913432 | 1.2071228317759897 | 0.0004997501249375 | 0.0019605581824472 | CAT1_ENRICHMENT_CREDITED_B2.12 | 507 | 0.9564277469651766 | 0.7249735680318201 | 1.216572059935623 | 0.0004997501249375 | 0.0018205183122724 | CAT1_ENRICHMENT_CREDITED_B2.12 | 0.4699858237214589 | underpowered_sensitivity | SENSITIVITY_ONLY_OPPOSITE_DIRECTION_VETO_NOT_REPLICATED |

All six legacy pointwise categories are unchanged between primary and no-purity. Under the complete B2.12 universal gene-LOO gate, GATA2, SOX9, HOXA9, and WT1 are uncredited/descriptive directional positives in both configurations because at least one mapped-gene deletion falls below |d|=0.50. PAX8 and LHX1 remain credited CAT1 in both. These are not new category labels and the four failed rows are not CAT5. The exact per-criterion record is in `results/CATEGORY_CHANGE_REPORT.tsv`.

The frozen TASK-029 numeric verdicts are retained as bytes and values. GATA2 and SOX9 are target-level directional concordances only, not confirmatory external replication of a B2.12-credited TCGA discovery. HOXA9 and WT1 are evaluable-not-replicated. PAX8 and LHX1 were frozen sensitivity-only in CPTAC.

## Purity coefficient/residual decomposition

GATA2 and SOX9 correlations, coefficients, residual SDs, standardized d values, and the exact symmetric Shapley decomposition are in `results/GATA2_SOX9_DECOMPOSITION.tsv` and `.md`. Every row verifies d=coefficient/residual_SD and the two attribution components sum to delta d within floating-point error.

## Design-matrix and estimand comparison

Overall verdict: COMPATIBLE-ONLY. The no-purity term omission and fixed C1/C2 contrast vectors are identical at the design-rule level. Cohort structure, subtype classifier implementation, within-dataset score scaling/filtering, site handling, inferential family, and the CPTAC two-stratum meta-analysis are not identical. Field-level evidence is in `results/COVARIATE_COMPATIBILITY.tsv`.

Required wording: model-adapted external replication, not identical covariate replication.

## Integrity checks

- All sealed TASK-028 and TASK-029 inputs verified with zero mismatch.
- Full PRIMARY and no-purity F1/F2 families regenerated from sealed scores/covariates and matched the sealed numerical model bytes at tolerance 1e-12.
- The full gene-LOO computation was rerun. Reconstructed scores matched at tolerance 1e-12; full family substitution was used, not BH-of-6.
- Every selected scientific result file was byte-identical across two complete runs: True.
- All 18 scientific files were byte-identical to cycle 1: True. The sole allowed prior-output mismatch is the declared base manuscript-candidate TSV narrative/credit correction: True.
- DOCX/PDF inventory, sizes, mtimes, and hashes were unchanged across the second complete run: True. The delivered package contains no DOCX/PDF.
- No threshold, weight, category rule, target, contrast, family, score definition, edge set, model term, or CPTAC verdict changed.

This package reports raw results only and does not assign a reviewer verdict.
