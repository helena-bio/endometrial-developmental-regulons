# Presentation and cross-format audit

Status: PASS.

- Long TSV and CSV are exactly equal after parsing (90 rows x 59 fields). JSON contains the same 90 records and values; JSON key sorting changes column order and decimal round-trip differences are at most 1.78e-15.
- `SIX_TARGET_COMPLETENESS.tsv` is exactly the long table. `GATA2_SOX9_HETEROGENEITY.tsv` is exactly the GATA2/SOX9 subset. `INTERPRETATION_TAXONOMY.tsv` matches the independently reconstructed 30 target/model states. The presentation table contains all 90 planned rows and all required columns.
- Both forest plots use identical model/contrast order, sign convention, zero line, labels, and color mapping. They show POLE-NSMP, MMRd-NSMP, and direct POLE-MMRd for four direct fits plus a separately labeled CPTAC fixed-effect row; no significance stars appear. Plotted values and intervals originate from the audited long table.
- The manuscript-ready caption begins `Post-hoc explanatory sensitivity` and explicitly says the results cannot change frozen categories or replication verdicts.
- GATA2: TCGA primary/no-purity and CPTAC Confirmatory/meta are same-direction compatible; CPTAC Discovery is same-direction unresolved. SOX9: TCGA primary, CPTAC Discovery, and CPTAC meta are same-direction compatible; TCGA no-purity is distinguishable below the inherited materiality floor; CPTAC Confirmatory is unresolved. Direct POLE-minus-MMRd intervals, not comparisons of nominal p values, drive these labels.
- The four-target appendix is complete. PAX8 Discovery and HOXA9 Confirmatory are point-heterogeneous but unresolved; PAX8/LHX1 CPTAC meta are POLE-dominant under the frozen materiality taxonomy; all other labels match the frozen rule order.
- Prose consistently treats non-significance as compatibility/unresolved uncertainty, never equality. It is explicitly post-hoc, descriptive, non-causal, subtype-level bulk expression, and not an individual-patient biomarker, treatment-predictive result, purity-independence result, or proof of a class-mix artifact.

