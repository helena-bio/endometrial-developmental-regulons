# TASK A analytical report

Status: COMPLETE producer execution; independent reviewer verdict pending.

This is a post-hoc explanatory sensitivity. It is descriptive and non-causal. It cannot change any frozen TCGA/CPTAC category or replication verdict, and no manuscript byte was edited.

## Guard and design result

All 21 pinned inputs, all nested seals, and the exact four cohort/count guards passed: TCGA primary 506 (49/148/146/163), TCGA no-purity 507 (49/148/147/163), CPTAC Discovery 95 (7/25/43/20), and CPTAC Confirmatory 135 (6/47/66/16), ordered POLE/MMRd/NSMP/p53abn. Each target/model used one unchanged full fit; no pairwise subset was refit.

## GATA2 and SOX9 per-class results

Values are coefficient (d), followed by the mechanically assigned direct same-fit state. The direct POLE-minus-MMRd model interval, not a comparison of nominal p values, determines class differentiation.

### GATA2

- TCGA_PRIMARY_CPE_N506: POLE-NSMP -5.36643807193979 (d=-0.682168557963428); MMRd-NSMP -3.33465442766363 (d=-0.423893161857245); direct POLE-MMRd -2.03178364427617 (d=-0.258275396106183, model 95% CI [-4.61677271858552, +0.553205430033192]); `SAME_DIRECTION_COMPATIBLE_BOTH_RESOLVED`.
- TCGA_NOPURITY_N507: POLE-NSMP -5.36079269280758 (d=-0.680326163708229); MMRd-NSMP -3.13505080664651 (d=-0.397862258538331); direct POLE-MMRd -2.22574188616107 (d=-0.282463905169898, model 95% CI [-4.80720630904612, +0.355722536723968]); `SAME_DIRECTION_COMPATIBLE_BOTH_RESOLVED`.
- CPTAC_DISCOVERY_N95: POLE-NSMP -4.77628445802661 (d=-0.627151723403765); MMRd-NSMP -5.71842661390091 (d=-0.750860033061692); direct POLE-MMRd +0.942142155874295 (d=+0.123708309657928, model 95% CI [-5.54313320278124, +7.42741751452983]); `SAME_DIRECTION_POINT_ESTIMATES_UNRESOLVED`.
- CPTAC_CONFIRMATORY_N135: POLE-NSMP -10.5080786235455 (d=-1.42267086857013); MMRd-NSMP -4.24454665626675 (d=-0.574661942919495); direct POLE-MMRd -6.26353196727874 (d=-0.848008925650636, model 95% CI [-12.6248481032991, +0.0977841687416303]); `SAME_DIRECTION_COMPATIBLE_BOTH_RESOLVED`.
- CPTAC fixed-effect meta: d_PN=-1.07116022884774, d_MN=-0.619906652517905, direct d_PM=-0.437389562696431 (95% normal CI [-1.07550814196303, +0.200729016570164]); `CPTAC_META_SAME_DIRECTION_COMPATIBLE_BOTH_RESOLVED`.

### SOX9

- TCGA_PRIMARY_CPE_N506: POLE-NSMP -5.55737261680548 (d=-0.690489527534985); MMRd-NSMP -2.92748355909365 (d=-0.363732446781112); direct POLE-MMRd -2.62988905771183 (d=-0.326757080753873, model 95% CI [-5.27459084572545, +0.0148127303017862]); `SAME_DIRECTION_COMPATIBLE_BOTH_RESOLVED`.
- TCGA_NOPURITY_N507: POLE-NSMP -5.58778537425003 (d=-0.689797682987325); MMRd-NSMP -2.65828394727454 (d=-0.32815829255047); direct POLE-MMRd -2.92950142697549 (d=-0.361639390436855, model 95% CI [-5.58332648012474, -0.275676373826232]); `SAME_DIRECTION_DISTINGUISHABLE_BELOW_MATERIALITY_FLOOR`.
- CPTAC_DISCOVERY_N95: POLE-NSMP -9.02863602816071 (d=-1.2888976105126); MMRd-NSMP -4.88479044277704 (d=-0.69733619894662); direct POLE-MMRd -4.14384558538367 (d=-0.591561411565977, model 95% CI [-10.1089023903846, +1.82121121961727]); `SAME_DIRECTION_COMPATIBLE_BOTH_RESOLVED`.
- CPTAC_CONFIRMATORY_N135: POLE-NSMP -5.23305969597327 (d=-0.714427003825959); MMRd-NSMP -4.15483260792515 (d=-0.567225444372903); direct POLE-MMRd -1.07822708804811 (d=-0.147201559453056, model 95% CI [-7.38672453430847, +5.23027035821224]); `SAME_DIRECTION_POINT_ESTIMATES_UNRESOLVED`.
- CPTAC fixed-effect meta: d_PN=-1.05009973419249, d_MN=-0.607093578794747, direct d_PM=-0.433721919315446 (95% normal CI [-1.040827833962, +0.173383995331106]); `CPTAC_META_SAME_DIRECTION_COMPATIBLE_BOTH_RESOLVED`.

For both primary targets, pooled C2 is directionally concordant with negative POLE-NSMP and MMRd-NSMP point effects in every direct cohort/stratum and CPTAC meta. GATA2 is compatible in TCGA and CPTAC Confirmatory/meta, with Discovery unresolved. SOX9 is compatible in TCGA primary, CPTAC Discovery/meta; TCGA no-purity has a supported magnitude distinction below the inherited materiality floor, and CPTAC Confirmatory is unresolved. None of these labels means equality or equivalence.

The equal-weight C2 is not sample-size weighted. POLE:MMRd sample proportions alone are not a class-mix artifact.

## Completeness, algebra, and multiplicity

All six targets and all 90 planned direct/meta rows are present. The four completeness targets show a mixture of compatible, unresolved, and point-heterogeneous states; the complete mechanical table is `results/INTERPRETATION_TAXONOMY.tsv`. No result was selected or repaired.

All 24 direct target/model coefficient and d identities passed scaled tolerance 1e-12, including every usable bootstrap replicate. CPTAC fixed-effect C2 identity is not required because residual scales and contrast-specific inverse-bootstrap-variance weights differ; its discrepancy and all weights are reported.

Exactly five separate descriptive BH-18 families were computed, each with 18 evaluable rows (no missing placeholders were needed). These q values do not confer confirmatory credit.

## Reproducibility and inference boundary

Two complete fresh-directory runs are byte-identical for all 27 produced files, including TSV, CSV, JSON, SVG, PNG, and PDF scientific artifacts. Upstream worktrees and all pinned bytes were unchanged.

Cross-cohort differences can reflect biology, composition, acquisition, platform, classifier, or scoring context; this analysis cannot identify cause. Bulk subtype-level results are not individual-patient biomarkers. Frozen verdicts are preserved and no manuscript edit was made.
