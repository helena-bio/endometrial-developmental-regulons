# Core developmental transcription-factor panel provenance

The M3 panel contains 20 transcription factors selected from developmental and Mullerian literature before tumour analysis. The factor inventory was fixed independently of TCGA and CPTAC expression results.

Primary targets were then assigned from the frozen CollecTRI signed-edge ledger. Only consensus-resolved edges enter primary scoring. The public `core_tf_panel.tsv` table reports, for each factor, the number of unique primary targets, positive and negative primary edges, unresolved excluded edges, and total admitted edges.

The panel contains 761 consensus-resolved primary edges. Counts are derived directly from `definitions/m3_primary_edge_ledger.tsv`; they are not reconstructed from manuscript prose.

Individual expected-direction priors are included only for the five factors for which the frozen source record states them explicitly: PAX8, SOX17, HOXA10, HOXA11, and WT1. No individual prior has been inferred for the remaining factors.

Source-record checksums:

- Factor inventory source SHA-256: `cf77b7106a09215461e72b7bea820bc006bd5a1178fbd6c0803b40b6b4ca56c9`
- Sealed edge ledger SHA-256: `96d3dd9a84f354a5d3676ea60ab8680ac10444b88cfceb6306a685572e607347`
