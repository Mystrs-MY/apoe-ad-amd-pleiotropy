# Analysis-to-output map

## Clinical layer

| Analysis | Primary code | Committed verification data |
|---|---|---|
| Random-effects meta-analysis | `01_meta/scripts/meta_analysis.R` | `01_meta/source_data/Figure_2_source_data.csv` |
| Effect-metric audit and HR-only sensitivity | `01_meta/scripts/effect_metric_sensitivity.R` | `01_meta/results/` |
| GRADE | `01_meta/scripts/grade_assessment.R` | `01_meta/tables/Table4_GRADE_SoF.csv` and `TableS7_GRADE_Evidence_Profile.csv` |
| PRISMA screening counts and checklist locations | Layout code intentionally not distributed | `01_meta/PRISMA/Article1_PRISMA_content.json` |

## Genetic architecture layer

| Analysis | Primary code | Committed verification data |
|---|---|---|
| LDSC | `02_genetic_arch/LDSC/0.LDSC.R` | `LDSC_Results_Formatted.csv`; uses `N_EFFECTIVE` |
| HDL | `02_genetic_arch/HDL/0.HDL.R` | `Figure_3_HDL_results.csv`; uses `N_EFFECTIVE` |
| MiXeR and conditional Q-Q | `02_genetic_arch/MiXeR/prepare_mixer_inputs.R` and `run_MiXeR_corrected.sh` | JSON fits and summary CSV files; uses `N_EFFECTIVE` |
| LAVA | `02_genetic_arch/LAVA/0_prepare_LAVA_totalN_inputs.R` and `1.LAVA.R` | Three final local-correlation scans; uses `N_TOTAL` and endpoint case/control counts |
| HyPrColoc | `02_genetic_arch/HyPrColoc_v1/HyPrColoc.R` | APOE-region source and result tables |
| MAGMA | `02_genetic_arch/MAGMA/run_magma_gene_analysis.sh` and `summarise_magma.R` | Four complete `.genes.out` files and figure source data; uses `N_TOTAL` |

## APOE anchoring layer

| Analysis | Primary code | Committed verification data |
|---|---|---|
| rs429358 and rs7412 contrasts | `03_causal_lock/02_wald_ratio.R`, `P0_isoform/01_rs7412_wald.R` | Variant-effect tables under each module's `results/` |
| Genome-wide bidirectional MR | `03_causal_lock/03_gw_bidirectional.R` | Summary tables under `03_causal_lock/results/` |
| APOE-region exclusion | `03_causal_lock/04_apoe_exclusion.R` | `table2_with_vs_without_apoe.csv` |
| Approximate conditional analysis | `P0_finemap/04_gcta_cojo_apoe_conditional.R` | COJO outputs and `FigS6_GCTA_COJO_source_data.tsv`; uses `N_EFFECTIVE` |

The corrected genetic-architecture outputs are tagged `correctedN_20260903`.
The LAVA scan evaluates 2,495 prespecified loci per comparison and retains
estimable bivariate rows in the committed result files. Only locus 2351, which
contains APOE, passes the prespecified `0.05/2495` threshold in all three
AD-AMD comparisons.

## Circulating-protein layer

| Analysis | Primary code | Committed verification data |
|---|---|---|
| Literature search, screening, and provenance | `A1_protein_upgrade/scripts/01_*` to `06_*` | `literature/` and provenance tables |
| APOE variant-to-protein alpha | `07_extract_apoe_alpha.py` | `APOE_variant_to_literature_proteins_alpha.tsv` |
| Protein-to-outcome beta | `09_reestimate_literature_panel_beta.R` | `literature_panel_beta_results.tsv` |
| Two-step mediation | `11_two_step_mediation.R` | `APOE_linkable_two_step_mediation.tsv` and sensitivity variants |
| Panel finalization and comparison | `12_finalize_panel_tables.py` | `literature_vs_biology_guided_panel_comparison.tsv` and exclusion tables |
| deCODE same-platform extension | `A1_dual_scale_mediation/scripts/04_*` to `18_*` | deCODE alpha, beta, mediation, PAV, and shared-instrument tables |

The exploratory biology-guided panel remains under `04_protein_mr/`, `05_mediation/`, `06_biology_category_models/`, and `P0_biology_guided_panel/` as a sensitivity analysis. It is an author-defined candidate panel and does not imply a de novo proteome-wide scan.
