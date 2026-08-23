# A1 Circulating-Protein Upgrade

## Scope

This directory implements the submission-facing circulating-protein layer for Article 1. The study uses published proteome-wide MR/PWAS findings to build an enriched, literature-prioritized candidate panel; it does not claim a de novo proteome-wide discovery scan.

The active A1 narrative has four sequential layers:

1. clinical anchoring by systematic review and cohort meta-analysis;
2. genome-wide pleiotropic architecture;
3. APOE isoform localization;
4. a literature-prioritized circulating-protein mediation stress test.

Drug-target, perturbational, pharmacovigilance, and derivative C3 full-region analyses are outside the scope of the manuscript and are not distributed in this repository.

## Protein Evidence Layers

1. `literature_source_universe`: verified positive protein-outcome records from published AD/AMD proteome-wide studies.
2. `primary_literature_prioritized_panel`: high-confidence protein-evidence entries after provenance, multiplicity, mapping, and evidence-tier review.
3. `APOE_linkable_mediation_subset`: assays with both direct APOE variant-to-protein alpha and re-estimable protein-to-outcome beta.
4. `exploratory_biology_guided_sensitivity_panel`: an author-defined 30-protein candidate panel retained only for sensitivity and biological context.

The five-protein PWAS crosswalk (BCAM, CD55, LILRB1, LILRB5, and SCARA5) is a prespecified sensitivity extension to layer 3. It does not redefine the 25-protein primary panel.

## Stable Results

- The integrated provenance master contains 345 row-level positive evidence records: 52 Tier 1 and 293 Tier 2. The four newly verified priority studies contribute 61 of the Tier 2 records.
- Layer 2 contains 41 protein-evidence entries representing 33 genes.
- Twenty-seven name-matched Olink assays have direct rs429358-C and rs7412-T alpha estimates.
- Twenty-five assays have both alpha and re-estimated beta and define 200 formal mediation paths.
- No individual mediation path survives FDR or Bonferroni correction in the unfiltered primary instrument analysis. One PAV-filtered rs429358-to-TREM2-to-AD path passes the planned path correction but fails the joint heterogeneity, cis-only, and Kunkle outcome robustness standard.
- The small rs429358-C-to-AD aggregate estimate under independent errors is not robust to positive error-correlation or high-confidence mapping sensitivity analyses.
- The five-protein extension retains corrected support for CD55-to-AD but not for APOE-linkable CD55 mediation and does not change the aggregate boundary.
- Under the planned cross-outcome 108-test protein-beta correction, only TREM2-to-AD and CSF2-to-wet-AMD remain significant; neither establishes robust APOE mediation.
- Four APOE-exclusion windows give 24 nonsignificant bidirectional-MR sensitivity estimates.
- The shared-instrument audit identifies 27 exact repeated pQTLs and 31 distinct high-LD SNP pairs; leave-one-protein-out estimates are reported as influence diagnostics rather than new inferential intervals.
- A Kunkle 2019 clinically diagnosed AD sensitivity covered 25/25 proteins and 385/409 prespecified instrument pairs; no corrected protein beta or mediation path was retained, and both aggregate intervals crossed zero.
- A published PAV/cis-eQTL audit linked 77/409 instrument pairs and identified two prespecified high-confidence epitope-risk instruments. Filtering attenuated CSF2-to-wet-AMD; the single corrected filtered TREM2 mediation path was not robust to heterogeneity, cis-only, and Kunkle outcome checks.
- A Lu 2026 externally defined post-publication sensitivity mapped five of eight candidates supported in at least two dataset-outcome cells and two independent datasets; it yielded no corrected protein beta or mediation path in genome-wide or cis-only analyses.

These results support a limited and assumption-sensitive mediation boundary. They do not estimate total mediation across the complete circulating proteome.

## Directory Structure

```text
A1_protein_upgrade/
|-- audit/                  project and protein-layer inventories
|-- config/                 inclusion, outcome, complex-locus, and resource rules
|-- literature/             search strategy, screening, evidence registry, and full-text records
|-- data_raw/               literature and non-redistributed source metadata
|-- data_processed/         candidate instruments, harmonized data, and intermediate objects
|-- tables/                 provenance, alpha, beta, mediation, comparison, and exclusion tables
|-- figures/                protein figure source data, contracts, previews, and QA
|-- scripts/                reproducible primary-panel scripts
|-- extensions/
|   |-- PWAS2026_crosswalk_extension/  prespecified five-protein sensitivity extension
|   `-- submission_validity_20260717/  Kunkle, PAV, and Lu 2026 bounded sensitivities
|-- logs/                   decisions, run logs, session information, and unresolved items
|-- manuscript_patch/       revised sections, table index, and change log
|-- run_all.ps1             primary protein-panel workflow
|-- run_PWAS2026_extension.ps1
|-- run_submission_core.ps1 submission-facing orchestration entry
`-- run_submission_QA.ps1   submission asset and table consistency checks
```

Controlled UKB-PPP archives are stored outside the repository at:

`<external_resource_root>/UKB-PPP/syn51365303_European_discovery`

Downloads use HTTP Range requests, resume from `.part` files, and require final Synapse MD5 verification. Tokens are read only from `SYNAPSE_AUTH_TOKEN` and must not be written to scripts, configuration files, tables, or logs.

## One-Command Run

With downloaded source archives already available:

```powershell
Set-Location <project_root>
powershell -ExecutionPolicy Bypass -File .\A1_protein_upgrade\run_submission_core.ps1 -SkipDownload
```

To download or resume controlled source files:

```powershell
$env:SYNAPSE_AUTH_TOKEN = '<personal access token>'
powershell -ExecutionPolicy Bypass -File .\A1_protein_upgrade\run_submission_core.ps1
Remove-Item Env:SYNAPSE_AUTH_TOKEN
```

Literature retrieval is refreshed only when `-RefreshLiterature` is explicitly supplied.

The targeted cross-module reproducibility workflow for the final technical checks is:

```powershell
powershell -ExecutionPolicy Bypass -File .\run_presubmission_technical_repro.ps1
```

It regenerates the effect-measure audit, GCTA-COJO conditional source tables, HyPrColoc prior sensitivity, APOE-window/protein sensitivities, complete MAGMA gene outputs, and final numerical QA. Submission figure layout and assembly remain outside the public package. Use `-SkipMAGMA` only when complete current `.genes.out` files already exist.

## Submission-facing data assets

- Figure source data: `<project_root>/figures_submission/source_data/`
- Included standalone plotting and data-reduction code: `<project_root>/figures_submission/code/`
- Supplementary tables: `<project_root>/tables_submission/supplementary_tables/`

The public package is intentionally limited to analysis reproduction and quantitative source-data verification. Bespoke layout code and assets for flowcharts, graphical abstracts, and multi-panel figure assembly are not distributed; the corresponding quantitative results remain available in the structured analysis tables.

## Key Files

- Inventory: `audit/A1_protein_layer_inventory.md`
- Search strategy: `literature/search_strategy.md`
- Provenance: `tables/Table_Literature_Prioritized_Protein_Provenance.tsv`
- Cross-platform mapping: `tables/protein_cross_platform_mapping.tsv`
- APOE alpha: `tables/APOE_variant_to_literature_proteins_alpha.tsv`
- Protein beta: `tables/literature_panel_beta_results.tsv`
- Main mediation: `tables/APOE_linkable_two_step_mediation.tsv`
- Cis-only mediation: `tables/APOE_linkable_two_step_mediation_cis_sensitivity.tsv`
- Strict mapping sensitivity: `tables/APOE_linkable_two_step_mediation_strict_sensitivity.tsv`
- Literature-versus-biology-guided comparison: `tables/literature_vs_biology_guided_panel_comparison.tsv`
- Decision log: `logs/decision_log.tsv`
- Submission table builder: `scripts/39_prepare_submission_facing_supplement_tables.R`
- Targeted sensitivity builder: `scripts/40_targeted_submission_sensitivities.R`
- Cross-module reproducibility entry: `../run_presubmission_technical_repro.ps1`

## Interpretation Boundaries

- Literature selection enriches the panel for previously positive protein-outcome beta estimates.
- The aggregate estimate is an enriched-panel stress test, not an unbiased estimate over all circulating proteins.
- Bootstrap intervals do not model publication selection, winner's curse, mapping uncertainty as a probability model, or empirical between-protein covariance.
- Cross-platform gene-symbol agreement alone does not authorize assay merging.
- TREM2 is retained with moderate mapping confidence because its measured cross-platform correspondence is traceable but its epitope and soluble or cleaved form remain incompletely resolved.
- The author-defined biology-category aggregation is not standard multivariable MR; its ridge-penalized output is reported only as a penalized multivariable sensitivity model.
