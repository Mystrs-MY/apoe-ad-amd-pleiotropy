# A1 submission-validity extension

This bounded extension was opened after the primary 25-protein analysis had been frozen. It does not redefine the primary literature-prioritized panel.

## Rules fixed before extension outcome analysis

1. **Clinically diagnosed AD sensitivity**: use Kunkle et al. 2019 Stage 1 GWAS. Full mediation and aggregate estimates require both APOE variants, at least 20 of 25 proteins, and at least 80% coverage of the 409 frozen protein-instrument pairs. Lower coverage permits only a partial direction-comparison table. Below 60% instrument coverage, or if either APOE variant cannot be harmonized, the aggregate branch stops.
2. **PAV/epitope-risk audit**: annotate exact variant-assay pairs. A direct or high-LD protein-altering variant is not automatically invalid. The primary exclusion sensitivity removes only high-confidence epitope-risk instruments, defined as target-gene PAV/high-LD PAV evidence without supporting target-gene cis-eQTL evidence when those fields are available. Unresolved antibody/aptamer epitopes remain unresolved rather than being imputed.
3. **Lu 2026 candidate gate**: this is an externally defined post-publication sensitivity, not a prespecified primary analysis. Candidates are taken from the source study's Fig. 6b observational APOE4-related mediator data, must be reported in at least two dataset-outcome cells and at least two independent datasets, and must be exact-assay linkable. Analysis proceeds only if at least five candidates can be mapped without gene-symbol-only coercion. These rules were fixed before the extension AD/AMD genetic results were inspected.

No module MVMR, full-region coloc, tissue-QTL bridge, larger-LD fine mapping, or cross-ancestry analysis is triggered unless one of the bounded analyses leaves a corrected, directionally stable result that could change the manuscript conclusion.

## Reproduction

Run the bounded extension from PowerShell:

```powershell
Set-Location <project_root>\A1_protein_upgrade\extensions\submission_validity_20260717
.\run_extension.ps1
```

The Synapse token is read only from the process environment variable `SYNAPSE_AUTH_TOKEN`. It is never written to a script, table, manifest, or log. Completed archives are checked against Synapse MD5 metadata; interrupted downloads remain as `.part` files and resume from their verified byte length.

The one-click entry preserves the frozen primary analysis. Generated sensitivity outputs are confined to this extension directory until explicitly packaged into submission-facing supplementary tables.

## Completed results and submission status

- The Kunkle 2019 Stage 1 gate passed: 25/25 proteins were estimable, 385/409 prespecified assay-instrument pairs were present in the outcome file, and 335 remained after conservative harmonization. No protein-to-AD estimate or two-step mediation path survived the planned correction. The aggregate estimates were 0.29% (95% bootstrap CI -1.18% to 1.81%) for rs429358-C and 1.08% (-4.17% to 6.26%) for rs7412-T.
- The PAV/epitope audit linked 77/409 pairs; 332 remain unresolved. CSF2 rs25882 and TREM2 rs188904277 met the prespecified high-confidence risk rule. Filtering attenuated CSF2-to-wet-AMD to P = 0.0919. One filtered rs429358-to-TREM2-to-AD path passed the planned 200-path correction, but it was not robust to heterogeneity, cis-only, and Kunkle outcome checks.
- The Lu 2026 observational gate yielded five exact-assay candidates from eight candidates supported in at least two dataset-outcome cells and two independent datasets. All direct APOE alpha estimates and all 20 protein-outcome estimates per instrument layer were estimable. Neither the genome-wide nor cis-only layer yielded a cross-outcome corrected protein beta or a corrected mediation path.
- No trigger condition was met for module MVMR, full-region multi-signal colocalization, tissue-QTL bridging, larger-reference fine mapping, or cross-ancestry analysis.
- Submission-facing copies are packaged as Tables S33a-S36c under `<project_root>/tables_submission/supplementary_tables/`. The primary 25-protein panel remains unchanged.
