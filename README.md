# APOE isoform-defining variants, AD-AMD pleiotropy, and circulating-protein mediation

This repository contains the analysis code and compact derived source data for a manuscript examining opposing genetic effects of APOE isoform-defining variants in Alzheimer's disease (AD) and age-related macular degeneration (AMD), together with the detectable limits of circulating-protein mediation.

## Repository status

This is a **pre-submission private repository**. The intended release sequence is:

1. Keep the repository private while the manuscript and code are frozen.
2. Create a versioned GitHub release at submission.
3. Make the repository public when disclosure is appropriate for the target journal.
4. Archive the release in Zenodo or another DOI-granting repository after the final public snapshot is verified.

No manuscript text, participant-level data, controlled summary statistics, access credentials, or third-party software bundles are stored here.

## Scientific scope

The code supports four linked evidence layers:

1. Systematic review and meta-analysis of the association between AMD and subsequent AD or dementia.
2. Genome-wide and local genetic architecture analyses using LDSC, HDL, MiXeR, LAVA, HyPrColoc, and MAGMA.
3. APOE anchoring through rs429358-C and rs7412-T, genome-wide bidirectional MR, APOE-region exclusion, and approximate conditional analysis.
4. Targeted triangulation of a literature-prioritized circulating-protein panel, followed by a non-UKB deCODE SomaScan same-platform coverage sensitivity analysis for excluded genes.

The protein layer is **not** a de novo proteome-wide discovery scan. Published proteome-wide MR or PWAS results define an enriched candidate set. Only proteins with both an APOE variant-to-protein estimate (alpha) and a protein-to-outcome estimate (beta) enter two-step mediation.

## Reproducibility levels

| Level | What can be reproduced | Additional access required |
|---|---|---|
| A | Included standalone quantitative figures from committed source-data tables | R/Python packages only |
| B | Meta-analysis, literature provenance, MAGMA summaries, and selected result audits | Public input files listed in `docs/data_access.md` |
| C | Full protein MR, deCODE sensitivity, conditional analysis, and genome-wide architecture workflows | Controlled or provider-authorized pQTL data, GWAS files, LD reference files, and external executables |

The committed result tables permit numerical and graphical verification even when an external provider does not allow redistribution of its raw files. Bespoke layout code and assets for flowcharts, graphical abstracts, and multi-panel figure assembly are intentionally excluded; the corresponding quantitative results remain auditable in the analysis tables.

## Repository map

| Path | Purpose |
|---|---|
| `01_meta/` | Meta-analysis, effect-metric sensitivity, GRADE tables, and PRISMA content |
| `02_genetic_arch/` | LDSC, HDL, MiXeR, LAVA, HyPrColoc, and MAGMA |
| `03_causal_lock/` | Variant-level and genome-wide MR, including APOE-region exclusion |
| `P0_isoform/` | rs429358/rs7412 isoform-defining variant contrasts |
| `P0_finemap/` | GCTA-COJO approximate conditional analysis and source data |
| `A1_protein_upgrade/` | Literature-prioritized Olink panel, alpha/beta separation, mediation, provenance, and sensitivity analyses |
| `A1_dual_scale_mediation/` | deCODE SomaScan non-UKB, same-platform coverage sensitivity branch only; not an independent replication |
| `04_protein_mr/`, `05_mediation/`, `06_biology_category_models/`, `P0_biology_guided_panel/` | Exploratory biology-guided panel sensitivity analyses |
| `figures_submission/code/` | Standalone quantitative plotting scripts and compact panel data; flowchart and multi-panel assembly code is excluded, including the S4/S5 composite layouts |
| `figures_submission/source_data/` | Figure-level source-data tables |
| `workflow/` | Public-release validation and convenience entry points |
| `environment/` | Package lists and recorded R session information |

See `docs/repository_map.md` for the analysis-to-output crosswalk.

## Quick start

Clone the repository and validate the public package before running analyses:

```powershell
python workflow/validate_release.py
```

Install Python dependencies:

```powershell
python -m pip install -r environment/python-requirements.txt
```

Install the R packages listed in `environment/R-packages.txt`. Exact package versions for completed analyses are recorded under `environment/session_info/` where available.

To regenerate quantitative figures whose plotting code is included:

```powershell
powershell -ExecutionPolicy Bypass -File workflow/run_figures.ps1
```

Some figure scripts require licensed fonts or external system libraries for PDF/TIFF export. The SVG output is the primary vector artifact.

## External resources

Set local resources through environment variables instead of editing scripts:

```powershell
$env:A1_PROJECT_ROOT = (Get-Location).Path
$env:A1_RESOURCE_ROOT = 'D:\path\to\external_resources'
$env:A1_LD_PREFIX = 'D:\path\to\1000G_EUR\EUR'
$env:PLINK_BIN = 'D:\path\to\plink.exe'
$env:GCTA_BIN = 'D:\path\to\gcta64.exe'
```

For UKB-PPP downloads, use `SYNAPSE_AUTH_TOKEN` only as a process environment variable. For deCODE S3 downloads, configure an AWS profile outside this repository. Never commit tokens, AWS credentials, `.synapseConfig`, `.Renviron`, or provider-downloaded raw files.

Detailed access instructions and filenames are in `docs/data_access.md`.

## Main entry points

- Primary protein workflow: `A1_protein_upgrade/run_submission_core.ps1`
- Full authorized-data protein workflow: `A1_protein_upgrade/run_all.ps1`
- Submission sensitivity QA: `A1_protein_upgrade/run_submission_QA.ps1`
- deCODE same-platform sensitivity: `A1_dual_scale_mediation/run_decode_same_platform_extension.ps1`
- Quantitative main and supplementary figure export: `workflow/run_figures.ps1`
- Public package validation: `workflow/validate_release.py`

The workflows use deterministic seeds where stochastic procedures are required. Missing values are not treated as zero, and assay/proteoform mappings are not inferred from gene symbols alone.

## Data and licensing boundaries

Code is released under the MIT License. Derived source-data tables are provided for manuscript verification but remain subject to the terms of their underlying data providers; the MIT License does not relicense third-party data. See `NOTICE.md`.

## Citation

Use the metadata in `CITATION.cff`. A DOI and final manuscript citation will be added when the public release is archived.

## Contact

Zhenglin Yang and Jinlei Chen
First School of Clinical Medicine, Guangdong Medical University, Guangdong, China
