# External data access

Raw and controlled inputs are intentionally excluded from this repository. Recreate the local resource tree under `data/external/`, or set `A1_RESOURCE_ROOT` to an equivalent location.

## Outcome GWAS

| Resource | Analysis role | Access route |
|---|---|---|
| Wightman et al. AD GWAS | Primary AD outcome | NHGRI-EBI GWAS Catalog accession `GCST013196`; the analysed public European-ancestry subset excludes UK Biobank and 23andMe and has variant-specific sample sizes up to 398,058 participants |
| Kunkle et al. 2019 Stage 1 | AD sensitivity analysis | GWAS Catalog accession `GCST007511`; the resumable downloader is in `A1_protein_upgrade/extensions/submission_validity_20260717/scripts/01_download_kunkle_stage1.ps1` |
| FinnGen Release 12 | Any, dry, and wet AMD outcomes | [FinnGen data-download documentation](https://finngen.gitbook.io/documentation/data-download) and the R12 manifest |

Expected harmonized filenames are listed in `A1_protein_upgrade/config/outcomes.yml`. The repository does not redistribute these files.

## Protein GWAS

| Resource | Analysis role | Access route |
|---|---|---|
| UKB-PPP European discovery | Primary Olink alpha and beta estimation | [Synapse syn51365303](https://www.synapse.org/Synapse%3Asyn51365303); provider approval is required |
| deCODE SomaScan v4 | Non-UKB same-platform coverage sensitivity for excluded genes; not an independent replication | Provider-authorized summary-data service for Ferkingstad et al., 2021, DOI `10.1038/s41588-021-00978-w` |

UKB-PPP scripts read `SYNAPSE_AUTH_TOKEN` from the process environment. deCODE scripts use an AWS-compatible named profile (default `decode-download`) configured outside the repository. Neither workflow writes credentials into code or result tables.

## Sample-size semantics

Harmonized case-control GWAS files retain both fields when available:

- `N_TOTAL`: total analysed participants at the variant; used by LAVA and
  MAGMA, which require the analysed sample count rather than the case-control
  effective sample size.
- `N_EFFECTIVE`: variant-specific case-control effective sample size; used by
  LDSC, HDL, MiXeR, and GCTA-COJO in this analysis.

For GCST013196, `N_TOTAL` is the public subset's variant-specific `n`, while
`N_EFFECTIVE` is the source-provided variant-specific `N_effective`. For the
FinnGen R12 endpoints, `N_TOTAL` is the endpoint-specific case plus control
count and `N_EFFECTIVE = 4/(1/N_cases + 1/N_controls)`. Case and control counts
are supplied separately to LAVA. Case-control colocalization/SuSiE routines use
`N_TOTAL` together with the case fraction. These fields must not be silently
interchanged.

The final FinnGen R12 counts used here are 8,570 cases and 329,258 controls for
dry AMD, 6,699 and 331,070 for wet AMD, and 12,495 and 461,686 for any AMD.

## LD and annotation references

- 1000 Genomes Phase 3 European reference panel, represented locally by the PLINK prefix `EUR/EUR`.
- MAGMA NCBI37.3 gene-location file and corresponding annotation resources.
- GRCh38-to-GRCh37 chain file for explicitly logged coordinate conversion in the deCODE branch.

## Third-party software

Install third-party software from its official distribution rather than copying it into this repository:

- PLINK 1.9 or a compatible release
- GCTA 1.95 or a compatible release
- MAGMA 1.10
- MiXeR
- HDL
- LAVA

Set `PLINK_BIN`, `GCTA_BIN`, `A1_LD_PREFIX`, and other paths through environment variables. Recorded versions and R package sessions are in `environment/` and module-specific logs.

Genetic-architecture modules additionally accept `A1_LDSC_LD_DIR`,
`A1_LDSC_WEIGHT_DIR`, `A1_HDL_LD_PATH`, `A1_LAVA_REF_PREFIX`,
`A1_MAGMA_GENE_LOC`, `MAGMA_BIN`, and `MIXER_PYTHON` as documented in their
entry-point scripts.

## Suggested local layout

```text
data/external/
|-- GWAS/
|-- EUR/
|   |-- EUR.bed
|   |-- EUR.bim
|   `-- EUR.fam
|-- ukbppp_proteins/
|   |-- ukbppp_merged/
|   `-- rsid_maps/
`-- UKB-PPP/
    `-- syn51365303_European_discovery/
```

Provider terms take precedence over this suggested layout. Do not upload restricted inputs to GitHub, cloud logs, or public archives.
