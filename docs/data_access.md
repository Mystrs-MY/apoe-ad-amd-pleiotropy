# External data access

Raw and controlled inputs are intentionally excluded from this repository. Recreate the local resource tree under `data/external/`, or set `A1_RESOURCE_ROOT` to an equivalent location.

## Outcome GWAS

| Resource | Analysis role | Access route |
|---|---|---|
| Wightman et al. AD GWAS | Primary AD outcome | NHGRI-EBI GWAS Catalog accession `GCST013196` and the original study data-access terms |
| Kunkle et al. 2019 Stage 1 | AD sensitivity analysis | GWAS Catalog accession `GCST007511`; the resumable downloader is in `A1_protein_upgrade/extensions/submission_validity_20260717/scripts/01_download_kunkle_stage1.ps1` |
| FinnGen Release 12 | Any, dry, and wet AMD outcomes | [FinnGen data-download documentation](https://finngen.gitbook.io/documentation/data-download) and the R12 manifest |

Expected harmonized filenames are listed in `A1_protein_upgrade/config/outcomes.yml`. The repository does not redistribute these files.

## Protein GWAS

| Resource | Analysis role | Access route |
|---|---|---|
| UKB-PPP European discovery | Primary Olink alpha and beta estimation | [Synapse syn51365303](https://www.synapse.org/Synapse%3Asyn51365303); provider approval is required |
| deCODE SomaScan | Independent same-platform sensitivity for excluded genes | Provider-authorized S3 distribution associated with Eldjarn et al., Nature 2023, DOI `10.1038/s41586-023-06563-x` |

UKB-PPP scripts read `SYNAPSE_AUTH_TOKEN` from the process environment. deCODE scripts use an AWS-compatible named profile (default `decode-download`) configured outside the repository. Neither workflow writes credentials into code or result tables.

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
