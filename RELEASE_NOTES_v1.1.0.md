# v1.1.0 correction release

This release freezes the correction-aligned reproducibility snapshot used by
the European Journal of Human Genetics submission.

## Corrected

- Identified the analyzed AD resource as the public European-ancestry
  GCST013196 subset excluding UK Biobank and 23andMe.
- Separated total analyzed sample size from case-control effective sample size
  and regenerated affected genetic-architecture outputs.
- Updated the LAVA scan to the corrected input configuration; only the
  APOE-containing block passes the prespecified threshold across all three
  AD-AMD comparisons.
- Reports MiXeR polygenic overlap with the symmetric Dice coefficient and
  labels asymmetric trait-specific shared fractions separately.
- Removes the unsupported `apoe_is_sole_causal_hub` field and uses neutral
  APOE-region attenuation diagnostics.

## Unchanged

The release does not alter the clinical meta-analysis pooled estimate, direct
rs429358/rs7412 effects, conventional MR estimates, or the primary Olink and
deCODE two-step mediation estimates.

## Data boundary

No participant-level UK Biobank data, controlled raw pQTL/GWAS files,
credentials, flowchart layout code, multi-panel assembly code, or Figure 5a
assets are included.
