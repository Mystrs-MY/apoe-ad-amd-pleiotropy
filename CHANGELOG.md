# Changelog

## v1.1.0 - 2026-09-06

This correction release freezes the code and derived outputs used by the EJHG
submission. It does not replace the immutable `v1.0.0` release or its Zenodo
record.

### Corrected

- Identified the analysed AD resource as the publicly available
  European-ancestry GCST013196 subset that excludes UK Biobank and 23andMe,
  rather than the full Wightman consortium sample.
- Separated `N_TOTAL` from `N_EFFECTIVE` during GWAS preprocessing and assigned
  the appropriate sample-size field to each downstream method.
- Recalculated LDSC, HDL, MiXeR, LAVA, MAGMA, and GCTA-COJO outputs from the
  corrected inputs. In the corrected LAVA scan, only the APOE-containing locus
  2351 passes `0.05/2495` for each AD-AMD comparison; previous non-APOE local
  claims are not retained.
- Updated compact figure source data and provenance records to point to the
  corrected output snapshots.
- Replaced the previous asymmetric MiXeR reporting field with the symmetric
  Dice overlap, while retaining explicitly labelled trait-specific shared
  fractions and recording that single-fit uncertainty intervals are unavailable.
- Removed the inferential field `apoe_is_sole_causal_hub`; APOE-region exclusion
  is now reported through neutral attenuation and nominal-support diagnostics.
- Replaced machine-specific paths in public scripts and manifests with
  environment variables, logical resource identifiers, or project-relative
  paths.

### Unchanged analytical layers

Direct rs429358/rs7412 effect estimates, conventional MR estimates, and the
primary protein alpha/beta and two-step mediation calculations use effect
estimates and standard errors rather than the corrected sample-size field.
They were audited for input identity but were not rerun solely because of this
sample-size correction.

### Public-data boundary

No participant-level UK Biobank data, controlled raw pQTL/GWAS files,
credentials, flowchart layout code, multi-panel assembly code, or Figure 5a
assets are included.
