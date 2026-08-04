# Independent deCODE SomaScan sensitivity

This module contains the same-platform deCODE SomaScan sensitivity analysis for eight genes excluded from the primary UKB-PPP Olink mediation panel.

The frozen branch evaluates nine exact aptamers and 72 two-step mediation paths. It remains separate from the 25-protein Olink panel and is not included in the primary aggregate mediation estimate.

Raw provider files, AWS credentials, and full pQTL objects are not distributed. Configure access outside the repository and run `run_decode_same_platform_extension.ps1` only after the provider-authorized inputs and external GWAS/LD resources are available.
