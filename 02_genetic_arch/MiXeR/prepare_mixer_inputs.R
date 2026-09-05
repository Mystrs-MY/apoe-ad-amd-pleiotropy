#!/usr/bin/env Rscript

required <- "data.table"
if (!requireNamespace(required, quietly = TRUE)) stop("Missing R package: data.table")
suppressPackageStartupMessages(library(data.table))

script_arg <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)][1])
SCRIPT_DIR <- dirname(normalizePath(script_file, winslash = "/", mustWork = TRUE))
PROJECT_ROOT <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), winslash = "/", mustWork = TRUE)
RESOURCE_ROOT <- Sys.getenv("A1_RESOURCE_ROOT", unset = file.path(PROJECT_ROOT, "data", "external"))
RUN_TAG <- Sys.getenv("A1_MIXER_RUN_TAG", unset = "correctedN_20260903")
GWAS_DIR <- file.path(RESOURCE_ROOT, "GWAS")
OUT_DIR <- file.path(SCRIPT_DIR, paste0("inputs_", RUN_TAG))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

files <- c(
  AD_Wightman = "AD_Wightman_cleaned_hg19.tsv.gz",
  AMD_Dry = "AMD_Dry_R12_cleaned_hg19.tsv.gz",
  AMD_Wet = "AMD_Wet_R12_cleaned_hg19.tsv.gz",
  AMD_Any = "AMD_H7_R12_cleaned_hg19.tsv.gz"
)

manifest <- list()
for (trait in names(files)) {
  input <- file.path(GWAS_DIR, files[[trait]])
  output <- file.path(OUT_DIR, paste0(trait, ".mixer.gz"))
  if (!file.exists(input)) stop("Missing GWAS file: ", input)

  dt <- fread(input, select = c("SNP", "A1", "A2", "BETA", "SE", "N"))
  dt <- dt[is.finite(BETA) & is.finite(SE) & SE > 0 & is.finite(N) & N > 0]
  mixer <- dt[, .(snp = SNP, a1 = A1, a2 = A2, z = BETA / SE, n = N)]
  mixer <- mixer[is.finite(z)]
  fwrite(mixer, output, sep = "\t", quote = FALSE)

  manifest[[trait]] <- data.table(
    trait = trait,
    source_file = normalizePath(input, winslash = "/"),
    output_file = normalizePath(output, winslash = "/"),
    rows = nrow(mixer),
    n_min = min(mixer$n),
    n_median = median(mixer$n),
    n_max = max(mixer$n),
    sample_size_field = "N (explicit copy of N_EFFECTIVE)",
    sample_size_definition = ifelse(
      trait == "AD_Wightman",
      "variant-specific_N_effective supplied by GCST013196",
      "4/(1/N_cases+1/N_controls) from final FinnGen R12 counts"
    ),
    analysis_run_tag = RUN_TAG
  )
  rm(dt, mixer)
  gc()
}

fwrite(rbindlist(manifest), file.path(OUT_DIR, "MiXeR_input_manifest.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
message("MiXeR inputs prepared: ", OUT_DIR)
