#!/usr/bin/env Rscript

# Global genetic correlation using HDL.
# Uses corrected variant-specific N_effective for AD and case-control effective
# sample sizes for FinnGen R12 AMD. Legacy outputs are not overwritten.

rm(list = ls())

required <- c("HDL", "data.table")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R package(s): ", paste(missing, collapse = ", "))

suppressPackageStartupMessages({
  library(HDL)
  library(data.table)
})

set.seed(20260902)
script_arg <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)][1])
SCRIPT_DIR <- dirname(normalizePath(script_file, winslash = "/", mustWork = TRUE))
PROJECT_ROOT <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), winslash = "/", mustWork = TRUE)
RESOURCE_ROOT <- Sys.getenv("A1_RESOURCE_ROOT", unset = file.path(PROJECT_ROOT, "data", "external"))
RUN_TAG <- Sys.getenv("A1_HDL_RUN_TAG", unset = "correctedN_20260903")
GWAS_DIR <- file.path(RESOURCE_ROOT, "GWAS")
LD_PATH <- Sys.getenv("A1_HDL_LD_PATH", unset = file.path(SCRIPT_DIR, "UKB_imputed_SVD_eigen99_extraction"))
OUT_DIR <- file.path(SCRIPT_DIR, paste0("results_", RUN_TAG))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(LD_PATH)) stop("Missing HDL reference panel: ", LD_PATH)

trait_files <- c(
  AD = "AD_Wightman_cleaned_hg19.tsv.gz",
  Dry_AMD = "AMD_Dry_R12_cleaned_hg19.tsv.gz",
  Wet_AMD = "AMD_Wet_R12_cleaned_hg19.tsv.gz",
  Any_AMD = "AMD_H7_R12_cleaned_hg19.tsv.gz"
)
input_files <- setNames(file.path(GWAS_DIR, unname(trait_files)), names(trait_files))
if (!identical(names(input_files), names(trait_files)) || length(input_files) != 4L) {
  stop("HDL input discovery did not preserve the four prespecified trait names.")
}
if (any(!file.exists(input_files))) {
  stop("Missing GWAS file(s):\n", paste(input_files[!file.exists(input_files)], collapse = "\n"))
}

prep_for_hdl <- function(file_path) {
  dt <- fread(file_path, select = c("SNP", "A1", "A2", "N", "BETA", "SE"))
  dt <- dt[is.finite(BETA) & is.finite(SE) & SE > 0 & is.finite(N) & N > 0]
  dt[, Z := BETA / SE]
  dt <- dt[is.finite(Z), .(SNP, A1, A2, N, Z)]
  if (!nrow(dt)) stop("No valid HDL rows in ", file_path)
  dt
}

gwas <- lapply(input_files, prep_for_hdl)
if (!identical(names(gwas), names(trait_files)) || length(gwas) != 4L) {
  stop("HDL prepared-data list must contain named AD, Dry_AMD, Wet_AMD, and Any_AMD elements.")
}

run_hdl <- function(trait1, trait2) {
  message("HDL: ", trait1, " vs ", trait2)
  res <- HDL.rg(
    gwas1.df = gwas[[trait1]],
    gwas2.df = gwas[[trait2]],
    LD.path = LD_PATH
  )
  p_value <- res$rg.p
  if (is.null(p_value) || !is.finite(p_value)) {
    p_value <- 2 * pnorm(abs(res$rg / res$rg.se), lower.tail = FALSE)
  }
  data.table(
    Trait1 = trait1,
    Trait2 = trait2,
    rg = res$rg,
    se = res$rg.se,
    p = p_value,
    ci_lower = res$rg - 1.96 * res$rg.se,
    ci_upper = res$rg + 1.96 * res$rg.se,
    analysis_run_tag = RUN_TAG
  )
}

pairs <- t(combn(names(gwas), 2L))
hdl_results <- rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
  run_hdl(pairs[i, 1], pairs[i, 2])
}))

fwrite(hdl_results, file.path(OUT_DIR, "HDL_Results_Formatted.csv"))
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
writeLines(
  c(
    paste0("analysis_run_tag=", RUN_TAG),
    "sample_size_field=N (explicit copy of N_EFFECTIVE)",
    "AD_N_definition=variant-specific_N_effective supplied by GCST013196",
    "AMD_N_definition=4/(1/N_cases+1/N_controls) from final FinnGen R12 counts",
    paste0("GWAS_directory=", normalizePath(GWAS_DIR, winslash = "/")),
    paste0("LD_reference=", normalizePath(LD_PATH, winslash = "/"))
  ),
  file.path(OUT_DIR, "run_metadata.txt")
)

print(hdl_results)
message("HDL completed: ", file.path(OUT_DIR, "HDL_Results_Formatted.csv"))
