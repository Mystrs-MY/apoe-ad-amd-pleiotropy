#!/usr/bin/env Rscript

# Global genetic correlation using GenomicSEM LDSC.
# The corrected GWAS files contain variant-specific N_effective for AD and
# case-control effective sample sizes for FinnGen R12 AMD. Results are written to a tagged
# directory so legacy outputs cannot be silently reused.

rm(list = ls())

required <- c("GenomicSEM", "data.table", "gdata")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R package(s): ", paste(missing, collapse = ", "))

suppressPackageStartupMessages({
  library(GenomicSEM)
  library(data.table)
})

set.seed(20260902)
script_arg <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)][1])
SCRIPT_DIR <- dirname(normalizePath(script_file, winslash = "/", mustWork = TRUE))
PROJECT_ROOT <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), winslash = "/", mustWork = TRUE)
RESOURCE_ROOT <- Sys.getenv("A1_RESOURCE_ROOT", unset = file.path(PROJECT_ROOT, "data", "external"))
RUN_TAG <- Sys.getenv("A1_LDSC_RUN_TAG", unset = "correctedN_20260903")
GWAS_DIR <- file.path(RESOURCE_ROOT, "GWAS")
OUT_DIR <- file.path(SCRIPT_DIR, paste0("results_", RUN_TAG))
LD_DIR <- Sys.getenv("A1_LDSC_LD_DIR", unset = file.path(SCRIPT_DIR, "LDscore"))
WEIGHT_DIR <- Sys.getenv("A1_LDSC_WEIGHT_DIR", unset = file.path(SCRIPT_DIR, "1000G_Phase3_weights_hm3_no_MHC"))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

input_files <- file.path(GWAS_DIR, c(
  "AD_Wightman_cleaned_hg19.tsv.gz",
  "AMD_Dry_R12_cleaned_hg19.tsv.gz",
  "AMD_Wet_R12_cleaned_hg19.tsv.gz",
  "AMD_H7_R12_cleaned_hg19.tsv.gz"
))
trait_names <- c("AD", "Dry_AMD", "Wet_AMD", "Any_AMD")
required_files <- c(input_files, file.path(LD_DIR, "w_hm3.snplist"))
if (any(!file.exists(required_files))) {
  stop("Missing required file(s):\n", paste(required_files[!file.exists(required_files)], collapse = "\n"))
}
if (!dir.exists(WEIGHT_DIR)) stop("Missing LDSC weight directory: ", WEIGHT_DIR)

old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)
setwd(OUT_DIR)

log_file <- file.path(OUT_DIR, "LDSC_run_console.log")
log_connection <- file(log_file, open = "at", encoding = "UTF-8")
sink(log_connection, split = TRUE)
sink(log_connection, type = "message")
on.exit({
  try(sink(type = "message"), silent = TRUE)
  try(sink(), silent = TRUE)
  try(close(log_connection), silent = TRUE)
}, add = TRUE)

munge(
  files = normalizePath(input_files, winslash = "/"),
  hm3 = normalizePath(file.path(LD_DIR, "w_hm3.snplist"), winslash = "/"),
  trait.names = trait_names,
  info.filter = 0.9,
  maf.filter = 0.01,
  column.names = list(
    SNP = "SNP", A1 = "A1", A2 = "A2", effect = "BETA",
    P = "P", N = "N"
  )
)

sumstats_files <- file.path(OUT_DIR, paste0(trait_names, ".sumstats.gz"))
if (any(!file.exists(sumstats_files))) stop("GenomicSEM munge did not produce all expected files.")

ldsc_results <- ldsc(
  traits = normalizePath(sumstats_files, winslash = "/"),
  sample.prev = rep(NA_real_, length(trait_names)),
  population.prev = rep(NA_real_, length(trait_names)),
  ld = normalizePath(LD_DIR, winslash = "/"),
  wld = normalizePath(WEIGHT_DIR, winslash = "/"),
  trait.names = trait_names,
  stand = TRUE
)
if (is.null(ldsc_results$S_Stand) || is.null(ldsc_results$V_Stand)) {
  stop("Standardized LDSC output is unavailable; inspect the LDSC log.")
}

r <- length(trait_names)
se_stand <- matrix(NA_real_, r, r)
se_stand[lower.tri(se_stand, diag = TRUE)] <- sqrt(diag(ldsc_results$V_Stand))
se_stand[upper.tri(se_stand)] <- t(se_stand)[upper.tri(se_stand)]
rownames(se_stand) <- colnames(se_stand) <- trait_names
rg_mat <- ldsc_results$S_Stand
rownames(rg_mat) <- colnames(rg_mat) <- trait_names

pairs <- t(combn(seq_along(trait_names), 2L))
formatted <- rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
  a <- pairs[i, 1]
  b <- pairs[i, 2]
  rg <- rg_mat[b, a]
  se <- se_stand[b, a]
  data.table(
    p1 = trait_names[a], p2 = trait_names[b], rg = rg, se = se,
    z_score = rg / se,
    p_value = 2 * pnorm(abs(rg / se), lower.tail = FALSE),
    ci_lower = rg - 1.96 * se,
    ci_upper = rg + 1.96 * se,
    analysis_run_tag = RUN_TAG
  )
}))

fwrite(formatted, file.path(OUT_DIR, "LDSC_Results_Formatted.csv"))
saveRDS(ldsc_results, file.path(OUT_DIR, "LDSC_full_result.rds"))
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
writeLines(
  c(
    paste0("analysis_run_tag=", RUN_TAG),
    "sample_size_field=N (explicit copy of N_EFFECTIVE)",
    "AD_N_definition=variant-specific_N_effective supplied by GCST013196",
    "AMD_N_definition=4/(1/N_cases+1/N_controls) from final FinnGen R12 counts",
    paste0("GWAS_directory=", normalizePath(GWAS_DIR, winslash = "/"))
  ),
  file.path(OUT_DIR, "run_metadata.txt")
)

print(formatted)
message("LDSC completed: ", file.path(OUT_DIR, "LDSC_Results_Formatted.csv"))
