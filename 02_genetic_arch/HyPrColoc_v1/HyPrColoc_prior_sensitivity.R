#!/usr/bin/env Rscript

rm(list = ls())

required <- c("data.table", "hyprcoloc")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R package(s): ", paste(missing, collapse = ", "))

suppressPackageStartupMessages({
  library(data.table)
  library(hyprcoloc)
})

script_arg <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)][1])
HERE <- dirname(normalizePath(script_file, winslash = "/", mustWork = TRUE))
ROOT <- normalizePath(file.path(HERE, "..", ".."), winslash = "/", mustWork = TRUE)
INPUT <- file.path(HERE, "HyPrColoc_SNP_Level_Data_FullLAVA.csv")
OUT <- file.path(
  ROOT, "tables_submission", "supplementary_tables",
  "TableS5f_HyPrColoc_Prior_Sensitivity.tsv"
)

if (!file.exists(INPUT)) stop("Missing harmonised HyPrColoc input: ", INPUT)
x <- fread(INPUT)
required_cols <- c(
  "SNP", "BETA_AD", "SE_AD", "BETA_Dry", "SE_Dry",
  "BETA_Wet", "SE_Wet", "BETA_Any", "SE_Any"
)
if (!all(required_cols %in% names(x))) stop("Unexpected HyPrColoc input schema.")
if (nrow(x) != 2582L) stop("Expected 2,582 shared APOE-region SNPs, found ", nrow(x))

betas <- as.matrix(x[, .(BETA_AD, BETA_Dry, BETA_Wet, BETA_Any)])
ses <- as.matrix(x[, .(SE_AD, SE_Dry, SE_Wet, SE_Any)])
rownames(betas) <- rownames(ses) <- x$SNP
colnames(betas) <- colnames(ses) <- c("AD", "Dry_AMD", "Wet_AMD", "Any_AMD")

grid <- CJ(prior_1 = c(1e-5, 1e-4, 1e-3), prior_c = c(0.01, 0.02, 0.05))
rows <- vector("list", nrow(grid))
for (i in seq_len(nrow(grid))) {
  fit <- hyprcoloc(
    effect.est = betas,
    effect.se = ses,
    trait.names = colnames(betas),
    snp.id = x$SNP,
    prior.1 = grid$prior_1[i],
    prior.c = grid$prior_c[i]
  )
  result <- as.data.table(fit$results)
  if (nrow(result) != 1L) stop("Expected one four-trait cluster at grid row ", i)
  rows[[i]] <- data.table(
    prior_1 = grid$prior_1[i],
    prior_c = grid$prior_c[i],
    n_shared_snps = nrow(x),
    traits = result$traits,
    posterior_probability = result$posterior_prob,
    regional_probability = result$regional_prob,
    candidate_snp = result$candidate_snp,
    posterior_explained_by_candidate_snp = result$posterior_explained_by_snp,
    model_boundary = paste(
      "Prior-grid sensitivity within the HyPrColoc single-shared-variant model;",
      "does not resolve multiple causal variants or prove biological causality."
    )
  )
}
out <- rbindlist(rows)
fwrite(out, OUT, sep = "\t", na = "NA")
writeLines(capture.output(sessionInfo()), file.path(HERE, "HyPrColoc_prior_sensitivity_sessionInfo.txt"))
print(out)
