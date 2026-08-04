#!/usr/bin/env Rscript

# FDR sensitivity analysis for the fixed, author-defined biology-guided panel.
# This script does not discover or screen proteins outside that candidate set.

suppressPackageStartupMessages({
  library(data.table)
})

set.seed(20240604)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this file with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = TRUE)
module_dir <- dirname(script_path)
project_root <- normalizePath(file.path(module_dir, ".."), winslash = "/", mustWork = TRUE)
result_dir <- file.path(module_dir, "results")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

mr_file <- file.path(project_root, "04_protein_mr", "C6_ukbppp_mr_all.csv")
alpha_file <- file.path(project_root, "04_protein_mr", "C6_rs429358_effects.csv")
if (!file.exists(mr_file)) stop("Missing biology-guided panel MR results: ", mr_file)
if (!file.exists(alpha_file)) stop("Missing biology-guided panel definition: ", alpha_file)

mr <- fread(mr_file)
panel <- fread(alpha_file)[Gene != "APOE", .(exposure = Gene, biology_category = Category)]

ivw <- mr[
  method == "Inverse variance weighted" & exposure %in% panel$exposure,
  .(exposure, outcome, b, se, pval)
]
if (!nrow(ivw)) stop("No eligible IVW estimates were found for the biology-guided panel.")
ivw <- merge(ivw, panel, by = "exposure", all.x = TRUE)
ivw[, `:=`(
  p_fdr_within_outcome = p.adjust(pval, method = "BH"),
  p_bonferroni_within_outcome = pmin(pval * .N, 1),
  log10p = -log10(pval)
), by = outcome]
ivw[, significance_class := fifelse(
  p_bonferroni_within_outcome < 0.05,
  "Bonferroni",
  fifelse(p_fdr_within_outcome < 0.05, "FDR", fifelse(pval < 0.05, "Nominal", "Not significant"))
)]

category_summary <- ivw[, .(
  n_proteins = uniqueN(exposure),
  n_bonferroni = sum(p_bonferroni_within_outcome < 0.05, na.rm = TRUE),
  n_fdr = sum(p_fdr_within_outcome < 0.05, na.rm = TRUE),
  n_nominal_only = sum(pval < 0.05 & p_fdr_within_outcome >= 0.05, na.rm = TRUE)
), by = .(biology_category, outcome)]

fwrite(
  ivw,
  file.path(result_dir, "biology_guided_panel_mr_fdr.csv")
)
fwrite(
  category_summary,
  file.path(result_dir, "biology_category_fdr_summary.csv")
)

cat(sprintf(
  "Biology-guided panel FDR sensitivity complete: %d proteins, %d protein-outcome estimates.\n",
  uniqueN(ivw$exposure),
  nrow(ivw)
))
