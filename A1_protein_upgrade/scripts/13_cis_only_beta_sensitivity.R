#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
  library(yaml)
})

set.seed(20260711)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
resources <- read_yaml(file.path(root, "config", "resources.yml"))
outcomes_config <- read_yaml(file.path(root, "config", "outcomes.yml"))

candidates <- fread(file.path(root, "data_processed", "literature_panel_pqtl_candidates.tsv"),
                    na.strings = c("NA", ""))
candidates <- candidates[cis_to_encoding_gene == TRUE & rsid_mapping_status == "mapped"]
provenance <- fread(file.path(root, "tables", "Table_Literature_Prioritized_Protein_Provenance.tsv"))
mapping <- fread(file.path(root, "tables", "APOE_linkable_target_assay_mapping.tsv"))
planned_family_size <- nrow(unique(mapping[tolower(eligible_for_primary) == "true",
                                           .(standardized_gene_symbol, UKB_PPP_OID)]))
output_file <- file.path(root, "tables", "literature_panel_beta_cis_sensitivity.tsv")
main_beta_file <- file.path(root, "tables", "literature_panel_beta_results.tsv")
iv_file <- file.path(root, "data_processed", "literature_panel_cis_only_IVs.tsv")
harm_file <- file.path(root, "data_processed", "literature_panel_cis_only_harmonized.tsv")
log_dir <- file.path(root, "logs", "cis_only_beta")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

plink_bin <- resources$paths$plink_binary
ld_ref <- resources$paths$ld_reference
outcome_specs <- list(
  AD = outcomes_config$outcomes$AD$local_file,
  any_AMD = outcomes_config$outcomes$any_AMD$local_file,
  dry_AMD = outcomes_config$outcomes$dry_AMD$local_file,
  wet_AMD = outcomes_config$outcomes$wet_AMD$local_file
)

clump_gene <- function(dat, gene) {
  dat <- dat[order(P_value)][, .SD[1], by = SNP]
  input <- file.path(log_dir, paste0(gene, "_cis_clump.tsv"))
  prefix <- file.path(log_dir, paste0(gene, "_cis_plink"))
  fwrite(dat[, .(SNP, P = P_value)], input, sep = "\t")
  status <- system2(plink_bin, c("--bfile", ld_ref, "--clump", input, "--clump-p1", "5e-8",
                                 "--clump-p2", "5e-8", "--clump-r2", "0.001", "--clump-kb", "10000",
                                 "--out", prefix),
                    stdout = paste0(prefix, ".stdout.txt"), stderr = paste0(prefix, ".stderr.txt"))
  if (status != 0 || !file.exists(paste0(prefix, ".clumped"))) stop("PLINK cis clump failed: ", gene)
  kept <- fread(paste0(prefix, ".clumped"), fill = TRUE)
  if (!"SNP" %in% names(kept)) return(dat[0])
  dat[SNP %in% kept$SNP]
}

iv_list <- lapply(sort(unique(candidates$gene_symbol)), function(gene) {
  ans <- clump_gene(candidates[gene_symbol == gene], gene)
  ans[, `:=`(F_statistic = (beta / SE)^2, analysis_set = "cis_only_sensitivity")]
  ans
})
names(iv_list) <- sort(unique(candidates$gene_symbol))
iv_list <- iv_list[vapply(iv_list, nrow, integer(1)) > 0]
if (!length(iv_list)) stop("No cis-only instruments remained after clumping.")
fwrite(rbindlist(iv_list, fill = TRUE), iv_file, sep = "\t", na = "NA")

format_exposure <- function(dat, gene) {
  dat <- copy(dat)
  exp <- format_data(as.data.frame(dat), type = "exposure", snp_col = "SNP",
                     beta_col = "beta", se_col = "SE", effect_allele_col = "effect_allele",
                     other_allele_col = "other_allele", eaf_col = "effect_allele_frequency",
                     pval_col = "P_value", samplesize_col = "sample_size",
                     chr_col = "chromosome", pos_col = "position_hg19_from_ID")
  exp$exposure <- gene
  exp$id.exposure <- gene
  exp
}
exposures <- lapply(names(iv_list), function(gene) format_exposure(iv_list[[gene]], gene))
names(exposures) <- names(iv_list)

results <- list()
harmonized_rows <- list()
for (outcome_name in names(outcome_specs)) {
  message("Loading ", outcome_name)
  out <- fread(outcome_specs[[outcome_name]], select = c("SNP", "CHR", "BP", "A1", "A2", "FREQ", "BETA", "SE", "P", "N"),
               showProgress = FALSE)
  needed <- unique(unlist(lapply(exposures, function(x) x$SNP)))
  out <- out[SNP %in% needed]
  for (gene in names(exposures)) {
    exposure <- exposures[[gene]]
    subset <- out[SNP %in% exposure$SNP]
    if (!nrow(subset)) next
    outcome <- format_data(as.data.frame(subset), type = "outcome", snp_col = "SNP",
                           beta_col = "BETA", se_col = "SE", effect_allele_col = "A1",
                           other_allele_col = "A2", eaf_col = "FREQ", pval_col = "P",
                           samplesize_col = "N", chr_col = "CHR", pos_col = "BP")
    outcome$outcome <- outcome_name
    outcome$id.outcome <- outcome_name
    harm <- harmonise_data(exposure, outcome, action = 2)
    harm <- harm[harm$mr_keep %in% TRUE, ]
    if (!nrow(harm)) next
    harm$F_statistic <- (harm$beta.exposure / harm$se.exposure)^2
    harmonized_rows[[paste(gene, outcome_name, sep = "_")]] <- as.data.table(harm)
    methods <- if (nrow(harm) == 1) "mr_wald_ratio" else c("mr_ivw", "mr_weighted_median", "mr_egger_regression")
    estimate <- as.data.table(mr(harm, method_list = methods))
    heterogeneity <- tryCatch(as.data.table(mr_heterogeneity(harm)), error = function(e) data.table())
    pleiotropy <- tryCatch(as.data.table(mr_pleiotropy_test(harm)), error = function(e) data.table())
    steiger <- tryCatch(as.data.table(directionality_test(harm)), error = function(e) data.table())
    if (nrow(heterogeneity)) estimate <- merge(estimate, heterogeneity[, .(method, Q, Q_df, Q_pval)], by = "method", all.x = TRUE)
    if (!"Q" %in% names(estimate)) estimate[, `:=`(Q = NA_real_, Q_df = NA_real_, Q_pval = NA_real_)]
    estimate[, `:=`(
      gene_symbol = gene, protein_name = gene, protein_form_or_isoform = "not_reported",
      outcome = outcome_name,
      method_role = fifelse(method %in% c("Inverse variance weighted", "Wald ratio"), "primary", "sensitivity"),
      CI_lower = b - 1.96 * se, CI_upper = b + 1.96 * se,
      OR = exp(b), OR_CI_lower = exp(b - 1.96 * se), OR_CI_upper = exp(b + 1.96 * se),
      Cochran_Q = Q, Q_P = Q_pval,
      Egger_intercept = if (nrow(pleiotropy)) pleiotropy$egger_intercept[1] else NA_real_,
      Egger_intercept_SE = if (nrow(pleiotropy)) pleiotropy$se[1] else NA_real_,
      Egger_intercept_P = if (nrow(pleiotropy)) pleiotropy$pval[1] else NA_real_,
      Steiger_direction = if (nrow(steiger)) steiger$correct_causal_direction[1] else NA,
      Steiger_P = if (nrow(steiger)) steiger$steiger_pval[1] else NA_real_,
      Steiger_method = "approximate_quantitative_trait_correlation_sensitivity",
      mean_F = mean(harm$F_statistic), min_F = min(harm$F_statistic),
      planned_family_size = planned_family_size, observed_primary_tests_in_family = NA_integer_,
      P_FDR_observed_reestimated = NA_real_, P_Bonferroni_planned_family = NA_real_,
      corrected_significance_FDR_observed = NA, corrected_significance_Bonferroni_planned = NA,
      beta_status = "reestimated_cis_only", exclusion_reason = "NA",
      beta_source = "primary_analysis_reestimated", availability_scope = "reestimated_subset",
      analysis_set = "cis_only_sensitivity"
    )]
    setnames(estimate, c("b", "se", "pval"), c("beta", "SE", "P_value"))
    results[[paste(gene, outcome_name, sep = "_")]] <- estimate
  }
  rm(out); gc()
}

result <- rbindlist(results, fill = TRUE)
assay_meta <- unique(candidates[, .(gene_symbol, assay_target_ID, UniProt_ID)])
stopifnot(!anyDuplicated(assay_meta$gene_symbol))
result <- merge(result, assay_meta, by = "gene_symbol", all.x = TRUE)
for (outcome_name in unique(result$outcome)) {
  idx <- which(result$outcome == outcome_name & result$method_role == "primary")
  result[idx, observed_primary_tests_in_family := length(idx)]
  result[idx, P_FDR_observed_reestimated := p.adjust(P_value, "BH")]
  result[idx, P_Bonferroni_planned_family := pmin(1, P_value * planned_family_size)]
  result[idx, corrected_significance_FDR_observed := P_FDR_observed_reestimated < 0.05]
  result[idx, corrected_significance_Bonferroni_planned := P_Bonferroni_planned_family < 0.05]
}
fwrite(result, output_file, sep = "\t", na = "NA")
fwrite(rbindlist(harmonized_rows, fill = TRUE, idcol = "analysis_id"), harm_file, sep = "\t", na = "NA")

main <- fread(main_beta_file, na.strings = c("NA", ""), fill = TRUE)
if (!"analysis_set" %in% names(main)) {
  main[, analysis_set := fifelse(beta_status == "reestimated", "genome_wide_instruments_primary", "planned_not_reestimable")]
}
main <- main[analysis_set != "cis_only_sensitivity"]
combined <- rbindlist(list(main, result), fill = TRUE)
fwrite(combined, main_beta_file, sep = "\t", na = "NA")
writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo.txt"))
message("Cis-only primary tests: ", nrow(result[method_role == "primary"]))
message("Output: ", output_file)
