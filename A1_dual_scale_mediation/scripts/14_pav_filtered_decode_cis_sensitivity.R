#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
})

set.seed(20260719)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)

harmonized_file <- file.path(root, "data_processed", "decode_smp_harmonized_data_cis_only.tsv.gz")
pav_file <- file.path(root, "tables", "decode_smp_PAV_epitope_instrument_audit.tsv")
output_file <- file.path(root, "tables", "decode_smp_beta_results_cis_only_PAV_filtered.tsv")
instrument_output <- file.path(root, "data_processed", "decode_smp_harmonized_instruments_cis_only_PAV_filtered.tsv")
log_dir <- file.path(root, "logs", "decode_smp_cis_PAV_filtered")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(harmonized_file), file.exists(pav_file))
harmonized <- fread(harmonized_file, na.strings = c("NA", ""))
pav <- fread(pav_file, na.strings = c("NA", ""))
flagged <- unique(pav[PAV_epitope_risk_classification == "target_gene_PAV_in_source_LD_class", .(assay_target_ID, SNP)])

planned_assays <- 9L
planned_global_beta_family <- 36L
assays <- sort(unique(harmonized$assay_target_ID))
outcomes <- c("AD", "any_AMD", "dry_AMD", "wet_AMD")
metadata <- unique(harmonized[, .(gene_symbol, assay_target_ID)])
universe <- fread(file.path(root, "tables", "APOE_variant_to_decode_somascan_alpha.tsv"),
                  select = c("gene_symbol", "assay_target_ID", "UniProt_ID"))
universe <- unique(universe)
stopifnot(nrow(universe) == planned_assays)

empty_result <- function(meta, outcome, reason, n_before, n_after) {
  data.table(
    gene_symbol = meta$gene_symbol, assay_target_ID = meta$assay_target_ID, UniProt_ID = meta$UniProt_ID,
    normalization = "smp", outcome = outcome, method = "not_estimable", method_role = "primary",
    nsnp = n_after, n_harmonized_before_PAV_filter = n_before, beta = NA_real_, SE = NA_real_,
    P_value = NA_real_, P_value_display = "NA", minus_log10_P = NA_real_,
    CI_lower = NA_real_, CI_upper = NA_real_, OR = NA_real_, OR_CI_lower = NA_real_, OR_CI_upper = NA_real_,
    Cochran_Q = NA_real_, Q_df = NA_real_, Q_P = NA_real_, Egger_intercept = NA_real_,
    Egger_intercept_SE = NA_real_, Egger_intercept_P = NA_real_, mean_F = NA_real_, min_F = NA_real_,
    beta_status = "not_estimable", exclusion_reason = reason,
    beta_source = "deCODE_SomaScan_SMP_cis_only_target_gene_PAV_filtered_to_frozen_A1_outcomes",
    instrument_scope = "cis_only_target_gene_PAV_filtered",
    planned_tests_per_outcome = planned_assays, planned_global_beta_tests = planned_global_beta_family,
    P_FDR_within_outcome = NA_real_, P_Bonferroni_within_outcome = NA_real_,
    P_FDR_global_36 = NA_real_, P_Bonferroni_global_36 = NA_real_
  )
}

result_rows <- list()
kept_rows <- list()
for (assay_id in universe$assay_target_ID) {
  meta <- universe[assay_target_ID == assay_id][1]
  for (outcome_name in outcomes) {
    dat <- harmonized[assay_target_ID == assay_id & outcome == outcome_name & mr_keep %in% TRUE]
    n_before <- nrow(dat)
    if (n_before) {
      dat <- merge(dat, flagged[, .(assay_target_ID, SNP, target_gene_PAV_flag = TRUE)],
                   by = c("assay_target_ID", "SNP"), all.x = TRUE)
      dat <- dat[is.na(target_gene_PAV_flag)]
    }
    if (!nrow(dat)) {
      result_rows[[length(result_rows) + 1L]] <- empty_result(
        meta, outcome_name,
        if (n_before) "All harmonized cis instruments removed by target-gene PAV filter" else "No harmonized cis instruments",
        n_before, 0L
      )
      next
    }
    dat[, F_statistic := (beta.exposure / se.exposure)^2]
    kept_rows[[paste(assay_id, outcome_name, sep = "_")]] <- dat
    methods <- if (nrow(dat) == 1L) "mr_wald_ratio" else c("mr_ivw", "mr_weighted_median", "mr_egger_regression")
    estimates <- as.data.table(mr(as.data.frame(dat), method_list = methods))
    heterogeneity <- tryCatch(as.data.table(mr_heterogeneity(as.data.frame(dat))), error = function(e) data.table())
    pleiotropy <- tryCatch(as.data.table(mr_pleiotropy_test(as.data.frame(dat))), error = function(e) data.table())
    if (nrow(heterogeneity)) {
      estimates <- merge(estimates, heterogeneity[, .(method, Q, Q_df, Q_pval)], by = "method", all.x = TRUE)
    } else {
      estimates[, `:=`(Q = NA_real_, Q_df = NA_real_, Q_pval = NA_real_)]
    }
    estimates[, `:=`(
      gene_symbol = meta$gene_symbol, assay_target_ID = assay_id, UniProt_ID = meta$UniProt_ID,
      normalization = "smp", outcome = outcome_name,
      method_role = fifelse(method %in% c("Inverse variance weighted", "Wald ratio"), "primary", "sensitivity"),
      n_harmonized_before_PAV_filter = n_before,
      CI_lower = b - 1.96 * se, CI_upper = b + 1.96 * se,
      OR = exp(b), OR_CI_lower = exp(b - 1.96 * se), OR_CI_upper = exp(b + 1.96 * se),
      Egger_intercept = if (nrow(pleiotropy)) pleiotropy$egger_intercept[1] else NA_real_,
      Egger_intercept_SE = if (nrow(pleiotropy)) pleiotropy$se[1] else NA_real_,
      Egger_intercept_P = if (nrow(pleiotropy)) pleiotropy$pval[1] else NA_real_,
      mean_F = mean(dat$F_statistic), min_F = min(dat$F_statistic),
      beta_status = "reestimated", exclusion_reason = "NA",
      beta_source = "deCODE_SomaScan_SMP_cis_only_target_gene_PAV_filtered_to_frozen_A1_outcomes",
      instrument_scope = "cis_only_target_gene_PAV_filtered",
      planned_tests_per_outcome = planned_assays, planned_global_beta_tests = planned_global_beta_family,
      P_FDR_within_outcome = NA_real_, P_Bonferroni_within_outcome = NA_real_,
      P_FDR_global_36 = NA_real_, P_Bonferroni_global_36 = NA_real_
    )]
    setnames(estimates, c("b", "se", "pval", "Q", "Q_pval"),
             c("beta", "SE", "P_value", "Cochran_Q", "Q_P"), skip_absent = TRUE)
    result_rows[[length(result_rows) + 1L]] <- estimates
  }
}

results <- rbindlist(result_rows, fill = TRUE)
results[, `:=`(minus_log10_P = NA_real_, P_value_display = "NA")]
estimable <- which(results$beta_status == "reestimated" & is.finite(results$beta) & is.finite(results$SE) & results$SE > 0)
if (length(estimable)) {
  log_p <- log(2) + pnorm(abs(results$beta[estimable] / results$SE[estimable]), lower.tail = FALSE, log.p = TRUE)
  results[estimable, minus_log10_P := -log_p / log(10)]
  results[estimable, P_value := ifelse(log_p < log(.Machine$double.xmin), 0, exp(log_p))]
  results[estimable, P_value_display := ifelse(minus_log10_P > 300, "<1e-300", formatC(P_value, format = "e", digits = 3))]
}
primary_idx <- which(results$method_role == "primary" & results$beta_status == "reestimated")
for (outcome_name in outcomes) {
  idx <- which(results$outcome == outcome_name & results$method_role == "primary" & results$beta_status == "reestimated")
  if (length(idx)) {
    results[idx, P_FDR_within_outcome := p.adjust(P_value, method = "BH", n = planned_assays)]
    results[idx, P_Bonferroni_within_outcome := pmin(1, P_value * planned_assays)]
  }
}
if (length(primary_idx)) {
  results[primary_idx, P_FDR_global_36 := p.adjust(P_value, method = "BH", n = planned_global_beta_family)]
  results[primary_idx, P_Bonferroni_global_36 := pmin(1, P_value * planned_global_beta_family)]
}
setorder(results, assay_target_ID, outcome, method_role, method)
fwrite(results, output_file, sep = "\t", na = "NA")
if (length(kept_rows)) fwrite(rbindlist(kept_rows, fill = TRUE), instrument_output, sep = "\t", na = "NA")
writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo.txt"))
writeLines(c(
  paste0("Target-gene PAV-linked instruments removed: ", paste(flagged$assay_target_ID, flagged$SNP, sep = ":", collapse = ";")),
  "PAV-unresolved instruments were retained and not treated as no risk.",
  "This is a sensitivity analysis and does not redefine the primary extension."
), file.path(log_dir, "analysis_rules.txt"))
message("PAV-filtered primary beta estimates: ", length(primary_idx))
message("Output: ", output_file)
