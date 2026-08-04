#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
  library(yaml)
})

set.seed(20260714)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = TRUE)
extension_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
upgrade_root <- normalizePath(file.path(extension_root, "..", ".."), winslash = "/", mustWork = TRUE)
resources <- read_yaml(file.path(upgrade_root, "config", "resources.yml"))
outcomes_config <- read_yaml(file.path(upgrade_root, "config", "outcomes.yml"))

candidate_file <- file.path(extension_root, "data_processed", "PWAS5_pqtl_candidates.tsv")
members_file <- file.path(extension_root, "config", "PWAS5_frozen_members.tsv")
output_file <- file.path(extension_root, "tables", "PWAS5_beta_main.tsv")
iv_output <- file.path(extension_root, "data_processed", "PWAS5_instruments_main.tsv")
harm_output <- file.path(extension_root, "tables", "PWAS5_harmonized_instruments_main.tsv")
log_dir <- file.path(extension_root, "logs", "beta_main")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

plink_bin <- resources$paths$plink_binary
ld_ref <- resources$paths$ld_reference
stopifnot(file.exists(plink_bin), file.exists(paste0(ld_ref, ".bed")))

outcome_specs <- list(
  AD = list(file = outcomes_config$outcomes$AD$local_file, label = "AD"),
  any_AMD = list(file = outcomes_config$outcomes$any_AMD$local_file, label = "any_AMD"),
  dry_AMD = list(file = outcomes_config$outcomes$dry_AMD$local_file, label = "dry_AMD"),
  wet_AMD = list(file = outcomes_config$outcomes$wet_AMD$local_file, label = "wet_AMD")
)

candidates <- fread(candidate_file, na.strings = c("NA", ""))
members <- fread(members_file, na.strings = c("NA", ""))
planned_assays <- unique(members[, .(gene_symbol, assay_target_ID = UKB_PPP_OID, UniProt_ID)])
planned_family_size <- nrow(planned_assays)
planned_family_total <- planned_family_size * length(outcome_specs)

candidates <- candidates[rsid_mapping_status == "mapped" & !is.na(SNP)]
setorder(candidates, gene_symbol, P_value)
candidates <- candidates[, .SD[1], by = .(gene_symbol, SNP)]

load_outcome <- function(path, label, keys) {
  stopifnot(file.exists(path))
  dat <- fread(path, select = c("SNP", "CHR", "BP", "A1", "A2", "FREQ", "BETA", "SE", "P", "N"),
               showProgress = FALSE)
  dat <- dat[SNP %in% keys]
  dat[, `:=`(outcome = label, id.outcome = label)]
  dat
}

plink_clump <- function(dat, gene) {
  clump_input <- file.path(log_dir, paste0(gene, "_clump_input.tsv"))
  prefix <- file.path(log_dir, paste0(gene, "_plink"))
  fwrite(dat[, .(SNP, P = P_value)], clump_input, sep = "\t")
  args <- c(
    "--bfile", ld_ref,
    "--clump", clump_input,
    "--clump-p1", "5e-8",
    "--clump-p2", "5e-8",
    "--clump-r2", "0.001",
    "--clump-kb", "10000",
    "--out", prefix
  )
  status <- system2(plink_bin, args = args, stdout = paste0(prefix, ".stdout.txt"),
                    stderr = paste0(prefix, ".stderr.txt"))
  clumped_file <- paste0(prefix, ".clumped")
  if (status != 0 || !file.exists(clumped_file)) {
    stop(sprintf("PLINK clumping failed for %s; inspect %s", gene, prefix))
  }
  clumped <- fread(clumped_file, fill = TRUE)
  if (!"SNP" %in% names(clumped)) return(dat[0])
  dat[SNP %in% clumped$SNP]
}

format_exposure <- function(dat, gene) {
  dat <- copy(dat)
  exposure <- format_data(
    as.data.frame(dat), type = "exposure",
    snp_col = "SNP", beta_col = "beta", se_col = "SE",
    effect_allele_col = "effect_allele", other_allele_col = "other_allele",
    eaf_col = "effect_allele_frequency", pval_col = "P_value",
    samplesize_col = "sample_size", chr_col = "chromosome", pos_col = "position_hg19_from_ID"
  )
  exposure$exposure <- gene
  exposure$id.exposure <- gene
  exposure
}

empty_result <- function(gene, outcome, status, reason, n_iv = 0L) {
  data.table(
    gene_symbol = gene, outcome = outcome, method = "not_estimable", method_role = "primary",
    nsnp = n_iv, beta = NA_real_, SE = NA_real_, P_value = NA_real_, CI_lower = NA_real_,
    CI_upper = NA_real_, OR = NA_real_, OR_CI_lower = NA_real_, OR_CI_upper = NA_real_,
    Cochran_Q = NA_real_, Q_df = NA_real_, Q_P = NA_real_, Egger_intercept = NA_real_,
    Egger_intercept_SE = NA_real_, Egger_intercept_P = NA_real_, Steiger_direction = NA,
    Steiger_P = NA_real_, Steiger_method = "not_estimable", mean_F = NA_real_, min_F = NA_real_,
    planned_family_size = planned_family_size, observed_primary_tests_in_family = NA_integer_,
    P_FDR_observed_reestimated = NA_real_, P_Bonferroni_planned_family = NA_real_,
    corrected_significance_FDR_observed = NA, corrected_significance_Bonferroni_planned = NA,
    beta_status = status, exclusion_reason = reason,
    beta_source = "PWAS2026_crosswalk_extension_same_assay",
    analysis_set = "PWAS5_genome_wide_instruments",
    planned_family_total = planned_family_total,
    P_FDR_all_estimable = NA_real_, P_Bonferroni_planned_20 = NA_real_,
    corrected_significance_FDR_all_estimable = NA,
    corrected_significance_Bonferroni_planned_20 = NA
  )
}

iv_rows <- list()
harm_rows <- list()
result_rows <- list()
clumped_by_gene <- list()

for (gene in sort(unique(candidates$gene_symbol))) {
  message("Clumping ", gene)
  gene_candidates <- candidates[gene_symbol == gene]
  clumped <- plink_clump(gene_candidates, gene)
  if (nrow(clumped)) {
    clumped[, F_statistic := (beta / SE)^2]
    clumped[, clump_r2 := 0.001]
    clumped[, clump_kb := 10000]
    iv_rows[[gene]] <- clumped
  }
  clumped_by_gene[[gene]] <- clumped
}

for (gene in setdiff(planned_assays$gene_symbol, names(clumped_by_gene))) {
  for (outcome in names(outcome_specs)) {
    result_rows[[length(result_rows) + 1L]] <- empty_result(
      gene, outcome, "not_estimable",
      "No genome-wide significant non-APOE-region pQTL candidate was available."
    )
  }
}

if (!length(iv_rows)) stop("No pQTL instruments remained after clumping.")
all_iv_keys <- unique(rbindlist(iv_rows, fill = TRUE)$SNP)
if (!length(all_iv_keys)) stop("No pQTL instruments remained after clumping.")
message("Loading each outcome GWAS once for ", length(all_iv_keys), " unique IV rsIDs")
outcomes_loaded <- lapply(names(outcome_specs), function(outcome_name) {
  spec <- outcome_specs[[outcome_name]]
  load_outcome(spec$file, outcome_name, all_iv_keys)
})
names(outcomes_loaded) <- names(outcome_specs)

for (gene in names(clumped_by_gene)) {
  clumped <- clumped_by_gene[[gene]]
  if (nrow(clumped) == 0) {
    for (outcome in names(outcome_specs)) {
      result_rows[[length(result_rows) + 1L]] <- empty_result(
        gene, outcome, "not_estimable", "No variants remained after PLINK clumping."
      )
    }
    next
  }
  exposure <- format_exposure(clumped, gene)

  for (outcome_name in names(outcome_specs)) {
    outcome_raw <- outcomes_loaded[[outcome_name]][SNP %in% exposure$SNP]
    if (nrow(outcome_raw) == 0) {
      result_rows[[length(result_rows) + 1L]] <- empty_result(
        gene, outcome_name, "not_estimable", "No clumped pQTL variants were present in the outcome GWAS.", nrow(clumped)
      )
      next
    }
    outcome <- format_data(
      as.data.frame(outcome_raw), type = "outcome",
      snp_col = "SNP", beta_col = "BETA", se_col = "SE",
      effect_allele_col = "A1", other_allele_col = "A2", eaf_col = "FREQ",
      pval_col = "P", samplesize_col = "N", chr_col = "CHR", pos_col = "BP"
    )
    outcome$outcome <- outcome_name
    outcome$id.outcome <- outcome_name
    harmonized <- harmonise_data(exposure, outcome, action = 2)
    harmonized <- harmonized[harmonized$mr_keep %in% TRUE, ]
    if (nrow(harmonized) == 0) {
      result_rows[[length(result_rows) + 1L]] <- empty_result(
        gene, outcome_name, "not_estimable", "No variants remained after allele harmonization.", nrow(clumped)
      )
      next
    }
    harmonized$F_statistic <- (harmonized$beta.exposure / harmonized$se.exposure)^2
    harmonized$steiger_limitation <- "Approximate quantitative-trait correlation used for binary outcome because calibrated liability prevalence was unavailable."
    harm_rows[[paste(gene, outcome_name, sep = "_")]] <- as.data.table(harmonized)

    methods <- if (nrow(harmonized) == 1) {
      "mr_wald_ratio"
    } else {
      c("mr_ivw", "mr_weighted_median", "mr_egger_regression")
    }
    estimates <- as.data.table(mr(harmonized, method_list = methods))
    if (nrow(estimates) == 0) {
      result_rows[[length(result_rows) + 1L]] <- empty_result(
        gene, outcome_name, "not_estimable", "TwoSampleMR returned no estimate.", nrow(harmonized)
      )
      next
    }
    heterogeneity <- tryCatch(as.data.table(mr_heterogeneity(harmonized)), error = function(e) data.table())
    pleiotropy <- tryCatch(as.data.table(mr_pleiotropy_test(harmonized)), error = function(e) data.table())
    steiger <- tryCatch(as.data.table(directionality_test(harmonized)), error = function(e) data.table())
    if (nrow(heterogeneity)) {
      estimates <- merge(estimates, heterogeneity[, .(method, Q, Q_df, Q_pval)], by = "method", all.x = TRUE)
    } else {
      estimates[, `:=`(Q = NA_real_, Q_df = NA_real_, Q_pval = NA_real_)]
    }
    estimates[, `:=`(
      gene_symbol = gene,
      outcome = outcome_name,
      method_role = fifelse(method %in% c("Inverse variance weighted", "Wald ratio"), "primary", "sensitivity"),
      CI_lower = b - 1.96 * se,
      CI_upper = b + 1.96 * se,
      OR = exp(b),
      OR_CI_lower = exp(b - 1.96 * se),
      OR_CI_upper = exp(b + 1.96 * se),
      Egger_intercept = if (nrow(pleiotropy)) pleiotropy$egger_intercept[1] else NA_real_,
      Egger_intercept_SE = if (nrow(pleiotropy)) pleiotropy$se[1] else NA_real_,
      Egger_intercept_P = if (nrow(pleiotropy)) pleiotropy$pval[1] else NA_real_,
      Steiger_direction = if (nrow(steiger)) steiger$correct_causal_direction[1] else NA,
      Steiger_P = if (nrow(steiger)) steiger$steiger_pval[1] else NA_real_,
      Steiger_method = "approximate_quantitative_trait_correlation_sensitivity",
      mean_F = mean(harmonized$F_statistic),
      min_F = min(harmonized$F_statistic),
      planned_family_size = planned_family_size,
      observed_primary_tests_in_family = NA_integer_,
      P_FDR_observed_reestimated = NA_real_,
      P_Bonferroni_planned_family = NA_real_,
      corrected_significance_FDR_observed = NA,
      corrected_significance_Bonferroni_planned = NA,
      beta_status = "reestimated",
      exclusion_reason = "NA",
      beta_source = "PWAS2026_crosswalk_extension_same_assay",
      analysis_set = "PWAS5_genome_wide_instruments",
      planned_family_total = planned_family_total,
      P_FDR_all_estimable = NA_real_, P_Bonferroni_planned_20 = NA_real_,
      corrected_significance_FDR_all_estimable = NA,
      corrected_significance_Bonferroni_planned_20 = NA
    )]
    setnames(estimates, c("b", "se", "pval", "Q", "Q_pval"),
             c("beta", "SE", "P_value", "Cochran_Q", "Q_P"), skip_absent = TRUE)
    result_rows[[length(result_rows) + 1L]] <- estimates
  }
}

results <- rbindlist(result_rows, fill = TRUE)
assay_meta <- planned_assays
stopifnot(!anyDuplicated(assay_meta$gene_symbol))
results <- merge(results, assay_meta, by = "gene_symbol", all.x = TRUE)
for (outcome_name in unique(results$outcome)) {
  idx <- which(results$outcome == outcome_name & results$method_role == "primary" & results$beta_status == "reestimated")
  if (length(idx)) {
    results[idx, observed_primary_tests_in_family := length(idx)]
    results[idx, P_FDR_observed_reestimated := p.adjust(P_value, method = "BH")]
    results[idx, P_Bonferroni_planned_family := pmin(1, P_value * planned_family_size)]
    results[idx, corrected_significance_FDR_observed := P_FDR_observed_reestimated < 0.05]
    results[idx, corrected_significance_Bonferroni_planned := P_Bonferroni_planned_family < 0.05]
  }
}
all_primary <- which(results$method_role == "primary" & results$beta_status == "reestimated")
if (length(all_primary)) {
  results[all_primary, P_FDR_all_estimable := p.adjust(P_value, method = "BH")]
  results[all_primary, P_Bonferroni_planned_20 := pmin(1, P_value * planned_family_total)]
  results[all_primary, corrected_significance_FDR_all_estimable := P_FDR_all_estimable < 0.05]
  results[all_primary, corrected_significance_Bonferroni_planned_20 := P_Bonferroni_planned_20 < 0.05]
}

ordered_columns <- c(
  "gene_symbol", "assay_target_ID", "UniProt_ID", "outcome", "method", "method_role", "nsnp", "beta", "SE", "P_value",
  "CI_lower", "CI_upper", "OR", "OR_CI_lower", "OR_CI_upper", "Cochran_Q", "Q_df", "Q_P",
  "Egger_intercept", "Egger_intercept_SE", "Egger_intercept_P", "Steiger_direction", "Steiger_P",
  "Steiger_method", "mean_F", "min_F", "planned_family_size", "observed_primary_tests_in_family",
  "P_FDR_observed_reestimated", "P_Bonferroni_planned_family", "corrected_significance_FDR_observed",
  "corrected_significance_Bonferroni_planned", "planned_family_total", "P_FDR_all_estimable",
  "P_Bonferroni_planned_20", "corrected_significance_FDR_all_estimable",
  "corrected_significance_Bonferroni_planned_20", "beta_status", "exclusion_reason", "beta_source", "analysis_set"
)
setcolorder(results, intersect(ordered_columns, names(results)))
fwrite(results, output_file, sep = "\t", na = "NA")
if (length(iv_rows)) fwrite(rbindlist(iv_rows, fill = TRUE), iv_output, sep = "\t", na = "NA")
if (length(harm_rows)) fwrite(rbindlist(harm_rows, fill = TRUE, idcol = "analysis_id"), harm_output, sep = "\t", na = "NA")
writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo.txt"))

message("Primary panel targets planned: ", planned_family_size)
message("Re-estimated primary tests: ", nrow(results[method_role == "primary" & beta_status == "reestimated"]))
message("Output: ", output_file)
