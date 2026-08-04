#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
  library(yaml)
})

set.seed(20260719)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = TRUE)
mode_arg <- grep("^--mode=", commandArgs(trailingOnly = TRUE), value = TRUE)
analysis_mode <- if (length(mode_arg)) sub("^--mode=", "", mode_arg[1]) else "genome_wide"
if (!analysis_mode %in% c("genome_wide", "cis_only")) stop("--mode must be genome_wide or cis_only")
suffix <- if (analysis_mode == "cis_only") "_cis_only" else ""
trailing_args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), trailing_args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[1]) else default
}
output_prefix <- arg_value("output-prefix", "decode_smp")
normalization_label <- arg_value("normalization-label", "smp")
inference_role <- arg_value("inference-role", "independent_same_platform_sensitivity")
planned_assays <- as.integer(arg_value("planned-assays", "9"))
root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
project_root <- normalizePath(file.path(root, ".."), winslash = "/", mustWork = TRUE)
resources <- read_yaml(file.path(project_root, "A1_protein_upgrade", "config", "resources.yml"))
outcomes_config <- read_yaml(file.path(project_root, "A1_protein_upgrade", "config", "outcomes.yml"))

candidate_file <- arg_value("candidate-file", file.path(root, "data_processed", "decode_smp_pqtl_candidates_for_clumping.tsv"))
result_file <- file.path(root, "tables", paste0(output_prefix, "_beta_results", suffix, ".tsv"))
iv_file <- file.path(root, "data_processed", paste0(output_prefix, "_clumped_instruments", suffix, ".tsv"))
harmonized_file <- file.path(root, "data_processed", paste0(output_prefix, "_harmonized_data", suffix, ".tsv.gz"))
clump_summary_file <- file.path(root, "tables", paste0(output_prefix, "_clumping_summary", suffix, ".tsv"))
log_dir <- file.path(root, "logs", paste0(output_prefix, "_beta_reestimation", suffix))
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

plink_bin <- resources$paths$plink_binary
ld_ref <- resources$paths$ld_reference
stopifnot(file.exists(candidate_file), file.exists(plink_bin), file.exists(paste0(ld_ref, ".bed")))

outcome_specs <- list(
  AD = outcomes_config$outcomes$AD$local_file,
  any_AMD = outcomes_config$outcomes$any_AMD$local_file,
  dry_AMD = outcomes_config$outcomes$dry_AMD$local_file,
  wet_AMD = outcomes_config$outcomes$wet_AMD$local_file
)
stopifnot(all(vapply(outcome_specs, file.exists, logical(1))))

planned_outcomes <- 4L
planned_global_beta_family <- planned_assays * planned_outcomes

all_candidates <- fread(candidate_file, na.strings = c("NA", ""), showProgress = FALSE)
candidates <- copy(all_candidates)
required_columns <- c("gene_symbol", "assay_target_ID", "UniProt_ID", "SNP", "beta", "SE", "P_value",
                      "effect_allele", "other_allele", "sample_size", "imputation_MAF",
                      "chromosome", "position_hg19_from_liftover", "cis_to_encoding_gene")
stopifnot(all(required_columns %in% names(candidates)))
stopifnot(uniqueN(candidates$assay_target_ID) == planned_assays)
assay_metadata <- unique(all_candidates[, .(gene_symbol, assay_target_ID, UniProt_ID)])
stopifnot(nrow(assay_metadata) == planned_assays)
if (analysis_mode == "cis_only") candidates <- candidates[cis_to_encoding_gene %in% TRUE]
setorder(candidates, assay_target_ID, SNP, P_value)
candidates <- candidates[, .SD[1], by = .(assay_target_ID, SNP)]

safe_id <- function(value) gsub("[^A-Za-z0-9_.-]", "_", value)

plink_clump <- function(dat, assay_id) {
  stem <- safe_id(assay_id)
  input <- file.path(log_dir, paste0(stem, "_clump_input.tsv"))
  prefix <- file.path(log_dir, paste0(stem, "_plink"))
  fwrite(dat[, .(SNP, P = P_value)], input, sep = "\t")
  args <- c(
    "--bfile", ld_ref,
    "--clump", input,
    "--clump-p1", "5e-8",
    "--clump-p2", "5e-8",
    "--clump-r2", "0.001",
    "--clump-kb", "10000",
    "--out", prefix
  )
  status <- system2(plink_bin, args = args, stdout = paste0(prefix, ".stdout.txt"),
                    stderr = paste0(prefix, ".stderr.txt"))
  clumped_path <- paste0(prefix, ".clumped")
  if (status != 0 || !file.exists(clumped_path)) {
    stop(sprintf("PLINK clumping failed for assay %s", assay_id))
  }
  clumped <- fread(clumped_path, fill = TRUE)
  if (!"SNP" %in% names(clumped)) return(dat[0])
  dat[SNP %in% clumped$SNP]
}

format_exposure <- function(dat, assay_id) {
  exposure <- format_data(
    as.data.frame(dat), type = "exposure",
    snp_col = "SNP", beta_col = "beta", se_col = "SE",
    effect_allele_col = "effect_allele", other_allele_col = "other_allele",
    pval_col = "P_value", samplesize_col = "sample_size",
    chr_col = "chromosome", pos_col = "position_hg19_from_liftover"
  )
  # deCODE provides ImpMAF, not effect-allele frequency. It is deliberately not passed as EAF.
  exposure$exposure <- assay_id
  exposure$id.exposure <- assay_id
  exposure
}

load_outcome <- function(path, label, keys) {
  dat <- fread(path, select = c("SNP", "CHR", "BP", "A1", "A2", "FREQ", "BETA", "SE", "P", "N"),
               showProgress = FALSE)
  dat <- dat[SNP %in% keys]
  dat[, `:=`(outcome = label, id.outcome = label)]
  dat
}

empty_result <- function(meta, outcome, status, reason, n_clumped = 0L, n_harmonized = 0L) {
  data.table(
    gene_symbol = meta$gene_symbol[1], assay_target_ID = meta$assay_target_ID[1], UniProt_ID = meta$UniProt_ID[1],
    normalization = normalization_label, outcome = outcome, method = "not_estimable", method_role = "primary",
    nsnp = n_harmonized, n_clumped = n_clumped, beta = NA_real_, SE = NA_real_, P_value = NA_real_,
    CI_lower = NA_real_, CI_upper = NA_real_, OR = NA_real_, OR_CI_lower = NA_real_, OR_CI_upper = NA_real_,
    Cochran_Q = NA_real_, Q_df = NA_real_, Q_P = NA_real_, Egger_intercept = NA_real_,
    Egger_intercept_SE = NA_real_, Egger_intercept_P = NA_real_, mean_F = NA_real_, min_F = NA_real_,
    beta_status = status, exclusion_reason = reason,
    beta_source = paste0("deCODE_SomaScan_", normalization_label, "_", analysis_mode, "_to_frozen_A1_outcomes"),
    instrument_scope = analysis_mode,
    inference_role = inference_role,
    planned_tests_per_outcome = planned_assays, planned_global_beta_tests = planned_global_beta_family,
    P_FDR_within_outcome = NA_real_, P_Bonferroni_within_outcome = NA_real_,
    P_FDR_global_36 = NA_real_, P_Bonferroni_global_36 = NA_real_
  )
}

iv_rows <- list()
clumped_by_assay <- list()
clump_summary <- list()
for (assay_id in sort(unique(candidates$assay_target_ID))) {
  message("Clumping ", assay_id)
  assay_candidates <- candidates[assay_target_ID == assay_id]
  clumped <- plink_clump(assay_candidates, assay_id)
  if (nrow(clumped)) {
    clumped[, `:=`(
      F_statistic = (beta / SE)^2,
      clump_r2 = 0.001,
      clump_kb = 10000L,
      exposure_EAF_status = "unavailable_ImpMAF_not_substituted"
    )]
    iv_rows[[assay_id]] <- clumped
  }
  clumped_by_assay[[assay_id]] <- clumped
  clump_summary[[assay_id]] <- data.table(
    gene_symbol = assay_candidates$gene_symbol[1], assay_target_ID = assay_id,
    n_candidates_in_EUR_reference = nrow(assay_candidates), n_clumped_instruments = nrow(clumped),
    n_clumped_cis = if (nrow(clumped)) sum(clumped$cis_to_encoding_gene %in% TRUE) else 0L,
    minimum_F = if (nrow(clumped)) min((clumped$beta / clumped$SE)^2) else NA_real_,
    median_F = if (nrow(clumped)) median((clumped$beta / clumped$SE)^2) else NA_real_,
    n_ImpMAF_lt_0_01 = if (nrow(clumped)) sum(clumped$imputation_MAF < 0.01, na.rm = TRUE) else 0L
  )
}

missing_assays <- setdiff(assay_metadata$assay_target_ID, names(clumped_by_assay))
if (length(missing_assays)) {
  for (assay_id in missing_assays) {
    meta <- assay_metadata[assay_target_ID == assay_id][1]
    clumped_by_assay[[assay_id]] <- candidates[0]
    clump_summary[[assay_id]] <- data.table(
      gene_symbol = meta$gene_symbol, assay_target_ID = assay_id,
      n_candidates_in_EUR_reference = 0L, n_clumped_instruments = 0L, n_clumped_cis = 0L,
      minimum_F = NA_real_, median_F = NA_real_, n_ImpMAF_lt_0_01 = 0L
    )
  }
}

all_ivs <- if (length(iv_rows)) rbindlist(iv_rows, fill = TRUE) else data.table()
if (!nrow(all_ivs)) stop("No instruments remained after clumping")
all_keys <- unique(all_ivs$SNP)
message("Loading frozen outcomes for ", length(all_keys), " unique clumped SNPs")
outcomes <- lapply(names(outcome_specs), function(name) load_outcome(outcome_specs[[name]], name, all_keys))
names(outcomes) <- names(outcome_specs)

result_rows <- list()
harmonized_rows <- list()
for (assay_id in names(clumped_by_assay)) {
  clumped <- clumped_by_assay[[assay_id]]
  meta <- assay_metadata[assay_target_ID == assay_id][1]
  if (!nrow(clumped)) {
    for (outcome_name in names(outcomes)) {
      result_rows[[length(result_rows) + 1L]] <- empty_result(meta, outcome_name, "not_estimable",
                                                              "No variants remained after PLINK clumping")
    }
    next
  }
  exposure <- format_exposure(clumped, assay_id)
  for (outcome_name in names(outcomes)) {
    outcome_raw <- outcomes[[outcome_name]][SNP %in% exposure$SNP]
    if (!nrow(outcome_raw)) {
      result_rows[[length(result_rows) + 1L]] <- empty_result(meta, outcome_name, "not_estimable",
                                                              "No clumped instrument was present in the outcome GWAS", nrow(clumped))
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
    harmonized_all <- harmonise_data(exposure, outcome, action = 2)
    harmonized_all$assay_target_ID <- assay_id
    harmonized_all$gene_symbol <- meta$gene_symbol[1]
    harmonized_all$normalization <- normalization_label
    harmonized_all$exposure_EAF_status <- "unavailable_ImpMAF_not_substituted"
    harmonized_rows[[paste(assay_id, outcome_name, sep = "_")]] <- as.data.table(harmonized_all)
    harmonized <- harmonized_all[harmonized_all$mr_keep %in% TRUE, ]
    if (!nrow(harmonized)) {
      result_rows[[length(result_rows) + 1L]] <- empty_result(meta, outcome_name, "not_estimable",
                                                              "No variants remained after pairwise allele harmonization", nrow(clumped), 0L)
      next
    }
    harmonized$F_statistic <- (harmonized$beta.exposure / harmonized$se.exposure)^2
    methods <- if (nrow(harmonized) == 1L) "mr_wald_ratio" else c("mr_ivw", "mr_weighted_median", "mr_egger_regression")
    estimates <- as.data.table(mr(harmonized, method_list = methods))
    if (!nrow(estimates)) {
      result_rows[[length(result_rows) + 1L]] <- empty_result(meta, outcome_name, "not_estimable",
                                                              "TwoSampleMR returned no estimate", nrow(clumped), nrow(harmonized))
      next
    }
    heterogeneity <- tryCatch(as.data.table(mr_heterogeneity(harmonized)), error = function(e) data.table())
    pleiotropy <- tryCatch(as.data.table(mr_pleiotropy_test(harmonized)), error = function(e) data.table())
    if (nrow(heterogeneity)) {
      estimates <- merge(estimates, heterogeneity[, .(method, Q, Q_df, Q_pval)], by = "method", all.x = TRUE)
    } else {
      estimates[, `:=`(Q = NA_real_, Q_df = NA_real_, Q_pval = NA_real_)]
    }
    estimates[, `:=`(
      gene_symbol = meta$gene_symbol[1], assay_target_ID = assay_id, UniProt_ID = meta$UniProt_ID[1],
      normalization = normalization_label, outcome = outcome_name,
      method_role = fifelse(method %in% c("Inverse variance weighted", "Wald ratio"), "primary", "sensitivity"),
      n_clumped = nrow(clumped), CI_lower = b - 1.96 * se, CI_upper = b + 1.96 * se,
      OR = exp(b), OR_CI_lower = exp(b - 1.96 * se), OR_CI_upper = exp(b + 1.96 * se),
      Egger_intercept = if (nrow(pleiotropy)) pleiotropy$egger_intercept[1] else NA_real_,
      Egger_intercept_SE = if (nrow(pleiotropy)) pleiotropy$se[1] else NA_real_,
      Egger_intercept_P = if (nrow(pleiotropy)) pleiotropy$pval[1] else NA_real_,
      mean_F = mean(harmonized$F_statistic), min_F = min(harmonized$F_statistic),
      beta_status = "reestimated", exclusion_reason = "NA",
      beta_source = paste0("deCODE_SomaScan_", normalization_label, "_", analysis_mode, "_to_frozen_A1_outcomes"),
      instrument_scope = analysis_mode,
      inference_role = inference_role,
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
estimable <- which(results$beta_status == "reestimated" & is.finite(results$beta) & is.finite(results$SE) & results$SE > 0)
results[, `:=`(minus_log10_P = NA_real_, P_value_display = NA_character_)]
if (length(estimable)) {
  log_p <- log(2) + pnorm(abs(results$beta[estimable] / results$SE[estimable]), lower.tail = FALSE, log.p = TRUE)
  results[estimable, minus_log10_P := -log_p / log(10)]
  results[estimable, P_value := ifelse(log_p < log(.Machine$double.xmin), 0, exp(log_p))]
  results[estimable, P_value_display := ifelse(
    minus_log10_P > 300, "<1e-300", formatC(P_value, format = "e", digits = 3)
  )]
}
primary_idx <- which(results$method_role == "primary" & results$beta_status == "reestimated")
for (outcome_name in names(outcome_specs)) {
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
results[, `:=`(
  P_FDR_global_planned = P_FDR_global_36,
  P_Bonferroni_global_planned = P_Bonferroni_global_36
)]

setorder(results, assay_target_ID, outcome, method_role, method)
fwrite(results, result_file, sep = "\t", na = "NA")
fwrite(all_ivs, iv_file, sep = "\t", na = "NA")
fwrite(rbindlist(clump_summary, fill = TRUE), clump_summary_file, sep = "\t", na = "NA")
if (length(harmonized_rows)) {
  fwrite(rbindlist(harmonized_rows, fill = TRUE, idcol = "analysis_id"), harmonized_file, sep = "\t", na = "NA", compress = "gzip")
}
writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo.txt"))
writeLines(c(
  paste0("Planned assays per outcome: ", planned_assays, "."),
  paste0("Planned global beta rows: ", planned_global_beta_family, "."),
  "ImpMAF was not substituted for exposure effect-allele frequency.",
  "Pairwise harmonization used TwoSampleMR action=2; unresolved palindromic variants were excluded.",
  paste0("Instrument scope: ", analysis_mode, "."),
  paste0("Inference role: ", inference_role, ".")
), file.path(log_dir, "analysis_rules.txt"))

message("Clumped assay-SNP instruments: ", nrow(all_ivs))
message("Primary beta estimates: ", length(primary_idx))
message("Output: ", result_file)
