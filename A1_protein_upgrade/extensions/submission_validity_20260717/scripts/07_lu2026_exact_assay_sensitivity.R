#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
  library(yaml)
})

set.seed(20260717)
n_boot <- 10000L

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = TRUE)
extension_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
upgrade_root <- normalizePath(file.path(extension_root, "..", ".."), winslash = "/", mustWork = TRUE)

resources <- read_yaml(file.path(upgrade_root, "config", "resources.yml"))
outcomes_config <- read_yaml(file.path(upgrade_root, "config", "outcomes.yml"))

candidate_file <- file.path(extension_root, "data_processed", "Lu2026_exact_assay_pqtl_candidates.tsv")
alpha_file <- file.path(extension_root, "tables", "Lu2026_APOE_variant_to_protein_alpha.tsv")
selection_file <- file.path(extension_root, "config", "lu2026_exact_assay_download_selection.tsv")
total_file <- file.path(upgrade_root, "tables", "APOE_variant_total_effects_primary_analysis.tsv")
table_dir <- file.path(extension_root, "tables")
processed_dir <- file.path(extension_root, "data_processed")
log_dir <- file.path(extension_root, "logs", "lu2026_exact_assay")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(candidate_file), file.exists(alpha_file), file.exists(selection_file), file.exists(total_file))
plink_bin <- resources$paths$plink_binary
ld_ref <- resources$paths$ld_reference
stopifnot(file.exists(plink_bin), file.exists(paste0(ld_ref, ".bed")))

outcome_specs <- list(
  AD = outcomes_config$outcomes$AD$local_file,
  any_AMD = outcomes_config$outcomes$any_AMD$local_file,
  dry_AMD = outcomes_config$outcomes$dry_AMD$local_file,
  wet_AMD = outcomes_config$outcomes$wet_AMD$local_file
)

selection <- fread(selection_file, na.strings = c("NA", ""))
selection <- selection[selection_status == "selected_for_download"]
if (nrow(selection) != 5L || uniqueN(selection$gene_symbol) != 5L) {
  stop("Lu 2026 feasibility gate assertion failed: expected exactly five selected genes.")
}
planned_genes <- sort(selection$gene_symbol)
planned_beta_family <- length(planned_genes)
planned_global_beta_family <- length(planned_genes) * length(outcome_specs)
planned_mediation_family <- length(planned_genes) * 2L * length(outcome_specs)

candidates <- fread(candidate_file, na.strings = c("NA", ""))
candidates <- candidates[gene_symbol %chin% planned_genes]
if (!all(planned_genes %chin% candidates$gene_symbol)) {
  stop("At least one gate-passing Lu 2026 gene has no extracted pQTL candidate row.")
}
candidates <- candidates[rsid_mapping_status == "mapped" & !is.na(SNP)]
setorder(candidates, gene_symbol, P_value)
candidates <- candidates[, .SD[1], by = .(gene_symbol, SNP)]

plink_clump <- function(dat, gene, analysis_set) {
  if (!nrow(dat)) return(dat)
  prefix_name <- paste(gene, analysis_set, sep = "_")
  clump_input <- file.path(log_dir, paste0(prefix_name, "_clump_input.tsv"))
  prefix <- file.path(log_dir, paste0(prefix_name, "_plink"))
  fwrite(dat[, .(SNP, P = P_value)], clump_input, sep = "\t")
  status <- system2(
    plink_bin,
    args = c(
      "--bfile", ld_ref, "--clump", clump_input,
      "--clump-p1", "5e-8", "--clump-p2", "5e-8",
      "--clump-r2", "0.001", "--clump-kb", "10000", "--out", prefix
    ),
    stdout = paste0(prefix, ".stdout.txt"), stderr = paste0(prefix, ".stderr.txt")
  )
  clumped_file <- paste0(prefix, ".clumped")
  if (status != 0 || !file.exists(clumped_file)) {
    stop(sprintf("PLINK clumping failed for %s (%s).", gene, analysis_set))
  }
  clumped <- fread(clumped_file, fill = TRUE)
  if (!"SNP" %in% names(clumped)) return(dat[0])
  dat[SNP %chin% clumped$SNP]
}

format_exposure <- function(dat, gene) {
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

empty_result <- function(gene, outcome, analysis_set, reason, n_iv = 0L) {
  meta <- unique(selection[gene_symbol == gene, .(UniProt_ID, Olink_target_ID)])
  data.table(
    gene_symbol = gene, assay_target_ID = meta$Olink_target_ID[1], UniProt_ID = meta$UniProt_ID[1],
    outcome = outcome, analysis_set = analysis_set, method = "not_estimable", method_role = "primary",
    nsnp = n_iv, beta = NA_real_, SE = NA_real_, P_value = NA_real_, CI_lower = NA_real_,
    CI_upper = NA_real_, OR = NA_real_, OR_CI_lower = NA_real_, OR_CI_upper = NA_real_,
    Cochran_Q = NA_real_, Q_df = NA_real_, Q_P = NA_real_, Egger_intercept = NA_real_,
    Egger_intercept_SE = NA_real_, Egger_intercept_P = NA_real_, mean_F = NA_real_, min_F = NA_real_,
    planned_family_size_per_outcome = planned_beta_family,
    planned_global_family_size = planned_global_beta_family,
    P_FDR_observed_per_outcome = NA_real_, P_Bonferroni_planned_per_outcome = NA_real_,
    P_FDR_global = NA_real_, P_Bonferroni_global = NA_real_,
    beta_status = "not_estimable", exclusion_reason = reason,
    evidence_layer = "externally_defined_post_publication_candidate_sensitivity"
  )
}

iv_rows <- list()
clumped_sets <- list()
for (analysis_set in c("genome_wide_instruments", "cis_only")) {
  for (gene in planned_genes) {
    gene_candidates <- candidates[gene_symbol == gene]
    if (analysis_set == "cis_only") gene_candidates <- gene_candidates[cis_to_encoding_gene %in% TRUE]
    clumped <- plink_clump(gene_candidates, gene, analysis_set)
    if (nrow(clumped)) {
      clumped[, `:=`(
        F_statistic = (beta / SE)^2,
        analysis_set = analysis_set,
        clump_r2 = 0.001,
        clump_kb = 10000L
      )]
      iv_rows[[paste(gene, analysis_set, sep = "_")]] <- clumped
    }
    clumped_sets[[paste(gene, analysis_set, sep = "_")]] <- clumped
  }
}
if (!length(iv_rows)) stop("No Lu 2026 pQTL instrument remained after clumping.")

all_iv_keys <- unique(rbindlist(iv_rows, fill = TRUE)$SNP)
outcomes_loaded <- lapply(names(outcome_specs), function(outcome_name) {
  path <- outcome_specs[[outcome_name]]
  stopifnot(file.exists(path))
  dat <- fread(path, select = c("SNP", "CHR", "BP", "A1", "A2", "FREQ", "BETA", "SE", "P", "N"), showProgress = FALSE)
  dat <- dat[SNP %chin% all_iv_keys]
  dat[, `:=`(outcome = outcome_name, id.outcome = outcome_name)]
  dat
})
names(outcomes_loaded) <- names(outcome_specs)

result_rows <- list()
harm_rows <- list()
for (analysis_set in c("genome_wide_instruments", "cis_only")) {
  for (gene in planned_genes) {
    clumped <- clumped_sets[[paste(gene, analysis_set, sep = "_")]]
    if (!nrow(clumped)) {
      for (outcome_name in names(outcome_specs)) {
        result_rows[[length(result_rows) + 1L]] <- empty_result(
          gene, outcome_name, analysis_set, "No mapped genome-wide significant pQTL remained after the prespecified cis restriction and/or clumping."
        )
      }
      next
    }
    exposure <- format_exposure(clumped, gene)
    for (outcome_name in names(outcome_specs)) {
      outcome_raw <- outcomes_loaded[[outcome_name]][SNP %chin% exposure$SNP]
      if (!nrow(outcome_raw)) {
        result_rows[[length(result_rows) + 1L]] <- empty_result(
          gene, outcome_name, analysis_set, "No clumped pQTL was present in the outcome GWAS.", nrow(clumped)
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
      if (!nrow(harmonized)) {
        result_rows[[length(result_rows) + 1L]] <- empty_result(
          gene, outcome_name, analysis_set, "No pQTL remained after allele harmonization.", nrow(clumped)
        )
        next
      }
      harmonized$F_statistic <- (harmonized$beta.exposure / harmonized$se.exposure)^2
      harm_rows[[paste(gene, outcome_name, analysis_set, sep = "_")]] <- as.data.table(harmonized)
      methods <- if (nrow(harmonized) == 1L) "mr_wald_ratio" else c("mr_ivw", "mr_weighted_median", "mr_egger_regression")
      estimates <- as.data.table(mr(harmonized, method_list = methods))
      heterogeneity <- tryCatch(as.data.table(mr_heterogeneity(harmonized)), error = function(e) data.table())
      pleiotropy <- tryCatch(as.data.table(mr_pleiotropy_test(harmonized)), error = function(e) data.table())
      if (nrow(heterogeneity)) {
        estimates <- merge(estimates, heterogeneity[, .(method, Q, Q_df, Q_pval)], by = "method", all.x = TRUE)
      } else {
        estimates[, `:=`(Q = NA_real_, Q_df = NA_real_, Q_pval = NA_real_)]
      }
      estimates[, `:=`(
        gene_symbol = gene,
        outcome = outcome_name,
        analysis_set = analysis_set,
        method_role = fifelse(method %chin% c("Inverse variance weighted", "Wald ratio"), "primary", "sensitivity"),
        CI_lower = b - 1.96 * se,
        CI_upper = b + 1.96 * se,
        OR = exp(b),
        OR_CI_lower = exp(b - 1.96 * se),
        OR_CI_upper = exp(b + 1.96 * se),
        Egger_intercept = if (nrow(pleiotropy)) pleiotropy$egger_intercept[1] else NA_real_,
        Egger_intercept_SE = if (nrow(pleiotropy)) pleiotropy$se[1] else NA_real_,
        Egger_intercept_P = if (nrow(pleiotropy)) pleiotropy$pval[1] else NA_real_,
        mean_F = mean(harmonized$F_statistic),
        min_F = min(harmonized$F_statistic),
        planned_family_size_per_outcome = planned_beta_family,
        planned_global_family_size = planned_global_beta_family,
        P_FDR_observed_per_outcome = NA_real_,
        P_Bonferroni_planned_per_outcome = NA_real_,
        P_FDR_global = NA_real_,
        P_Bonferroni_global = NA_real_,
        beta_status = "reestimated",
        exclusion_reason = "NA",
        evidence_layer = "externally_defined_post_publication_candidate_sensitivity"
      )]
      setnames(estimates, c("b", "se", "pval", "Q", "Q_pval"), c("beta", "SE", "P_value", "Cochran_Q", "Q_P"), skip_absent = TRUE)
      meta <- unique(selection[gene_symbol == gene, .(UniProt_ID, Olink_target_ID)])
      estimates[, `:=`(assay_target_ID = meta$Olink_target_ID[1], UniProt_ID = meta$UniProt_ID[1])]
      result_rows[[length(result_rows) + 1L]] <- estimates
    }
  }
}

results <- rbindlist(result_rows, fill = TRUE)
for (analysis_set_i in unique(results$analysis_set)) {
  primary_all <- which(results$analysis_set == analysis_set_i & results$method_role == "primary" & results$beta_status == "reestimated")
  if (length(primary_all)) {
    results[primary_all, P_FDR_global := p.adjust(P_value, method = "BH")]
    results[primary_all, P_Bonferroni_global := pmin(1, P_value * planned_global_beta_family)]
  }
  for (outcome_i in names(outcome_specs)) {
    idx <- which(results$analysis_set == analysis_set_i & results$outcome == outcome_i &
                   results$method_role == "primary" & results$beta_status == "reestimated")
    if (length(idx)) {
      results[idx, P_FDR_observed_per_outcome := p.adjust(P_value, method = "BH")]
      results[idx, P_Bonferroni_planned_per_outcome := pmin(1, P_value * planned_beta_family)]
    }
  }
}
setorder(results, analysis_set, outcome, gene_symbol, method_role, method)
fwrite(results, file.path(table_dir, "Lu2026_exact_assay_beta_sensitivity.tsv"), sep = "\t", na = "NA")
fwrite(rbindlist(iv_rows, fill = TRUE), file.path(processed_dir, "Lu2026_exact_assay_clumped_IVs.tsv"), sep = "\t", na = "NA")
if (length(harm_rows)) {
  fwrite(rbindlist(harm_rows, fill = TRUE, idcol = "analysis_id"),
         file.path(processed_dir, "Lu2026_exact_assay_harmonized_data.tsv"), sep = "\t", na = "NA")
}

alpha <- fread(alpha_file, na.strings = c("NA", ""))
alpha_use <- unique(alpha[
  gene_symbol %chin% planned_genes & availability_status == "direct_variant_available",
  .(gene_symbol, variant, alpha = as.numeric(beta), SE_alpha = as.numeric(SE), alpha_P = as.numeric(P_value),
    alpha_effect_allele = requested_effect_allele, alpha_assay_target_ID = assay_target_ID,
    alpha_source, alpha_sample_size = sample_size)
], by = c("gene_symbol", "variant"))

totals <- fread(total_file, na.strings = c("NA", ""))
total_use <- totals[availability_status == "direct_variant_available", .(
  variant, outcome, total_effect = as.numeric(beta), SE_total = as.numeric(SE),
  total_effect_P = as.numeric(P_value), total_effect_source = outcome_GWAS
)]

beta_use <- results[
  method_role == "primary" & beta_status == "reestimated",
  .(gene_symbol, assay_target_ID, outcome, analysis_set, beta = as.numeric(beta), SE_beta = as.numeric(SE),
    beta_P = as.numeric(P_value), beta_method = method, beta_nsnp = nsnp,
    beta_Q_P = Q_P, beta_Egger_intercept_P = Egger_intercept_P,
    beta_P_FDR_per_outcome = P_FDR_observed_per_outcome,
    beta_P_Bonferroni_per_outcome = P_Bonferroni_planned_per_outcome,
    beta_P_FDR_global = P_FDR_global, beta_P_Bonferroni_global = P_Bonferroni_global)
]

med <- merge(alpha_use, beta_use, by = "gene_symbol", allow.cartesian = TRUE)
med <- med[alpha_assay_target_ID == assay_target_ID]
med <- merge(med, total_use, by = c("variant", "outcome"), all.x = TRUE)
if (anyNA(med$total_effect)) stop("APOE total effect missing from Lu 2026 mediation rows.")

med[, `:=`(
  indirect_effect = alpha * beta,
  SE_indirect_delta = sqrt(beta^2 * SE_alpha^2 + alpha^2 * SE_beta^2),
  extension_variant_role = fifelse(
    variant == "rs429358", "source_aligned_APOE4_analysis", "cross_isoform_rs7412_secondary_analysis"
  )
)]
med[, indirect_P_delta := 2 * pnorm(-abs(indirect_effect / SE_indirect_delta))]
med[, mediated_proportion := indirect_effect / total_effect]
med[, SE_mediated_proportion_delta := sqrt(
  (beta / total_effect)^2 * SE_alpha^2 +
    (alpha / total_effect)^2 * SE_beta^2 +
    (alpha * beta / total_effect^2)^2 * SE_total^2
)]
med[, direction_classification := fifelse(
  sign(indirect_effect) == sign(total_effect), "concordant_mediation", "opposing_or_suppressing_mediation"
)]

for (analysis_set_i in unique(med$analysis_set)) {
  idx <- which(med$analysis_set == analysis_set_i)
  med[idx, `:=`(
    planned_mediation_family_size = planned_mediation_family,
    observed_mediation_tests = length(idx),
    indirect_P_FDR_observed = p.adjust(indirect_P_delta, method = "BH"),
    indirect_P_Bonferroni_planned = pmin(1, indirect_P_delta * planned_mediation_family)
  )]
}
for (i in seq_len(nrow(med))) {
  alpha_draw <- rnorm(n_boot, med$alpha[i], med$SE_alpha[i])
  beta_draw <- rnorm(n_boot, med$beta[i], med$SE_beta[i])
  total_draw <- rnorm(n_boot, med$total_effect[i], med$SE_total[i])
  indirect_draw <- alpha_draw * beta_draw
  proportion_draw <- indirect_draw / total_draw
  med$indirect_CI_lower_bootstrap[i] <- quantile(indirect_draw, 0.025, na.rm = TRUE)
  med$indirect_CI_upper_bootstrap[i] <- quantile(indirect_draw, 0.975, na.rm = TRUE)
  med$mediated_proportion_CI_lower_bootstrap[i] <- quantile(proportion_draw, 0.025, na.rm = TRUE)
  med$mediated_proportion_CI_upper_bootstrap[i] <- quantile(proportion_draw, 0.975, na.rm = TRUE)
}
med[, `:=`(
  bootstrap_n = n_boot,
  evidence_layer = "externally_defined_post_publication_candidate_sensitivity",
  interpretation_boundary = paste(
    "Lu et al. 2026 candidates were defined by observational cross-dataset APOE4-associated proteomics, not proteome-wide MR.",
    "This exact-assay extension is sensitivity evidence and does not modify the prespecified primary panel."
  )
)]
setorder(med, analysis_set, variant, outcome, gene_symbol)
fwrite(med, file.path(table_dir, "Lu2026_exact_assay_two_step_mediation_sensitivity.tsv"), sep = "\t", na = "NA")

summary <- rbindlist(lapply(c("genome_wide_instruments", "cis_only"), function(set_i) {
  beta_set <- results[analysis_set == set_i & method_role == "primary" & beta_status == "reestimated"]
  med_set <- med[analysis_set == set_i]
  data.table(
    analysis_set = set_i,
    planned_genes = planned_beta_family,
    estimable_beta_tests = nrow(beta_set),
    planned_beta_tests = planned_global_beta_family,
    beta_FDR_signals = sum(beta_set$P_FDR_global < 0.05, na.rm = TRUE),
    beta_Bonferroni_signals = sum(beta_set$P_Bonferroni_global < 0.05, na.rm = TRUE),
    direct_alpha_rows = nrow(alpha_use),
    mediation_tests = nrow(med_set),
    planned_mediation_tests = planned_mediation_family,
    mediation_FDR_signals = sum(med_set$indirect_P_FDR_observed < 0.05, na.rm = TRUE),
    mediation_Bonferroni_signals = sum(med_set$indirect_P_Bonferroni_planned < 0.05, na.rm = TRUE),
    primary_panel_modified = FALSE
  )
}))
fwrite(summary, file.path(table_dir, "Lu2026_exact_assay_run_summary.tsv"), sep = "\t", na = "NA")
writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo.txt"))

message("Lu 2026 exact-assay sensitivity complete: ", nrow(results), " beta rows and ", nrow(med), " mediation rows.")
