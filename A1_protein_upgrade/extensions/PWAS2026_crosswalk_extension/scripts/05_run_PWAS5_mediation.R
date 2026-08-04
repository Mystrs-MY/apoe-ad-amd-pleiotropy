#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

set.seed(20260714)
n_boot <- 10000L

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = TRUE)
extension_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
upgrade_root <- normalizePath(file.path(extension_root, "..", ".."), winslash = "/", mustWork = TRUE)

alpha <- fread(file.path(extension_root, "tables", "PWAS5_APOE_alpha.tsv"), na.strings = c("NA", ""))
members <- fread(file.path(extension_root, "config", "PWAS5_frozen_members.tsv"), na.strings = c("NA", ""))
totals <- fread(file.path(upgrade_root, "tables", "APOE_variant_total_effects_current_A1.tsv"), na.strings = c("NA", ""))
planned_family <- 5L * 2L * 4L

total_use <- totals[availability_status == "direct_variant_available",
  .(variant, outcome, total_effect = as.numeric(beta), SE_total = as.numeric(SE),
    total_effect_P = as.numeric(P_value), total_effect_allele = requested_effect_allele,
    total_effect_source = outcome_GWAS)]

run_set <- function(label, beta_file, beta_status_required, strict_mapping = FALSE) {
  beta <- fread(beta_file, na.strings = c("NA", ""), fill = TRUE)
  beta_primary <- beta[method_role == "primary"]
  beta_primary <- beta_primary[order(gene_symbol, outcome, P_value)]
  beta_primary <- beta_primary[, .SD[1], by = .(gene_symbol, outcome)]

  grid <- CJ(
    gene_symbol = members$gene_symbol,
    variant = c("rs429358", "rs7412"),
    outcome = c("AD", "any_AMD", "dry_AMD", "wet_AMD"),
    unique = TRUE
  )
  mapping <- members[, .(
    gene_symbol,
    UniProt_ID,
    assay_target_ID = UKB_PPP_OID,
    mapping_confidence,
    strict_mapping_eligible = tolower(as.character(strict_mapping_eligible)) == "true"
  )]
  out <- merge(grid, mapping, by = "gene_symbol", all.x = TRUE)
  alpha_use <- alpha[, .(
    gene_symbol, variant,
    alpha = as.numeric(beta), SE_alpha = as.numeric(SE), alpha_P = as.numeric(P_value),
    alpha_effect_allele = requested_effect_allele,
    alpha_assay_target_ID = assay_target_ID,
    alpha_sample_size = as.numeric(sample_size),
    alpha_available = availability_status == "direct_variant_available" & variant_match_count == 1
  )]
  beta_use <- beta_primary[, .(
    gene_symbol, outcome,
    beta = as.numeric(beta), SE_beta = as.numeric(SE), beta_P = as.numeric(P_value),
    beta_method = method, beta_nsnp = nsnp, beta_mean_F = mean_F, beta_min_F = min_F,
    beta_Q_P = Q_P, beta_Egger_intercept_P = Egger_intercept_P,
    beta_P_FDR_per_outcome = P_FDR_observed_reestimated,
    beta_P_Bonferroni_planned_5 = P_Bonferroni_planned_family,
    beta_P_FDR_all_estimable = P_FDR_all_estimable,
    beta_P_Bonferroni_planned_20 = P_Bonferroni_planned_20,
    beta_assay_target_ID = assay_target_ID,
    beta_status
  )]
  out <- merge(out, alpha_use, by = c("gene_symbol", "variant"), all.x = TRUE)
  out <- merge(out, beta_use, by = c("gene_symbol", "outcome"), all.x = TRUE)
  out <- merge(out, total_use, by = c("variant", "outcome"), all.x = TRUE)
  out[, same_assay_alpha_beta := !is.na(alpha_assay_target_ID) & !is.na(beta_assay_target_ID) &
        alpha_assay_target_ID == beta_assay_target_ID]
  out[, mapping_pass := if (strict_mapping) strict_mapping_eligible else TRUE]
  out[, mediation_status := fifelse(
    !mapping_pass, "not_estimable_mapping_not_high_confidence",
    fifelse(!(alpha_available %in% TRUE), "not_estimable_alpha_missing_or_nonunique",
      fifelse(is.na(beta) | !(beta_status %in% beta_status_required), "not_estimable_beta_missing",
        fifelse(!(same_assay_alpha_beta %in% TRUE), "not_estimable_assay_mismatch",
          fifelse(is.na(total_effect), "not_estimable_total_effect_missing", "estimable"))
      )
    )
  )]
  out[, `:=`(
    indirect_effect = NA_real_, SE_indirect_delta = NA_real_, indirect_P_delta = NA_real_,
    mediated_proportion = NA_real_, SE_mediated_proportion_delta = NA_real_,
    indirect_CI_lower_bootstrap = NA_real_, indirect_CI_upper_bootstrap = NA_real_,
    mediated_proportion_CI_lower_bootstrap = NA_real_, mediated_proportion_CI_upper_bootstrap = NA_real_,
    total_effect_stable = !is.na(total_effect) & !is.na(SE_total) & abs(total_effect / SE_total) >= 2,
    direction_classification = "not_estimable",
    planned_mediation_family_size = planned_family,
    actual_estimable_family_size = sum(mediation_status == "estimable"),
    indirect_P_FDR_actual = NA_real_, indirect_P_Bonferroni_planned_40 = NA_real_,
    corrected_significance_FDR_actual = NA,
    corrected_significance_Bonferroni_planned_40 = NA,
    analysis_set = label, bootstrap_n = n_boot
  )]

  estimable <- which(out$mediation_status == "estimable")
  for (i in estimable) {
    out$indirect_effect[i] <- out$alpha[i] * out$beta[i]
    out$SE_indirect_delta[i] <- sqrt(out$beta[i]^2 * out$SE_alpha[i]^2 + out$alpha[i]^2 * out$SE_beta[i]^2)
    out$indirect_P_delta[i] <- 2 * pnorm(-abs(out$indirect_effect[i] / out$SE_indirect_delta[i]))
    out$mediated_proportion[i] <- out$indirect_effect[i] / out$total_effect[i]
    out$SE_mediated_proportion_delta[i] <- sqrt(
      (out$beta[i] / out$total_effect[i])^2 * out$SE_alpha[i]^2 +
      (out$alpha[i] / out$total_effect[i])^2 * out$SE_beta[i]^2 +
      (out$alpha[i] * out$beta[i] / out$total_effect[i]^2)^2 * out$SE_total[i]^2
    )
    alpha_draw <- rnorm(n_boot, out$alpha[i], out$SE_alpha[i])
    beta_draw <- rnorm(n_boot, out$beta[i], out$SE_beta[i])
    total_draw <- rnorm(n_boot, out$total_effect[i], out$SE_total[i])
    indirect_draw <- alpha_draw * beta_draw
    proportion_draw <- indirect_draw / total_draw
    out$indirect_CI_lower_bootstrap[i] <- quantile(indirect_draw, 0.025, na.rm = TRUE)
    out$indirect_CI_upper_bootstrap[i] <- quantile(indirect_draw, 0.975, na.rm = TRUE)
    out$mediated_proportion_CI_lower_bootstrap[i] <- quantile(proportion_draw, 0.025, na.rm = TRUE)
    out$mediated_proportion_CI_upper_bootstrap[i] <- quantile(proportion_draw, 0.975, na.rm = TRUE)
    out$direction_classification[i] <- if (!out$total_effect_stable[i]) {
      "unstable_total_effect_ratio_not_interpreted"
    } else if (sign(out$indirect_effect[i]) == sign(out$total_effect[i])) {
      "concordant_mediation"
    } else {
      "opposing_or_suppressing_mediation"
    }
  }
  if (length(estimable)) {
    out[estimable, indirect_P_FDR_actual := p.adjust(indirect_P_delta, "BH")]
    out[estimable, indirect_P_Bonferroni_planned_40 := pmin(1, indirect_P_delta * planned_family)]
    out[estimable, corrected_significance_FDR_actual := indirect_P_FDR_actual < 0.05]
    out[estimable, corrected_significance_Bonferroni_planned_40 := indirect_P_Bonferroni_planned_40 < 0.05]
  }
  out[, interpretation_boundary := paste(
    "Pre-specified five-protein cross-platform extension sensitivity; not a replacement for the frozen 25-protein primary analysis.",
    "Independent normal alpha/beta bootstrap excludes selection uncertainty, winner's curse, outcome reuse and empirical cross-protein covariance."
  )]
  setorder(out, variant, outcome, gene_symbol)
  out
}

sets <- list(
  main = list(
    beta_file = file.path(extension_root, "tables", "PWAS5_beta_main.tsv"),
    required = "reestimated", strict = FALSE,
    output = file.path(extension_root, "tables", "PWAS5_two_step_mediation_main.tsv")
  ),
  cis = list(
    beta_file = file.path(extension_root, "tables", "PWAS5_beta_cis.tsv"),
    required = "reestimated_cis_only", strict = FALSE,
    output = file.path(extension_root, "tables", "PWAS5_two_step_mediation_cis.tsv")
  ),
  strict = list(
    beta_file = file.path(extension_root, "tables", "PWAS5_beta_main.tsv"),
    required = "reestimated", strict = TRUE,
    output = file.path(extension_root, "tables", "PWAS5_two_step_mediation_strict.tsv")
  )
)

for (name in names(sets)) {
  spec <- sets[[name]]
  result <- run_set(name, spec$beta_file, spec$required, spec$strict)
  fwrite(result, spec$output, sep = "\t", na = "NA")
  message(name, ": estimable paths = ", nrow(result[mediation_status == "estimable"]), "/", planned_family)
}

writeLines(capture.output(sessionInfo()), file.path(extension_root, "logs", "mediation_sessionInfo.txt"))
