#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

set.seed(20260711)
n_boot <- 10000L
analysis_args <- commandArgs(trailingOnly = TRUE)
analysis_set_arg <- grep("^--analysis-set=", analysis_args, value = TRUE)
analysis_set <- if (length(analysis_set_arg)) sub("^--analysis-set=", "", analysis_set_arg[1]) else "genome_wide_instruments_primary"
if (!analysis_set %in% c("genome_wide_instruments_primary", "cis_only_sensitivity", "strict_annotation_sensitivity")) {
  stop("Unsupported --analysis-set: ", analysis_set)
}
selected_analysis_set <- if (analysis_set == "strict_annotation_sensitivity") "genome_wide_instruments_primary" else analysis_set

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = TRUE)
upgrade_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)

alpha_file <- file.path(upgrade_root, "tables", "APOE_variant_to_literature_proteins_alpha.tsv")
beta_file <- file.path(upgrade_root, "tables", "literature_panel_beta_results.tsv")
total_file <- file.path(upgrade_root, "tables", "APOE_variant_total_effects_current_A1.tsv")
provenance_file <- file.path(upgrade_root, "tables", "Table_Literature_Prioritized_Protein_Provenance.tsv")
target_mapping_file <- file.path(upgrade_root, "tables", "APOE_linkable_target_assay_mapping.tsv")
output_file <- file.path(
  upgrade_root, "tables",
  if (analysis_set == "cis_only_sensitivity") "APOE_linkable_two_step_mediation_cis_sensitivity.tsv" else
    if (analysis_set == "strict_annotation_sensitivity") "APOE_linkable_two_step_mediation_strict_sensitivity.tsv" else
      "APOE_linkable_two_step_mediation.tsv"
)
eligibility_file <- file.path(upgrade_root, "tables", "APOE_linkable_eligibility.tsv")
flow_file <- file.path(upgrade_root, "tables", "APOE_linkable_subset_flow.tsv")
bootstrap_file <- file.path(
  upgrade_root, "data_processed",
  if (analysis_set == "cis_only_sensitivity") "APOE_linkable_cis_only_bootstrap_summary.tsv" else
    if (analysis_set == "strict_annotation_sensitivity") "APOE_linkable_strict_bootstrap_summary.tsv" else
      "APOE_linkable_bootstrap_summary.tsv"
)
log_file <- file.path(
  upgrade_root, "logs",
  if (analysis_set == "cis_only_sensitivity") "two_step_mediation_cis_only_sessionInfo.txt" else
    if (analysis_set == "strict_annotation_sensitivity") "two_step_mediation_strict_sessionInfo.txt" else
      "two_step_mediation_sessionInfo.txt"
)

alpha <- fread(alpha_file, na.strings = c("NA", ""))
alpha[, eligible_for_two_step_mapping := tolower(as.character(eligible_for_two_step_mapping)) == "true"]
beta <- fread(beta_file, na.strings = c("NA", ""))
if (!"analysis_set" %in% names(beta)) {
  beta[, analysis_set := fifelse(beta_status == "reestimated", "genome_wide_instruments_primary", "planned_not_reestimable")]
}
totals <- fread(total_file, na.strings = c("NA", ""))
provenance <- fread(provenance_file, na.strings = c("NA", ""))
target_mapping <- fread(target_mapping_file, na.strings = c("NA", ""))
target_mapping[, mapping_eligible_bool := if (analysis_set == "strict_annotation_sensitivity") {
  tolower(as.character(eligible_for_strict_sensitivity)) == "true"
} else {
  tolower(as.character(eligible_for_primary)) == "true"
}]

primary_targets <- unique(provenance[inclusion_status == "primary_high_confidence_panel",
                                     .(gene_symbol, protein_name, protein_form_or_isoform)])
mapping_genes <- unique(target_mapping[mapping_eligible_bool == TRUE, standardized_gene_symbol])
alpha_counts <- alpha[
  availability_status == "direct_variant_available" & eligible_for_two_step_mapping == TRUE &
    gene_symbol %in% mapping_genes,
  .(n_variants = uniqueN(variant)), by = gene_symbol
]
alpha_available <- alpha_counts[n_variants == 2L, gene_symbol]
beta_counts <- beta[method_role == "primary" & analysis_set == selected_analysis_set &
                      beta_status %in% c("reestimated", "reestimated_cis_only") &
                      gene_symbol %in% mapping_genes,
                    .(n_outcomes = uniqueN(outcome)), by = gene_symbol]
beta_available <- beta_counts[n_outcomes == 4L, gene_symbol]
eligible_genes <- intersect(alpha_available, beta_available)

eligibility <- copy(primary_targets)
eligibility <- merge(
  eligibility,
  target_mapping[, .(
    gene_symbol, protein_name, protein_form_or_isoform, local_Olink_target_ID,
    target_assay_mapping_status = mapping_status,
    target_assay_mapping_confidence = mapping_confidence,
    target_assay_mapping_eligible = mapping_eligible_bool
  )],
  by = c("gene_symbol", "protein_name", "protein_form_or_isoform"),
  all.x = TRUE
)
eligibility[, local_Olink_assay_available := !is.na(local_Olink_target_ID) & local_Olink_target_ID != "mapping_unresolved"]
eligibility[, rs429358_alpha_available := target_assay_mapping_eligible %in% TRUE & gene_symbol %in% alpha[
  variant == "rs429358" & availability_status == "direct_variant_available", gene_symbol
]]
eligibility[, rs7412_alpha_available := target_assay_mapping_eligible %in% TRUE & gene_symbol %in% alpha[
  variant == "rs7412" & availability_status == "direct_variant_available", gene_symbol
]]
eligibility[, beta_reestimated_current_A1 := target_assay_mapping_eligible %in% TRUE & gene_symbol %in% beta_available]
eligibility[, eligible_for_two_step_MR := target_assay_mapping_eligible %in% TRUE & gene_symbol %in% eligible_genes]
eligibility[, exclusion_reason := fifelse(
  eligible_for_two_step_MR, "NA",
  fifelse(!local_Olink_assay_available,
          "local_Olink_assay_unavailable",
          fifelse(!(target_assay_mapping_eligible %in% TRUE),
                  "target_to_Olink_assay_mapping_unresolved",
          fifelse(!rs429358_alpha_available | !rs7412_alpha_available,
                  "direct_APOE_variant_alpha_unavailable", "beta_not_reestimable")))
)]
if (analysis_set == "genome_wide_instruments_primary") {
  fwrite(eligibility, eligibility_file, sep = "\t", na = "NA")
}

matched <- target_mapping[tolower(as.character(eligible_for_primary)) == "true"]
matched_gene_assays <- unique(matched[, .(standardized_gene_symbol, UKB_PPP_OID)])
alpha_429 <- unique(alpha[availability_status == "direct_variant_available" & variant == "rs429358" &
                            gene_symbol %in% matched_gene_assays$standardized_gene_symbol, gene_symbol])
alpha_7412 <- unique(alpha[availability_status == "direct_variant_available" & variant == "rs7412" &
                             gene_symbol %in% matched_gene_assays$standardized_gene_symbol, gene_symbol])
both_alpha <- intersect(alpha_429, alpha_7412)
beta_any <- beta_counts[n_outcomes >= 1L, gene_symbol]
beta_four <- beta_counts[n_outcomes == 4L, gene_symbol]
flow <- data.table(
  stage = c(
    "primary_high_confidence_target_entries", "primary_high_confidence_unique_genes",
    "name_matched_unique_genes", "name_matched_unique_Olink_assays",
    "genes_with_rs429358_alpha", "genes_with_rs7412_alpha", "genes_with_both_alpha",
    "genes_with_at_least_one_estimable_beta", "genes_with_four_outcome_beta",
    "final_mediation_genes", "final_mediation_assays", "mediation_path_count"
  ),
  n = c(
    nrow(primary_targets), uniqueN(primary_targets$gene_symbol),
    uniqueN(matched_gene_assays$standardized_gene_symbol), uniqueN(matched_gene_assays$UKB_PPP_OID),
    length(alpha_429), length(alpha_7412), length(both_alpha), length(beta_any), length(beta_four),
    length(eligible_genes), uniqueN(matched_gene_assays[standardized_gene_symbol %in% eligible_genes, UKB_PPP_OID]),
    length(eligible_genes) * 2L * 4L
  ),
  denominator = c("target entries", "unique genes", "unique genes", "unique assays",
                  rep("unique genes", 6), "unique assays", "protein-variant-outcome paths"),
  notes = c(
    "Literature target entries retain distinct reported protein forms.",
    "Gene-level deduplication does not imply exact originating-assay replication.",
    "High- and moderate-confidence approved-name/gene-symbol mappings are included.",
    "One selected Olink assay per included gene; one-to-many mappings are not auto-selected.",
    "Direct rs429358-C alpha; no proxy.", "Direct rs7412-T alpha; no proxy.",
    "Both direct APOE alpha estimates are required for the eight-path design.",
    "At least one current-A1 outcome beta is estimable.",
    "All four current-A1 outcome betas are estimable.",
    "Same-assay alpha and beta with both APOE variants and all four outcomes.",
    "Unique assays represented by final mediation genes.",
    "Asserted as final genes x 2 APOE variants x 4 outcomes."
  )
)
if (analysis_set == "genome_wide_instruments_primary") {
  fwrite(flow, flow_file, sep = "\t")
}

alpha_use <- alpha[
  gene_symbol %in% eligible_genes & availability_status == "direct_variant_available" &
    eligible_for_two_step_mapping == TRUE,
  .(gene_symbol, variant, alpha = as.numeric(beta), SE_alpha = as.numeric(SE), alpha_P = as.numeric(P_value),
    alpha_effect_allele = requested_effect_allele, alpha_assay_target_ID = assay_target_ID,
    alpha_sample_size = sample_size, alpha_source)
]
alpha_use <- unique(alpha_use, by = c("gene_symbol", "variant"))
beta_use <- beta[
  gene_symbol %in% eligible_genes & method_role == "primary" & analysis_set == selected_analysis_set &
    beta_status %in% c("reestimated", "reestimated_cis_only"),
  .(gene_symbol, outcome, beta = as.numeric(beta), SE_beta = as.numeric(SE), beta_P = as.numeric(P_value),
    beta_method = method, beta_nsnp = nsnp, beta_mean_F = mean_F, beta_min_F = min_F,
    beta_Q_P = Q_P, beta_Egger_intercept_P = Egger_intercept_P,
    beta_P_FDR_observed = P_FDR_observed_reestimated,
    beta_P_Bonferroni_planned = P_Bonferroni_planned_family,
    beta_assay_target_ID = assay_target_ID)
]
total_use <- totals[availability_status == "direct_variant_available",
                    .(variant, outcome, total_effect = as.numeric(beta), SE_total = as.numeric(SE),
                      total_effect_P = as.numeric(P_value), total_effect_allele = requested_effect_allele,
                      total_effect_source = outcome_GWAS)]

med <- merge(alpha_use, beta_use, by = "gene_symbol", allow.cartesian = TRUE)
med <- merge(med, total_use, by = c("variant", "outcome"), all.x = TRUE)
stopifnot(nrow(med) == length(eligible_genes) * 2L * 4L, !anyNA(med$total_effect))
stopifnot(all(med$alpha_assay_target_ID == med$beta_assay_target_ID))
gene_confidence <- target_mapping[mapping_eligible_bool == TRUE,
  .(mapping_confidence = if (any(mapping_confidence == "high")) "high" else "moderate"),
  by = .(gene_symbol = standardized_gene_symbol)]
med <- merge(med, gene_confidence, by = "gene_symbol", all.x = TRUE)
stopifnot(!anyNA(med$mapping_confidence))

med[, indirect_effect := alpha * beta]
med[, SE_indirect_delta := sqrt(beta^2 * SE_alpha^2 + alpha^2 * SE_beta^2)]
med[, indirect_P_delta := 2 * pnorm(-abs(indirect_effect / SE_indirect_delta))]
med[, mediated_proportion := indirect_effect / total_effect]
med[, SE_mediated_proportion_delta := sqrt(
  (beta / total_effect)^2 * SE_alpha^2 +
    (alpha / total_effect)^2 * SE_beta^2 +
    (alpha * beta / total_effect^2)^2 * SE_total^2
)]
med[, total_effect_stable := abs(total_effect / SE_total) >= 2]
med[, direction_classification := fifelse(
  !total_effect_stable, "unstable_total_effect_ratio_not_interpreted",
  fifelse(sign(indirect_effect) == sign(total_effect), "concordant_mediation", "opposing_or_suppressing_mediation")
)]

bootstrap_rows <- list()
for (i in seq_len(nrow(med))) {
  a <- rnorm(n_boot, med$alpha[i], med$SE_alpha[i])
  b <- rnorm(n_boot, med$beta[i], med$SE_beta[i])
  t <- rnorm(n_boot, med$total_effect[i], med$SE_total[i])
  indirect <- a * b
  proportion <- indirect / t
  med$indirect_CI_lower_bootstrap[i] <- quantile(indirect, 0.025, na.rm = TRUE)
  med$indirect_CI_upper_bootstrap[i] <- quantile(indirect, 0.975, na.rm = TRUE)
  med$mediated_proportion_CI_lower_bootstrap[i] <- quantile(proportion, 0.025, na.rm = TRUE)
  med$mediated_proportion_CI_upper_bootstrap[i] <- quantile(proportion, 0.975, na.rm = TRUE)
  bootstrap_rows[[length(bootstrap_rows) + 1L]] <- data.table(
    row_type = "protein", gene_symbol = med$gene_symbol[i], variant = med$variant[i], outcome = med$outcome[i],
    n_boot = n_boot, indirect_median = median(indirect), proportion_median = median(proportion),
    assumption = "Independent normal coefficient draws; cross-protein covariance and selection uncertainty not included."
  )
}

med[, mediation_family_size := .N]
med[, indirect_P_FDR := p.adjust(indirect_P_delta, method = "BH")]
med[, indirect_P_Bonferroni := pmin(1, indirect_P_delta * mediation_family_size)]
med[, row_type := "protein"]
med[, analysis_set := analysis_set]

total_rows <- list()
mediation_groups <- unique(med[, .(variant, outcome)])
for (group_index in seq_len(nrow(mediation_groups))) {
  group_variant <- mediation_groups$variant[group_index]
  group_outcome <- mediation_groups$outcome[group_index]
  subset <- med[variant == group_variant & outcome == group_outcome]
  total_draw <- rnorm(n_boot, subset$total_effect[1], subset$SE_total[1])
  indirect_draws <- matrix(NA_real_, nrow = n_boot, ncol = nrow(subset))
  for (j in seq_len(nrow(subset))) {
    indirect_draws[, j] <- rnorm(n_boot, subset$alpha[j], subset$SE_alpha[j]) *
      rnorm(n_boot, subset$beta[j], subset$SE_beta[j])
  }
  summed_indirect <- rowSums(indirect_draws)
  summed_proportion <- summed_indirect / total_draw
  total_indirect <- sum(subset$indirect_effect)
  total_se_delta <- sqrt(sum(subset$SE_indirect_delta^2))
  total_proportion <- total_indirect / subset$total_effect[1]
  total_rows[[length(total_rows) + 1L]] <- data.table(
    row_type = "total", gene_symbol = "TOTAL_ELIGIBLE_PROTEINS",
    variant = group_variant, outcome = group_outcome,
    alpha = NA_real_, SE_alpha = NA_real_, alpha_P = NA_real_, beta = NA_real_, SE_beta = NA_real_, beta_P = NA_real_,
    indirect_effect = total_indirect, SE_indirect_delta = total_se_delta,
    indirect_P_delta = 2 * pnorm(-abs(total_indirect / total_se_delta)),
    indirect_CI_lower_bootstrap = quantile(summed_indirect, 0.025),
    indirect_CI_upper_bootstrap = quantile(summed_indirect, 0.975),
    total_effect = subset$total_effect[1], SE_total = subset$SE_total[1],
    total_effect_P = subset$total_effect_P[1], total_effect_allele = subset$total_effect_allele[1],
    total_effect_source = subset$total_effect_source[1],
    mediated_proportion = total_proportion,
    SE_mediated_proportion_delta = NA_real_,
    mediated_proportion_CI_lower_bootstrap = quantile(summed_proportion, 0.025),
    mediated_proportion_CI_upper_bootstrap = quantile(summed_proportion, 0.975),
    total_effect_stable = subset$total_effect_stable[1],
    direction_classification = fifelse(
      !subset$total_effect_stable[1], "unstable_total_effect_ratio_not_interpreted",
      fifelse(sign(total_indirect) == sign(subset$total_effect[1]), "concordant_total_mediation", "opposing_or_suppressing_total_mediation")
    ),
    mediation_family_size = nrow(med), indirect_P_FDR = NA_real_, indirect_P_Bonferroni = NA_real_,
    number_of_eligible_proteins = nrow(subset),
    bootstrap_assumption = "Independent normal alpha/beta draws and common total-effect draw; cross-protein covariance, literature selection, publication bias and winner's curse not included."
  )
  total_rows[[length(total_rows)]][, analysis_set := analysis_set]
  bootstrap_rows[[length(bootstrap_rows) + 1L]] <- data.table(
    row_type = "total", gene_symbol = "TOTAL_ELIGIBLE_PROTEINS", variant = group_variant,
    outcome = group_outcome, n_boot = n_boot, indirect_median = median(summed_indirect),
    proportion_median = median(summed_proportion),
    assumption = "Cross-protein covariance and literature-selection uncertainty not included."
  )
}

output <- rbindlist(list(med, rbindlist(total_rows, fill = TRUE)), fill = TRUE)
output[, bootstrap_n := n_boot]
output[, interpretation_boundary := paste(
  "Enriched-panel mediation stress test; not an unbiased estimate across the circulating proteome.",
  "Bootstrap intervals exclude literature selection, publication bias, winner's curse, outcome reuse and cross-protein covariance."
)]
setorder(output, variant, outcome, row_type, gene_symbol)
fwrite(output, output_file, sep = "\t", na = "NA")
fwrite(rbindlist(bootstrap_rows, fill = TRUE), bootstrap_file, sep = "\t", na = "NA")
writeLines(capture.output(sessionInfo()), log_file)

message("Eligible genes: ", paste(eligible_genes, collapse = ", "))
message("Analysis set: ", analysis_set)
message("Protein-level mediation rows: ", nrow(med))
message("Output: ", output_file)
