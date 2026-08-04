#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

set.seed(20260711)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
source_dir <- file.path(root, "figures", "source_data")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

provenance <- fread(file.path(root, "tables", "Table_Literature_Prioritized_Protein_Provenance.tsv"), na.strings = c("NA", ""))
alpha <- fread(file.path(root, "tables", "APOE_variant_to_literature_proteins_alpha.tsv"), na.strings = c("NA", ""))
beta <- fread(file.path(root, "tables", "literature_panel_beta_results.tsv"), na.strings = c("NA", ""))
med_main <- fread(file.path(root, "tables", "APOE_linkable_two_step_mediation.tsv"), na.strings = c("NA", ""))
med_cis <- fread(file.path(root, "tables", "APOE_linkable_two_step_mediation_cis_sensitivity.tsv"), na.strings = c("NA", ""))
med_strict <- fread(file.path(root, "tables", "APOE_linkable_two_step_mediation_strict_sensitivity.tsv"), na.strings = c("NA", ""))
comparison <- fread(file.path(root, "tables", "new_vs_legacy_panel_comparison.tsv"), na.strings = c("NA", ""))

final_n <- uniqueN(med_main[row_type == "protein", gene_symbol])
stopifnot(final_n == 25L, med_main[row_type == "protein", .N] == final_n * 2L * 4L)

alpha[, mapping_eligible := tolower(as.character(eligible_for_two_step_mapping)) == "true"]
alpha_use <- unique(alpha[
  availability_status == "direct_variant_available" & mapping_eligible == TRUE,
  .(gene_symbol, variant, alpha = as.numeric(beta), SE_alpha = as.numeric(SE),
    alpha_P = as.numeric(P_value), assay_target_ID, target_mapping_confidence)
], by = c("gene_symbol", "variant"))

beta_use <- beta[
  method_role == "primary" & beta_status %in% c("reestimated", "reestimated_cis_only"),
  .(gene_symbol, outcome, analysis_set, beta = as.numeric(beta), SE_beta = as.numeric(SE),
    beta_P = as.numeric(P_value), beta_P_FDR = as.numeric(P_FDR_observed_reestimated),
    beta_P_Bonferroni = as.numeric(P_Bonferroni_planned_family), nsnp = as.integer(nsnp),
    Q_P = as.numeric(Q_P), Egger_intercept_P = as.numeric(Egger_intercept_P),
    beta_assay_target_ID = assay_target_ID)
]
beta_use[, instrument_set := fifelse(
  analysis_set == "cis_only_sensitivity", "Cis-only", "Genome-wide pQTL"
)]
linkage <- merge(alpha_use, beta_use, by = "gene_symbol", allow.cartesian = TRUE)
stopifnot(all(linkage$assay_target_ID == linkage$beta_assay_target_ID))
linkage[, pleiotropy_flag := !is.na(Egger_intercept_P) & Egger_intercept_P < 0.05]
linkage[, significance := fifelse(
  beta_P_Bonferroni < 0.05, "Planned Bonferroni",
  fifelse(beta_P_FDR < 0.05, "Observed-set FDR", "Not corrected-significant")
)]
provenance[, PP_H4_num := suppressWarnings(as.numeric(PP_H4))]
coloc_genes <- unique(provenance[inclusion_status == "primary_high_confidence_panel" & PP_H4_num >= 0.8, gene_symbol])
linkage[, coloc_supported := gene_symbol %in% coloc_genes]
fwrite(linkage, file.path(source_dir, "Figure_5b_APOE_linkage_scatter.csv"))

forest <- med_main[row_type == "protein", .(
  gene_symbol, variant, outcome, indirect_effect, SE_indirect_delta,
  indirect_CI_lower_bootstrap, indirect_CI_upper_bootstrap,
  mediated_proportion, mediated_proportion_CI_lower_bootstrap,
  mediated_proportion_CI_upper_bootstrap, indirect_P_delta, indirect_P_FDR,
  indirect_P_Bonferroni, direction_classification, beta_nsnp, beta_Q_P,
  beta_Egger_intercept_P, mapping_confidence, alpha_assay_target_ID, beta_assay_target_ID
)]
forest[, `:=`(
  mediated_percent = 100 * mediated_proportion,
  CI_lower_percent = 100 * mediated_proportion_CI_lower_bootstrap,
  CI_upper_percent = 100 * mediated_proportion_CI_upper_bootstrap
)]
gene_order <- forest[, .(max_abs = max(abs(mediated_percent), na.rm = TRUE)), by = gene_symbol][order(max_abs), gene_symbol]
forest[, gene_order_rank := match(gene_symbol, gene_order)]
setorder(forest, gene_order_rank, variant, outcome)
fwrite(forest, file.path(source_dir, "Figure_5c_protein_mediation_forest.csv"))

total_from_mediation <- function(dat, label) {
  dat[row_type == "total", .(
    panel = "Literature-prioritized", instrument_set = label, variant, outcome,
    n = as.integer(number_of_eligible_proteins),
    mediated_percent = 100 * as.numeric(mediated_proportion),
    CI_lower_percent = 100 * as.numeric(mediated_proportion_CI_lower_bootstrap),
    CI_upper_percent = 100 * as.numeric(mediated_proportion_CI_upper_bootstrap),
    CI_available = TRUE
  )]
}
boundary <- rbindlist(list(
  total_from_mediation(med_main, "Genome-wide pQTL"),
  total_from_mediation(med_cis, "Cis-only"),
  total_from_mediation(med_strict, "Strict annotation"),
  comparison[record_type == "total_mediation" & panel == "legacy", .(
    panel = "Legacy biology-guided", instrument_set = "Legacy rs429358",
    variant, outcome, n = as.integer(n), mediated_percent = 100 * as.numeric(value),
    CI_lower_percent = NA_real_, CI_upper_percent = NA_real_, CI_available = FALSE
  )]
), fill = TRUE)
membership_summary <- comparison[record_type == "panel_summary", .(metric, value = as.numeric(value))]
boundary[, overlap_genes := membership_summary[metric == "overlap_genes", value]]
boundary[, jaccard_index := membership_summary[metric == "Jaccard_index", value]]
fwrite(boundary, file.path(source_dir, "Figure_5d_panel_boundary_comparison.csv"))
fwrite(comparison[record_type == "membership"], file.path(source_dir, "Figure_5d_panel_membership.csv"))
message("Independent figure source data written to: ", source_dir)
