#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(TwoSampleMR)
  library(yaml)
})

set.seed(20260717)
n_boot <- 10000L

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = TRUE)
extension_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
upgrade_root <- normalizePath(file.path(extension_root, "..", ".."), winslash = "/", mustWork = TRUE)

eldjarn_file <- file.path(
  extension_root, "data_raw", "eldjarn_2023", "extracted", "41586_2023_6563_MOESM3_ESM.xlsx"
)
harmonized_file <- file.path(extension_root, "freeze_before_extension", "literature_panel_harmonized_data.tsv")
beta_file <- file.path(extension_root, "freeze_before_extension", "literature_panel_beta_results.tsv")
mediation_file <- file.path(extension_root, "freeze_before_extension", "APOE_linkable_two_step_mediation.tsv")
alpha_file <- file.path(extension_root, "freeze_before_extension", "APOE_variant_to_literature_proteins_alpha.tsv")
total_file <- file.path(upgrade_root, "tables", "APOE_variant_total_effects_primary_analysis.tsv")
resources <- read_yaml(file.path(upgrade_root, "config", "resources.yml"))

table_dir <- file.path(extension_root, "tables")
processed_dir <- file.path(extension_root, "data_processed")
log_dir <- file.path(extension_root, "logs", "pav_epitope_audit")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(
  file.exists(eldjarn_file), file.exists(harmonized_file), file.exists(beta_file),
  file.exists(mediation_file), file.exists(alpha_file), file.exists(total_file),
  file.exists(resources$paths$plink_binary), file.exists(paste0(resources$paths$ld_reference, ".bed"))
)

harm <- fread(harmonized_file, na.strings = c("NA", ""), showProgress = FALSE)
iv <- unique(harm[, .(
  gene_symbol = exposure,
  SNP,
  chromosome = as.integer(chr.exposure),
  position_hg19 = as.integer(pos.exposure)
)])
if (nrow(iv) != 409L || uniqueN(iv$gene_symbol) != 25L) {
  stop("Frozen instrument assertion failed before the PAV audit.")
}

frozen_beta <- fread(beta_file, na.strings = c("NA", ""))
assay_meta <- unique(frozen_beta[
  outcome == "AD" & method_role == "primary" & analysis_set == "genome_wide_instruments_primary",
  .(gene_symbol, assay_target_ID, UniProt_ID)
])
if (nrow(assay_meta) != 25L || anyDuplicated(assay_meta$gene_symbol)) {
  stop("Expected one frozen Olink assay for each of 25 proteins.")
}
iv <- merge(iv, assay_meta, by = "gene_symbol", all.x = TRUE)
stopifnot(!anyNA(iv$assay_target_ID))

message("Reading Eldjarn et al. Supplementary Table 21.")
st21 <- as.data.table(read_excel(eldjarn_file, sheet = "ST21_olink_ukb_bi_pQTLs", skip = 3))
st21[, st21_row_id := .I]
st21_core <- st21[, .(
  st21_row_id,
  assay_target_ID = as.character(OlinkID),
  affected_gene = as.character(`Affected Protein Gene name`),
  st21_pQTL_ID = as.character(`pQTL ID`),
  st21_cis_trans = tolower(as.character(CisOrTrans)),
  st21_marker = as.character(marker),
  st21_rsName = as.character(rsName),
  st21_rank_within_locus = suppressWarnings(as.integer(`Rank Within Locus`)),
  st21_r2_to_sentinel = suppressWarnings(as.numeric(`r2 sentinel`)),
  st21_cis_eQTL_genes = as.character(`cis eQTL genes`),
  st21_PAV_genes = as.character(`PAV genes`),
  st21_cis_eQTL_same_gene = toupper(as.character(`cis eQTL same gene`)),
  st21_PAV_same_gene = toupper(as.character(`PAV same gene`)),
  st21_UniProt = as.character(UniProt),
  st21_Olink_version = as.character(`Olink version`),
  st21_protein_on_SomaScan = toupper(as.character(`protein on soma`))
)]
st21_core <- st21_core[assay_target_ID %chin% iv$assay_target_ID]
st21_rsid <- st21_core[!is.na(st21_rsName) & st21_rsName != "."]
st21_rsid <- st21_rsid[, .(st21_SNP = trimws(unlist(strsplit(st21_rsName, ",", fixed = TRUE)))), by = setdiff(names(st21_rsid), "st21_rsName")]
st21_rsid <- st21_rsid[grepl("^rs[0-9]+$", st21_SNP)]

variant_list <- sort(unique(c(iv$SNP, st21_rsid$st21_SNP)))
extract_file <- file.path(processed_dir, "PAV_audit_LD_variant_list.txt")
writeLines(variant_list, extract_file)

ld_prefix <- file.path(log_dir, "PAV_audit_1000G_EUR_r2_0.8")
plink_args <- c(
  "--bfile", resources$paths$ld_reference,
  "--extract", extract_file,
  "--r2",
  "--ld-window", "999999",
  "--ld-window-kb", "10000",
  "--ld-window-r2", "0.8",
  "--allow-no-sex",
  "--out", ld_prefix
)
status <- system2(
  resources$paths$plink_binary,
  args = plink_args,
  stdout = paste0(ld_prefix, ".stdout.txt"),
  stderr = paste0(ld_prefix, ".stderr.txt")
)
if (status != 0L || !file.exists(paste0(ld_prefix, ".ld"))) {
  stop("PLINK LD calculation failed; inspect the PAV audit log directory.")
}

ld <- fread(paste0(ld_prefix, ".ld"), fill = TRUE)
if (nrow(ld)) {
  ld_edges <- unique(rbindlist(list(
    ld[, .(current_SNP = SNP_A, st21_SNP = SNP_B, LD_R2 = as.numeric(R2))],
    ld[, .(current_SNP = SNP_B, st21_SNP = SNP_A, LD_R2 = as.numeric(R2))]
  )))
} else {
  ld_edges <- data.table(current_SNP = character(), st21_SNP = character(), LD_R2 = numeric())
}

exact_matches <- merge(
  iv,
  st21_rsid,
  by.x = c("assay_target_ID", "SNP"),
  by.y = c("assay_target_ID", "st21_SNP"),
  allow.cartesian = TRUE
)
if (nrow(exact_matches)) {
  exact_matches[, `:=`(current_SNP = SNP, st21_SNP = SNP, LD_R2 = 1, link_type = "exact_rsID_same_assay")]
}

ld_matches <- merge(
  iv,
  ld_edges,
  by.x = "SNP",
  by.y = "current_SNP",
  allow.cartesian = TRUE
)
ld_matches <- merge(
  ld_matches,
  st21_rsid,
  by = c("assay_target_ID", "st21_SNP"),
  allow.cartesian = TRUE
)
if (nrow(ld_matches)) {
  ld_matches[, link_type := "r2_ge_0.80_same_assay"]
}

match_columns <- union(names(exact_matches), names(ld_matches))
for (column in setdiff(match_columns, names(exact_matches))) exact_matches[, (column) := NA]
for (column in setdiff(match_columns, names(ld_matches))) ld_matches[, (column) := NA]
matches <- unique(rbindlist(list(exact_matches[, ..match_columns], ld_matches[, ..match_columns]), fill = TRUE),
                  by = c("gene_symbol", "assay_target_ID", "SNP", "st21_row_id", "st21_SNP"))
if (nrow(matches)) {
  matches[, high_confidence_epitope_risk :=
            st21_cis_trans == "cis" & st21_PAV_same_gene == "Y" & st21_cis_eQTL_same_gene != "Y"]
  matches[, target_gene_PAV_with_cis_eQTL_support :=
            st21_cis_trans == "cis" & st21_PAV_same_gene == "Y" & st21_cis_eQTL_same_gene == "Y"]
  setorder(matches, gene_symbol, SNP, -high_confidence_epitope_risk, -LD_R2, st21_row_id)
}
fwrite(matches, file.path(table_dir, "PAV_epitope_annotation_matches.tsv"), sep = "\t", na = "NA")

match_summary <- if (nrow(matches)) {
  matches[, .(
    n_linked_Eldjarn_records = .N,
    best_LD_R2 = max(LD_R2, na.rm = TRUE),
    exact_same_assay_rsID_match = any(link_type == "exact_rsID_same_assay"),
    any_target_gene_PAV = any(st21_cis_trans == "cis" & st21_PAV_same_gene == "Y", na.rm = TRUE),
    any_target_gene_cis_eQTL = any(st21_cis_trans == "cis" & st21_cis_eQTL_same_gene == "Y", na.rm = TRUE),
    high_confidence_epitope_risk = any(high_confidence_epitope_risk, na.rm = TRUE),
    PAV_with_cis_eQTL_support = any(target_gene_PAV_with_cis_eQTL_support, na.rm = TRUE),
    linked_st21_pQTL_IDs = paste(sort(unique(st21_pQTL_ID)), collapse = ";"),
    linked_st21_SNPs = paste(sort(unique(st21_SNP)), collapse = ";"),
    linked_PAV_genes = paste(sort(unique(na.omit(st21_PAV_genes))), collapse = ";"),
    linked_cis_eQTL_genes = paste(sort(unique(na.omit(st21_cis_eQTL_genes))), collapse = ";")
  ), by = .(gene_symbol, assay_target_ID, UniProt_ID, SNP, chromosome, position_hg19)]
} else {
  data.table()
}

audit <- merge(
  iv,
  match_summary,
  by = c("gene_symbol", "assay_target_ID", "UniProt_ID", "SNP", "chromosome", "position_hg19"),
  all.x = TRUE
)
audit[, annotation_status := fifelse(
  is.na(n_linked_Eldjarn_records), "unresolved_no_exact_or_high_LD_same_assay_record",
  fifelse(high_confidence_epitope_risk %in% TRUE, "high_confidence_epitope_risk",
    fifelse(PAV_with_cis_eQTL_support %in% TRUE, "PAV_with_target_gene_cis_eQTL_support_not_excluded",
      fifelse(any_target_gene_PAV %in% TRUE, "target_gene_PAV_annotation_uncertain",
              "no_target_gene_PAV_in_linked_Eldjarn_record")))
)]
audit[, exclusion_in_primary_PAV_sensitivity := annotation_status == "high_confidence_epitope_risk"]
audit[, interpretation_boundary := fifelse(
  annotation_status == "unresolved_no_exact_or_high_LD_same_assay_record",
  "Unresolved is not equivalent to no epitope risk.",
  "Eldjarn ST21 annotation uses protein-altering variants and cis-eQTLs in high LD (r2 > 0.8) with the pQTL."
)]
setorder(audit, gene_symbol, SNP)
fwrite(audit, file.path(table_dir, "PAV_epitope_instrument_audit.tsv"), sep = "\t", na = "NA")

risk_pairs <- audit[exclusion_in_primary_PAV_sensitivity == TRUE, .(gene_symbol, SNP)]
risk_keys <- if (nrow(risk_pairs)) paste(risk_pairs$gene_symbol, risk_pairs$SNP, sep = "|") else character()
harm[, pair_key := paste(exposure, SNP, sep = "|")]
harm_filtered <- harm[!pair_key %chin% risk_keys & mr_keep %in% TRUE]

estimate_filtered_beta <- function(dat, gene, outcome_name) {
  if (!nrow(dat)) {
    return(data.table(
      gene_symbol = gene, outcome = outcome_name, method = "not_estimable", method_role = "primary",
      nsnp = 0L, beta = NA_real_, SE = NA_real_, P_value = NA_real_, Cochran_Q = NA_real_, Q_df = NA_real_,
      Q_P = NA_real_, Egger_intercept = NA_real_, Egger_intercept_SE = NA_real_, Egger_intercept_P = NA_real_,
      beta_status = "not_estimable_after_high_confidence_epitope_filter"
    ))
  }
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
    gene_symbol = gene,
    outcome = outcome_name,
    method_role = fifelse(method %chin% c("Inverse variance weighted", "Wald ratio"), "primary", "sensitivity"),
    CI_lower = b - 1.96 * se,
    CI_upper = b + 1.96 * se,
    OR = exp(b),
    OR_CI_lower = exp(b - 1.96 * se),
    OR_CI_upper = exp(b + 1.96 * se),
    Cochran_Q = Q,
    Q_P = Q_pval,
    Egger_intercept = if (nrow(pleiotropy)) pleiotropy$egger_intercept[1] else NA_real_,
    Egger_intercept_SE = if (nrow(pleiotropy)) pleiotropy$se[1] else NA_real_,
    Egger_intercept_P = if (nrow(pleiotropy)) pleiotropy$pval[1] else NA_real_,
    mean_F = mean((dat$beta.exposure / dat$se.exposure)^2),
    min_F = min((dat$beta.exposure / dat$se.exposure)^2),
    beta_status = "reestimated_after_high_confidence_epitope_filter",
    excluded_high_confidence_epitope_instruments = sum(risk_pairs$gene_symbol == gene),
    audit_source = "Eldjarn_2023_ST21_same_OID_exact_or_1000G_EUR_r2_ge_0.80"
  )]
  setnames(estimates, c("b", "se", "pval"), c("beta", "SE", "P_value"), skip_absent = TRUE)
  estimates
}

filtered_beta_rows <- list()
for (gene in sort(unique(iv$gene_symbol))) {
  for (outcome_name in c("AD", "any_AMD", "dry_AMD", "wet_AMD")) {
    dat <- harm_filtered[exposure == gene & outcome == outcome_name]
    filtered_beta_rows[[paste(gene, outcome_name, sep = "_")]] <- estimate_filtered_beta(dat, gene, outcome_name)
  }
}
filtered_beta <- rbindlist(filtered_beta_rows, fill = TRUE)
filtered_beta <- merge(filtered_beta, assay_meta, by = "gene_symbol", all.x = TRUE)
filtered_beta[, `:=`(
  planned_family_size = 25L,
  observed_primary_tests_in_family = NA_integer_,
  P_FDR_observed = NA_real_,
  P_Bonferroni_planned = NA_real_
)]
for (outcome_name in c("AD", "any_AMD", "dry_AMD", "wet_AMD")) {
  idx <- which(filtered_beta$outcome == outcome_name & filtered_beta$method_role == "primary" & !is.na(filtered_beta$P_value))
  filtered_beta[idx, observed_primary_tests_in_family := length(idx)]
  filtered_beta[idx, P_FDR_observed := p.adjust(P_value, method = "BH")]
  filtered_beta[idx, P_Bonferroni_planned := pmin(1, P_value * 25L)]
}
fwrite(filtered_beta, file.path(table_dir, "PAV_filtered_protein_beta_sensitivity.tsv"), sep = "\t", na = "NA")

filtered_primary <- filtered_beta[method_role == "primary" & !is.na(beta)]
eligible_genes <- filtered_primary[, .(n_outcomes = uniqueN(outcome)), by = gene_symbol][n_outcomes == 4L, gene_symbol]
alpha <- fread(alpha_file, na.strings = c("NA", ""))
alpha_use <- unique(alpha[
  availability_status == "direct_variant_available" &
    tolower(as.character(eligible_for_two_step_mapping)) == "true" &
    gene_symbol %chin% eligible_genes,
  .(gene_symbol, variant, alpha = as.numeric(beta), SE_alpha = as.numeric(SE), alpha_P = as.numeric(P_value),
    alpha_assay_target_ID = assay_target_ID)
], by = c("gene_symbol", "variant"))
beta_use <- filtered_primary[gene_symbol %chin% eligible_genes, .(
  gene_symbol, outcome, beta = as.numeric(beta), SE_beta = as.numeric(SE), beta_P = as.numeric(P_value),
  beta_method = method, beta_nsnp = nsnp, beta_Q_P = Q_P,
  beta_Egger_intercept_P = Egger_intercept_P, beta_assay_target_ID = assay_target_ID
)]
totals <- fread(total_file, na.strings = c("NA", ""))
total_use <- totals[availability_status == "direct_variant_available", .(
  variant, outcome, total_effect = as.numeric(beta), SE_total = as.numeric(SE),
  total_effect_P = as.numeric(P_value), total_effect_allele = requested_effect_allele,
  total_effect_source = outcome_GWAS
)]

med <- merge(alpha_use, beta_use, by = "gene_symbol", allow.cartesian = TRUE)
med <- merge(med, total_use, by = c("variant", "outcome"), all.x = TRUE)
if (nrow(med) != length(eligible_genes) * 2L * 4L || anyNA(med$total_effect)) {
  stop("PAV-filtered mediation row-count or total-effect assertion failed.")
}
if (!all(med$alpha_assay_target_ID == med$beta_assay_target_ID)) stop("Alpha/beta assay mismatch after PAV filtering.")

med[, `:=`(
  indirect_effect = alpha * beta,
  SE_indirect_delta = sqrt(beta^2 * SE_alpha^2 + alpha^2 * SE_beta^2)
)]
med[, indirect_P_delta := 2 * pnorm(-abs(indirect_effect / SE_indirect_delta))]
med[, mediated_proportion := indirect_effect / total_effect]
med[, direction_classification := fifelse(
  sign(indirect_effect) == sign(total_effect), "concordant_mediation", "opposing_or_suppressing_mediation"
)]
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
}
med[, `:=`(
  row_type = "protein",
  planned_mediation_family_size = 200L,
  observed_mediation_paths = .N,
  indirect_P_FDR_observed = p.adjust(indirect_P_delta, method = "BH"),
  indirect_P_Bonferroni_planned = pmin(1, indirect_P_delta * 200L),
  bootstrap_n = n_boot,
  analysis_label = "high_confidence_epitope_risk_filtered_sensitivity"
)]

aggregate_rows <- list()
groups <- unique(med[, .(variant, outcome)])
for (i in seq_len(nrow(groups))) {
  variant_i <- groups$variant[i]
  outcome_i <- groups$outcome[i]
  subset <- med[variant == variant_i & outcome == outcome_i]
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
  aggregate_rows[[paste(variant_i, outcome_i)]] <- data.table(
    row_type = "total",
    gene_symbol = "TOTAL_ELIGIBLE_PROTEINS",
    variant = variant_i,
    outcome = outcome_i,
    indirect_effect = total_indirect,
    SE_indirect_delta = total_se_delta,
    indirect_P_delta = 2 * pnorm(-abs(total_indirect / total_se_delta)),
    indirect_CI_lower_bootstrap = quantile(summed_indirect, 0.025),
    indirect_CI_upper_bootstrap = quantile(summed_indirect, 0.975),
    total_effect = subset$total_effect[1],
    SE_total = subset$SE_total[1],
    mediated_proportion = total_indirect / subset$total_effect[1],
    mediated_proportion_CI_lower_bootstrap = quantile(summed_proportion, 0.025),
    mediated_proportion_CI_upper_bootstrap = quantile(summed_proportion, 0.975),
    number_of_eligible_proteins = nrow(subset),
    planned_mediation_family_size = 200L,
    observed_mediation_paths = nrow(med),
    bootstrap_n = n_boot,
    analysis_label = "high_confidence_epitope_risk_filtered_sensitivity",
    bootstrap_assumption = "Independent normal coefficient draws; empirical cross-protein covariance, unresolved epitope risk and mapping uncertainty are not included."
  )
}

filtered_mediation <- rbindlist(list(med, rbindlist(aggregate_rows, fill = TRUE)), fill = TRUE)
setorder(filtered_mediation, variant, outcome, row_type, gene_symbol)
fwrite(filtered_mediation, file.path(table_dir, "PAV_filtered_two_step_mediation_sensitivity.tsv"), sep = "\t", na = "NA")

current_primary <- frozen_beta[
  method_role == "primary" & analysis_set == "genome_wide_instruments_primary",
  .(gene_symbol, outcome, current_nsnp = nsnp, current_beta = as.numeric(beta), current_SE = as.numeric(SE), current_P = as.numeric(P_value))
]
beta_comparison <- merge(
  current_primary,
  filtered_primary[, .(gene_symbol, outcome, filtered_nsnp = nsnp, filtered_beta = as.numeric(beta), filtered_SE = as.numeric(SE), filtered_P = as.numeric(P_value))],
  by = c("gene_symbol", "outcome"), all = TRUE
)
beta_comparison[, `:=`(
  instruments_removed = current_nsnp - filtered_nsnp,
  beta_absolute_difference = abs(current_beta - filtered_beta),
  direction_changed = sign(current_beta) != sign(filtered_beta)
)]
fwrite(beta_comparison, file.path(table_dir, "PAV_filtered_beta_comparison.tsv"), sep = "\t", na = "NA")

current_mediation <- fread(mediation_file, na.strings = c("NA", ""))
aggregate_comparison <- merge(
  current_mediation[row_type == "total", .(
    variant, outcome,
    unfiltered_number_of_proteins = number_of_eligible_proteins,
    unfiltered_mediated_proportion = mediated_proportion,
    unfiltered_CI_lower = mediated_proportion_CI_lower_bootstrap,
    unfiltered_CI_upper = mediated_proportion_CI_upper_bootstrap
  )],
  filtered_mediation[row_type == "total", .(
    variant, outcome,
    filtered_number_of_proteins = number_of_eligible_proteins,
    filtered_mediated_proportion = mediated_proportion,
    filtered_CI_lower = mediated_proportion_CI_lower_bootstrap,
    filtered_CI_upper = mediated_proportion_CI_upper_bootstrap
  )],
  by = c("variant", "outcome"), all = TRUE
)
aggregate_comparison[, mediated_proportion_absolute_difference :=
                       abs(unfiltered_mediated_proportion - filtered_mediated_proportion)]
fwrite(aggregate_comparison, file.path(table_dir, "PAV_filtered_aggregate_comparison.tsv"), sep = "\t", na = "NA")

status_summary <- audit[, .(n = .N), by = .(category = annotation_status)]
status_summary[, summary_type := "instrument_annotation_status"]
setcolorder(status_summary, c("summary_type", "category", "n"))
audit_summary <- rbindlist(list(
  status_summary,
  data.table(
    summary_type = "analysis_counts",
    category = c(
      "frozen_protein_SNP_pairs", "exact_or_high_LD_same_assay_annotated_pairs",
      "unresolved_pairs", "high_confidence_epitope_risk_pairs", "eligible_proteins_after_filter",
      "protein_mediation_paths_after_filter"
    ),
    n = c(
      nrow(audit), audit[!grepl("^unresolved", annotation_status), .N],
      audit[grepl("^unresolved", annotation_status), .N], nrow(risk_pairs),
      length(eligible_genes), nrow(med)
    )
  )
), fill = TRUE)
fwrite(audit_summary, file.path(table_dir, "PAV_epitope_audit_summary.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo.txt"))

message("PAV/epitope audit completed.")
message("Annotated exact-or-high-LD same-assay pairs: ", audit[!grepl("^unresolved", annotation_status), .N], "/", nrow(audit))
message("High-confidence epitope-risk pairs excluded in sensitivity: ", nrow(risk_pairs))
message("Eligible proteins after filtering: ", length(eligible_genes))
