#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
})

set.seed(20260717)
n_boot <- 10000L

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = TRUE)
extension_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
upgrade_root <- normalizePath(file.path(extension_root, "..", ".."), winslash = "/", mustWork = TRUE)

raw_file <- file.path(extension_root, "data_raw", "kunkle_2019", "Kunkle_etal_Stage1_results.txt")
frozen_harmonized_file <- file.path(extension_root, "freeze_before_extension", "literature_panel_harmonized_data.tsv")
frozen_beta_file <- file.path(extension_root, "freeze_before_extension", "literature_panel_beta_results.tsv")
frozen_mediation_file <- file.path(extension_root, "freeze_before_extension", "APOE_linkable_two_step_mediation.tsv")
alpha_file <- file.path(extension_root, "freeze_before_extension", "APOE_variant_to_literature_proteins_alpha.tsv")

table_dir <- file.path(extension_root, "tables")
processed_dir <- file.path(extension_root, "data_processed")
log_dir <- file.path(extension_root, "logs")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(raw_file), file.exists(frozen_harmonized_file), file.exists(alpha_file))

harm <- fread(frozen_harmonized_file, na.strings = c("NA", ""), showProgress = FALSE)
iv <- unique(harm[, .(
  gene_symbol = exposure, SNP,
  effect_allele = effect_allele.exposure,
  other_allele = other_allele.exposure,
  beta_exposure = as.numeric(beta.exposure),
  se_exposure = as.numeric(se.exposure),
  eaf_exposure = as.numeric(eaf.exposure),
  pval_exposure = as.numeric(pval.exposure),
  samplesize_exposure = as.numeric(samplesize.exposure),
  chr_exposure = as.integer(chr.exposure),
  pos_exposure = as.integer(pos.exposure)
)])

if (nrow(iv) != 409L || uniqueN(iv$gene_symbol) != 25L) {
  stop(sprintf("Frozen instrument assertion failed: expected 409 pairs and 25 proteins; observed %d and %d.",
               nrow(iv), uniqueN(iv$gene_symbol)))
}
if (iv[, .N, by = .(gene_symbol, SNP)][N > 1L, .N] > 0L) {
  stop("Frozen protein-SNP pairs are not unique.")
}

requested_snps <- unique(c(iv$SNP, "rs429358", "rs7412"))
requested_cache <- file.path(processed_dir, "Kunkle_2019_requested_variants.tsv")
if (file.exists(requested_cache)) {
  message("Reading the verified Kunkle requested-variant cache.")
  kunkle <- fread(requested_cache, na.strings = c("NA", ""))
} else {
  message("Reading Kunkle Stage 1 summary statistics and retaining ", length(requested_snps), " requested rsIDs.")
  kunkle <- fread(
    raw_file,
    select = c("Chromosome", "Position", "MarkerName", "Effect_allele", "Non_Effect_allele", "Beta", "SE", "Pvalue"),
    showProgress = TRUE
  )
  kunkle <- kunkle[MarkerName %chin% requested_snps]
  kunkle[, `:=`(
    Chromosome = as.integer(Chromosome),
    Position = as.integer(Position),
    Beta = as.numeric(Beta),
    SE = as.numeric(SE),
    Pvalue = as.numeric(Pvalue)
  )]
  setorder(kunkle, MarkerName, Pvalue)
  kunkle <- kunkle[, .SD[1], by = MarkerName]
  setnames(kunkle, "MarkerName", "SNP")
  fwrite(kunkle, requested_cache, sep = "\t", na = "NA")
}
stopifnot(all(c("SNP", "Chromosome", "Position", "Effect_allele", "Non_Effect_allele", "Beta", "SE", "Pvalue") %in% names(kunkle)))

format_exposure <- function(dat, gene) {
  out <- format_data(
    as.data.frame(dat), type = "exposure",
    snp_col = "SNP", beta_col = "beta_exposure", se_col = "se_exposure",
    effect_allele_col = "effect_allele", other_allele_col = "other_allele",
    eaf_col = "eaf_exposure", pval_col = "pval_exposure",
    samplesize_col = "samplesize_exposure", chr_col = "chr_exposure", pos_col = "pos_exposure"
  )
  out$exposure <- gene
  out$id.exposure <- gene
  out
}

format_outcome <- function(dat) {
  dat <- copy(dat)
  dat[, samplesize := 21982 + 41944]
  out <- format_data(
    as.data.frame(dat), type = "outcome",
    snp_col = "SNP", beta_col = "Beta", se_col = "SE",
    effect_allele_col = "Effect_allele", other_allele_col = "Non_Effect_allele",
    pval_col = "Pvalue", samplesize_col = "samplesize",
    chr_col = "Chromosome", pos_col = "Position"
  )
  out$outcome <- "AD_Kunkle_2019_clinically_diagnosed"
  out$id.outcome <- "AD_Kunkle_2019_clinically_diagnosed"
  out
}

result_rows <- list()
harm_rows <- list()
coverage_rows <- list()

for (gene in sort(unique(iv$gene_symbol))) {
  gene_iv <- iv[gene_symbol == gene]
  present <- gene_iv[SNP %chin% kunkle$SNP]
  if (!nrow(present)) {
    coverage_rows[[gene]] <- data.table(
      gene_symbol = gene, prespecified_instrument_pairs = nrow(gene_iv), present_in_Kunkle = 0L,
      retained_after_harmonization = 0L, coverage_fraction = 0,
      status = "not_estimable", exclusion_reason = "No frozen instruments were present in Kunkle Stage 1."
    )
    next
  }

  exposure <- format_exposure(present, gene)
  outcome <- format_outcome(kunkle[SNP %chin% present$SNP])
  harmonized <- harmonise_data(exposure, outcome, action = 3)
  harmonized <- harmonized[harmonized$mr_keep %in% TRUE, ]

  coverage_rows[[gene]] <- data.table(
    gene_symbol = gene,
    prespecified_instrument_pairs = nrow(gene_iv),
    present_in_Kunkle = nrow(present),
    retained_after_harmonization = nrow(harmonized),
    coverage_fraction = nrow(present) / nrow(gene_iv),
    status = if (nrow(harmonized)) "estimable" else "not_estimable",
    exclusion_reason = if (nrow(harmonized)) "NA" else "No instruments remained after conservative allele harmonization (palindromic variants removed)."
  )

  if (!nrow(harmonized)) next
  harmonized$F_statistic <- (harmonized$beta.exposure / harmonized$se.exposure)^2
  harm_rows[[gene]] <- as.data.table(harmonized)

  methods <- if (nrow(harmonized) == 1L) {
    "mr_wald_ratio"
  } else {
    c("mr_ivw", "mr_weighted_median", "mr_egger_regression")
  }
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
    outcome = "AD_Kunkle_2019_clinically_diagnosed",
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
    mean_F = mean(harmonized$F_statistic),
    min_F = min(harmonized$F_statistic),
    beta_status = "reestimated_clinically_diagnosed_AD_sensitivity",
    beta_source = "Kunkle_2019_Stage1",
    harmonization_rule = "TwoSampleMR action=3; palindromic variants removed because outcome EAF was unavailable"
  )]
  setnames(estimates, c("b", "se", "pval"), c("beta", "SE", "P_value"), skip_absent = TRUE)
  result_rows[[gene]] <- estimates
}

coverage <- rbindlist(coverage_rows, fill = TRUE)
fwrite(coverage, file.path(table_dir, "Kunkle_AD_instrument_coverage.tsv"), sep = "\t", na = "NA")
if (length(harm_rows)) {
  fwrite(rbindlist(harm_rows, fill = TRUE, idcol = "gene_symbol"),
         file.path(processed_dir, "Kunkle_AD_harmonized_protein_instruments.tsv"), sep = "\t", na = "NA")
}

beta <- rbindlist(result_rows, fill = TRUE)
primary_idx <- which(beta$method_role == "primary")
beta[, `:=`(
  planned_family_size = 25L,
  observed_primary_tests = length(primary_idx),
  P_FDR_observed = NA_real_,
  P_Bonferroni_planned = NA_real_
)]
if (length(primary_idx)) {
  beta[primary_idx, P_FDR_observed := p.adjust(P_value, method = "BH")]
  beta[primary_idx, P_Bonferroni_planned := pmin(1, P_value * 25)]
}

frozen_beta <- fread(frozen_beta_file, na.strings = c("NA", ""))
assay_meta <- unique(frozen_beta[
  outcome == "AD" & method_role == "primary" & analysis_set == "genome_wide_instruments_primary",
  .(gene_symbol, assay_target_ID, UniProt_ID)
])
beta <- merge(beta, assay_meta, by = "gene_symbol", all.x = TRUE)
setcolorder(beta, c("gene_symbol", "assay_target_ID", "UniProt_ID", setdiff(names(beta), c("gene_symbol", "assay_target_ID", "UniProt_ID"))))
fwrite(beta, file.path(table_dir, "Kunkle_AD_protein_beta_sensitivity.tsv"), sep = "\t", na = "NA")

complement <- function(x) chartr("ACGT", "TGCA", toupper(x))
align_total <- function(variant, requested_effect, requested_other) {
  row <- kunkle[SNP == variant]
  if (!nrow(row)) {
    return(data.table(variant = variant, availability_status = "variant_absent", exclusion_reason = "Variant not present in Kunkle Stage 1."))
  }
  ea <- toupper(row$Effect_allele[1])
  oa <- toupper(row$Non_Effect_allele[1])
  req_ea <- toupper(requested_effect)
  req_oa <- toupper(requested_other)
  orientation <- if (ea == req_ea && oa == req_oa) {
    1
  } else if (ea == req_oa && oa == req_ea) {
    -1
  } else if (complement(ea) == req_ea && complement(oa) == req_oa) {
    1
  } else if (complement(ea) == req_oa && complement(oa) == req_ea) {
    -1
  } else {
    NA_real_
  }
  if (is.na(orientation)) {
    return(data.table(variant = variant, availability_status = "alleles_unresolved", exclusion_reason = "Kunkle alleles could not be aligned to the prespecified APOE effect allele."))
  }
  data.table(
    variant = variant,
    outcome = "AD_Kunkle_2019_clinically_diagnosed",
    requested_effect_allele = req_ea,
    requested_other_allele = req_oa,
    original_effect_allele = ea,
    original_other_allele = oa,
    allele_flipped = orientation == -1,
    beta = orientation * row$Beta[1],
    SE = row$SE[1],
    P_value = row$Pvalue[1],
    sample_size = 21982 + 41944,
    outcome_GWAS = "Kunkle_2019_Stage1_clinically_diagnosed_AD",
    availability_status = "direct_variant_available",
    exclusion_reason = "NA"
  )
}

totals <- rbindlist(list(
  align_total("rs429358", "C", "T"),
  align_total("rs7412", "T", "C")
), fill = TRUE)
fwrite(totals, file.path(table_dir, "Kunkle_AD_APOE_total_effects.tsv"), sep = "\t", na = "NA")

pair_coverage <- sum(iv$SNP %chin% kunkle$SNP) / nrow(iv)
n_estimable_proteins <- uniqueN(beta[method_role == "primary", gene_symbol])
both_totals <- all(c("rs429358", "rs7412") %chin% totals[availability_status == "direct_variant_available", variant])
full_gate_passed <- n_estimable_proteins >= 20L && pair_coverage >= 0.80 && both_totals

gate <- data.table(
  criterion = c("estimable_proteins", "frozen_protein_SNP_pair_coverage", "both_direct_APOE_total_effects", "full_mediation_gate"),
  observed = c(as.character(n_estimable_proteins), sprintf("%.6f", pair_coverage), as.character(both_totals), as.character(full_gate_passed)),
  threshold = c(">=20 of 25", ">=0.80 of 409 pairs", "TRUE", "all criteria TRUE"),
  passed = c(n_estimable_proteins >= 20L, pair_coverage >= 0.80, both_totals, full_gate_passed)
)
fwrite(gate, file.path(table_dir, "Kunkle_AD_analysis_gate.tsv"), sep = "\t")

if (!full_gate_passed) {
  writeLines(capture.output(sessionInfo()), file.path(log_dir, "Kunkle_AD_sensitivity_sessionInfo.txt"))
  stop("Kunkle full-mediation gate failed. Coverage and beta outputs were retained; mediation and aggregate were not computed.")
}

alpha <- fread(alpha_file, na.strings = c("NA", ""))
estimable_genes <- beta[method_role == "primary", gene_symbol]
alpha_use <- unique(alpha[
  availability_status == "direct_variant_available" &
    tolower(as.character(eligible_for_two_step_mapping)) == "true" &
    gene_symbol %chin% estimable_genes,
  .(gene_symbol, variant, alpha = as.numeric(beta), SE_alpha = as.numeric(SE), alpha_P = as.numeric(P_value),
    alpha_effect_allele = requested_effect_allele, alpha_assay_target_ID = assay_target_ID)
], by = c("gene_symbol", "variant"))
beta_use <- beta[method_role == "primary", .(
  gene_symbol, beta = as.numeric(beta), SE_beta = as.numeric(SE), beta_P = as.numeric(P_value),
  beta_method = method, beta_nsnp = nsnp, beta_Q_P = Q_P,
  beta_Egger_intercept_P = Egger_intercept_P, beta_assay_target_ID = assay_target_ID
)]
total_use <- totals[availability_status == "direct_variant_available", .(
  variant, total_effect = as.numeric(beta), SE_total = as.numeric(SE), total_effect_P = as.numeric(P_value),
  total_effect_allele = requested_effect_allele, total_effect_source = outcome_GWAS
)]

med <- merge(alpha_use, beta_use, by = "gene_symbol", allow.cartesian = TRUE)
med <- merge(med, total_use, by = "variant", all.x = TRUE)
if (anyNA(med$total_effect)) stop("Missing Kunkle APOE total effect after the full-analysis gate passed.")
if (!all(med$alpha_assay_target_ID == med$beta_assay_target_ID)) stop("Alpha/beta assay mismatch in Kunkle mediation.")

med[, `:=`(
  outcome = "AD_Kunkle_2019_clinically_diagnosed",
  indirect_effect = alpha * beta,
  SE_indirect_delta = sqrt(beta^2 * SE_alpha^2 + alpha^2 * SE_beta^2)
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

for (i in seq_len(nrow(med))) {
  a <- rnorm(n_boot, med$alpha[i], med$SE_alpha[i])
  b <- rnorm(n_boot, med$beta[i], med$SE_beta[i])
  t <- rnorm(n_boot, med$total_effect[i], med$SE_total[i])
  ind <- a * b
  prop <- ind / t
  med$indirect_CI_lower_bootstrap[i] <- quantile(ind, 0.025, na.rm = TRUE)
  med$indirect_CI_upper_bootstrap[i] <- quantile(ind, 0.975, na.rm = TRUE)
  med$mediated_proportion_CI_lower_bootstrap[i] <- quantile(prop, 0.025, na.rm = TRUE)
  med$mediated_proportion_CI_upper_bootstrap[i] <- quantile(prop, 0.975, na.rm = TRUE)
}
med[, `:=`(
  row_type = "protein",
  mediation_family_size = 50L,
  indirect_P_FDR = p.adjust(indirect_P_delta, method = "BH"),
  indirect_P_Bonferroni = pmin(1, indirect_P_delta * 50L),
  bootstrap_n = n_boot,
  analysis_label = "clinically_diagnosed_AD_outcome_sensitivity"
)]

aggregate_rows <- list()
for (variant_i in c("rs429358", "rs7412")) {
  subset <- med[variant == variant_i]
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
  aggregate_rows[[variant_i]] <- data.table(
    row_type = "total",
    gene_symbol = "TOTAL_ELIGIBLE_PROTEINS",
    variant = variant_i,
    outcome = "AD_Kunkle_2019_clinically_diagnosed",
    indirect_effect = total_indirect,
    SE_indirect_delta = total_se_delta,
    indirect_P_delta = 2 * pnorm(-abs(total_indirect / total_se_delta)),
    indirect_CI_lower_bootstrap = quantile(summed_indirect, 0.025),
    indirect_CI_upper_bootstrap = quantile(summed_indirect, 0.975),
    total_effect = subset$total_effect[1],
    SE_total = subset$SE_total[1],
    total_effect_P = subset$total_effect_P[1],
    total_effect_allele = subset$total_effect_allele[1],
    total_effect_source = subset$total_effect_source[1],
    mediated_proportion = total_indirect / subset$total_effect[1],
    mediated_proportion_CI_lower_bootstrap = quantile(summed_proportion, 0.025),
    mediated_proportion_CI_upper_bootstrap = quantile(summed_proportion, 0.975),
    direction_classification = fifelse(
      sign(total_indirect) == sign(subset$total_effect[1]),
      "concordant_total_mediation", "opposing_or_suppressing_total_mediation"
    ),
    number_of_eligible_proteins = nrow(subset),
    mediation_family_size = 50L,
    bootstrap_n = n_boot,
    analysis_label = "clinically_diagnosed_AD_outcome_sensitivity",
    bootstrap_assumption = "Independent normal alpha/beta draws and a common total-effect draw; cross-protein covariance, mapping uncertainty and literature-selection uncertainty not included."
  )
}

mediation <- rbindlist(list(med, rbindlist(aggregate_rows, fill = TRUE)), fill = TRUE)
setorder(mediation, variant, row_type, gene_symbol)
fwrite(mediation, file.path(table_dir, "Kunkle_AD_two_step_mediation_sensitivity.tsv"), sep = "\t", na = "NA")

current_beta <- frozen_beta[
  outcome == "AD" & method_role == "primary" & analysis_set == "genome_wide_instruments_primary",
  .(gene_symbol, Wightman_beta = as.numeric(beta), Wightman_SE = as.numeric(SE), Wightman_P = as.numeric(P_value), Wightman_nsnp = nsnp)
]
comparison <- merge(
  current_beta,
  beta[method_role == "primary", .(gene_symbol, Kunkle_beta = as.numeric(beta), Kunkle_SE = as.numeric(SE), Kunkle_P = as.numeric(P_value), Kunkle_nsnp = nsnp)],
  by = "gene_symbol", all = TRUE
)
comparison[, `:=`(
  same_direction = sign(Wightman_beta) == sign(Kunkle_beta),
  absolute_beta_difference = abs(Wightman_beta - Kunkle_beta)
)]
fwrite(comparison, file.path(table_dir, "Kunkle_AD_comparison_with_Wightman.tsv"), sep = "\t", na = "NA")

current_med <- fread(frozen_mediation_file, na.strings = c("NA", ""))
aggregate_comparison <- merge(
  current_med[outcome == "AD" & row_type == "total", .(
    variant,
    Wightman_number_of_proteins = number_of_eligible_proteins,
    Wightman_mediated_proportion = mediated_proportion,
    Wightman_CI_lower = mediated_proportion_CI_lower_bootstrap,
    Wightman_CI_upper = mediated_proportion_CI_upper_bootstrap
  )],
  mediation[row_type == "total", .(
    variant,
    Kunkle_number_of_proteins = number_of_eligible_proteins,
    Kunkle_mediated_proportion = mediated_proportion,
    Kunkle_CI_lower = mediated_proportion_CI_lower_bootstrap,
    Kunkle_CI_upper = mediated_proportion_CI_upper_bootstrap
  )],
  by = "variant", all = TRUE
)
fwrite(aggregate_comparison, file.path(table_dir, "Kunkle_AD_aggregate_comparison_with_Wightman.tsv"), sep = "\t", na = "NA")

run_summary <- data.table(
  metric = c(
    "frozen_protein_SNP_pairs", "unique_frozen_SNPs", "requested_SNPs_found_in_Kunkle",
    "pair_coverage_fraction", "estimable_proteins", "primary_beta_tests",
    "protein_mediation_paths", "aggregate_rows", "beta_direction_concordance_fraction"
  ),
  value = c(
    nrow(iv), uniqueN(iv$SNP), nrow(kunkle), pair_coverage, n_estimable_proteins,
    nrow(beta[method_role == "primary"]), nrow(med), nrow(mediation[row_type == "total"]),
    comparison[!is.na(same_direction), mean(same_direction)]
  )
)
fwrite(run_summary, file.path(table_dir, "Kunkle_AD_run_summary.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(log_dir, "Kunkle_AD_sensitivity_sessionInfo.txt"))

message("Kunkle clinically diagnosed AD sensitivity completed.")
message("Pair coverage: ", sprintf("%.1f%%", 100 * pair_coverage))
message("Estimable proteins: ", n_estimable_proteins)
message("Output directory: ", table_dir)
