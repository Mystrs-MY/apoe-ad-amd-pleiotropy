#!/usr/bin/env Rscript

rm(list = ls())

required <- c("data.table", "TwoSampleMR", "plinkbinr")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R package(s): ", paste(missing, collapse = ", "))

suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
})

set.seed(20260717)

script_arg <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)][1])
SCRIPT_DIR <- dirname(normalizePath(script_file, winslash = "/", mustWork = TRUE))
ROOT <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), winslash = "/", mustWork = TRUE)
UPGRADE <- file.path(ROOT, "A1_protein_upgrade")
OUT_DIR <- file.path(UPGRADE, "tables", "targeted_submission_sensitivities")
LOG_DIR <- file.path(UPGRADE, "logs")
SUBMISSION_DIR <- file.path(ROOT, "tables_submission", "supplementary_tables")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SUBMISSION_DIR, recursive = TRUE, showWarnings = FALSE)

copy_submission <- function(source, target_name) {
  target <- file.path(SUBMISSION_DIR, target_name)
  ok <- file.copy(source, target, overwrite = TRUE)
  if (!ok) stop("Failed to copy submission table: ", target)
  target
}

# -----------------------------------------------------------------------------
# 1. Cross-outcome global correction for the primary protein-beta family.
# -----------------------------------------------------------------------------
beta_file <- file.path(SUBMISSION_DIR, "TableS23_Same_Assay_Protein_Beta.tsv")
if (!file.exists(beta_file)) stop("Missing Table S23: ", beta_file)
beta <- fread(beta_file)

beta[, `:=`(
  global_planned_family_size = NA_integer_,
  global_observed_primary_tests = NA_integer_,
  P_FDR_global_observed_primary = NA_real_,
  P_Bonferroni_global_planned_108 = NA_real_,
  corrected_significance_FDR_global_observed_primary = NA,
  corrected_significance_Bonferroni_global_planned_108 = NA
)]

primary_idx <- which(
  beta$analysis_set == "genome_wide_instruments_primary" &
    beta$method_role == "primary" &
    beta$beta_status == "reestimated" &
    is.finite(beta$P_value)
)
if (length(primary_idx) != 100L) {
  stop("Expected 100 estimable primary protein-outcome tests, found ", length(primary_idx))
}
beta[primary_idx, global_planned_family_size := 27L * 4L]
beta[primary_idx, global_observed_primary_tests := length(primary_idx)]
beta[primary_idx, P_FDR_global_observed_primary := p.adjust(P_value, method = "BH")]
beta[primary_idx, P_Bonferroni_global_planned_108 := pmin(1, P_value * 108)]
beta[primary_idx, corrected_significance_FDR_global_observed_primary :=
       P_FDR_global_observed_primary < 0.05]
beta[primary_idx, corrected_significance_Bonferroni_global_planned_108 :=
       P_Bonferroni_global_planned_108 < 0.05]

fwrite(beta, beta_file, sep = "\t", na = "NA")

global_summary <- beta[primary_idx, .(
  gene_symbol, assay_target_ID, outcome, method, nsnp, beta, SE, P_value,
  P_FDR_per_outcome_observed = P_FDR_observed_reestimated,
  P_Bonferroni_per_outcome_planned_27 = P_Bonferroni_planned_family,
  P_FDR_global_observed_primary,
  P_Bonferroni_global_planned_108,
  corrected_significance_FDR_global_observed_primary,
  corrected_significance_Bonferroni_global_planned_108
)][order(P_value)]
global_file <- file.path(OUT_DIR, "protein_beta_global_108_sensitivity.tsv")
fwrite(global_summary, global_file, sep = "\t", na = "NA")

# -----------------------------------------------------------------------------
# 2. APOE exclusion-window sensitivity using the frozen harmonised IV sets.
# -----------------------------------------------------------------------------
mr_rds <- file.path(ROOT, "03_causal_lock", "gw_mr_with_apoe.rds")
if (!file.exists(mr_rds)) stop("Missing harmonised genome-wide MR object: ", mr_rds)
mr_objects <- readRDS(mr_rds)

APOE_GENE_START <- 45409039L
APOE_GENE_END <- 45412650L
windows <- data.table(
  window_id = c("none", "APOE_gene_plus_minus_500kb", "APOE_gene_plus_minus_1Mb",
                "prespecified_chr19_44_46_5Mb", "APOE_gene_plus_minus_2Mb"),
  start_bp = c(NA, APOE_GENE_START - 500000L, APOE_GENE_START - 1000000L,
               44000000L, APOE_GENE_START - 2000000L),
  end_bp = c(NA, APOE_GENE_END + 500000L, APOE_GENE_END + 1000000L,
             46500000L, APOE_GENE_END + 2000000L),
  role = c("baseline", rep("window_sensitivity", 4))
)

window_rows <- list()
for (pair_name in names(mr_objects)) {
  dat0 <- as.data.table(mr_objects[[pair_name]]$dat)
  if (!nrow(dat0)) next
  for (j in seq_len(nrow(windows))) {
    w <- windows[j]
    dat <- copy(dat0)
    if (w$window_id != "none") {
      dat <- dat[!(chr.exposure == 19L & pos.exposure >= w$start_bp & pos.exposure <= w$end_bp)]
    }
    if (nrow(dat) >= 3L) {
      fit <- TwoSampleMR::mr(dat, method_list = "mr_ivw")
      fit <- fit[fit$method == "Inverse variance weighted", ]
      beta_ivw <- fit$b[1]
      se_ivw <- fit$se[1]
      p_ivw <- fit$pval[1]
    } else {
      beta_ivw <- se_ivw <- p_ivw <- NA_real_
    }
    window_rows[[length(window_rows) + 1L]] <- data.table(
      pair_id = pair_name,
      exposure = unique(dat0$exposure)[1],
      outcome = unique(dat0$outcome)[1],
      window_id = w$window_id,
      window_role = w$role,
      chromosome = ifelse(w$window_id == "none", NA_integer_, 19L),
      exclusion_start_bp = w$start_bp,
      exclusion_end_bp = w$end_bp,
      n_iv_before_exclusion = nrow(dat0),
      n_iv_after_exclusion = nrow(dat),
      n_iv_removed = nrow(dat0) - nrow(dat),
      mean_F_after_exclusion = if (nrow(dat)) mean((dat$beta.exposure / dat$se.exposure)^2) else NA_real_,
      IVW_beta = beta_ivw,
      IVW_SE = se_ivw,
      IVW_P = p_ivw,
      nominal_P_lt_0_05 = is.finite(p_ivw) && p_ivw < 0.05,
      interpretation = ifelse(
        w$window_id == "none",
        "Baseline harmonised and clumped instrument set.",
        "Sensitivity to the APOE-region exclusion boundary; not a new discovery analysis."
      )
    )
  }
}
window_result <- rbindlist(window_rows, fill = TRUE)
window_file <- file.path(OUT_DIR, "TableS2c_APOE_Exclusion_Window_Sensitivity.tsv")
fwrite(window_result, window_file, sep = "\t", na = "NA")
copy_submission(window_file, basename(window_file))

# -----------------------------------------------------------------------------
# 3. Shared and high-LD pQTL instrument audit across proteins.
# -----------------------------------------------------------------------------
harmonized_file <- file.path(UPGRADE, "data_processed", "literature_panel_harmonized_data.tsv")
if (!file.exists(harmonized_file)) stop("Missing harmonised protein MR data: ", harmonized_file)
harm <- fread(harmonized_file)
iv <- unique(harm[mr_keep == TRUE, .(
  protein = exposure,
  SNP,
  chromosome = chr.exposure,
  position_bp = pos.exposure,
  cis_trans_status = NA_character_
)], by = c("protein", "SNP"))

exact <- iv[, .(
  proteins = paste(sort(unique(protein)), collapse = ";"),
  protein_count = uniqueN(protein),
  chromosome = chromosome[1],
  position_bp = position_bp[1]
), by = SNP][protein_count > 1L]
exact[, `:=`(
  row_type = "exact_shared_SNP",
  SNP_A = SNP,
  SNP_B = SNP,
  R2 = 1,
  proteins_A = proteins,
  proteins_B = proteins,
  audit_note = "The same pQTL contributes to more than one protein model; aggregate estimates may reuse a genetic perturbation."
)]
exact_out <- exact[, .(
  row_type, SNP_A, SNP_B, chromosome_A = chromosome, position_A_bp = position_bp,
  chromosome_B = chromosome, position_B_bp = position_bp, R2,
  proteins_A, proteins_B, protein_count_A = protein_count,
  protein_count_B = protein_count, audit_note
)]

plink <- plinkbinr::get_plink_exe()
resource_root <- Sys.getenv("A1_RESOURCE_ROOT", unset = file.path(ROOT, "data", "external"))
ld_prefix <- Sys.getenv("A1_LD_PREFIX", unset = file.path(resource_root, "EUR", "EUR"))
ld_work <- file.path(UPGRADE, "data_processed", "shared_instrument_ld_audit")
dir.create(ld_work, recursive = TRUE, showWarnings = FALSE)
extract_file <- file.path(ld_work, "protein_instrument_snps.txt")
writeLines(sort(unique(iv$SNP)), extract_file)
plink_prefix <- file.path(ld_work, "protein_instrument_pairs_r2_0_8")
ld_file <- paste0(plink_prefix, ".ld.gz")
if (!file.exists(ld_file)) {
  args <- c(
    "--bfile", ld_prefix,
    "--extract", extract_file,
    "--r2", "gz",
    "--ld-window", "99999",
    "--ld-window-kb", "1000000",
    "--ld-window-r2", "0.8",
    "--out", plink_prefix
  )
  status <- system2(plink, args = args, stdout = TRUE, stderr = TRUE)
  writeLines(status, file.path(ld_work, "plink_r2_console.log"))
  if (!file.exists(ld_file)) stop("PLINK LD audit failed; inspect ", ld_work)
}

ld <- fread(ld_file)
if (nrow(ld)) {
  gene_map <- iv[, .(proteins = paste(sort(unique(protein)), collapse = ";"),
                     protein_count = uniqueN(protein)), by = SNP]
  ld <- merge(ld, gene_map, by.x = "SNP_A", by.y = "SNP", all.x = TRUE)
  setnames(ld, c("proteins", "protein_count"), c("proteins_A", "protein_count_A"))
  ld <- merge(ld, gene_map, by.x = "SNP_B", by.y = "SNP", all.x = TRUE)
  setnames(ld, c("proteins", "protein_count"), c("proteins_B", "protein_count_B"))
  ld <- ld[!is.na(proteins_A) & !is.na(proteins_B)]
  ld[, row_type := "distinct_SNP_pair_in_high_LD"]
  ld[, audit_note := "Distinct protein instruments have R2 >= 0.80 in the 1000 Genomes EUR reference; aggregate reuse is possible."]
  ld_out <- ld[, .(
    row_type, SNP_A, SNP_B,
    chromosome_A = CHR_A, position_A_bp = BP_A,
    chromosome_B = CHR_B, position_B_bp = BP_B,
    R2, proteins_A, proteins_B, protein_count_A, protein_count_B, audit_note
  )]
} else {
  ld_out <- exact_out[0]
}

shared_audit <- rbindlist(list(exact_out, ld_out), fill = TRUE)
setorder(shared_audit, row_type, -R2, SNP_A, SNP_B)
shared_file <- file.path(OUT_DIR, "TableS11b_Shared_Protein_Instrument_Audit.tsv")
fwrite(shared_audit, shared_file, sep = "\t", na = "NA")
copy_submission(shared_file, basename(shared_file))

# -----------------------------------------------------------------------------
# 4. Leave-one-protein-out influence analysis for the eight aggregate estimates.
# -----------------------------------------------------------------------------
med_file <- file.path(
  UPGRADE, "tables", "supplementary_name_match_revision",
  "TableS25_Expanded_Primary_Two_Step_Mediation.tsv"
)
if (!file.exists(med_file)) stop("Missing primary mediation table: ", med_file)
med <- fread(med_file)
proteins <- med[row_type == "protein" & is.finite(indirect_effect)]
if (proteins[, uniqueN(gene_symbol)] != 25L) stop("Expected 25 eligible proteins for LOO analysis.")

loo_rows <- list()
loo_groups <- unique(proteins[, .(variant, outcome)])
for (group_index in seq_len(nrow(loo_groups))) {
  key <- loo_groups[group_index]
  dat <- proteins[variant == key$variant & outcome == key$outcome]
  full_indirect <- sum(dat$indirect_effect)
  total_effect <- unique(dat$total_effect)
  if (length(total_effect) != 1L || !is.finite(total_effect)) stop("Invalid total effect in LOO group.")
  full_prop <- full_indirect / total_effect
  for (gene in sort(dat$gene_symbol)) {
    loo_indirect <- sum(dat[gene_symbol != gene]$indirect_effect)
    loo_prop <- loo_indirect / total_effect
    loo_rows[[length(loo_rows) + 1L]] <- data.table(
      row_type = "leave_one_protein_out",
      variant = key$variant,
      outcome = key$outcome,
      omitted_protein = gene,
      n_proteins_full = nrow(dat),
      n_proteins_after_omission = nrow(dat) - 1L,
      full_indirect_effect = full_indirect,
      leave_one_out_indirect_effect = loo_indirect,
      absolute_change_indirect_effect = loo_indirect - full_indirect,
      full_mediated_proportion = full_prop,
      leave_one_out_mediated_proportion = loo_prop,
      absolute_change_mediated_proportion = loo_prop - full_prop,
      sign_changed = sign(loo_prop) != sign(full_prop),
      influence_rank = NA_integer_,
      interpretation = "Influence diagnostic only; no confidence interval and no new inferential test."
    )
  }
}
loo <- rbindlist(loo_rows)
loo[, influence_rank := frank(-abs(absolute_change_mediated_proportion), ties.method = "min"),
    by = .(variant, outcome)]
setorder(loo, variant, outcome, influence_rank, omitted_protein)
loo_file <- file.path(OUT_DIR, "TableS11c_Leave_One_Protein_Out_Aggregate.tsv")
fwrite(loo, loo_file, sep = "\t", na = "NA")
copy_submission(loo_file, basename(loo_file))

summary_file <- file.path(OUT_DIR, "targeted_sensitivity_summary.tsv")
summary <- rbindlist(list(
  data.table(
    analysis = "global_108_protein_beta_correction",
    n_rows = nrow(global_summary),
    key_result = paste0(
      sum(global_summary$corrected_significance_Bonferroni_global_planned_108),
      " of 100 estimable tests passed the planned 108-test Bonferroni sensitivity."
    )
  ),
  data.table(
    analysis = "APOE_exclusion_window_sensitivity",
    n_rows = nrow(window_result),
    key_result = paste0(
      sum(window_result[window_role == "window_sensitivity"]$nominal_P_lt_0_05, na.rm = TRUE),
      " nominally significant estimates across 24 exclusion-window sensitivity rows."
    )
  ),
  data.table(
    analysis = "shared_protein_instrument_audit",
    n_rows = nrow(shared_audit),
    key_result = paste0(nrow(exact_out), " exact shared SNPs and ", nrow(ld_out), " distinct high-LD SNP pairs recorded.")
  ),
  data.table(
    analysis = "leave_one_protein_out_aggregate",
    n_rows = nrow(loo),
    key_result = paste0(sum(loo$sign_changed), " of ", nrow(loo), " omissions changed the aggregate sign.")
  )
), fill = TRUE)
fwrite(summary, summary_file, sep = "\t", na = "NA")
writeLines(capture.output(sessionInfo()), file.path(LOG_DIR, "targeted_submission_sensitivities_sessionInfo.txt"))

cat("Targeted submission sensitivity analyses completed.\n")
print(summary)
