# ============================================================
# 08_biology_category_models.R
# B: Biology-category aggregated genetic score analysis
# D: Penalized multivariable sensitivity model
# ============================================================

library(data.table)
library(TwoSampleMR)
library(dplyr)
library(glmnet)

set.seed(20240604)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this file with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = TRUE)
SCRIPT_DIR <- dirname(script_path)
PROJECT_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), winslash = "/", mustWork = TRUE)
RESOURCE_ROOT <- Sys.getenv("A1_RESOURCE_ROOT", unset = file.path(PROJECT_ROOT, "data", "external"))
RESULTS_DIR <- SCRIPT_DIR

# ---- Load exposures & outcomes ----
exposure_file <- file.path(PROJECT_ROOT, "04_protein_mr", "results", "ukbppp_exposures.rds")
if (!file.exists(exposure_file)) stop("Missing protein instrument object: ", exposure_file)
exposures <- readRDS(exposure_file)
cat(sprintf("Loaded %d proteins\n", length(exposures)))

GWAS_DIR <- file.path(RESOURCE_ROOT, "GWAS")
load_outcome <- function(path, trait_name) {
  dat <- fread(path, select=c("SNP","CHR","BP","A1","A2","FREQ","BETA","SE","P","N"))
  setnames(dat, c("CHR","BP","A1","A2","FREQ","BETA","SE","P","N"),
           c("chr.outcome","pos.outcome","effect_allele.outcome","other_allele.outcome",
             "eaf.outcome","beta.outcome","se.outcome","pval.outcome","samplesize.outcome"))
  dat$outcome <- trait_name; dat$id.outcome <- trait_name
  return(dat)
}

outcomes <- list(
  AD      = load_outcome(file.path(GWAS_DIR, "AD_Wightman_cleaned_hg19.tsv.gz"), "AD"),
  Dry_AMD = load_outcome(file.path(GWAS_DIR, "AMD_Dry_R12_cleaned_hg19.tsv.gz"), "Dry_AMD"),
  Wet_AMD = load_outcome(file.path(GWAS_DIR, "AMD_Wet_R12_cleaned_hg19.tsv.gz"), "Wet_AMD"),
  Any_AMD = load_outcome(file.path(GWAS_DIR, "AMD_H7_R12_cleaned_hg19.tsv.gz"), "Any_AMD")
)
cat("Outcomes loaded\n")

# ---- Biology category definitions (excluding APOE) ----
biology_categories <- list(
  Complement   = c("C1QA","C2","C3","C5","CFB","CFD","CFH","CFI","CFP","SERPING1"),
  Inflammatory = c("CSF1","IFNG","IL10","IL18","IL1B","IL6","TGFB1","TNF"),
  Chemokine    = c("CCL2","CCL5","CX3CL1","CXCL10","CXCL12"),
  Immune       = c("CD14","CD40","TLR4","TREM2"),
  Lipid        = c("LPA","PON1")
)

# ============================================================
# B: BIOLOGY-CATEGORY AGGREGATED GENETIC SCORE ANALYSIS
# ============================================================
cat("\n", paste(rep("=",60), collapse=""), "\n")
cat("B: BIOLOGY-CATEGORY AGGREGATED GENETIC SCORE ANALYSIS\n")
cat(paste(rep("=",60), collapse=""), "\n")

# For each biology category, build a single aggregated exposure from all member proteins.
build_biology_category_exposure <- function(pw_genes, exposures_list) {
  pw_exp <- exposures_list[names(exposures_list) %in% pw_genes]
  if (length(pw_exp) == 0) return(NULL)

  # Collect all unique SNPs across the biology category.
  all_snps <- unique(unlist(lapply(pw_exp, function(x) x$SNP)))

  # For each SNP, compute IVW meta-analysis across proteins that have it
  snp_data <- data.frame(SNP = all_snps, stringsAsFactors = FALSE)
  snp_data$chr.exposure <- NA_integer_
  snp_data$pos.exposure <- NA_integer_
  snp_data$effect_allele.exposure <- NA_character_
  snp_data$other_allele.exposure <- NA_character_
  snp_data$beta_agg <- 0
  snp_data$w_agg <- 0
  snp_data$n_prot <- 0

  for (gene in names(pw_exp)) {
    exp <- pw_exp[[gene]]
    idx <- match(exp$SNP, all_snps)
    w <- 1 / exp$se.exposure^2
    snp_data$beta_agg[idx] <- snp_data$beta_agg[idx] + exp$beta.exposure * w
    snp_data$w_agg[idx]    <- snp_data$w_agg[idx] + w
    snp_data$n_prot[idx]   <- snp_data$n_prot[idx] + 1

    # Fill metadata from first occurrence
    na_idx <- which(is.na(snp_data$chr.exposure[idx]))
    snp_data$chr.exposure[idx[na_idx]] <- exp$chr.exposure[na_idx]
    snp_data$pos.exposure[idx[na_idx]] <- exp$pos.exposure[na_idx]
    snp_data$effect_allele.exposure[idx[na_idx]] <- exp$effect_allele.exposure[na_idx]
    snp_data$other_allele.exposure[idx[na_idx]] <- exp$other_allele.exposure[na_idx]
  }

  # IVW pooled beta = sum(beta*w) / sum(w)
  snp_data$beta.exposure <- snp_data$beta_agg / snp_data$w_agg
  snp_data$se.exposure   <- 1 / sqrt(snp_data$w_agg)
  snp_data$pval.exposure <- 2 * pnorm(-abs(snp_data$beta.exposure / snp_data$se.exposure))
  snp_data$SNP <- snp_data$SNP  # keep rsID
  snp_data$id.exposure <- paste0("BC_", paste(pw_genes, collapse="_"))
  snp_data$exposure <- paste0("BC_", paste(pw_genes, collapse="_"))

  # Retain variants represented in at least one member-protein instrument set.
  snp_data <- snp_data[snp_data$n_prot >= 1 & snp_data$se.exposure < 1, ]

  return(snp_data)
}

# Run the aggregated genetic score analysis for each biology category and outcome
b_results <- list()

for (pw_name in names(biology_categories)) {
  pw_genes <- biology_categories[[pw_name]]
  pw_genes_avail <- pw_genes[pw_genes %in% names(exposures)]
  if (length(pw_genes_avail) < 2) next

  cat(sprintf("\n--- %s (%d proteins) ---\n", pw_name, length(pw_genes_avail)))

  pw_exp <- build_biology_category_exposure(pw_genes_avail, exposures)
  if (is.null(pw_exp) || nrow(pw_exp) < 3) { cat("  Failed to build exposure\n"); next }

  # Ensure chr/pos are integer
  pw_exp$chr.exposure <- as.integer(pw_exp$chr.exposure)
  pw_exp$pos.exposure <- as.integer(pw_exp$pos.exposure)

  for (out_name in names(outcomes)) {
    harm <- tryCatch(
      harmonise_data(pw_exp, outcomes[[out_name]], action = 2),
      error = function(e) NULL)
    if (is.null(harm) || nrow(harm) < 3) next

    # MR
    mr_res <- tryCatch(
      mr(harm, method_list = c("mr_ivw", "mr_weighted_median", "mr_egger_regression")),
      error = function(e) NULL)
    if (is.null(mr_res)) next

    ivw_row <- mr_res[mr_res$method == "Inverse variance weighted", ]
    if (nrow(ivw_row) == 0) next

    # F-statistic
    mean_f <- mean((harm$beta.exposure / harm$se.exposure)^2, na.rm = TRUE)

    b_results[[length(b_results) + 1]] <- data.frame(
      biology_category       = pw_name,
      outcome       = out_name,
      n_proteins    = length(pw_genes_avail),
      n_ivs         = nrow(harm),
      mean_f        = mean_f,
      beta          = ivw_row$b,
      se            = ivw_row$se,
      pval          = ivw_row$pval,
      stringsAsFactors = FALSE
    )

    sig_mark <- if(ivw_row$pval < 0.05/5) "***" else if(ivw_row$pval < 0.05) "*" else ""
    cat(sprintf("  %-8s: %2d IVs  F=%.0f  b=%.4f  se=%.4f  P=%.2e %s\n",
        out_name, nrow(harm), mean_f, ivw_row$b, ivw_row$se, ivw_row$pval, sig_mark))
  }
}

b_all <- rbindlist(b_results)
fwrite(b_all, file.path(RESULTS_DIR, "E2_biology_category_aggregated_score.csv"))

# Bonferroni summary
cat("\n=== B: Bonferroni-significant biology categories (P<0.05/5=0.01) ===\n")
b_sig <- b_all[pval < 0.01, ][order(pval), ]
if (nrow(b_sig) > 0) {
  for (i in 1:nrow(b_sig))
    cat(sprintf("  %-14s %-8s b=%.4f P=%.2e F=%.0f\n",
        b_sig$biology_category[i], b_sig$outcome[i], b_sig$beta[i], b_sig$pval[i], b_sig$mean_f[i]))
} else {
  cat("  None at Bonferroni threshold. Nominally significant:\n")
  b_nom <- b_all[pval < 0.05, ][order(pval), ]
  for (i in 1:nrow(b_nom))
    cat(sprintf("  %-14s %-8s b=%.4f P=%.2e\n",
        b_nom$biology_category[i], b_nom$outcome[i], b_nom$beta[i], b_nom$pval[i]))
}

# ============================================================
# D: PENALIZED MULTIVARIABLE SENSITIVITY MODEL
# ============================================================
cat("\n", paste(rep("=",60), collapse=""), "\n")
cat("D: PENALIZED MULTIVARIABLE SENSITIVITY MODEL\n")
cat(paste(rep("=",60), collapse=""), "\n")

ridge_results <- list()

for (pw_name in names(biology_categories)) {
  pw_genes <- biology_categories[[pw_name]]
  pw_exp <- exposures[names(exposures) %in% pw_genes]
  if (length(pw_exp) < 2) next

  # Collect all unique SNPs
  all_snps <- unique(unlist(lapply(pw_exp, function(x) x$SNP)))
  n_prot <- length(pw_exp)
  gene_names <- names(pw_exp)

  # Build exposure matrix: rows=SNPs, cols=proteins
  bx <- matrix(0, nrow = length(all_snps), ncol = n_prot)
  colnames(bx) <- gene_names
  for (j in seq_along(gene_names)) {
    exp <- pw_exp[[gene_names[j]]]
    idx <- match(all_snps, exp$SNP)
    bx[, j] <- ifelse(is.na(idx), 0, exp$beta.exposure[idx])
  }

  for (out_name in names(outcomes)) {
    template <- pw_exp[[1]]
    harm <- tryCatch(
      harmonise_data(template[template$SNP %in% all_snps, ], outcomes[[out_name]], action = 2),
      error = function(e) NULL)
    if (is.null(harm)) next

    by <- harm$beta.outcome[match(all_snps, harm$SNP)]
    se_by <- harm$se.outcome[match(all_snps, harm$SNP)]

    valid <- !is.na(by) & !is.na(se_by) & se_by > 0 & rowSums(abs(bx)) > 0
    if (sum(valid) < n_prot + 1) next

    bx_v <- bx[valid, , drop = FALSE]
    by_v <- by[valid]
    se_by_v <- se_by[valid]
    w <- 1 / se_by_v^2
    bx_w <- bx_v * sqrt(w)
    by_w <- by_v * sqrt(w)

    # Ridge with cross-validation
    fit <- tryCatch(
      cv.glmnet(bx_w, by_w, alpha = 0, standardize = TRUE, nfolds = min(10, nrow(bx_w))),
      error = function(e) NULL)
    if (is.null(fit)) next

    # Coefficients at lambda.min
    coefs <- coef(fit, s = "lambda.min")[-1, 1]  # drop intercept

    # Bootstrap SE
    boot_betas <- matrix(NA, nrow = 200, ncol = n_prot)
    colnames(boot_betas) <- gene_names
    n <- nrow(bx_w)
    for (b in 1:200) {
      boot_idx <- sample(n, replace = TRUE)
      bfit <- tryCatch(
        glmnet(bx_w[boot_idx, ], by_w[boot_idx], alpha = 0, lambda = fit$lambda.min, standardize = TRUE),
        error = function(e) NULL)
      if (!is.null(bfit)) {
        bcoef <- coef(bfit)[-1, 1]
        boot_betas[b, ] <- bcoef
      }
    }
    se_ridge <- apply(boot_betas, 2, sd, na.rm = TRUE)

    for (j in seq_along(gene_names)) {
      ridge_results[[length(ridge_results) + 1]] <- data.frame(
        biology_category   = pw_name,
        outcome   = out_name,
        gene      = gene_names[j],
        n_snps    = nrow(bx_v),
        n_proteins = n_prot,
        beta_ridge = coefs[j],
        se_ridge   = se_ridge[j],
        pval_ridge = 2 * pnorm(-abs(coefs[j] / se_ridge[j])),
        lambda     = fit$lambda.min,
        stringsAsFactors = FALSE
      )
    }

    sig_n <- sum(2 * pnorm(-abs(coefs / se_ridge)) < 0.05, na.rm = TRUE)
    cat(sprintf("  %-14s %-8s: %d SNPs, %d/%d sig, lambda=%.3f\n",
        pw_name, out_name, nrow(bx_v), sig_n, n_prot, fit$lambda.min))
  }
}

r_all <- rbindlist(ridge_results)
fwrite(r_all, file.path(RESULTS_DIR, "E2_penalized_multivariable_sensitivity.csv"))

# Ridge summary
cat("\n=== D: Nominally significant penalized coefficients (P<0.05) ===\n")
r_sig <- r_all[pval_ridge < 0.05, ][order(pval_ridge), ]
if (nrow(r_sig) > 0) {
  for (i in 1:min(20, nrow(r_sig)))
    cat(sprintf("  %-12s %-8s %-14s b=%.4f se=%.4f P=%.2e\n",
        r_sig$gene[i], r_sig$outcome[i], r_sig$biology_category[i],
        r_sig$beta_ridge[i], r_sig$se_ridge[i], r_sig$pval_ridge[i]))
} else {
  cat("  None; robust to regularization\n")
}

cat("\nE2 complete: biology-category aggregated score analysis and penalized multivariable sensitivity\n")
