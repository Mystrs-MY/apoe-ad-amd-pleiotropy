# ============================================================
# Biology-guided two-step mediation sensitivity analysis
# rs429358 -> protein -> outcome
# α = rs429358->protein (single-SNP), β = protein->outcome (IVW)
# Mediation = α × β, SE via Delta method
# ============================================================

library(data.table)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this file with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = TRUE)
SCRIPT_DIR <- dirname(script_path)
PROJECT_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), winslash = "/", mustWork = TRUE)
INPUT_DIR <- file.path(PROJECT_ROOT, "04_protein_mr", "results")
if (!file.exists(file.path(INPUT_DIR, "C6_rs429358_effects.csv"))) {
  INPUT_DIR <- file.path(PROJECT_ROOT, "04_protein_mr")
}

# ---- Load data ----
alpha <- fread(file.path(INPUT_DIR, "C6_rs429358_effects.csv"))  # rs429358 -> protein
beta  <- fread(file.path(INPUT_DIR, "C6_ukbppp_mr_all.csv"))     # protein -> outcome

# Only IVW for beta
beta <- beta[method == "Inverse variance weighted", ]

# Exclude APOE protein (cis self-mediation, tautological)
alpha <- alpha[Gene != "APOE", ]
beta  <- beta[exposure != "APOE", ]
cat(sprintf("Excluded APOE: %d proteins remain\n", nrow(alpha)))

# ---- Two-step product-of-coefficients analysis ----
cat("=== Biology-guided two-step mediation sensitivity ===\n\n")

# Author-defined biology categories
biology_categories <- list(
  Complement   = c("C1QA","C2","C3","C5","CFB","CFD","CFH","CFI","CFP","SERPING1"),
  Inflammatory = c("CSF1","IFNG","IL10","IL18","IL1B","IL6","TGFB1","TNF"),
  Chemokine    = c("CCL2","CCL5","CX3CL1","CXCL10","CXCL12"),
  Immune       = c("CD14","CD40","TLR4","TREM2"),
  Lipid        = c("APOE","LPA","PON1")
)

results <- list()
outcomes <- unique(beta$outcome)

for (out in outcomes) {
  cat(sprintf("\n--- %s ---\n", out))
  beta_o <- beta[outcome == out, ]

  for (i in 1:nrow(alpha)) {
    gene <- alpha$Gene[i]
    a  <- alpha$BETA[i];  se_a <- alpha$SE[i]
    b_row <- beta_o[exposure == gene, ]
    if (nrow(b_row) == 0) next

    b  <- b_row$b;        se_b <- b_row$se
    category <- alpha$Category[i]

    # Mediation effect = α × β
    ind_eff <- a * b
    # Delta method SE = sqrt(α²·SE_β² + β²·SE_α²)
    se_ind  <- sqrt(a^2 * se_b^2 + b^2 * se_a^2)
    z_ind   <- ind_eff / se_ind
    p_ind   <- 2 * pnorm(-abs(z_ind))
    ci_lo   <- ind_eff - 1.96 * se_ind
    ci_hi   <- ind_eff + 1.96 * se_ind

    results[[length(results) + 1]] <- data.frame(
      outcome   = out,
      gene      = gene,
      category  = category,
      alpha     = a,
      se_alpha  = se_a,
      beta      = b,
      se_beta   = se_b,
      mediation = ind_eff,
      se_med    = se_ind,
      z         = z_ind,
      pval      = p_ind,
      ci_lo     = ci_lo,
      ci_hi     = ci_hi,
      stringsAsFactors = FALSE
    )
  }
}

nv <- rbindlist(results)
fwrite(nv, file.path(SCRIPT_DIR, "D_two_step_mediation_all.csv"))
cat(sprintf("\nDone: %d mediation paths\n", nrow(nv)))

# ---- Significant mediators ----
cat("\n=== Significant mediators (nominal P<0.05) ===\n")
for (out in outcomes) {
  sub <- nv[outcome == out & pval < 0.05, ][order(pval), ]
  if (nrow(sub) == 0) next
  cat(sprintf("\n--- %s (%d significant) ---\n", out, nrow(sub)))
  for (i in 1:min(nrow(sub), 10)) {
    cat(sprintf("  %-12s α=%.4f β=%.4f mediation=%.4f [%.4f, %.4f] P=%.2e\n",
        sub$gene[i], sub$alpha[i], sub$beta[i], sub$mediation[i], sub$ci_lo[i], sub$ci_hi[i], sub$pval[i]))
  }
}

# ---- Biology-category descriptive summary ----
cat("\n\n=== Biology-category mediation summary ===\n")

# Total rs429358-C effects on the log-odds scale.
total_effect <- c(AD = 1.1275, Dry_AMD = -0.2109, Wet_AMD = -0.2216, Any_AMD = -0.2031)

for (out in outcomes) {
  cat(sprintf("\n--- %s (total effect ≈ %.4f) ---\n", out, total_effect[out]))

  for (pw_name in names(biology_categories)) {
    pw_genes <- biology_categories[[pw_name]]
    pw <- nv[outcome == out & gene %in% pw_genes, ]
    if (nrow(pw) == 0) next

    # Descriptive sum within an author-defined biology category.
    pw_sum <- sum(pw$mediation, na.rm = TRUE)
    # SE via pooled variance (sum of squared SEs)
    pw_se <- sqrt(sum(pw$se_med^2, na.rm = TRUE))
    pw_pct <- pw_sum / total_effect[out] * 100
    pw_ci_lo <- pw_sum - 1.96 * pw_se
    pw_ci_hi <- pw_sum + 1.96 * pw_se

    cat(sprintf("  %-14s sum=%.4f (%.1f%%) SE=%.4f CI=[%.4f, %.4f]\n",
        paste0(pw_name, ":"), pw_sum, pw_pct, pw_se, pw_ci_lo, pw_ci_hi))

    # List significant members
    pw_sig <- pw[pval < 0.05, ][order(pval), ]
    if (nrow(pw_sig) > 0) {
      cat(sprintf("    Significant: %s\n", paste(pw_sig$gene, collapse=", ")))
    }
  }

  # Residual (unmediated)
  all_med <- sum(nv[outcome == out, ]$mediation, na.rm = TRUE)
  residual <- total_effect[out] - all_med
  cat(sprintf("  %-14s sum=%.4f (%.1f%%)  RESIDUAL=%.4f (%.1f%%)\n",
      "ALL MEDIATED:", all_med, all_med/total_effect[out]*100,
      residual, residual/total_effect[out]*100))
}

cat("\nBiology-guided two-step mediation sensitivity complete.\n")
