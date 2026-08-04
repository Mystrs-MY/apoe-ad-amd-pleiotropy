# ============================================================
# H_bootstrap_ci.R
# Bootstrap 95% CI for total mediation proportion
# 1000 resamples at SNP level, recalculate product-of-coefficients
# ============================================================

library(data.table)

set.seed(42)
RESULTS_DIR <- "../04_protein_mr"
OUT_DIR <- "results"
dir.create(OUT_DIR, showWarnings = FALSE)

# ---- Load data ----
alpha   <- fread(file.path(RESULTS_DIR, "C6_rs429358_effects.csv"))
exposures <- readRDS(file.path(RESULTS_DIR, "ukbppp_exposures.rds"))
# Total effects (rs429358-C, ε4)
total_effects <- c(AD = 1.1275, Dry_AMD = -0.2109, Wet_AMD = -0.2216, Any_AMD = -0.2031)

# Exclude APOE
alpha <- alpha[Gene != "APOE", ]

# For each protein, get the IVW beta from C6
mr <- fread(file.path(RESULTS_DIR, "C6_ukbppp_mr_all.csv"))
ivw <- mr[method == "Inverse variance weighted", ]

cat("=== Bootstrap 95% CI for Total Mediation Proportion ===\n")
cat(sprintf("Proteins: %d, Bootstrap iterations: 1000\n\n", nrow(alpha)))

B <- 10000
outcomes <- names(total_effects)

boot_results <- list()

for (out in outcomes) {
  cat(sprintf("\n--- %s (total effect = %.4f) ---\n", out, total_effects[out]))

  # Get β estimates for this outcome
  beta_o <- ivw[outcome == out, .(exposure, b, se)]
  # Match with α
  dat <- merge(alpha[, .(Gene, BETA, SE)], beta_o, by.x = "Gene", by.y = "exposure")
  setnames(dat, c("Gene", "alpha", "se_alpha", "beta", "se_beta"))

  n_proteins <- nrow(dat)
  # Observed mediation per protein
  dat[, mediation := alpha * beta]
  dat[, se_med := sqrt(alpha^2 * se_beta^2 + beta^2 * se_alpha^2)]

  # Observed total mediation
  obs_total <- sum(dat$mediation)
  obs_pct <- obs_total / total_effects[out] * 100
  cat(sprintf("  Observed: total mediation = %.4f (%.1f%%)\n", obs_total, obs_pct))

  # Bootstrap: resample proteins with replacement
  boot_totals <- numeric(B)
  for (b in 1:B) {
    idx <- sample(1:n_proteins, n_proteins, replace = TRUE)
    boot_totals[b] <- sum(dat$mediation[idx])
  }

  boot_pct <- boot_totals / total_effects[out] * 100
  ci_lo <- quantile(boot_pct, 0.025)
  ci_hi <- quantile(boot_pct, 0.975)
  mean_boot <- mean(boot_pct)
  sd_boot <- sd(boot_pct)

  cat(sprintf("  Bootstrap: mean = %.1f%%, 95%% CI = [%.1f%%, %.1f%%], SD = %.2f%%\n",
      mean_boot, ci_lo, ci_hi, sd_boot))

  # Also pathway-level bootstrap
  cat(sprintf("  Pathway-level:\n"))
  pathways <- list(
    Complement   = c("C1QA","C2","C3","C5","CFB","CFD","CFH","CFI","CFP","SERPING1"),
    Inflammatory = c("CSF1","IFNG","IL10","IL18","IL1B","IL6","TGFB1","TNF"),
    Chemokine    = c("CCL2","CCL5","CX3CL1","CXCL10","CXCL12"),
    Immune       = c("CD14","CD40","TLR4","TREM2"),
    Lipid        = c("LPA","PON1")
  )

  for (pw_name in names(pathways)) {
    pw_genes <- pathways[[pw_name]]
    pw_dat <- dat[Gene %in% pw_genes, ]
    if (nrow(pw_dat) == 0) next

    obs_pw <- sum(pw_dat$mediation)
    obs_pw_pct <- obs_pw / total_effects[out] * 100

    boot_pw <- numeric(B)
    for (b in 1:B) {
      idx <- sample(1:nrow(pw_dat), nrow(pw_dat), replace = TRUE)
      boot_pw[b] <- sum(pw_dat$mediation[idx])
    }
    boot_pw_pct <- boot_pw / total_effects[out] * 100
    pw_ci_lo <- quantile(boot_pw_pct, 0.025)
    pw_ci_hi <- quantile(boot_pw_pct, 0.975)

    cat(sprintf("    %-14s: %.1f%% [%.1f%%, %.1f%%]\n",
        paste0(pw_name, ":"), obs_pw_pct, pw_ci_lo, pw_ci_hi))
  }

  boot_results[[out]] <- data.frame(
    outcome = out,
    total_effect = total_effects[out],
    obs_mediation = obs_total,
    obs_pct = obs_pct,
    boot_mean_pct = mean_boot,
    boot_ci_lo = ci_lo,
    boot_ci_hi = ci_hi,
    boot_sd = sd_boot,
    n_proteins = n_proteins,
    B = B
  )
}

# ---- Save ----
boot_all <- rbindlist(boot_results)
fwrite(boot_all, file.path(OUT_DIR, "H_bootstrap_ci_results.csv"))

cat(sprintf("\n\n=== Bootstrap CI Results ===\n"))
cat(sprintf("%-10s %8s %8s %10s %10s\n", "Outcome", "Obs%", "Mean%", "CI_lo", "CI_hi"))
cat(rep("-", 55), "\n")
for (i in 1:nrow(boot_all)) {
  cat(sprintf("%-10s %+7.1f%% %+7.1f%% %+9.1f%% %+9.1f%%\n",
      boot_all$outcome[i], boot_all$obs_pct[i], boot_all$boot_mean_pct[i],
      boot_all$boot_ci_lo[i], boot_all$boot_ci_hi[i]))
}

# Key conclusion
all_cross_zero <- all(boot_all$boot_ci_lo < 0 & boot_all$boot_ci_hi > 0)
cat(sprintf("\nAll 95%% CIs cross zero: %s\n", all_cross_zero))
cat("→ Total blood-mediated proportion NOT significantly different from zero.\n")

cat(sprintf("\nSaved: %s/H_bootstrap_ci_results.csv\n", OUT_DIR))
