#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(MASS)
})

set.seed(20260712)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
input <- file.path(root, "tables", "supplementary_name_match_revision",
                   "TableS25_Expanded_Primary_Two_Step_Mediation.tsv")
output <- file.path(root, "tables", "covariance_mapping_bootstrap_sensitivity.tsv")
session_path <- file.path(root, "logs", "covariance_mapping_bootstrap_sessionInfo.txt")

n_boot <- 10000L
rho_grid <- c(0, 0.25, 0.50, 0.75)
x <- fread(input)
x <- x[row_type == "protein" & is.finite(alpha) & is.finite(SE_alpha) &
         is.finite(beta) & is.finite(SE_beta) & is.finite(total_effect) & is.finite(SE_total)]

draw_correlated <- function(mu, se, rho, n) {
  k <- length(mu)
  correlation <- matrix(rho, k, k)
  diag(correlation) <- 1
  covariance <- correlation * tcrossprod(se)
  MASS::mvrnorm(n = n, mu = mu, Sigma = covariance, tol = 1e-8)
}

rows <- list()
for (variant_name in unique(x$variant)) {
  for (outcome_name in unique(x$outcome)) {
    base <- x[variant == variant_name & outcome == outcome_name]
    if (!nrow(base)) next
    for (mapping_scenario in c("expanded_name_matched", "high_confidence_only")) {
      dat <- if (mapping_scenario == "high_confidence_only") base[mapping_confidence == "high"] else base
      if (!nrow(dat)) next
      for (rho in rho_grid) {
        alpha_draw <- draw_correlated(dat$alpha, dat$SE_alpha, rho, n_boot)
        beta_draw <- draw_correlated(dat$beta, dat$SE_beta, rho, n_boot)
        indirect_draw <- rowSums(alpha_draw * beta_draw)
        total_draw <- rnorm(n_boot, mean = dat$total_effect[1], sd = dat$SE_total[1])
        stable <- abs(total_draw) > 3 * dat$SE_total[1]
        proportion_draw <- rep(NA_real_, n_boot)
        proportion_draw[stable] <- indirect_draw[stable] / total_draw[stable]
        rows[[length(rows) + 1L]] <- data.table(
          variant = variant_name,
          outcome = outcome_name,
          mapping_scenario = mapping_scenario,
          n_proteins = nrow(dat),
          assumed_common_error_correlation = rho,
          bootstrap_n = n_boot,
          stable_total_effect_draw_fraction = mean(stable),
          indirect_effect_median = median(indirect_draw),
          indirect_CI_lower = quantile(indirect_draw, 0.025),
          indirect_CI_upper = quantile(indirect_draw, 0.975),
          mediated_proportion_median = median(proportion_draw, na.rm = TRUE),
          mediated_proportion_CI_lower = quantile(proportion_draw, 0.025, na.rm = TRUE),
          mediated_proportion_CI_upper = quantile(proportion_draw, 0.975, na.rm = TRUE),
          covariance_model = "Equicorrelation sensitivity applied separately to alpha and beta errors; empirical cross-protein covariance unavailable.",
          mapping_uncertainty_model = "Scenario envelope: expanded name-matched set versus high-confidence-only set; no arbitrary assay-validity probabilities assigned."
        )
      }
    }
  }
}

result <- rbindlist(rows)
fwrite(result, output, sep = "\t")
writeLines(capture.output(sessionInfo()), session_path)
stopifnot(nrow(result) == 64L, all(is.finite(result$mediated_proportion_CI_lower)))
message("Wrote covariance/mapping bootstrap sensitivity rows: ", nrow(result))
