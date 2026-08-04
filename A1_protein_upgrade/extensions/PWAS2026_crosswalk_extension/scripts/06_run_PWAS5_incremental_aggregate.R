#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

set.seed(20260714)
n_boot <- 10000L
rho_grid <- c(0, 0.25, 0.50, 0.75)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = TRUE)
extension_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
upgrade_root <- normalizePath(file.path(extension_root, "..", ".."), winslash = "/", mustWork = TRUE)

draw_correlated <- function(mu, se, rho, n) {
  k <- length(mu)
  if (!k) return(matrix(numeric(), nrow = n, ncol = 0))
  common <- rnorm(n)
  independent <- matrix(rnorm(n * k), nrow = n, ncol = k)
  z <- sqrt(rho) * common + sqrt(1 - rho) * independent
  sweep(sweep(z, 2, se, `*`), 2, mu, `+`)
}

source_specs <- list(
  main = list(
    primary = file.path(upgrade_root, "tables", "APOE_linkable_two_step_mediation.tsv"),
    extension = file.path(extension_root, "tables", "PWAS5_two_step_mediation_main.tsv")
  ),
  cis = list(
    primary = file.path(upgrade_root, "tables", "APOE_linkable_two_step_mediation_cis_sensitivity.tsv"),
    extension = file.path(extension_root, "tables", "PWAS5_two_step_mediation_cis.tsv")
  )
)

results <- list()
qa_rows <- list()
for (analysis_set in names(source_specs)) {
  primary_all <- fread(source_specs[[analysis_set]]$primary, na.strings = c("NA", ""), fill = TRUE)
  primary_total_reported <- primary_all[row_type == "total"]
  primary <- primary_all[row_type == "protein"]
  extension <- fread(source_specs[[analysis_set]]$extension, na.strings = c("NA", ""), fill = TRUE)
  extension <- extension[mediation_status == "estimable"]
  primary[, source_panel := "frozen_primary_25"]
  extension[, source_panel := "PWAS5_extension"]

  for (mapping_scope in c("expanded", "high_confidence")) {
    primary_scope <- if (mapping_scope == "high_confidence") primary[mapping_confidence == "high"] else primary
    extension_scope <- if (mapping_scope == "high_confidence") extension[mapping_confidence == "high"] else extension
    combined <- rbindlist(list(primary_scope, extension_scope), fill = TRUE)
    groups <- unique(rbindlist(list(
      primary_scope[, .(variant, outcome)], extension_scope[, .(variant, outcome)]
    ), fill = TRUE))
    for (g in seq_len(nrow(groups))) {
      variant_i <- groups$variant[g]
      outcome_i <- groups$outcome[g]
      dat <- combined[variant == variant_i & outcome == outcome_i]
      if (!nrow(dat)) next
      primary_index <- which(dat$source_panel == "frozen_primary_25")
      extension_index <- which(dat$source_panel == "PWAS5_extension")
      total_effect <- dat$total_effect[1]
      se_total <- dat$SE_total[1]
      for (rho in rho_grid) {
        alpha_draw <- draw_correlated(dat$alpha, dat$SE_alpha, rho, n_boot)
        beta_draw <- draw_correlated(dat$beta, dat$SE_beta, rho, n_boot)
        indirect_matrix <- alpha_draw * beta_draw
        primary_draw <- if (length(primary_index)) rowSums(indirect_matrix[, primary_index, drop = FALSE]) else rep(0, n_boot)
        extension_draw <- if (length(extension_index)) rowSums(indirect_matrix[, extension_index, drop = FALSE]) else rep(0, n_boot)
        combined_draw <- primary_draw + extension_draw
        total_draw <- rnorm(n_boot, total_effect, se_total)

        panels <- list(
          frozen_primary = list(draw = primary_draw, idx = primary_index),
          PWAS5_extension = list(draw = extension_draw, idx = extension_index),
          combined_sensitivity = list(draw = combined_draw, idx = seq_len(nrow(dat))),
          combined_minus_primary_increment = list(draw = extension_draw, idx = extension_index)
        )
        for (panel_name in names(panels)) {
          item <- panels[[panel_name]]
          point <- if (length(item$idx)) sum(dat$indirect_effect[item$idx]) else 0
          proportion_draw <- item$draw / total_draw
          results[[length(results) + 1L]] <- data.table(
            analysis_set = analysis_set,
            mapping_scope = mapping_scope,
            assumed_common_error_correlation = rho,
            variant = variant_i,
            outcome = outcome_i,
            panel = panel_name,
            n_primary_proteins = length(primary_index),
            n_extension_proteins = length(extension_index),
            n_total_proteins = length(item$idx),
            summed_indirect_effect = point,
            indirect_CI_lower = quantile(item$draw, 0.025, na.rm = TRUE),
            indirect_CI_upper = quantile(item$draw, 0.975, na.rm = TRUE),
            total_effect = total_effect,
            SE_total = se_total,
            mediated_proportion = point / total_effect,
            mediated_proportion_CI_lower = quantile(proportion_draw, 0.025, na.rm = TRUE),
            mediated_proportion_CI_upper = quantile(proportion_draw, 0.975, na.rm = TRUE),
            bootstrap_n = n_boot,
            paired_draw_design = TRUE,
            interpretation_boundary = paste(
              "Sensitivity analysis for the pre-specified PWAS5 extension; primary 25-protein results remain frozen.",
              "Equicorrelated alpha and beta error draws are assumptions, not empirically estimated covariance."
            )
          )
        }
      }
      if (mapping_scope == "expanded") {
        reported <- primary_total_reported[variant == variant_i & outcome == outcome_i]
        if (nrow(reported)) {
          recomputed <- sum(primary_scope[variant == variant_i & outcome == outcome_i]$indirect_effect)
          qa_rows[[length(qa_rows) + 1L]] <- data.table(
            analysis_set = analysis_set, variant = variant_i, outcome = outcome_i,
            frozen_reported_total = reported$indirect_effect[1],
            recomputed_from_frozen_protein_rows = recomputed,
            absolute_difference = abs(reported$indirect_effect[1] - recomputed),
            exact_within_1e_12 = abs(reported$indirect_effect[1] - recomputed) < 1e-12
          )
        }
      }
    }
  }
}

result <- rbindlist(results, fill = TRUE)
qa <- rbindlist(qa_rows, fill = TRUE)
fwrite(result, file.path(extension_root, "tables", "PWAS5_incremental_aggregate_mediation.tsv"), sep = "\t", na = "NA")
fwrite(result[assumed_common_error_correlation %in% rho_grid],
       file.path(extension_root, "tables", "PWAS5_covariance_mapping_sensitivity.tsv"), sep = "\t", na = "NA")
fwrite(qa, file.path(extension_root, "tables", "PWAS5_primary_reproduction_QA.tsv"), sep = "\t", na = "NA")
if (nrow(qa) != 16L || !all(qa$exact_within_1e_12)) {
  stop("Frozen primary aggregate reproduction failed; inspect PWAS5_primary_reproduction_QA.tsv")
}
writeLines(capture.output(sessionInfo()), file.path(extension_root, "logs", "incremental_aggregate_sessionInfo.txt"))
message("Aggregate sensitivity rows: ", nrow(result))
