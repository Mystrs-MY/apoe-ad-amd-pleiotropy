#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = TRUE)
project_root <- if (length(args) >= 1L) {
  normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
} else {
  Sys.getenv("A1_PROJECT_ROOT", unset = normalizePath(file.path(dirname(script_path), "..", "..", "..", ".."), winslash = "/", mustWork = TRUE))
}
ext_root <- file.path(project_root, "A1_protein_upgrade/extensions/PWAS2026_crosswalk_extension")
table_dir <- file.path(ext_root, "tables")
figure_dir <- file.path(project_root, "figures_submission/code/grouped_supplementary_panel_assets")
preview_dir <- file.path(project_root, "figures_submission/code/grouped_supplementary_panel_assets/previews")
source_dir <- file.path(project_root, "figures_submission/source_data")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(preview_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

palette <- c(BCAM = "#3C5488", CD55 = "#E64B35", LILRB1 = "#00A087", LILRB5 = "#4DBBD5", SCARA5 = "#F39B7F")
theme_a1 <- function() {
  theme_classic(base_size = 8.5, base_family = "Arial") +
    theme(
      plot.title = element_blank(), axis.text = element_text(color = "#333333"),
      axis.title = element_text(color = "#222222"),
      axis.line = element_line(linewidth = 0.45, color = "#333333"),
      axis.ticks = element_line(linewidth = 0.4, color = "#333333"),
      legend.title = element_blank(), strip.background = element_blank(),
      strip.text = element_text(face = "bold", color = "#222222"),
      panel.spacing = grid::unit(3, "mm"), plot.margin = margin(5, 7, 5, 5)
    )
}
export_plot <- function(plot, stem, width, height) {
  ggsave(file.path(figure_dir, paste0(stem, ".pdf")), plot, width = width, height = height,
         units = "in", device = cairo_pdf)
  ggsave(file.path(preview_dir, paste0(stem, ".png")), plot, width = width, height = height,
         units = "in", dpi = 300, bg = "white")
}

crosswalk <- fread(file.path(table_dir, "PWAS5_crosswalk_mapping.tsv"))
flow <- crosswalk[, .(
  gene_symbol,
  stage1 = "2026 SomaScan PWAS",
  stage2 = paste0(assay_target_ID, " Olink"),
  stage3 = ifelse(strict_mapping_eligible, "Exact assay", "Gene-level crosswalk"),
  mapping_confidence,
  frozen_before_results
)]
flow[, gene_symbol := factor(gene_symbol, levels = rev(c("BCAM", "CD55", "LILRB1", "LILRB5", "SCARA5")))]
flow_long <- melt(flow, id.vars = c("gene_symbol", "mapping_confidence", "frozen_before_results"),
                  measure.vars = c("stage1", "stage2", "stage3"), variable.name = "stage", value.name = "label")
flow_long[, stage := factor(stage, levels = c("stage1", "stage2", "stage3"),
                            labels = c("Literature candidate", "UKB-PPP assay", "Mapping boundary"))]
fwrite(flow_long, file.path(source_dir, "FigS7a_PWAS5_crosswalk_source_data.tsv"), sep = "\t")

p12a <- ggplot(flow_long, aes(stage, gene_symbol, group = gene_symbol, color = gene_symbol)) +
  geom_line(linewidth = 0.75, alpha = 0.75) +
  geom_point(size = 2.7) +
  geom_text(aes(label = label), color = "#282828", size = 2.35, nudge_y = 0.18, check_overlap = TRUE) +
  scale_color_manual(values = palette) +
  labs(x = NULL, y = NULL) +
  theme_a1() + theme(legend.position = "none", axis.ticks.x = element_blank())
export_plot(p12a, "FigS7a_PWAS5_frozen_crosswalk", 7.0, 3.4)

alpha <- fread(file.path(table_dir, "PWAS5_APOE_alpha.tsv"))
alpha <- alpha[, .(gene_symbol, variant, beta, SE, P_value)]
alpha[, `:=`(lower = beta - 1.96 * SE, upper = beta + 1.96 * SE)]
alpha[, gene_symbol := factor(gene_symbol, levels = rev(c("BCAM", "CD55", "LILRB1", "LILRB5", "SCARA5")))]
alpha[, variant := factor(variant, levels = c("rs429358", "rs7412"), labels = c("APOE rs429358", "APOE rs7412"))]
fwrite(alpha, file.path(source_dir, "FigS7b_PWAS5_APOE_alpha_source_data.tsv"), sep = "\t")

p12b <- ggplot(alpha, aes(beta, gene_symbol, color = gene_symbol)) +
  geom_vline(xintercept = 0, linewidth = 0.45, color = "#777777") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0, linewidth = 0.62) +
  geom_point(size = 2.05) +
  facet_wrap(~variant, ncol = 2, scales = "free_x") +
  scale_color_manual(values = palette) +
  labs(x = "APOE allele effect on Olink protein (SD units)", y = NULL) +
  theme_a1() + theme(legend.position = "none")
export_plot(p12b, "FigS7b_PWAS5_APOE_alpha_forest", 6.8, 3.15)

main <- fread(file.path(table_dir, "PWAS5_beta_main.tsv"))
cis <- fread(file.path(table_dir, "PWAS5_beta_cis.tsv"))
main_primary <- main[method_role == "primary", .(gene_symbol, outcome, beta, SE, P_value,
                                                  analysis_set = "Genome-wide instruments", beta_status)]
cis_primary <- cis[method_role == "primary", .(gene_symbol, outcome, beta, SE, P_value,
                                                analysis_set = "Cis instruments", beta_status)]
beta <- rbindlist(list(main_primary, cis_primary), fill = TRUE)
beta <- beta[grepl("^reestimated", beta_status) & is.finite(beta) & is.finite(SE)]
beta[, `:=`(lower = beta - 1.96 * SE, upper = beta + 1.96 * SE)]
beta[, gene_symbol := factor(gene_symbol, levels = rev(c("BCAM", "CD55", "LILRB1", "LILRB5", "SCARA5")))]
beta[, outcome := factor(outcome, levels = c("AD", "dry_AMD", "wet_AMD", "any_AMD"),
                         labels = c("AD", "Dry AMD", "Wet AMD", "Any AMD"))]
fwrite(beta, file.path(source_dir, "FigS7c_PWAS5_beta_source_data.tsv"), sep = "\t")

p12c <- ggplot(beta, aes(beta, gene_symbol, color = gene_symbol, shape = analysis_set)) +
  geom_vline(xintercept = 0, linewidth = 0.45, color = "#777777") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0, linewidth = 0.55,
                 position = position_dodge(width = 0.34)) +
  geom_point(size = 1.75, position = position_dodge(width = 0.34)) +
  facet_wrap(~outcome, ncol = 2, scales = "free_x") +
  scale_color_manual(values = palette) +
  labs(x = "Protein effect on outcome (log-odds scale)", y = NULL) +
  theme_a1() + theme(legend.position = "top")
export_plot(p12c, "FigS7c_PWAS5_main_and_cis_beta_forest", 7.0, 5.6)

agg <- fread(file.path(table_dir, "PWAS5_incremental_aggregate_mediation.tsv"))
agg <- agg[analysis_set == "main" & mapping_scope == "expanded" &
             assumed_common_error_correlation %chin% c(0, 0.25, 0.5, 0.75) &
             panel %chin% c("frozen_primary", "combined_sensitivity")]
agg[, panel_label := factor(panel, levels = c("frozen_primary", "combined_sensitivity"),
                            labels = c("Frozen 25-protein panel", "Panel plus PWAS5 crosswalk"))]
agg[, outcome := factor(outcome, levels = c("AD", "dry_AMD", "wet_AMD", "any_AMD"),
                        labels = c("AD", "Dry AMD", "Wet AMD", "Any AMD"))]
agg[, rho_label := factor(sprintf("rho = %.2f", assumed_common_error_correlation),
                          levels = sprintf("rho = %.2f", c(0, 0.25, 0.5, 0.75)))]
fwrite(agg, file.path(source_dir, "FigS7d_PWAS5_aggregate_mediation_source_data.tsv"), sep = "\t")

p12d <- ggplot(agg, aes(mediated_proportion, rho_label, color = panel_label, shape = panel_label)) +
  geom_vline(xintercept = 0, linewidth = 0.45, color = "#777777") +
  geom_errorbarh(aes(xmin = mediated_proportion_CI_lower, xmax = mediated_proportion_CI_upper),
                 height = 0, linewidth = 0.55, position = position_dodge(width = 0.36)) +
  geom_point(size = 2, position = position_dodge(width = 0.36)) +
  facet_grid(variant ~ outcome, scales = "free_x") +
  scale_color_manual(values = c("Frozen 25-protein panel" = "#3C5488", "Panel plus PWAS5 crosswalk" = "#E64B35")) +
  labs(x = "Aggregate mediated proportion (sensitivity estimate)", y = "Assumed common error correlation") +
  theme_a1() + theme(legend.position = "top", axis.text.y = element_text(size = 7.2))
export_plot(p12d, "FigS7d_PWAS5_incremental_aggregate_sensitivity", 7.4, 6.2)

message("PWAS5 figure assets exported: FigS7a-FigS7d")
