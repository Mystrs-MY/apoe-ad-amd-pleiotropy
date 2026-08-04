#!/usr/bin/env Rscript

rm(list = ls())

required <- c("data.table", "ggplot2", "svglite", "ragg")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R package(s): ", paste(missing, collapse = ", "))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

script_arg <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)][1])
MAGMA_ROOT <- dirname(normalizePath(script_file, winslash = "/", mustWork = TRUE))
ROOT <- normalizePath(file.path(MAGMA_ROOT, "..", ".."), winslash = "/", mustWork = TRUE)
RESULT_DIR <- file.path(MAGMA_ROOT, "results")
TABLE_DIR <- file.path(ROOT, "tables_submission", "supplementary_tables")
FIG_DIR <- file.path(ROOT, "figures_submission", "supplementary_figures")
PREVIEW_DIR <- file.path(ROOT, "figures_submission", "previews")
SOURCE_DIR <- file.path(ROOT, "figures_submission", "source_data")
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PREVIEW_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SOURCE_DIR, recursive = TRUE, showWarnings = FALSE)

traits <- data.table(
  trait_id = c("AD", "Dry_AMD", "Wet_AMD", "Any_AMD"),
  trait = c("AD", "Dry AMD", "Wet AMD", "Any AMD")
)

loci <- data.table(
  locus = c("APOE", "IGH", "chr2", "chr20"),
  chromosome = c(19L, 14L, 2L, 20L),
  locus_start_bp = c(45040933L, 105091341L, 44104111L, 51533751L),
  locus_end_bp = c(45893307L, 106347143L, 45189468L, 52411532L),
  locus_source = c(
    "LAVA locus 2351", "LAVA locus 2031",
    "LAVA locus 249", "LAVA locus 2415"
  )
)

gene_results <- list()
for (i in seq_len(nrow(traits))) {
  file <- file.path(RESULT_DIR, paste0(traits$trait_id[i], ".genes.out"))
  if (!file.exists(file) || file.info(file)$size == 0) {
    stop("Missing complete MAGMA gene output: ", file)
  }
  dat <- fread(file)
  expected <- c("GENE", "CHR", "START", "STOP", "NSNPS", "NPARAM", "N", "ZSTAT", "P")
  if (!identical(names(dat), expected)) {
    stop("Unexpected MAGMA schema in ", file, ": ", paste(names(dat), collapse = ", "))
  }
  dat[, `:=`(
    trait_id = traits$trait_id[i],
    trait = traits$trait[i],
    genes_tested = nrow(dat),
    bonferroni_threshold = 0.05 / nrow(dat),
    bonferroni_significant = P < 0.05 / nrow(dat),
    magma_version = "1.10",
    gene_location = "Rev.NCBI37.3.gene.loc",
    gene_window = "upstream_35kb_downstream_10kb",
    LD_reference = "1000_Genomes_EUR_n503"
  )]
  gene_results[[i]] <- dat
}
all_genes <- rbindlist(gene_results, use.names = TRUE)
if (all_genes[, uniqueN(trait)] != 4L) stop("Not all four traits were loaded.")

full_file <- file.path(TABLE_DIR, "TableS5g_MAGMA_Full_Gene_Results.tsv.gz")
fwrite(all_genes, full_file, sep = "\t", na = "NA", compress = "gzip")

summary_rows <- list()
for (trait_name in traits$trait) {
  trait_data <- all_genes[trait == trait_name]
  threshold <- unique(trait_data$bonferroni_threshold)
  if (length(threshold) != 1L) stop("Non-unique MAGMA threshold for ", trait_name)
  for (j in seq_len(nrow(loci))) {
    loc <- loci[j]
    genes <- trait_data[
      CHR == loc$chromosome & STOP >= loc$locus_start_bp & START <= loc$locus_end_bp
    ]
    if (!nrow(genes)) stop("No genes overlap ", loc$locus, " in ", trait_name)
    best <- genes[order(P)][1]
    summary_rows[[length(summary_rows) + 1L]] <- data.table(
      trait = trait_name,
      locus = loc$locus,
      chromosome = loc$chromosome,
      locus_start_bp = loc$locus_start_bp,
      locus_end_bp = loc$locus_end_bp,
      locus_source = loc$locus_source,
      genes_tested_genome_wide = nrow(trait_data),
      bonferroni_threshold = threshold,
      genes_overlapping_locus = nrow(genes),
      bonferroni_significant_genes_in_locus = sum(genes$P < threshold),
      best_gene = best$GENE,
      best_Z = best$ZSTAT,
      best_P = best$P,
      max_log10P = -log10(max(best$P, 1e-300)),
      interpretation = ifelse(
        any(genes$P < threshold),
        "At least one gene in this prespecified LAVA locus passed the phenotype-specific genome-wide gene Bonferroni threshold.",
        "No gene in this prespecified LAVA locus passed the phenotype-specific genome-wide gene Bonferroni threshold."
      )
    )
  }
}
summary <- rbindlist(summary_rows)
summary[, trait := factor(trait, levels = traits$trait)]
summary[, locus := factor(locus, levels = loci$locus)]
setorder(summary, trait, locus)

summary_file <- file.path(TABLE_DIR, "TableS5e_MAGMA_Locus_Summary.csv")
fwrite(summary, summary_file, bom = TRUE, na = "NA")

source_file <- file.path(SOURCE_DIR, "FigS2_MAGMA_source_data.tsv")
fwrite(summary, source_file, sep = "\t", na = "NA")

palette <- c(
  "APOE" = "#C82433",
  "IGH" = "#4C78A8",
  "chr2" = "#59A14F",
  "chr20" = "#F28E2B"
)
strict_threshold <- min(summary$bonferroni_threshold)
plot_data <- copy(summary)
plot_data[, sig_label := fifelse(
  bonferroni_significant_genes_in_locus > 0,
  paste0("n=", bonferroni_significant_genes_in_locus),
  ""
)]

p <- ggplot(plot_data, aes(x = trait, y = max_log10P, fill = locus)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.72) +
  geom_hline(
    yintercept = -log10(strict_threshold),
    color = "#A61B29", linetype = "dashed", linewidth = 0.55
  ) +
  geom_text(
    aes(label = sig_label),
    position = position_dodge(width = 0.78),
    vjust = -0.35, size = 3.1, fontface = "bold", color = "#333333"
  ) +
  annotate(
    "text", x = 4.42, y = -log10(strict_threshold) + 0.25,
    label = "Bonferroni", color = "#A61B29", hjust = 1,
    size = 3.2
  ) +
  scale_fill_manual(values = palette, drop = FALSE) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.08)),
    limits = c(0, max(plot_data$max_log10P) * 1.10)
  ) +
  labs(
    x = NULL,
    y = expression(-log[10](P[max])),
    fill = "Locus"
  ) +
  theme_bw(base_family = "Arial") +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "#E4E4E4", linewidth = 0.28),
    panel.border = element_rect(color = "#333333", fill = NA, linewidth = 0.75),
    axis.text = element_text(color = "#333333", size = 10),
    axis.title.y = element_text(size = 11),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 10),
    legend.key.height = grid::unit(0.34, "cm"),
    legend.key.width = grid::unit(0.52, "cm"),
    plot.margin = margin(8, 12, 4, 10)
  )

pdf_file <- file.path(FIG_DIR, "FigS2_MAGMA_Bar.pdf")
svg_file <- file.path(FIG_DIR, "FigS2_MAGMA_Bar.svg")
tiff_file <- file.path(FIG_DIR, "FigS2_MAGMA_Bar.tiff")
preview_file <- file.path(PREVIEW_DIR, "FigS2_MAGMA_Bar.png")

ggsave(pdf_file, p, width = 9.2, height = 5.6, device = cairo_pdf, family = "Arial")
ggsave(svg_file, p, width = 9.2, height = 5.6, device = svglite::svglite)
ggsave(
  tiff_file, p, width = 9.2, height = 5.6, units = "in", dpi = 600,
  device = ragg::agg_tiff, compression = "lzw"
)
ggsave(
  preview_file, p, width = 9.2, height = 5.6, units = "in", dpi = 180,
  device = ragg::agg_png
)

writeLines(capture.output(sessionInfo()), file.path(MAGMA_ROOT, "logs", "summary_plot_sessionInfo.txt"))
cat("MAGMA summary, full table, source data and Fig. S2 exported.\n")
print(summary)
