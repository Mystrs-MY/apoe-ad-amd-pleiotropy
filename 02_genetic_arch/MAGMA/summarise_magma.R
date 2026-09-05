#!/usr/bin/env Rscript

rm(list = ls())

required <- c("data.table")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R package(s): ", paste(missing, collapse = ", "))

suppressPackageStartupMessages({
  library(data.table)
})

script_arg <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)][1])
MAGMA_ROOT <- dirname(normalizePath(script_file, winslash = "/", mustWork = TRUE))
ROOT <- normalizePath(file.path(MAGMA_ROOT, "..", ".."), winslash = "/", mustWork = TRUE)
RUN_TAG <- Sys.getenv("A1_MAGMA_RUN_TAG", unset = "correctedN_20260903")
RESULT_DIR <- file.path(MAGMA_ROOT, paste0("results_", RUN_TAG))
LAVA_RESULT_DIR <- file.path(ROOT, "02_genetic_arch", "LAVA", paste0("results_", RUN_TAG))
TABLE_DIR <- file.path(ROOT, "tables_submission", "supplementary_tables")
SOURCE_DIR <- file.path(ROOT, "figures_submission", "source_data")
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SOURCE_DIR, recursive = TRUE, showWarnings = FALSE)

traits <- data.table(
  trait_id = c("AD", "Dry_AMD", "Wet_AMD", "Any_AMD"),
  trait = c("AD", "Dry AMD", "Wet AMD", "Any AMD")
)

lava_files <- file.path(
  LAVA_RESULT_DIR,
  paste0(
    "LAVA_FullScan_AD_vs_", c("DryAMD", "WetAMD", "AnyAMD"),
    "_", RUN_TAG, ".csv"
  )
)
if (!all(file.exists(lava_files))) {
  stop("Missing corrected LAVA full-scan output(s): ",
       paste(lava_files[!file.exists(lava_files)], collapse = ", "))
}
lava <- rbindlist(lapply(lava_files, fread), fill = TRUE)
lava[, p := suppressWarnings(as.numeric(p))]
lava_significant <- unique(
  lava[is.finite(p) & p < 0.05 / 2495,
       .(locus_id = as.integer(locus), chromosome = as.integer(CHR),
         locus_start_bp = as.integer(START), locus_end_bp = as.integer(STOP))]
)
if (!nrow(lava_significant)) stop("No corrected LAVA locus passed 0.05/2495.")
loci <- copy(lava_significant)
loci[, locus := fifelse(
  chromosome == 19L & locus_start_bp <= 45411941L & locus_end_bp >= 45411941L,
  "APOE", paste0("LAVA_", locus_id)
)]
loci[, locus_source := paste0("Corrected LAVA locus ", locus_id)]
setorder(loci, chromosome, locus_start_bp, locus_id)
if (anyDuplicated(loci$locus)) stop("Non-unique corrected LAVA locus labels.")

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

tagged_log_dir <- file.path(MAGMA_ROOT, paste0("logs_", RUN_TAG))
dir.create(tagged_log_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(capture.output(sessionInfo()), file.path(tagged_log_dir, "summary_sessionInfo.txt"))
cat("MAGMA full table, locus summary, and source data exported.\n")
print(summary)
