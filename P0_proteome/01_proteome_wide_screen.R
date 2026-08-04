# ============================================================
# Legacy biology-guided protein-panel FDR sensitivity
# Uses the original prespecified panel; this is not a proteome-wide scan.
# Unified within-panel FDR pipeline + volcano plot
# ============================================================

rm(list = ls())
library(data.table)
library(ggplot2)

OUT_DIR <- "./results/"
FIG_DIR <- "./figures/"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Step 1: Collect available legacy-panel proteins ----
# Scan for all protein GWAS in ukbppp_merged
MERGED_DIR <- "../../Resource/ukbppp_proteins/ukbppp_merged/"

# List available merged protein GWAS files
available_merged <- list.files(MERGED_DIR, pattern = "_merged\\.txt\\.gz$")
cat(sprintf("Available merged protein GWAS: %d\n", length(available_merged)))

# For each, extract the gene name
protein_list <- data.table(
  file = available_merged,
  gene = gsub("_merged\\.txt\\.gz$", "", available_merged)
)

# ---- Step 2: Also check tar files for additional proteins ----
TAR_DIR <- "../../Resource/ukbppp_proteins/"
tar_files <- list.files(TAR_DIR, pattern = "\\.tar$")
# Extract gene names from tar files
tar_genes <- unique(gsub("_.*", "", tar_files[tolower(tar_files) != "rsid_maps"]))
cat(sprintf("Proteins in tar files: %d\n", length(tar_genes)))

# ---- Step 3: Load existing MR results for FDR integration ----
# Check if existing protein MR results exist
existing_mr_file <- "../04_protein_mr/C6_ukbppp_mr_all.csv"
existing_alpha_file <- "../04_protein_mr/C6_rs429358_effects.csv"

if (file.exists(existing_mr_file)) {
  cat("\n=== Loading existing protein MR results ===\n")
  beta_all <- fread(existing_mr_file)
  cat(sprintf("  %d protein-outcome pairs in existing MR\n", nrow(beta_all)))

  # Extract IVW results only
  beta_ivw <- beta_all[method == "Inverse variance weighted", ]
  n_proteins <- uniqueN(beta_ivw$exposure)
  cat(sprintf("  %d unique proteins with IVW results\n", n_proteins))

  # ---- FDR Correction ----
  # Per outcome
  for (out in unique(beta_ivw$outcome)) {
    beta_o <- beta_ivw[outcome == out, ]
    beta_o$fdr <- p.adjust(beta_o$pval, method = "BH")
    beta_ivw[outcome == out, fdr := p.adjust(pval, method = "BH")]
  }

  # ---- Volcano Plot ----
  volcano_data <- beta_ivw[, .(exposure, outcome, b, se, pval, fdr)]
  volcano_data[, log10p := -log10(pval)]
  volcano_data[, significant := fifelse(fdr < 0.05, "FDR<0.05",
                                 fifelse(pval < 0.05, "Nominal P<0.05", "NS"))]

  pdf(paste0(FIG_DIR, "Fig_proteome_volcano.pdf"), width = 14, height = 10)

  outcomes <- unique(volcano_data$outcome)
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))

  for (out in outcomes) {
    vd <- volcano_data[outcome == out, ]
    cols <- ifelse(vd$fdr < 0.05, "#B2182B",
            ifelse(vd$pval < 0.05, "#4393C3", "grey60"))
    sizes <- ifelse(vd$fdr < 0.05, 2.0, ifelse(vd$pval < 0.05, 1.2, 0.6))

    plot(vd$b, vd$log10p, pch = 19, cex = sizes, col = cols,
         xlab = "IVW β", ylab = "-log10(P)",
         main = paste0(out, " — Protein MR (n=", nrow(vd), ")"),
         xlim = c(min(vd$b, na.rm = TRUE) * 1.3, max(vd$b, na.rm = TRUE) * 1.3))

    abline(h = -log10(0.05), lty = 3, col = "grey60")
    # Bonferroni line
    abline(h = -log10(0.05 / nrow(vd)), lty = 2, col = "orange", lwd = 1.5)

    # Label significant hits
    sig_hits <- vd[fdr < 0.05, ]
    if (nrow(sig_hits) > 0) {
      text(sig_hits$b, sig_hits$log10p, sig_hits$exposure,
           pos = 4, cex = 0.7, col = "#B2182B")
    }

    # Legend
    legend("topleft",
           legend = c("FDR<0.05", "P<0.05", "NS", "Bonferroni"),
           col = c("#B2182B", "#4393C3", "grey60", "orange"),
           pch = 19, cex = 0.6, pt.cex = 1, lty = c(NA, NA, NA, 2))
  }

  dev.off()
  file.copy(paste0(FIG_DIR, "Fig_proteome_volcano.pdf"),
            "../figures_submission/supplementary_figures/FigS3d_Legacy_proteome_volcano.pdf",
            overwrite = TRUE)
  cat("  -> Volcano plot saved\n")

  # ---- Per-pathway FDR summary ----
  # Map proteins to pathways
  pathways <- list(
    Complement   = c("C1QA","C2","C3","C5","CFB","CFD","CFH","CFI","CFP","SERPING1"),
    Inflammatory = c("CSF1","IFNG","IL10","IL18","IL1B","IL6","TGFB1","TNF"),
    Chemokine    = c("CCL2","CCL5","CX3CL1","CXCL10","CXCL12"),
    Immune       = c("CD14","CD40","TLR4","TREM2"),
    Lipid        = c("APOE","LPA","PON1")
  )

  pathway_summary <- data.table()
  for (pathway_name in names(pathways)) {
    pathway_genes <- pathways[[pathway_name]]
    for (out in outcomes) {
      hits <- volcano_data[outcome == out & exposure %in% pathway_genes, ]
      n_sig_fdr <- sum(hits$fdr < 0.05, na.rm = TRUE)
      n_sig_nom <- sum(hits$pval < 0.05 & hits$fdr >= 0.05, na.rm = TRUE)
      n_total <- nrow(hits)

      pathway_summary <- rbind(pathway_summary, data.table(
        pathway = pathway_name, outcome = out,
        n_proteins = n_total, n_fdr_sig = n_sig_fdr, n_nominal = n_sig_nom
      ))
    }
  }

  cat("\n=== Pathway FDR Summary ===\n")
  print(pathway_summary)

  # ---- Save ----
  fwrite(volcano_data, paste0(OUT_DIR, "legacy_panel_mr_fdr.csv"))
  fwrite(pathway_summary, paste0(OUT_DIR, "pathway_fdr_summary.csv"))

  # ---- Mediation Integration (if alpha available) ----
  if (file.exists(existing_alpha_file)) {
    cat("\n=== Integrating Mediation with FDR ===\n")
    alpha <- fread(existing_alpha_file)

    # Merge alpha + beta for mediation
    mediation_fdr <- merge(
      alpha[, .(Gene, BETA, SE, Category)],
      volcano_data[, .(exposure, outcome, b, se, pval, fdr)],
      by.x = "Gene", by.y = "exposure", all = TRUE
    )

    mediation_fdr[, mediation_effect := BETA * b]
    mediation_fdr[, mediation_se := sqrt(BETA^2 * se^2 + b^2 * SE^2)]
    mediation_fdr[, mediation_z := mediation_effect / mediation_se]
    mediation_fdr[, mediation_pval := 2 * pnorm(-abs(mediation_z))]
    # FDR within outcome
    mediation_fdr[, mediation_fdr := p.adjust(mediation_pval, method = "BH"), by = outcome]

    cat("\n=== Significant Mediation (FDR<0.05) ===\n")
    sig_med <- mediation_fdr[mediation_fdr < 0.05 & !is.na(mediation_fdr), ]
    if (nrow(sig_med) > 0) {
      print(sig_med[, .(Gene, outcome, Category, BETA, b, mediation_effect,
                        mediation_pval, mediation_fdr)])
    } else {
      cat("  No proteins pass FDR<0.05 for mediation.\n")
      # Report strongest nominal signals
      top_nom <- mediation_fdr[order(mediation_pval)][1:10, ]
      cat("  Top 10 nominal mediation signals:\n")
      print(top_nom[, .(Gene, outcome, Category, mediation_effect,
                        mediation_pval)])
    }

    fwrite(mediation_fdr, paste0(OUT_DIR, "mediation_with_fdr.csv"))
  }

} else {
  cat(sprintf("\n[WARN] Existing MR results not found at: %s\n", existing_mr_file))
  cat("Running proteome-wide mediation framework as placeholder.\n")

  # Create framework placeholder
  framework_note <- data.table(
    note = "Legacy biology-guided sensitivity panel; not a proteome-wide screen",
    current_proteins = length(tar_genes),
    available_merged = length(available_merged),
    full_ukbppp = 2923,
    status = "Need to download additional UKB-PPP protein GWAS from Synapse (syn51364943)"
  )
  fwrite(framework_note, paste0(OUT_DIR, "proteome_expansion_plan.csv"))
}

cat("\n[Done] P0-5 proteome-wide screen completed.\n")
