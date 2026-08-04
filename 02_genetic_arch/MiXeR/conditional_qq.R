# ============================================================
# Conditional QQ plot — 分箱版 (轻量矢量 PDF)
# Reads GWAS summary stats, stratifies by conditional trait,
# bins expected -log10(P) and plots mean observed per bin
# ============================================================

library(ggplot2)
library(data.table)

OUT_FIG <- "../../figures"
dir.create(OUT_FIG, showWarnings=FALSE, recursive=TRUE)

GWAS_DIR <- "../../../Resource/GWAS"

# ---- 1. Load GWAS --------------------------------------------------------
message("Loading GWAS data...")
ad   <- fread(file.path(GWAS_DIR, "AD_Wightman_cleaned_hg19.tsv.gz"),
              select=c("SNP","CHR","BP","A1","A2","BETA","SE","P","N"))
dry  <- fread(file.path(GWAS_DIR, "AMD_Dry_R12_cleaned_hg19.tsv.gz"),
              select=c("SNP","CHR","BP","A1","A2","BETA","SE","P","N"))

message(sprintf("AD SNPs: %s  |  Dry AMD SNPs: %s",
                format(nrow(ad), big.mark=","),
                format(nrow(dry), big.mark=",")))

# ---- 2. Merge on common SNPs ---------------------------------------------
setkey(ad, SNP); setkey(dry, SNP)
common <- merge(ad[, .(SNP, P_AD=P)], dry[, .(SNP, P_Dry=P)], by="SNP")
common <- common[!is.na(P_AD) & !is.na(P_Dry)]
message(sprintf("Common SNPs: %s", format(nrow(common), big.mark=",")))

# ---- 3. Binned QQ helper -------------------------------------------------
# Bin expected -log10(P) into nbins, compute mean ± SE observed
binned_qq <- function(p_primary, p_conditional, strata, nbins=200) {
  result <- data.table()
  for (s in unique(strata)) {
    idx <- which(strata == s)
    if (length(idx) < 100) next  # skip empty / tiny strata

    p_sub <- p_primary[idx]  # AD P

    # Sort by P and compute expected
    p_sorted  <- sort(p_sub)
    n         <- length(p_sorted)
    expected  <- -log10(ppoints(n))   # uniform quantiles
    observed  <- -log10(p_sorted)

    # Bin on expected axis
    bin_edges <- seq(min(expected), max(expected), length.out=nbins+1)
    bin_mid   <- (bin_edges[-1] + bin_edges[-(nbins+1)]) / 2
    bin_idx   <- findInterval(expected, bin_edges, rightmost.closed=TRUE)

    bin_mean  <- tapply(observed, bin_idx, mean, na.rm=TRUE)
    bin_sd    <- tapply(observed, bin_idx, sd, na.rm=TRUE)
    bin_n     <- tapply(observed, bin_idx, length)
    bin_se    <- bin_sd / sqrt(bin_n)

    has_bin   <- as.integer(names(bin_mean))
    result <- rbind(result, data.table(
      stratum        = s,
      n_snps         = n,
      expected       = bin_mid[has_bin],
      observed       = as.numeric(bin_mean),
      se             = as.numeric(bin_se),
      ci_lower       = as.numeric(bin_mean) - 1.96 * as.numeric(bin_se),
      ci_upper       = as.numeric(bin_mean) + 1.96 * as.numeric(bin_se)
    ))
  }
  result
}

# ---- 4. AD | Dry AMD stratification --------------------------------------
message("Computing binned QQ: AD given Dry AMD association...")

# Stratify by Dry AMD P-value
common[, stratum_dry := cut(P_Dry,
  breaks        = c(0, 1e-4, 1e-3, 1e-2, 1e-1, 1),
  labels        = c("P < 1e-4", "1e-4 < P < 0.001", "0.001 < P < 0.01",
                     "0.01 < P < 0.1", "P > 0.1"),
  include.lowest = TRUE)]

# Reverse order for legend (most significant on top)
common[, stratum_dry := factor(stratum_dry,
  levels = rev(c("P < 1e-4", "1e-4 < P < 0.001", "0.001 < P < 0.01",
                 "0.01 < P < 0.1", "P > 0.1")))]

# Generate binned data
qq_ad_given_dry <- binned_qq(common$P_AD, common$P_Dry, common$stratum_dry)

# Cap axes at max expected -log10(P) to avoid coord_fixed distortion from
# extreme APOE signal (P ~ 1e-300 → -log10 ≈ 300 vs max expected ≈ 7)
max_axis <- max(qq_ad_given_dry$expected, na.rm=TRUE) * 1.05
n_trunc <- sum(-log10(common$P_AD) > max_axis)
trunc_note <- sprintf("(%s SNPs with -log10(P) > %.0f truncated)", format(n_trunc, big.mark=","), max_axis)

message(sprintf("Axis cap: %.1f  |  Truncated: %s SNPs", max_axis, format(n_trunc, big.mark=",")))

# Colors: warm (red) = most significant Dry stratum
stratum_colors <- c(
  "P < 1e-4"               = "#D73027",
  "1e-4 < P < 0.001"       = "#FC8D59",
  "0.001 < P < 0.01"       = "#FEE090",
  "0.01 < P < 0.1"         = "#91BFDB",
  "P > 0.1"                = "#4575B4"
)

p1 <- ggplot(qq_ad_given_dry, aes(x=expected, y=observed, color=stratum, fill=stratum)) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="grey40", linewidth=0.5) +
  geom_ribbon(aes(ymin=ci_lower, ymax=ci_upper), alpha=0.10, color=NA) +
  geom_line(linewidth=0.7) +
  scale_color_manual(values=stratum_colors, name="Dry AMD P stratum") +
  scale_fill_manual(values=stratum_colors, guide="none") +
  coord_cartesian(xlim=c(0, max_axis), ylim=c(0, max_axis)) +
  labs(
    x        = expression("Expected " * -log[10](italic(P))),
    y        = expression("Observed " * -log[10](italic(P)) ~ " (AD)")
  ) +
  theme_bw(base_size=12) +
  theme(
    legend.position  = c(0.85, 0.15),
    legend.background = element_rect(fill="white", color="grey70", linewidth=0.3),
    legend.title     = element_text(size=9),
    legend.text      = element_text(size=8)
  )

ggsave(file.path(OUT_FIG, "Fig3a_Conditional_QQ_AD_Dry.pdf"), p1, width=7, height=6.5)
message("Fig3a saved (", format(file.size(file.path(OUT_FIG, "Fig3a_Conditional_QQ_AD_Dry.pdf"))/1024, digits=3), " KB)")

# ---- 5. Dry AMD | AD stratification --------------------------------------
message("Computing binned QQ: Dry AMD given AD association...")

common[, stratum_ad := cut(P_AD,
  breaks        = c(0, 1e-4, 1e-3, 1e-2, 1e-1, 1),
  labels        = c("P < 1e-4", "1e-4 < P < 0.001", "0.001 < P < 0.01",
                     "0.01 < P < 0.1", "P > 0.1"),
  include.lowest =TRUE)]

common[, stratum_ad := factor(stratum_ad,
  levels = rev(c("P < 1e-4", "1e-4 < P < 0.001", "0.001 < P < 0.01",
                 "0.01 < P < 0.1", "P > 0.1")))]

qq_dry_given_ad <- binned_qq(common$P_Dry, common$P_AD, common$stratum_ad)

max_axis2 <- max(qq_dry_given_ad$expected, na.rm=TRUE) * 1.05
n_trunc2 <- sum(-log10(common$P_Dry) > max_axis2)
trunc_note2 <- sprintf("(%s SNPs with -log10(P) > %.0f truncated)", format(n_trunc2, big.mark=","), max_axis2)

message(sprintf("Axis cap: %.1f  |  Truncated: %s SNPs", max_axis2, format(n_trunc2, big.mark=",")))

p2 <- ggplot(qq_dry_given_ad, aes(x=expected, y=observed, color=stratum, fill=stratum)) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="grey40", linewidth=0.5) +
  geom_ribbon(aes(ymin=ci_lower, ymax=ci_upper), alpha=0.10, color=NA) +
  geom_line(linewidth=0.7) +
  scale_color_manual(values=stratum_colors, name="AD P stratum") +
  scale_fill_manual(values=stratum_colors, guide="none") +
  coord_cartesian(xlim=c(0, max_axis2), ylim=c(0, max_axis2)) +
  labs(
    x        = expression("Expected " * -log[10](italic(P))),
    y        = expression("Observed " * -log[10](italic(P)) ~ " (Dry AMD)")
  ) +
  theme_bw(base_size=12) +
  theme(
    legend.position  = c(0.85, 0.15),
    legend.background = element_rect(fill="white", color="grey70", linewidth=0.3),
    legend.title     = element_text(size=9),
    legend.text      = element_text(size=8)
  )

ggsave(file.path(OUT_FIG, "Fig3b_Conditional_QQ_Dry_AD.pdf"), p2, width=7, height=6.5)
message("Fig3b saved (", format(file.size(file.path(OUT_FIG, "Fig3b_Conditional_QQ_Dry_AD.pdf"))/1024, digits=3), " KB)")

# ---- 6. Stratum SNP counts ------------------------------------------------
message("\n--- Stratum SNP counts: AD | Dry ---")
print(common[, .(N_SNPs = .N), by=stratum_dry])

message("\n--- Stratum SNP counts: Dry | AD ---")
print(common[, .(N_SNPs = .N), by=stratum_ad])

message("\nDone. Both conditional QQ plots saved.\n")
