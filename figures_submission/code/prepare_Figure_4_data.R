###############################################################################
# Prepare plotting data snapshots for Figure 4.
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = TRUE)
CODE_DIR <- dirname(script_path)
ARTICLE_ROOT <- normalizePath(file.path(CODE_DIR, "..", ".."), winslash = "/", mustWork = TRUE)

copy_required <- function(src, dest_name) {
  if (!file.exists(src)) stop("Missing required source file: ", src, call. = FALSE)
  dt <- fread(src)
  out <- file.path(CODE_DIR, dest_name)
  fwrite(dt, out)
  out
}

mixer_src <- file.path(ARTICLE_ROOT, "tables", "TableS3b_MiXeR_Bivariate.csv")
hypr_src <- file.path(ARTICLE_ROOT, "02_genetic_arch", "HyPrColoc_v1",
                      "HyPrColoc_SNP_Level_Data_FullLAVA.csv")
wald_src <- file.path(ARTICLE_ROOT, "03_causal_lock", "results", "table1_wald_ratio.csv")
iso_src <- file.path(ARTICLE_ROOT, "P0_isoform", "results", "rs429358_vs_rs7412_wald.csv")
lava_sources <- c(
  `Dry AMD` = file.path(ARTICLE_ROOT, "02_genetic_arch", "LAVA", "LAVA_FullScan_AD_vs_DryAMD_Final.csv"),
  `Wet AMD` = file.path(ARTICLE_ROOT, "02_genetic_arch", "LAVA", "LAVA_FullScan_AD_vs_WetAMD_Final.csv"),
  `Any AMD` = file.path(ARTICLE_ROOT, "02_genetic_arch", "LAVA", "LAVA_FullScan_AD_vs_AnyAMD_Final.csv")
)

for (path in c(mixer_src, hypr_src, wald_src, iso_src, lava_sources)) {
  if (!file.exists(path)) stop("Missing required source file: ", path, call. = FALSE)
}

mixer <- fread(mixer_src)
fwrite(mixer, file.path(CODE_DIR, "Figure_4a_MiXeR_bivariate.csv"))

lava <- rbindlist(lapply(names(lava_sources), function(type) {
  dt <- fread(lava_sources[[type]])
  dt[, subtype := type]
  dt
}), fill = TRUE)
lava <- lava[!is.na(p) & p > 0]
lava[, logP := -log10(p)]
chr_info <- lava[, .(chr_len = max(STOP, na.rm = TRUE)), by = CHR][order(CHR)]
chr_info[, tot := cumsum(as.numeric(chr_len)) - chr_len]
axis_df <- merge(
  lava[, .(min_bp = min(START), max_bp = max(STOP)), by = CHR],
  chr_info[, .(CHR, tot)],
  by = "CHR"
)
axis_df[, center := tot + (min_bp + max_bp) / 2]
lava <- merge(lava, chr_info[, .(CHR, tot)], by = "CHR", all.x = TRUE)
lava[, BP_cum := START + tot]
fwrite(lava[, .(subtype, locus, CHR, START, STOP, BP_cum, rho, rho.lower,
                rho.upper, p, logP)],
       file.path(CODE_DIR, "Figure_4b_LAVA_fullscan.csv"))
fwrite(axis_df[, .(CHR, center)], file.path(CODE_DIR, "Figure_4b_LAVA_axis.csv"))

hypr <- fread(hypr_src)
hypr_rs429358 <- hypr[SNP == "rs429358"]
if (nrow(hypr_rs429358) != 1) stop("Expected one rs429358 row in HyPrColoc data.", call. = FALSE)
hypr_long <- rbindlist(list(
  data.table(SNP = "rs429358", Trait = "AD", beta = hypr_rs429358$BETA_AD, se = hypr_rs429358$SE_AD),
  data.table(SNP = "rs429358", Trait = "Dry AMD", beta = hypr_rs429358$BETA_Dry, se = hypr_rs429358$SE_Dry),
  data.table(SNP = "rs429358", Trait = "Wet AMD", beta = hypr_rs429358$BETA_Wet, se = hypr_rs429358$SE_Wet),
  data.table(SNP = "rs429358", Trait = "Any AMD", beta = hypr_rs429358$BETA_Any, se = hypr_rs429358$SE_Any)
))
hypr_long[, `:=`(ci_lower = beta - 1.96 * se, ci_upper = beta + 1.96 * se)]
fwrite(hypr_long, file.path(CODE_DIR, "Figure_4c_HyPrColoc_rs429358.csv"))

wald <- fread(wald_src)
fwrite(wald, file.path(CODE_DIR, "Figure_4d_rs429358_wald_ratio.csv"))

iso <- fread(iso_src)
iso <- iso[exposure == "AD"]
iso[, outcome_label := fifelse(outcome == "Dry_AMD", "Dry AMD",
                               fifelse(outcome == "Wet_AMD", "Wet AMD", "Any AMD"))]
iso[, isoform := fifelse(snp == "rs429358", "epsilon4 proxy (rs429358)",
                         "epsilon2 proxy (rs7412)")]
iso[, `:=`(ci_lower = wald_beta - 1.96 * wald_se,
           ci_upper = wald_beta + 1.96 * wald_se)]
fwrite(iso, file.path(CODE_DIR, "Figure_4e_isoform_contrast.csv"))

provenance <- data.frame(
  data_file = c(
    "Figure_4a_MiXeR_bivariate.csv",
    "Figure_4b_LAVA_fullscan.csv",
    "Figure_4b_LAVA_axis.csv",
    "Figure_4c_HyPrColoc_rs429358.csv",
    "Figure_4d_rs429358_wald_ratio.csv",
    "Figure_4e_isoform_contrast.csv"
  ),
  source = c(
    "TableS3b_MiXeR_Bivariate.csv",
    "Three LAVA_FullScan_AD_vs_*AMD_Final.csv files",
    "Derived chromosome centers from LAVA full scan",
    "HyPrColoc_SNP_Level_Data_FullLAVA.csv; rs429358 only",
    "03_causal_lock/results/table1_wald_ratio.csv",
    "P0_isoform/results/rs429358_vs_rs7412_wald.csv; AD exposure rows"
  ),
  stringsAsFactors = FALSE
)
write.csv(provenance, file.path(CODE_DIR, "Figure_4_data_provenance.csv"), row.names = FALSE)

message("Figure 4 plotting data written to: ", normalizePath(CODE_DIR, winslash = "/"))
