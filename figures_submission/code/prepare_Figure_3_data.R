###############################################################################
# Prepare small plotting data snapshots for Figure 3.
# This script performs non-graphical data reduction only.
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = TRUE)
code_dir <- dirname(script_path)
ARTICLE_ROOT <- normalizePath(file.path(code_dir, "..", ".."), winslash = "/", mustWork = TRUE)

RESOURCE_ROOT <- Sys.getenv(
  "A1_RESOURCE_ROOT",
  unset = file.path(ARTICLE_ROOT, "data", "external")
)
GWAS_DIR <- file.path(RESOURCE_ROOT, "GWAS")

gwas_ad <- file.path(GWAS_DIR, "AD_Wightman_cleaned_hg19.tsv.gz")
gwas_dry <- file.path(GWAS_DIR, "AMD_Dry_R12_cleaned_hg19.tsv.gz")
ldsc_source <- file.path(ARTICLE_ROOT, "02_genetic_arch", "LDSC", "LDSC_Results_Formatted.csv")
hdl_source <- file.path(ARTICLE_ROOT, "tables", "TableS3d_HDL_Results.csv")

for (path in c(gwas_ad, gwas_dry, ldsc_source, hdl_source)) {
  if (!file.exists(path)) stop("Missing required source file: ", path, call. = FALSE)
}

p_from_z <- function(z) {
  pmax(2 * pnorm(abs(z), lower.tail = FALSE), .Machine$double.xmin)
}

binned_qq <- function(p_primary, strata, nbins = 200) {
  result <- data.table()
  for (s in levels(strata)) {
    idx <- which(strata == s)
    if (length(idx) < 100) next
    p_sorted <- sort(p_primary[idx])
    n <- length(p_sorted)
    expected <- -log10(ppoints(n))
    observed <- -log10(p_sorted)
    bin_edges <- seq(min(expected), max(expected), length.out = nbins + 1)
    bin_mid <- (bin_edges[-1] + bin_edges[-(nbins + 1)]) / 2
    bin_idx <- findInterval(expected, bin_edges, rightmost.closed = TRUE)
    bin_mean <- tapply(observed, bin_idx, mean, na.rm = TRUE)
    bin_sd <- tapply(observed, bin_idx, sd, na.rm = TRUE)
    bin_n <- tapply(observed, bin_idx, length)
    bin_se <- bin_sd / sqrt(bin_n)
    has_bin <- as.integer(names(bin_mean))
    result <- rbind(result, data.table(
      stratum = s,
      n_snps = n,
      expected = bin_mid[has_bin],
      observed = as.numeric(bin_mean),
      se = as.numeric(bin_se),
      ci_lower = as.numeric(bin_mean) - 1.96 * as.numeric(bin_se),
      ci_upper = as.numeric(bin_mean) + 1.96 * as.numeric(bin_se)
    ))
  }
  result
}

message("Reading original GWAS files for Figure 3 conditional QQ...")
ad <- fread(gwas_ad, select = c("SNP", "CHR", "BP", "P"))
dry <- fread(gwas_dry, select = c("SNP", "CHR", "BP", "P"))
setnames(ad, "P", "P_AD")
setnames(dry, "P", "P_Dry")
ad <- ad[!is.na(P_AD) & is.finite(P_AD) & P_AD > 0 & P_AD <= 1]
dry <- dry[!is.na(P_Dry) & is.finite(P_Dry) & P_Dry > 0 & P_Dry <= 1]
setkey(ad, SNP)
setkey(dry, SNP)
common <- merge(ad, dry, by = "SNP")
merge_key <- "SNP"

if (nrow(common) < 100000) {
  message("SNP ID overlap is limited (", nrow(common),
          "); retrying conditional QQ merge by CHR+BP.")
  setkey(ad, CHR, BP)
  setkey(dry, CHR, BP)
  common <- merge(ad, dry, by = c("CHR", "BP"), suffixes = c("_AD", "_Dry"))
  merge_key <- "CHR+BP"
}

if (nrow(common) < 100000) {
  stop("Insufficient overlapping variants for conditional QQ after SNP and CHR+BP merge: ",
       nrow(common), call. = FALSE)
}

stratum_levels <- c("P < 1e-4", "1e-4 to 1e-3", "1e-3 to 1e-2",
                    "1e-2 to 0.1", "P > 0.1")
common[, stratum_dry := cut(
  P_Dry,
  breaks = c(0, 1e-4, 1e-3, 1e-2, 1e-1, 1),
  labels = stratum_levels,
  include.lowest = TRUE
)]
common[, stratum_ad := cut(
  P_AD,
  breaks = c(0, 1e-4, 1e-3, 1e-2, 1e-1, 1),
  labels = stratum_levels,
  include.lowest = TRUE
)]
common[, stratum_dry := factor(stratum_dry, levels = rev(stratum_levels))]
common[, stratum_ad := factor(stratum_ad, levels = rev(stratum_levels))]

qq_ad_given_dry <- binned_qq(common$P_AD, common$stratum_dry)
qq_dry_given_ad <- binned_qq(common$P_Dry, common$stratum_ad)

fwrite(qq_ad_given_dry, file.path(code_dir, "Figure_3a_conditional_QQ_AD_given_Dry.csv"))
fwrite(qq_dry_given_ad, file.path(code_dir, "Figure_3b_conditional_QQ_Dry_given_AD.csv"))
fwrite(common[, .N, by = stratum_dry], file.path(code_dir, "Figure_3a_stratum_counts.csv"))
fwrite(common[, .N, by = stratum_ad], file.path(code_dir, "Figure_3b_stratum_counts.csv"))

ldsc <- fread(ldsc_source)
hdl <- fread(hdl_source)
fwrite(ldsc, file.path(code_dir, "Figure_3_LDSC_results.csv"))
fwrite(hdl, file.path(code_dir, "Figure_3_HDL_results.csv"))

provenance <- data.frame(
  data_file = c(
    "Figure_3a_conditional_QQ_AD_given_Dry.csv",
    "Figure_3b_conditional_QQ_Dry_given_AD.csv",
    "Figure_3_LDSC_results.csv",
    "Figure_3_HDL_results.csv"
  ),
  source = c(
    paste0("AD_Wightman_cleaned_hg19.tsv.gz + AMD_Dry_R12_cleaned_hg19.tsv.gz; original P column; merge key=", merge_key),
    paste0("AD_Wightman_cleaned_hg19.tsv.gz + AMD_Dry_R12_cleaned_hg19.tsv.gz; original P column; merge key=", merge_key),
    "LDSC_Results_Formatted.csv",
    "TableS3d_HDL_Results.csv"
  ),
  rows = c(nrow(qq_ad_given_dry), nrow(qq_dry_given_ad), nrow(ldsc), nrow(hdl)),
  stringsAsFactors = FALSE
)
write.csv(provenance, file.path(code_dir, "Figure_3_data_provenance.csv"), row.names = FALSE)

message("Figure 3 plotting data written to: ", normalizePath(code_dir, winslash = "/"))
