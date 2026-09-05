# ============================================================
# NetMR Article 1 — 00_setup_and_load.R
# 环境搭建 + 4 GWAS 加载 + 格式验证 + rs429358 QC
# ============================================================

rm(list = ls())

# ---- 安装包（首次运行取消注释）----
# packages <- c("TwoSampleMR", "MRPRESSO", "MendelianRandomization",
#               "RMediation", "MVMR", "data.table", "dplyr", "ggplot2",
#               "ggsci", "forestploter", "ieugwasr", "plinkbinr",
#               "stringr", "tidyr", "tibble", "ggpubr")
# install.packages(setdiff(packages, installed.packages()[,"Package"]))

# ---- 加载包 ----
library(TwoSampleMR)
library(MRPRESSO)
library(data.table)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ieugwasr)
library(plinkbinr)
library(ggplot2)
library(ggsci)

# ---- Reproducible paths ----
script_arg <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)][1])
SCRIPT_DIR <- dirname(normalizePath(script_file, winslash = "/", mustWork = TRUE))
PROJECT_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), winslash = "/", mustWork = TRUE)
RESOURCE_ROOT <- Sys.getenv("A1_RESOURCE_ROOT", unset = file.path(PROJECT_ROOT, "data", "external"))
GWAS_DIR <- file.path(RESOURCE_ROOT, "GWAS")
RES_DIR  <- RESOURCE_ROOT
OUT_DIR  <- file.path(SCRIPT_DIR, "results")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(SCRIPT_DIR, "figures"), showWarnings = FALSE)
dir.create(file.path(SCRIPT_DIR, "tables"), showWarnings = FALSE)

# PLINK
PLINK_BIN <- plinkbinr::get_plink_exe()
LD_REF <- file.path(RESOURCE_ROOT, "EUR", "EUR")

cat("PLINK binary:", PLINK_BIN, "\n")
cat("LD reference:", LD_REF, "\n")
cat("GWAS dir:", normalizePath(GWAS_DIR), "\n")
cat("Output dir:", normalizePath(OUT_DIR), "\n")

# ---- 加载工具函数 ----
source(file.path(SCRIPT_DIR, "utils", "mr_helpers.R"))

# ============================================================
# 加载 4 个 GWAS
# ============================================================

cat("\n==================== Loading GWAS ====================\n")

ad_gwas <- load_local_gwas(
  file.path(GWAS_DIR, "AD_Wightman_cleaned_hg19.tsv.gz"),
  "AD"
)

dry_gwas <- load_local_gwas(
  file.path(GWAS_DIR, "AMD_Dry_R12_cleaned_hg19.tsv.gz"),
  "Dry_AMD"
)

wet_gwas <- load_local_gwas(
  file.path(GWAS_DIR, "AMD_Wet_R12_cleaned_hg19.tsv.gz"),
  "Wet_AMD"
)

any_gwas <- load_local_gwas(
  file.path(GWAS_DIR, "AMD_H7_R12_cleaned_hg19.tsv.gz"),
  "Any_AMD"
)

# ---- 打包 ----
gwases <- list(
  AD      = ad_gwas,
  Dry_AMD = dry_gwas,
  Wet_AMD = wet_gwas,
  Any_AMD = any_gwas
)

# ============================================================
# QC 检查
# ============================================================

cat("\n==================== QC Checks ====================\n")

# 1. 基本统计
for (nm in names(gwases)) {
  g <- gwases[[nm]]
  cat(sprintf("  %-10s: %d SNPs, N=%.0f, BETA range=[%.3f, %.3f]\n",
    nm, nrow(g), mean(g$samplesize.exposure, na.rm = TRUE),
    min(g$beta.exposure, na.rm = TRUE), max(g$beta.exposure, na.rm = TRUE)))
}

# 2. rs429358 验证 ⚠️ 最关键 QC
verify_rs429358(gwases)

# 3. 交集 SNP 数
common_snps <- Reduce(intersect, lapply(gwases, function(x) x$SNP))
cat(sprintf("  Common SNPs across all 4 GWAS: %d\n", length(common_snps)))
cat(sprintf("  rs429358 in common: %s\n", "rs429358" %in% common_snps))

# 4. 保存
saveRDS(gwases, file.path(OUT_DIR, "gwas_formatted.rds"))
cat("\n[Done] Saved to", file.path(OUT_DIR, "gwas_formatted.rds"), "\n")

# ---- 环境变量（供后续脚本使用）----
cat("\n==================== Environment Ready ====================\n")
cat("  gwases: AD, Dry_AMD, Wet_AMD, Any_AMD\n")
cat("  PLINK_BIN:", PLINK_BIN, "\n")
cat("  LD_REF:", LD_REF, "\n")
cat("  OUT_DIR:", OUT_DIR, "\n")
cat("==========================================================\n")
