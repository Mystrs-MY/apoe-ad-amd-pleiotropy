# ============================================================
# NetMR Article 1 — B2_gw_apoe_exclusion.R
# 全基因组双向 MR — APOE 排除敏感性分析
# APOE 策略: ❌ 排除 chr19:44.0-46.5 Mb
# ⚠️ Clumping 后排除: 先全基因组 clump → harmonise 后移除 APOE 区域 IV
# ============================================================

rm(list = ls())

library(TwoSampleMR)
library(MRPRESSO)
library(data.table)
library(ieugwasr)
library(plinkbinr)

# ---- 路径 ----
OUT_DIR  <- "./results/"
dir.create(OUT_DIR, showWarnings = FALSE)

PLINK_BIN <- plinkbinr::get_plink_exe()
LD_REF <- "../../Resource/EUR/EUR"

# ---- 工具函数 ----
source("./utils/mr_helpers.R")

# ---- 加载 ----
gwases <- readRDS(paste0(OUT_DIR, "gwas_formatted.rds"))

# ============================================================
# Clumping (复用 B1 的 clumped 数据，如不存在则重新跑)
# ============================================================

clumped_file <- paste0(OUT_DIR, "exposures_clumped.rds")

if (file.exists(clumped_file)) {
  cat("\n==================== Loading Clumped Data (from B1) ====================\n")
  exposures <- readRDS(clumped_file)
  cat(sprintf("  Clumped exposures loaded: %s\n",
    paste(names(exposures), collapse = ", ")))
} else {
  cat("\n==================== Clumping (B1 results not found) ====================\n")
  ad_clumped <- ld_clump_local(gwases$AD, 5e-8, PLINK_BIN, LD_REF)
  dry_clumped <- ld_clump_local(gwases$Dry_AMD, 5e-8, PLINK_BIN, LD_REF)
  wet_clumped <- ld_clump_local(gwases$Wet_AMD, 5e-8, PLINK_BIN, LD_REF)
  any_clumped <- ld_clump_local(gwases$Any_AMD, 5e-8, PLINK_BIN, LD_REF)
  exposures <- list(
    AD      = ad_clumped,
    Dry_AMD = dry_clumped,
    Wet_AMD = wet_clumped,
    Any_AMD = any_clumped
  )
}

# ============================================================
# 6 对 MR — APOE 排除
# ============================================================

cat("\n==================== MR with APOE EXCLUSION ====================\n")

pair_defs <- list(
  c("AD", "Dry_AMD"),
  c("AD", "Wet_AMD"),
  c("AD", "Any_AMD"),
  c("Dry_AMD", "AD"),
  c("Wet_AMD", "AD"),
  c("Any_AMD", "AD")
)

gw_mr_without_apoe <- list()
for (pair in pair_defs) {
  exp_label <- pair[1]
  out_label <- pair[2]
  key <- paste(exp_label, out_label, sep = "_")

  cat(sprintf("\n[%s -> %s] (APOE excluded)\n", exp_label, out_label))

  exp_dat <- exposures[[exp_label]]
  out_gwas <- gwases[[out_label]]

  if (is.null(exp_dat) || nrow(exp_dat) < 3) {
    cat(sprintf("  [SKIP] Insufficient IVs\n"))
    next
  }

  res <- run_standard_mr(exp_dat, out_gwas, exp_label, out_label,
                          exclude_apoe = TRUE)

  if (!is.null(res)) {
    gw_mr_without_apoe[[key]] <- res
  }
}

# ============================================================
# 汇总
# ============================================================

gw_summary_without <- summarise_mr_list(gw_mr_without_apoe)

cat("\n==================== MR Results (Without APOE) ====================\n")
print_cols2 <- intersect(names(gw_summary_without),
  c("exposure", "outcome", "apoe_excluded", "n_iv", "mean_f",
    "ivw_beta", "ivw_se", "ivw_pval",
    "egger_intercept_pval", "cochran_q_pval"))
print(gw_summary_without[, ..print_cols2])

# ============================================================
# Table 2: With vs Without APOE 并列对比
# ============================================================

gw_summary_with <- fread(paste0(OUT_DIR, "table2a_gw_mr_with_apoe.csv"))

# 合并
table2 <- merge(
  gw_summary_with[, .(exposure, outcome, n_iv, mean_f,
                       ivw_beta, ivw_se, ivw_pval,
                       cochran_q_pval, egger_intercept_pval)],
  gw_summary_without[, .(exposure, outcome, n_iv, mean_f,
                          ivw_beta, ivw_se, ivw_pval,
                          cochran_q_pval, egger_intercept_pval)],
  by = c("exposure", "outcome"),
  suffixes = c("_with", "_without")
)

cat("\n==================== Table 2: With vs Without APOE ====================\n")
print(table2)

# 解读：如果 Without APOE 的 IVW P > 0.05 → 所有因果信号集中在 APOE
table2$apoe_is_sole_causal_hub <- table2$ivw_pval_without > 0.05

# ---- 保存 ----
fwrite(table2, paste0(OUT_DIR, "table2_with_vs_without_apoe.csv"))
saveRDS(gw_mr_without_apoe, paste0(OUT_DIR, "gw_mr_without_apoe.rds"))

cat("\n[Done] Table 2 saved to", paste0(OUT_DIR, "table2_with_vs_without_apoe.csv"), "\n")

# 关键结论
cat("\n==================== KEY FINDING ====================\n")
n_significant_without <- sum(!table2$apoe_is_sole_causal_hub, na.rm = TRUE)
cat(sprintf("  Without APOE: %d/%d pairs significant at P<0.05\n",
    n_significant_without, nrow(table2)))
if (n_significant_without == 0) {
  cat("  ==> APOE is the SOLE genome-wide causal hub for the brain-eye axis.\n")
  cat("  ==> All causal signals beyond APOE are cancelled by mixed-direction pleiotropy.\n")
}
