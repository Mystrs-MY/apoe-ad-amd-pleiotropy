# ============================================================
# NetMR Article 1 — B1_gw_bidirectional.R
# 全基因组双向 MR (含 APOE)
# 6 对方向: AD↔Dry/Wet/Any
# APOE 策略: 不排除（基线）
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

# ---- 加载格式化 GWAS ----
gwases <- readRDS(paste0(OUT_DIR, "gwas_formatted.rds"))

# ============================================================
# Clumping 每个暴露
# ============================================================

cat("\n==================== Clumping Exposures ====================\n")

ad_clumped <- ld_clump_local(gwases$AD, p_threshold = 5e-8,
                              plink_bin = PLINK_BIN, bfile = LD_REF)
dry_clumped <- ld_clump_local(gwases$Dry_AMD, p_threshold = 5e-8,
                               plink_bin = PLINK_BIN, bfile = LD_REF)
wet_clumped <- ld_clump_local(gwases$Wet_AMD, p_threshold = 5e-8,
                               plink_bin = PLINK_BIN, bfile = LD_REF)
any_clumped <- ld_clump_local(gwases$Any_AMD, p_threshold = 5e-8,
                               plink_bin = PLINK_BIN, bfile = LD_REF)

# ---- 打包 ----
exposures <- list(
  AD      = ad_clumped,
  Dry_AMD = dry_clumped,
  Wet_AMD = wet_clumped,
  Any_AMD = any_clumped
)

# 保存 clumped 数据（供 B2 复用，避免重复 clumping）
saveRDS(exposures, paste0(OUT_DIR, "exposures_clumped.rds"))
cat("[Saved] exposures_clumped.rds\n")

# ============================================================
# 6 对双向 MR (含 APOE)
# ============================================================

cat("\n==================== Running 6-pair Bidirectional MR ====================\n")

pair_defs <- list(
  c("AD", "Dry_AMD"),
  c("AD", "Wet_AMD"),
  c("AD", "Any_AMD"),
  c("Dry_AMD", "AD"),
  c("Wet_AMD", "AD"),
  c("Any_AMD", "AD")
)

gw_mr_with_apoe <- list()
for (pair in pair_defs) {
  exp_label <- pair[1]
  out_label <- pair[2]
  key <- paste(exp_label, out_label, sep = "_")

  cat(sprintf("\n[%s -> %s]\n", exp_label, out_label))

  exp_dat <- exposures[[exp_label]]
  out_gwas <- gwases[[out_label]]

  if (is.null(exp_dat) || nrow(exp_dat) < 3) {
    cat(sprintf("  [SKIP] Insufficient IVs for %s (%s)\n", exp_label,
        if(is.null(exp_dat)) "NULL" else nrow(exp_dat)))
    next
  }

  res <- run_standard_mr(exp_dat, out_gwas, exp_label, out_label,
                          exclude_apoe = FALSE)

  if (!is.null(res)) {
    gw_mr_with_apoe[[key]] <- res
  }
}

# ============================================================
# 汇总 Table 2 (With APOE 部分)
# ============================================================

gw_summary_with <- summarise_mr_list(gw_mr_with_apoe)

cat("\n==================== MR Results (With APOE) ====================\n")
print_cols <- intersect(names(gw_summary_with),
  c("exposure", "outcome", "n_iv", "mean_f",
    "ivw_beta", "ivw_se", "ivw_pval",
    "egger_intercept_pval", "cochran_q_pval",
    "steiger_correct_direction"))
print(gw_summary_with[, ..print_cols])

# ---- 保存 ----
fwrite(gw_summary_with, paste0(OUT_DIR, "table2a_gw_mr_with_apoe.csv"))
saveRDS(gw_mr_with_apoe, paste0(OUT_DIR, "gw_mr_with_apoe.rds"))

cat("\n[Done] Saved to", paste0(OUT_DIR, "table2a_gw_mr_with_apoe.csv"), "\n")
