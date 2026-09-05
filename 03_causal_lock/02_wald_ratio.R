# ============================================================
# NetMR Article 1 — A_apoe_wald_ratio.R
# rs429358 Wald Ratio (单 SNP MR)
# APOE 策略: ✅ 保留 — rs429358 就是研究对象
# ============================================================

rm(list = ls())

library(TwoSampleMR)
library(data.table)

# ---- 路径 ----
OUT_DIR  <- "./results/"
dir.create(OUT_DIR, showWarnings = FALSE)

# ---- 工具函数 ----
source("./utils/mr_helpers.R")

# ---- 加载格式化 GWAS ----
gwases <- readRDS(file.path(OUT_DIR, "gwas_formatted.rds"))

# ============================================================
# Wald Ratio: AD → 3 AMD
# ============================================================

cat("\n==================== Wald Ratio: AD -> AMD ====================\n")

wald_ad_dry <- run_wald_ratio(gwases$AD, gwases$Dry_AMD,
                               "rs429358", "AD", "Dry_AMD")
wald_ad_wet <- run_wald_ratio(gwases$AD, gwases$Wet_AMD,
                               "rs429358", "AD", "Wet_AMD")
wald_ad_any <- run_wald_ratio(gwases$AD, gwases$Any_AMD,
                               "rs429358", "AD", "Any_AMD")

# ============================================================
# 反向: 3 AMD → AD (补充验证)
# ============================================================

cat("\n==================== Wald Ratio: AMD -> AD (reverse) ====================\n")

wald_dry_ad <- run_wald_ratio(gwases$Dry_AMD, gwases$AD,
                               "rs429358", "Dry_AMD", "AD")
wald_wet_ad <- run_wald_ratio(gwases$Wet_AMD, gwases$AD,
                               "rs429358", "Wet_AMD", "AD")
wald_any_ad <- run_wald_ratio(gwases$Any_AMD, gwases$AD,
                               "rs429358", "Any_AMD", "AD")

# ============================================================
# 汇总 Table 1
# ============================================================

wald_all <- rbindlist(list(
  wald_ad_dry$summary, wald_ad_wet$summary, wald_ad_any$summary,
  wald_dry_ad$summary, wald_wet_ad$summary, wald_any_ad$summary
), fill = TRUE)

# 方向标注
wald_all[, direction := fifelse(wald_beta < 0, "Antagonistic", "Synergistic")]
wald_all[, interpretation := fifelse(
  wald_beta < 0,
  sprintf("T allele (ε3) decreases %s risk but INCREASES %s risk",
    exposure, outcome),
  ""
)]

cat("\n==================== Table 1: rs429358 Wald Ratio ====================\n")
print(wald_all[, .(exposure, outcome, wald_beta, wald_se, wald_pval,
                    alleles_flipped, direction)])

# 手动验算对齐
cat("\n  Manual vs TwoSampleMR check:\n")
for (i in 1:nrow(wald_all)) {
  delta <- abs(wald_all$wald_beta[i] - wald_all$wald_manual[i])
  status <- if (delta < 1e-4) "OK" else "MISMATCH!"
  cat(sprintf("  %s -> %s: TwoSampleMR=%.4f, Manual=%.4f, delta=%.2e [%s]\n",
    wald_all$exposure[i], wald_all$outcome[i],
    wald_all$wald_beta[i], wald_all$wald_manual[i], delta, status))
}

# ---- 保存 ----
fwrite(wald_all, file.path(OUT_DIR, "table1_wald_ratio.csv"))
saveRDS(list(wald_ad_dry = wald_ad_dry, wald_ad_wet = wald_ad_wet,
             wald_ad_any = wald_ad_any, wald_dry_ad = wald_dry_ad,
             wald_wet_ad = wald_wet_ad, wald_any_ad = wald_any_ad),
        file.path(OUT_DIR, "wald_ratio_results.rds"))

cat("\n[Done] Table 1 saved to", file.path(OUT_DIR, "table1_wald_ratio.csv"), "\n")
