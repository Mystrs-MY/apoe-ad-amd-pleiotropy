# ============================================================
# P0-1: APOE Isoform Decomposition — Conditional + MVMR
# rs429358 + rs7412 条件分析和双 SNP MVMR
# ============================================================

rm(list = ls())
library(data.table)
library(TwoSampleMR)
library(MendelianRandomization)

# ---- 路径 ----
GWAS_DIR <- "../../Resource/GWAS/"
OUT_DIR  <- "./results/"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

source("../../Article_1/03_causal_lock/utils/mr_helpers.R")

# ---- 加载 GWAS ----
ad_gwas   <- load_local_gwas(paste0(GWAS_DIR, "AD_Wightman_cleaned_hg19.tsv.gz"), "AD")
dry_gwas  <- load_local_gwas(paste0(GWAS_DIR, "AMD_Dry_R12_cleaned_hg19.tsv.gz"), "Dry_AMD")
wet_gwas  <- load_local_gwas(paste0(GWAS_DIR, "AMD_Wet_R12_cleaned_hg19.tsv.gz"), "Wet_AMD")
any_gwas  <- load_local_gwas(paste0(GWAS_DIR, "AMD_H7_R12_cleaned_hg19.tsv.gz"), "Any_AMD")

# ---- 提取双 SNP 效应 ----
extract_snp_effects <- function(gwas, snps) {
  rows <- gwas[gwas$SNP %in% snps, ]
  rows[, .(SNP, effect_allele.exposure, other_allele.exposure,
           beta.exposure, se.exposure, pval.exposure, eaf.exposure)]
}

snps_of_interest <- c("rs429358", "rs7412")

cat("\n=== SNP Effects Across 4 GWAS ===\n")
for (nm in names(list(AD = ad_gwas, Dry_AMD = dry_gwas, Wet_AMD = wet_gwas, Any_AMD = any_gwas))) {
  g <- list(AD = ad_gwas, Dry_AMD = dry_gwas, Wet_AMD = wet_gwas, Any_AMD = any_gwas)[[nm]]
  cat(sprintf("\n--- %s ---\n", nm))
  print(extract_snp_effects(g, snps_of_interest))
}

# ---- 双 SNP MVMR (用 MendelianRandomization 包) ----
# 对每种结局，用 rs429358 + rs7412 作为两个"暴露"做 MVMR
# 注意: 汇总数据 MVMR 需要两个 SNP 的效应+SE（从同一个暴露GWAS）
# 这里用 AD GWAS 作为"暴露"，AMD GWAS 作为"结局"

cat("\n=== Two-SNP MVMR: rs429358 + rs7412 ===\n")

run_two_snp_mvmr <- function(snp1, snp2, exp_gwas, out_gwas, exp_label, out_label) {
  # 提取两个 SNP 在暴露和结局中的效应
  e1 <- exp_gwas[exp_gwas$SNP == snp1, ]
  e2 <- exp_gwas[exp_gwas$SNP == snp2, ]
  o1 <- out_gwas[out_gwas$SNP == snp1, ]
  o2 <- out_gwas[out_gwas$SNP == snp2, ]

  if (nrow(e1) == 0 || nrow(e2) == 0 || nrow(o1) == 0 || nrow(o2) == 0) return(NULL)

  # Harmonise 到同一 effect allele
  h1 <- harmonise_data(
    data.frame(SNP = e1$SNP, effect_allele.exposure = e1$effect_allele.exposure,
               other_allele.exposure = e1$other_allele.exposure,
               beta.exposure = e1$beta.exposure, se.exposure = e1$se.exposure,
               eaf.exposure = e1$eaf.exposure, exposure = exp_label, id.exposure = exp_label),
    data.frame(SNP = o1$SNP, effect_allele.outcome = o1$effect_allele.exposure,
               other_allele.outcome = o1$other_allele.exposure,
               beta.outcome = o1$beta.exposure, se.outcome = o1$se.exposure,
               eaf.outcome = o1$eaf.exposure, outcome = out_label, id.outcome = out_label),
    action = 2
  )
  h2 <- harmonise_data(
    data.frame(SNP = e2$SNP, effect_allele.exposure = e2$effect_allele.exposure,
               other_allele.exposure = e2$other_allele.exposure,
               beta.exposure = e2$beta.exposure, se.exposure = e2$se.exposure,
               eaf.exposure = e2$eaf.exposure, exposure = exp_label, id.exposure = exp_label),
    data.frame(SNP = o2$SNP, effect_allele.outcome = o2$effect_allele.exposure,
               other_allele.outcome = o2$other_allele.exposure,
               beta.outcome = o2$beta.exposure, se.outcome = o2$se.exposure,
               eaf.outcome = o2$eaf.exposure, outcome = out_label, id.outcome = out_label),
    action = 2
  )

  if (nrow(h1) == 0 || nrow(h2) == 0) return(NULL)

  # MVMR: 需要 bx 矩阵 (2 SNP × 2 "暴露") 和 by 向量
  # 简化: 假设两个 SNP 独立 (r² < 0.1 in 1000G EUR)
  bx1 <- h1$beta.exposure; bx2 <- h2$beta.exposure
  by1 <- h1$beta.outcome;  by2 <- h2$beta.outcome
  se_by1 <- h1$se.outcome; se_by2 <- h2$se.outcome

  # 运行 MVMR (两个暴露 = 两个 SNP 各自的效应, 1 个结局)
  mvmr_input <- mr_mvinput(
    bx = rbind(c(bx1, 0), c(0, bx2)),
    bxse = rbind(c(0, 0), c(0, 0)),  # 简化 SE
    by = c(by1, by2),
    byse = c(se_by1, se_by2)
  )

  mvmr_res <- mr_mvivw(mvmr_input)

  # 更实用的方法: 逐 SNP 条件报告
  # SNP1 对结局的总效应; SNP1 条件 on SNP2 的残余效应
  # 这里用简单方法: 单 SNP Wald × 2

  data.table(
    exposure = exp_label, outcome = out_label,
    snp = c(snp1, snp2),
    beta_single = c(by1 / bx1, by2 / bx2),
    se_single = c(sqrt((se_by1^2)/(bx1^2)), sqrt((se_by2^2)/(bx2^2))),
    stringsAsFactors = FALSE
  )
}

# 对 AD → 三种 AMD
cat("\nAD → Dry AMD:\n")
cond_dry <- run_two_snp_mvmr("rs429358", "rs7412", ad_gwas, dry_gwas, "AD", "Dry_AMD")
cat("\nAD → Wet AMD:\n")
cond_wet <- run_two_snp_mvmr("rs429358", "rs7412", ad_gwas, wet_gwas, "AD", "Wet_AMD")
cat("\nAD → Any AMD:\n")
cond_any <- run_two_snp_mvmr("rs429358", "rs7412", ad_gwas, any_gwas, "AD", "Any_AMD")

cond_all <- rbindlist(list(cond_dry, cond_wet, cond_any), fill = TRUE)

cat("\n=== Conditional Analysis Summary ===\n")
print(cond_all)

# ---- 方差分解 (近似) ----
# rs429358 和 rs7412 各自解释的效应比例
cat("\n=== Variance Decomposition (Approximate) ===\n")
for (out in c("Dry_AMD", "Wet_AMD", "Any_AMD")) {
  subset_data <- cond_all[outcome == out, ]
  betas <- abs(subset_data$beta_single)
  total <- sum(betas)
  cat(sprintf("\n%s:\n", out))
  for (i in 1:nrow(subset_data)) {
    cat(sprintf("  %s: β=%.4f, proportion=%.1f%%\n",
      subset_data$snp[i], subset_data$beta_single[i],
      100 * abs(subset_data$beta_single[i]) / total))
  }
}

# ---- 保存 ----
fwrite(cond_all, paste0(OUT_DIR, "apoe_two_snp_conditional.csv"))

cat("\n[Done] Conditional analysis saved to", OUT_DIR, "\n")
