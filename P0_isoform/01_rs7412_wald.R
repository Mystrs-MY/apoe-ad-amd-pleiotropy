# ============================================================
# P0-1: APOE Isoform Decomposition — rs7412 Wald Ratio
# rs7412-T = APOE ε2 定义性等位基因
# 预期: ε2 AD↓ (保护), AMD↑ (风险) — 与 ε4 镜像相反
# ============================================================

rm(list = ls())
library(TwoSampleMR)
library(data.table)
library(ggplot2)

# ---- 路径 ----
GWAS_DIR <- "../../Resource/GWAS/"
OUT_DIR  <- "./results/"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- 工具函数 ----
source("../../Article_1/03_causal_lock/utils/mr_helpers.R")

# ---- 加载 GWAS ----
cat("\n=== Loading GWAS ===\n")
ad_gwas   <- load_local_gwas(paste0(GWAS_DIR, "AD_Wightman_cleaned_hg19.tsv.gz"), "AD")
dry_gwas  <- load_local_gwas(paste0(GWAS_DIR, "AMD_Dry_R12_cleaned_hg19.tsv.gz"), "Dry_AMD")
wet_gwas  <- load_local_gwas(paste0(GWAS_DIR, "AMD_Wet_R12_cleaned_hg19.tsv.gz"), "Wet_AMD")
any_gwas  <- load_local_gwas(paste0(GWAS_DIR, "AMD_H7_R12_cleaned_hg19.tsv.gz"), "Any_AMD")

gwases <- list(AD = ad_gwas, Dry_AMD = dry_gwas, Wet_AMD = wet_gwas, Any_AMD = any_gwas)

# ---- 验证 rs7412 存在 ----
cat("\n=== rs7412 Verification ===\n")
for (nm in names(gwases)) {
  row <- gwases[[nm]][gwases[[nm]]$SNP == "rs7412", ]
  if (nrow(row) == 0) {
    cat(sprintf("[WARN] %s: rs7412 NOT FOUND!\n", nm))
  } else {
    cat(sprintf("  %-10s: A1=%s, A2=%s, EAF=%.3f, BETA=%+.4f, SE=%.4f, P=%.2e\n",
      nm, row$effect_allele.exposure, row$other_allele.exposure,
      row$eaf.exposure, row$beta.exposure, row$se.exposure, row$pval.exposure))
  }
}

# ---- rs7412 Wald Ratio (6对) ----
cat("\n=== rs7412 Wald Ratio: AD -> AMD ===\n")

run_wald_snp <- function(exp_gwas, out_gwas, snp, exp_label, out_label) {
  exp_snp <- exp_gwas[exp_gwas$SNP == snp, ]
  out_snp <- out_gwas[out_gwas$SNP == snp, ]
  if (nrow(exp_snp) == 0 || nrow(out_snp) == 0) return(NULL)

  exp_snp$id.exposure <- exp_label; exp_snp$exposure <- exp_label
  out_formatted <- out_snp
  colnames(out_formatted) <- gsub("exposure", "outcome", colnames(out_formatted))
  out_formatted$outcome <- out_label; out_formatted$id.outcome <- out_label

  dat <- harmonise_data(exp_snp, out_formatted, action = 2)
  if (!dat$mr_keep) return(NULL)

  mr_res <- mr(dat, method_list = c("mr_wald_ratio"))
  wald_manual <- dat$beta.outcome / dat$beta.exposure

  data.table(
    snp = snp, exposure = exp_label, outcome = out_label,
    beta_exposure = dat$beta.exposure, se_exposure = dat$se.exposure,
    beta_outcome = dat$beta.outcome, se_outcome = dat$se.outcome,
    wald_beta = mr_res$b, wald_se = mr_res$se, wald_pval = mr_res$pval,
    wald_manual = wald_manual,
    exp_a1 = dat$effect_allele.exposure, out_a1 = dat$effect_allele.outcome,
    alleles_flipped = dat$effect_allele.exposure != dat$effect_allele.outcome
  )
}

# AD -> 3 AMD
rs7412_ad_dry <- run_wald_snp(ad_gwas, dry_gwas, "rs7412", "AD", "Dry_AMD")
rs7412_ad_wet <- run_wald_snp(ad_gwas, wet_gwas, "rs7412", "AD", "Wet_AMD")
rs7412_ad_any <- run_wald_snp(ad_gwas, any_gwas, "rs7412", "AD", "Any_AMD")

# 反向: 3 AMD -> AD
rs7412_dry_ad <- run_wald_snp(dry_gwas, ad_gwas, "rs7412", "Dry_AMD", "AD")
rs7412_wet_ad <- run_wald_snp(wet_gwas, ad_gwas, "rs7412", "Wet_AMD", "AD")
rs7412_any_ad <- run_wald_snp(any_gwas, ad_gwas, "rs7412", "Any_AMD", "AD")

# ---- 合并 rs429358 + rs7412 结果 ----
rs7412_all <- rbindlist(list(
  rs7412_ad_dry, rs7412_ad_wet, rs7412_ad_any,
  rs7412_dry_ad, rs7412_wet_ad, rs7412_any_ad
), fill = TRUE)

# 加载 rs429358 已有结果用于对比
rs429358_all <- fread("../03_causal_lock/results/table1_wald_ratio.csv")

# ---- 合并双 SNP 对比 ----
combined <- rbind(
  rs429358_all[, .(snp, exposure, outcome, wald_beta, wald_se, wald_pval)],
  rs7412_all[, .(snp, exposure, outcome, wald_beta, wald_se, wald_pval)]
)
combined[, direction := fifelse(wald_beta < 0, "Negative (Antagonistic)", "Positive (Synergistic)")]

cat("\n=== rs429358 vs rs7412 Wald Ratio Comparison ===\n")
print(combined)

# ---- 保存 ----
fwrite(rs7412_all, paste0(OUT_DIR, "rs7412_wald_ratio.csv"))
fwrite(combined, paste0(OUT_DIR, "rs429358_vs_rs7412_wald.csv"))

# ---- Isoform 方向对比图 ----
# 聚焦 AD→AMD 方向，展示 ε4 (rs429358-C) vs ε2 (rs7412-T) 的方向对比
plot_data <- combined[exposure == "AD", ]
plot_data[, isoform := fifelse(snp == "rs429358", "ε4 (rs429358-C)", "ε2 (rs7412-T)")]
plot_data[, outcome_label := fifelse(outcome == "Dry_AMD", "Dry AMD",
                              fifelse(outcome == "Wet_AMD", "Wet AMD", "Any AMD"))]

p <- ggplot(plot_data, aes(x = outcome_label, y = wald_beta, color = isoform, shape = isoform)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 4, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = wald_beta - 1.96 * wald_se, ymax = wald_beta + 1.96 * wald_se),
                width = 0.2, position = position_dodge(width = 0.5)) +
  scale_color_manual(values = c("ε4 (rs429358-C)" = "#B2182B", "ε2 (rs7412-T)" = "#2166AC")) +
  labs(
    title = "APOE Isoform-Specific Antagonistic Pleiotropy",
    subtitle = "rs429358-C (ε4) AD↑/AMD↓ vs rs7412-T (ε2) AD↓/AMD↑ — Mirror Opposite Directions",
    x = "AMD Subtype", y = "Wald Ratio β (AD → AMD)",
    color = "APOE Isoform", shape = "APOE Isoform"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

ggsave(paste0(OUT_DIR, "Fig_isoform_contrast.pdf"), p, width = 8, height = 6)

# ---- 关键定量报告 ----
cat("\n========================================\n")
cat("KEY TAKEAWAY: APOE ε2 vs ε4 Mirror Antagonism\n")
cat("========================================\n")
for (i in 1:nrow(rs7412_all)) {
  cat(sprintf("%s -> %s: Wald β=%.4f, SE=%.4f, P=%.2e\n",
    rs7412_all$exposure[i], rs7412_all$outcome[i],
    rs7412_all$wald_beta[i], rs7412_all$wald_se[i], rs7412_all$wald_pval[i]))
}
cat("\n[Done] Full results in", OUT_DIR, "\n")
