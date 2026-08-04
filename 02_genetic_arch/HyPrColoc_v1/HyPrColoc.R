# =====================================================================
# 脚本名称: 1.HyPrColoc_Full_LAVA_Module.R
# 目标：从标准化 TSV 提取数据，进行多疾病大区共定位 (AD vs 3xAMD)
# 核心升级：对接 Wightman 与 FinnGen R12 终极数据，采用 Any AMD 命名
# =====================================================================

library(data.table)
library(hyprcoloc)

# 1. 设定客观边界 (Chr 19: APOE 核心区域，基于 LAVA 显著模块)
chr_target <- 19
bp_min <- 45040933
bp_max <- 45893307

# 2. 定义【极简极速版】提取函数 (直接读取标准列)
extract_apoe_tsv <- function(file_path, trait_name) {
  message("正在提取 ", trait_name, " 数据...")

  # 直接读取我们刚才清洗好的标准化文件
  dt <- fread(file_path)

  # 锁定物理区域 (防止 CHR 类型不一，统一转换为字符或数字，这里按数字过滤)
  dt_apoe <- dt[CHR == chr_target & BP >= bp_min & BP <= bp_max]

  # 提取我们需要的 6 列核心数据
  res <- data.table(
    SNP  = dt_apoe$SNP,
    BP   = dt_apoe$BP,
    A1   = dt_apoe$A1,
    A2   = dt_apoe$A2,
    BETA = dt_apoe$BETA,
    SE   = dt_apoe$SE
  )

  # 质控：剔除缺失值
  return(res[!is.na(BETA) & !is.na(SE) & SNP != ""])
}

# 3. 运行提取 (读取刚刚数据工厂生成的最新文件)
message("\n--- 第一阶段：局部数据提取 ---")
# 🚨 注意：输入文件名已完美匹配数据工厂的输出
ad_apoe  <- extract_apoe_tsv("../GWAS数据/AD_Wightman_cleaned_hg19.tsv.gz", "AD")
dry_apoe <- extract_apoe_tsv("../GWAS数据/AMD_Dry_R12_cleaned_hg19.tsv.gz", "Dry")
wet_apoe <- extract_apoe_tsv("../GWAS数据/AMD_Wet_R12_cleaned_hg19.tsv.gz", "Wet")
any_apoe <- extract_apoe_tsv("../GWAS数据/AMD_H7_R12_cleaned_hg19.tsv.gz", "Any") # H7 已升级为 Any

# 4. 完美对齐引擎 (方向翻转与列拼接)
message("\n--- 第二阶段：执行等位基因严格对齐 ---")

# 以 AD 为基准尺
base_apoe <- ad_apoe[, .(SNP, BP = BP, A1_ref = A1, A2_ref = A2, BETA_AD = BETA, SE_AD = SE)]

add_trait_standalone <- function(dt_base, dt_trait, trait_name) {
  # 合并交集
  merged <- merge(dt_base, dt_trait[, .(SNP, A1, A2, BETA, SE)], by = "SNP")

  # 方向判断
  idx_keep <- merged$A1_ref == merged$A1 & merged$A2_ref == merged$A2
  idx_flip <- merged$A1_ref == merged$A2 & merged$A2_ref == merged$A1

  # 创建新列并处理方向
  merged[, paste0("BETA_", trait_name) := NA_real_]
  merged[, paste0("SE_", trait_name) := NA_real_]
  merged[idx_keep, paste0("BETA_", trait_name) := BETA]
  merged[idx_keep, paste0("SE_", trait_name) := SE]
  merged[idx_flip, paste0("BETA_", trait_name) := -BETA]
  merged[idx_flip, paste0("SE_", trait_name) := SE]

  # 删除 trait 原始列，保持工作区干净
  merged[, c("A1", "A2", "BETA", "SE") := NULL]

  return(merged[!is.na(get(paste0("BETA_", trait_name)))])
}

# 流水线拼接
m1 <- add_trait_standalone(base_apoe, dry_apoe, "Dry")
m2 <- add_trait_standalone(m1, wet_apoe, "Wet")
final_apoe <- add_trait_standalone(m2, any_apoe, "Any")

# 保存用于下游绘图 (LocusZoom) 的绝对对齐 SNP 级数据
write.csv(final_apoe, "HyPrColoc_SNP_Level_Data_FullLAVA.csv", row.names = FALSE)
message(">> 对齐完成！共保留 ", nrow(final_apoe), " 个重叠 SNP。")

# 5. 运行多维贝叶斯聚类 (HyPrColoc)
message("\n--- 第三阶段：启动 HyPrColoc 贝叶斯共定位 ---")

betas <- as.matrix(final_apoe[, .(BETA_AD, BETA_Dry, BETA_Wet, BETA_Any)])
ses   <- as.matrix(final_apoe[, .(SE_AD, SE_Dry, SE_Wet, SE_Any)])
rownames(betas) <- rownames(ses) <- final_apoe$SNP
# 标签命名更新为 Any_AMD
colnames(betas) <- colnames(ses) <- c("AD", "Dry_AMD", "Wet_AMD", "Any_AMD")

res <- hyprcoloc(effect.est = betas, effect.se = ses,
                 trait.names = colnames(betas),
                 snp.id = final_apoe$SNP)

# 导出汇总结果
write.csv(res$results, "HyPrColoc_4Traits_APOE_Result_Wightman_FullLAVA.csv", row.names = FALSE)

message("\n🎉 宏观疾病级共定位完成！请打开 CSV 查看聚类结果。")