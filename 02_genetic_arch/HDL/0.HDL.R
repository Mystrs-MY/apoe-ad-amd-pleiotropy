# =====================================================================
# 脚本名称: 0.HDL_Ultimate_Pipeline.R
# 核心功能: 执行高精度全局遗传相关性分析 (HDL)
# 架构优势:
#   1. 无缝对接 Ultimate GWAS Factory 产出的标准 hg19 数据
#   2. 内存优化：动态推导 Z-score 并释放冗余列
#   3. 视觉统一：输出顶刊级别的“下三角”高对比度热图
# =====================================================================

# ---------------------------------------------------------------------
# 模块一: 环境准备与参考面板配置
# ---------------------------------------------------------------------
message("==================================================")
message("🌐 模块一: 初始化 HDL 引擎...")

# 若未安装 HDL，请取消注释下一行
# devtools::install_github("zhenin/HDL")

library(HDL)
library(data.table)
library(dplyr)
library(ggplot2)

# 🚨🚨🚨 【请确认路径】 指向你下载解压好的 UKB 参考面板文件夹
LD_path <- "UKB_imputed_SVD_eigen99_extraction"

if(!dir.exists(LD_path)) {
  stop("❌ 找不到 HDL 参考面板！请去 HDL GitHub 下载 UKB_imputed_SVD_eigen99_extraction 并正确解压。")
}

# ---------------------------------------------------------------------
# 模块二: 极速读取与 Z-score 动态转化
# ---------------------------------------------------------------------
message("\n==================================================")
message("🚀 模块二: 读取工厂数据并适配 HDL 矩阵...")

# 定义 HDL 专属读取函数 (计算 Z 值，同时极限压缩内存占用)
prep_for_hdl <- function(file_path) {
  message("   >> 载入并转化: ", basename(file_path))
  dt <- fread(file_path)

  # 核心转化：HDL 需要 SNP, A1, A2, N, Z 这 5 列
  dt[, Z := BETA / SE]

  # 仅保留 HDL 所需列，防止 R 语言内存溢出 (OOM)
  dt_hdl <- dt[, .(SNP, A1, A2, N, Z)]

  # 清理极端异常值
  dt_hdl <- dt_hdl[!is.na(Z) & is.finite(Z)]
  return(dt_hdl)
}

# 读取神级工厂刚产出的 Wightman AD 与 FinnGen R12 数据
GWAS_DIR <- "../../../Resource/GWAS"
OUT_FIG <- "../../figures"
OUT_TAB <- "../../tables"
df_ad  <- prep_for_hdl(file.path(GWAS_DIR, "AD_Wightman_cleaned_hg19.tsv.gz"))
df_dry <- prep_for_hdl(file.path(GWAS_DIR, "AMD_Dry_R12_cleaned_hg19.tsv.gz"))
df_wet <- prep_for_hdl(file.path(GWAS_DIR, "AMD_Wet_R12_cleaned_hg19.tsv.gz"))
df_h7  <- prep_for_hdl(file.path(GWAS_DIR, "AMD_H7_R12_cleaned_hg19.tsv.gz"))

message("   ✅ 所有数据已成功转化为 HDL Standard 格式！")

# ---------------------------------------------------------------------
# 模块三: 核心 HDL 交叉扫描
# ---------------------------------------------------------------------
message("\n==================================================")
message("🧬 模块三: 启动高精度本征值分解遗传相关性推断...")

run_hdl <- function(trait1_name, trait1_df, trait2_name, trait2_df) {
  message("\n   [分析中] ", trait1_name, " vs ", trait2_name, " ...")

  # 抑制 HDL 包自带的冗长 log 刷屏，保持控制台整洁
  res <- tryCatch({
    suppressMessages(
      HDL.rg(gwas1.df = trait1_df, gwas2.df = trait2_df, LD.path = LD_path)
    )
  }, error = function(e) {
    message("   ❌ 计算出错: ", e$message)
    return(NULL)
  })

  if (is.null(res)) return(NULL)

  # 如果底层没有直接暴露 p 值，手动通过 Z 分数精确计算双侧 P 值
  rg_val <- res$rg
  se_val <- res$rg.se
  p_val  <- res$rg.p
  if (is.null(p_val) || is.na(p_val)) {
    p_val <- 2 * pnorm(abs(rg_val / se_val), lower.tail = FALSE)
  }

  res_df <- data.frame(
    Trait1 = trait1_name, Trait2 = trait2_name,
    rg = rg_val, se = se_val, p = p_val
  )

  message(sprintf("   >>> 结果: rg = %.4f (se: %.4f, p = %.2e)", rg_val, se_val, p_val))
  return(res_df)
}

# 执行 6 种两两组合扫描
res_1 <- run_hdl("AD", df_ad, "Dry_AMD", df_dry)
res_2 <- run_hdl("AD", df_ad, "Wet_AMD", df_wet)
res_3 <- run_hdl("AD", df_ad, "H7_AMD",  df_h7)
res_4 <- run_hdl("Dry_AMD", df_dry, "Wet_AMD", df_wet)
res_5 <- run_hdl("Dry_AMD", df_dry, "H7_AMD",  df_h7)
res_6 <- run_hdl("Wet_AMD", df_wet, "H7_AMD",  df_h7)

# 汇总并保存
hdl_results <- bind_rows(res_1, res_2, res_3, res_4, res_5, res_6)
dir.create(OUT_TAB, showWarnings=FALSE, recursive=TRUE)
write.csv(hdl_results, file.path(OUT_TAB, "TableS3d_HDL_Results.csv"), row.names = FALSE)
message("\n   ✅ HDL 结果: TableS3d_HDL_Results.csv")

# ---------------------------------------------------------------------
# 模块四: 顶刊级下三角热图渲染 (与 LDSC 完全统一)
# ---------------------------------------------------------------------
message("\n==================================================")
message("📊 模块四: 正在渲染高精度下三角热图...")

# 1. 重命名映射 — 将 HDL 内部名称统一为显示名称
name_map <- c(
  "AD"      = "AD",
  "Dry_AMD" = "Dry AMD",
  "Wet_AMD" = "Wet AMD",
  "H7_AMD"  = "Any AMD"
)

hdl_display <- hdl_results %>%
  mutate(
    Trait1 = name_map[Trait1],
    Trait2 = name_map[Trait2]
  )

# Cap rg > 1 (HDL known issue with highly correlated subtypes) and < -1
hdl_display <- hdl_display %>%
  mutate(rg = pmax(pmin(rg, 1), -1))

# 2. 构建全矩阵
base_pairs <- hdl_display %>% select(Trait1, Trait2, rg)
mirror_pairs <- base_pairs %>% select(Trait1 = Trait2, Trait2 = Trait1, rg)
traits <- unique(c(hdl_display$Trait1, hdl_display$Trait2))
self_cor <- data.frame(Trait1 = traits, Trait2 = traits, rg = 1)

heatmap_data <- bind_rows(base_pairs, mirror_pairs, self_cor) %>%
  distinct(Trait1, Trait2, .keep_all = TRUE)

# 3. 坐标轴排序
lvls_x <- c("AD", "Dry AMD", "Wet AMD", "Any AMD")
lvls_y <- rev(lvls_x)  # AD 在底部

heatmap_data$Trait1 <- factor(heatmap_data$Trait1, levels = lvls_x)
heatmap_data$Trait2 <- factor(heatmap_data$Trait2, levels = lvls_y)

# Drop rows where factor conversion failed (safety)
heatmap_data <- heatmap_data[!is.na(Trait1) & !is.na(Trait2), ]

# 4. 🔪 核心切割：只保留下三角与对角线
heatmap_data <- heatmap_data %>%
  filter(as.numeric(Trait1) + as.numeric(Trait2) <= 5)

# 4. 绘图
p1 <- ggplot(heatmap_data, aes(Trait1, Trait2, fill = rg)) +
  geom_tile(color = "white", size = 1) +
  # HDL 算出的同类亚型 rg 极高 (甚至有时微超1)，此处放宽上限至 1.2 防止溢出破图
  scale_fill_gradient2(low = "#0085B1", mid = "white", high = "#D94E49",
                       midpoint = 0, limit = c(-1, 1.2),
                       name = "HDL Genetic\nCorrelation (rg)") +
  geom_text(aes(label = sprintf("%.3f", rg)), size = 5, fontface = "bold", color = "black") +
  theme_minimal() +
  labs(x = "", y = "") +
  coord_fixed() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 11, color = "black"),
    axis.text.y = element_text(face = "bold", size = 11, color = "black"),
    panel.grid = element_blank(),
    panel.background = element_blank()
  )

print(p1)
dir.create(OUT_FIG, showWarnings=FALSE, recursive=TRUE)
ggsave(file.path(OUT_FIG, "Fig3e_HDL_Triangle_Heatmap.pdf"), p1, width = 6, height = 5)
message("Fig3c_HDL_Triangle_Heatmap saved")

# ---------------------------------------------------------------------
# 模块五: HDL_Results_Formatted
# ---------------------------------------------------------------------

message("--- 正在绘制 Figure 1B: HDL 核心背景森林图 ---")

# 2. 读取 HDL 格式化后的结果长表
hdl_res <- read.csv(file.path(OUT_TAB, "TableS3d_HDL_Results.csv"))

# 3. 🧠 核心揉捏：筛选出所有包含 AD 的组合，并动态计算 95% 置信区间 (CI)
ad_forest_data <- hdl_res %>%
  # 筛选 AD 与其他性状的对比
  filter(Trait1 == "AD" | Trait2 == "AD") %>%
  # 动态判定谁是对比目标 (Target)
  mutate(Target = ifelse(Trait1 == "AD", Trait2, Trait1)) %>%
  # 严格依据统计学原理计算 95% 置信区间
  mutate(
    ci_lower = rg - 1.96 * se,
    ci_upper = rg + 1.96 * se
  )

# 4. 🎨 顶刊级级视觉渲染 (100% 继承你的绝美冷暖配色体系)
p2 <- ggplot(ad_forest_data, aes(x = reorder(Target, rg), y = rg)) +
  # 绘制 0 线（全局无相关性的基准线，使用高级灰虚线）
  geom_hline(yintercept = 0, linetype = "dashed", color = "#7F8C8D", size = 0.8) +

  # 绘制高对比度误差棒 (置信区间)
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.15, size = 1, color = "#2C3E50") +

  # 绘制核心效应点 (通过填入的 rg 动态映射冷暖色调，与热图无缝呼应)
  geom_point(aes(color = rg), size = 6) +

  # 配色升级：标准的红蓝双极渐变色，以 0 为中界线
  scale_color_gradient2(low = "#0085B1", mid = "#E5E7E9", high = "#D94E49", midpoint = 0) +

  coord_flip() +
  theme_bw() +
  labs(x = "", y = "Genetic Correlation (rg) ± 95% CI") +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12, face = "bold", color = "black"),
    axis.title = element_text(size = 13, face = "bold"),
    panel.grid.minor = element_blank()
  )

# 5. 打印并保存为高清矢量 PDF
print(p2)
ggsave(file.path(OUT_FIG, "Fig3f_HDL_ForestPlot.pdf"), p2, width = 7, height = 4)
message("Fig3f_HDL_ForestPlot saved")
