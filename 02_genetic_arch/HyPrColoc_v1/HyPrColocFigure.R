# =====================================================================
# 脚本名称: 05_HyPrColoc_EffectSize_Antagonism_Plot.R
# 目标：可视化 APOE 双子星位点在 AD 与 3 种 AMD 中的拮抗多效性
# 升级：应用顶级红蓝渐变自定义调色板 (区分不同 AMD 亚型)
# =====================================================================

library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)

message("--- 正在加载 HyPrColoc 终极对齐数据 ---")
# 1. 读取绝对对齐数据
dt <- fread("HyPrColoc_SNP_Level_Data_FullLAVA.csv")

# 2. 提取 rs429358（APOE ε4 定义性变异）
target_snps <- c("rs429358")

dt_sub <- dt[SNP %in% target_snps]
if(nrow(dt_sub) < 1) stop("🚨 未找到 rs429358，请检查数据！")

# 3. 数据重构 (长表转换)
dt_beta <- dt_sub[, .(SNP, AD = BETA_AD, Dry_AMD = BETA_Dry, Wet_AMD = BETA_Wet, Any_AMD = BETA_Any)]
melt_beta <- melt(dt_beta, id.vars = "SNP", variable.name = "Trait", value.name = "BETA")

dt_se <- dt_sub[, .(SNP, AD = SE_AD, Dry_AMD = SE_Dry, Wet_AMD = SE_Wet, Any_AMD = SE_Any)]
melt_se <- melt(dt_se, id.vars = "SNP", variable.name = "Trait", value.name = "SE")

# 合并并计算 95% 置信区间
plot_df <- merge(melt_beta, melt_se, by = c("SNP", "Trait"))
plot_df <- plot_df %>%
  mutate(
    Lower = BETA - 1.96 * SE,
    Upper = BETA + 1.96 * SE
  )

# 固定疾病的显示顺序 (决定了 y 轴从上到下的顺序)
plot_df$Trait <- factor(plot_df$Trait, levels = c("AD", "Dry_AMD", "Wet_AMD", "Any_AMD"))

# =====================================================================
# 4. 定义你的专属高级调色板
# =====================================================================
# 严格按照你给定的十六进制代码进行映射
custom_colors <- c(
  "AD"      = "#B2182B",  # 深红色 (致病风险)
  "Dry_AMD" = "#4393C3",  # 浅蓝色
  "Wet_AMD" = "#2166AC",  # 中蓝色
  "Any_AMD" = "#053061"   # 深蓝色
)

# 为了让图例文字更优美，这里再做一个标签映射字典
trait_labels <- c(
  "AD"      = "Alzheimer's Disease",
  "Dry_AMD" = "Dry AMD",
  "Wet_AMD" = "Wet AMD",
  "Any_AMD" = "Any AMD"
)


# =====================================================================
# 5. 绘制终极对冲森林图
# =====================================================================
message("--- 正在渲染极化森林图 ---")
p_antagonism <- ggplot(plot_df, aes(x = BETA, y = reorder(Trait, desc(Trait)), color = Trait)) +

  # 添加零点基准线（虚线）
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey30", linewidth = 0.5) +

  # 绘制置信区间误差棒
  geom_errorbar(aes(xmin = Lower, xmax = Upper), width = 0.2, linewidth = 0.8) +

  # 绘制效应量点 (带有高亮边框和纯白填充)
  geom_point(size = 2.5, shape = 21, fill = "white", stroke = 1.2) +

  # 应用你的专属调色板，并修改图例标签
  scale_color_manual(values = custom_colors, labels = trait_labels) +

  theme_bw(base_size = 22) +
  labs(x = "Effect Size (Aligned BETA +/- 95% CI)", y = "") +
  theme(
    axis.text.x = element_text(size = 18, color = "black"),
    axis.title.x = element_text(size = 20, face = "bold", margin = margin(t = 16)),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 18, face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  # 替换 y 轴显示的刻度标签
  scale_y_discrete(labels = trait_labels)

OUT_FIG <- "../../figures"
dir.create(OUT_FIG, showWarnings=FALSE, recursive=TRUE)

print(p_antagonism)
ggsave(file.path(OUT_FIG, "Fig4c_HyPrColoc_APOE_Antagonism.pdf"), p_antagonism, width=4.8, height=4.8)

message("Fig3c saved to figures/")