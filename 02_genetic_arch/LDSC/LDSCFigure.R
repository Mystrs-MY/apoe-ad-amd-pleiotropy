# =====================================================================
# 脚本名称: 04_Figure1_LDSC_Plotting_Fixed.R
# 目标: 基于 LDSC_Results_Formatted.csv 生成 Figure 1A 与 1B
# =====================================================================

# 1. 加载必要的包
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("reshape2")) install.packages("reshape2")
if (!require("dplyr")) install.packages("dplyr")

library(ggplot2)
library(reshape2)
library(dplyr)

# 统一读取上一轮生成的标准结果文件
message("--- 读取数据并执行名称洗脱 ---")
ldsc_res <- read.csv("LDSC_Results_Formatted.csv")

# 【核心修复 1】：自动替换底层表型名称，告别下划线与黑话，彻底对齐你的 levels
ldsc_res <- ldsc_res %>%
  mutate(
    p1 = case_when(
      p1 == "Dry_AMD" ~ "Dry AMD",
      p1 == "Wet_AMD" ~ "Wet AMD",
      p1 == "H7_AMD"  ~ "Any AMD",
      TRUE ~ p1
    ),
    p2 = case_when(
      p2 == "Dry_AMD" ~ "Dry AMD",
      p2 == "Wet_AMD" ~ "Wet AMD",
      p2 == "H7_AMD"  ~ "Any AMD",
      TRUE ~ p2
    )
  )

# =====================================================================
# 图表 A：全景遗传相关性热图 (Figure 1A) - 顶刊下三角版本
# =====================================================================
message("--- 正在绘制 Figure 1A: 下三角热图 ---")

# 1. 构建全矩阵基础数据
base_pairs <- ldsc_res %>% select(Trait1 = p1, Trait2 = p2, rg)
mirror_pairs <- base_pairs %>% select(Trait1 = Trait2, Trait2 = Trait1, rg)
traits <- unique(c(ldsc_res$p1, ldsc_res$p2))
self_cor <- data.frame(Trait1 = traits, Trait2 = traits, rg = 1)

heatmap_data <- bind_rows(base_pairs, mirror_pairs, self_cor) %>%
  distinct(Trait1, Trait2, .keep_all = TRUE)

# 2. 坐标轴方向与排序 (现已与修复后的底层数据完美匹配)
lvls_x <- c("AD", "Dry AMD", "Wet AMD", "Any AMD")
lvls_y <- c("Any AMD", "Wet AMD", "Dry AMD", "AD")

heatmap_data$Trait1 <- factor(heatmap_data$Trait1, levels = lvls_x)
heatmap_data$Trait2 <- factor(heatmap_data$Trait2, levels = lvls_y)

# 3. 🔪 核心切割魔法：只保留下三角与对角线
heatmap_data <- heatmap_data %>%
  filter(as.numeric(Trait1) + as.numeric(Trait2) <= 5)

# 4. 绘图
p1 <- ggplot(heatmap_data, aes(Trait1, Trait2, fill = rg)) +
  geom_tile(color = "white", size = 1) +
  scale_fill_gradient2(low = "#0085B1", mid = "white", high = "#D94E49",
                       midpoint = 0, limit = c(-1, 1),
                       name = "Genetic\nCorrelation (rg)") +
  geom_text(aes(label = sprintf("%.3f", rg)), size = 5, fontface = "bold") +
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
ggsave("../../figures/Fig3c_LDSC_Triangle_Heatmap.pdf", p1, width = 6, height = 5)

# =====================================================================
# 图表 B：AD 核心拮抗背景森林图 (Figure 1B)
# =====================================================================
message("--- 正在绘制 Figure 1B: 森林图 ---")

# 【核心修复 2】：利用 rg 和 se 动态计算置信区间，彻底解决找不到对象的问题
ad_forest_data <- ldsc_res %>%
  filter(p1 == "AD" | p2 == "AD") %>%
  mutate(Target = ifelse(p1 == "AD", p2, p1)) %>%
  mutate(
    ci_lower = rg - 1.96 * se,
    ci_upper = rg + 1.96 * se
  )

# 绘图
p2 <- ggplot(ad_forest_data, aes(x = reorder(Target, rg), y = rg)) +
  # 绘制 0 线（基准线）
  geom_hline(yintercept = 0, linetype = "dashed", color = "#7F8C8D", size = 0.8) +
  # 绘制误差棒 (现在可以完美调用刚才计算的 ci_lower 和 ci_upper 了)
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.15, size = 1, color = "#2C3E50") +
  # 绘制点
  geom_point(aes(color = rg), size = 6) +
  # 配色：与热图统一
  scale_color_gradient2(low = "#0085B1", mid = "white", high = "#D94E49", midpoint = 0) +
  coord_flip() + # 横向展示
  theme_bw() +
  labs(x = "", y = "Genetic Correlation (rg) ± 95% CI") +
  theme(legend.position = "none",
        axis.text = element_text(size = 12, face = "bold", color = "black"),
        axis.title = element_text(size = 13, face = "bold"),
        panel.grid.minor = element_blank())

print(p2)
ggsave("../../figures/Fig3d_LDSC_ForestPlot.pdf", p2, width = 7, height = 4)

message("🎉 绘图完成！PDF 文件已保存在工作目录。")