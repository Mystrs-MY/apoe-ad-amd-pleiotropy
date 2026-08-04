# Fix HDL figures from existing TableS3d_HDL_Results.csv
# (skips expensive HDL recomputation)
library(ggplot2)
library(dplyr)

# Run from HDL/ directory (or set working directory accordingly)
# setwd("...02_genetic_arch/HDL")
OUT_FIG <- "../../figures"
OUT_TAB <- "../../tables"
dir.create(OUT_FIG, showWarnings=FALSE, recursive=TRUE)

# ---- Load existing results ----
hdl_results <- read.csv(file.path(OUT_TAB, "TableS3d_HDL_Results.csv"))
message("Loaded HDL results:")
print(hdl_results)

# ---- Name mapping ----
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

# Cap rg > 1 (HDL convergence issue for highly correlated subtypes)
hdl_display <- hdl_display %>%
  mutate(rg = pmax(pmin(rg, 1), -1))

message("\nAfter renaming & capping:")
print(hdl_display)

# ---- 1. Lower-triangle heatmap ----
base_pairs <- hdl_display %>% select(Trait1, Trait2, rg)
mirror_pairs <- base_pairs %>% select(Trait1 = Trait2, Trait2 = Trait1, rg)
traits <- unique(c(hdl_display$Trait1, hdl_display$Trait2))
self_cor <- data.frame(Trait1 = traits, Trait2 = traits, rg = 1)

heatmap_data <- bind_rows(base_pairs, mirror_pairs, self_cor) %>%
  distinct(Trait1, Trait2, .keep_all = TRUE)

lvls_x <- c("AD", "Dry AMD", "Wet AMD", "Any AMD")
lvls_y <- rev(lvls_x)

heatmap_data$Trait1 <- factor(heatmap_data$Trait1, levels = lvls_x)
heatmap_data$Trait2 <- factor(heatmap_data$Trait2, levels = lvls_y)

heatmap_data <- subset(heatmap_data, !is.na(Trait1) & !is.na(Trait2))

# Keep lower triangle + diagonal
heatmap_data <- heatmap_data %>%
  filter(as.numeric(Trait1) + as.numeric(Trait2) <= 5)

message("\nHeatmap data (", nrow(heatmap_data), " cells):")
print(heatmap_data[, c("Trait1","Trait2","rg")])

p1 <- ggplot(heatmap_data, aes(Trait1, Trait2, fill = rg)) +
  geom_tile(color = "white", size = 1) +
  scale_fill_gradient2(low = "#0085B1", mid = "white", high = "#D94E49",
                       midpoint = 0, limit = c(-1, 1),
                       name = "HDL Genetic\nCorrelation (rg)") +
  geom_text(aes(label = sprintf("%.3f", rg)), size = 5, fontface = "bold", color = "black") +
  theme_minimal(base_size = 12) +
  labs(title = "High-Definition Likelihood (HDL) Matrix",
       subtitle = "Global genetic correlation among AD and AMD subtypes",
       x = "", y = "") +
  coord_fixed() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 11, color = "black"),
    axis.text.y = element_text(face = "bold", size = 11, color = "black"),
    panel.grid = element_blank(),
    panel.background = element_blank()
  )

ggsave(file.path(OUT_FIG, "Fig3c_HDL_Triangle_Heatmap.pdf"), p1, width = 6, height = 5)
message("Fig3c saved (", round(file.size(file.path(OUT_FIG, "Fig3c_HDL_Triangle_Heatmap.pdf"))/1024, 1), " KB)")

# ---- 2. Forest plot: AD vs AMD subtypes ----
ad_forest <- hdl_display %>%
  filter(Trait1 == "AD" | Trait2 == "AD") %>%
  mutate(
    Target = ifelse(Trait1 == "AD", Trait2, Trait1),
    ci_lower = rg - 1.96 * se,
    ci_upper = rg + 1.96 * se
  )

message("\nForest plot data:")
print(ad_forest[, c("Target","rg","se","p","ci_lower","ci_upper")])

p2 <- ggplot(ad_forest, aes(x = reorder(Target, rg), y = rg)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#7F8C8D", linewidth = 0.8) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.15, linewidth = 1, color = "#2C3E50") +
  geom_point(aes(color = rg), size = 6) +
  scale_color_gradient2(low = "#0085B1", mid = "#E5E7E9", high = "#D94E49", midpoint = 0) +
  coord_flip() +
  theme_bw(base_size = 12) +
  labs(
    title = "HDL Genetic Correlation with AD",
    subtitle = "AD vs. AMD subtypes (rg ± 95% CI)",
    x = "",
    y = "Genetic Correlation (rg) ± 95% CI"
  ) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12, face = "bold", color = "black"),
    axis.title = element_text(size = 13, face = "bold"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 10, face = "italic", color = "grey30")
  )

ggsave(file.path(OUT_FIG, "Fig3d_HDL_ForestPlot.pdf"), p2, width = 7, height = 4)
message("Fig3d saved (", round(file.size(file.path(OUT_FIG, "Fig3d_HDL_ForestPlot.pdf"))/1024, 1), " KB)")

message("\nDone. Both HDL figures regenerated.\n")
