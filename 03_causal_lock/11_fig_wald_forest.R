# ============================================================
# Fig 3d: rs429358 Wald Ratio — Dual Forest Plot
# Left: AD→AMD (forward)  |  Right: AMD→AD (reverse)
# Style: forestploter — matches SMRForestFigure.R
# ============================================================

library(data.table)
library(forestploter)
library(grid)

dt <- fread("./results/table1_wald_ratio.csv")

# Forward: AD -> AMD
dt_fwd <- dt[exposure == "AD", ]
dt_fwd$ci_lo <- dt_fwd$wald_beta - 1.96 * dt_fwd$wald_se
dt_fwd$ci_hi <- dt_fwd$wald_beta + 1.96 * dt_fwd$wald_se
dt_fwd$effect_str <- sprintf("%.3f (%.3f to %.3f)", dt_fwd$wald_beta, dt_fwd$ci_lo, dt_fwd$ci_hi)
dt_fwd$pval_str  <- formatC(dt_fwd$wald_pval, format = "e", digits = 2)
dt_fwd$label <- paste0("AD → ", dt_fwd$outcome)

# Reverse: AMD -> AD
dt_rev <- dt[exposure != "AD", ]
dt_rev$ci_lo <- dt_rev$wald_beta - 1.96 * dt_rev$wald_se
dt_rev$ci_hi <- dt_rev$wald_beta + 1.96 * dt_rev$wald_se
dt_rev$effect_str <- sprintf("%.3f (%.3f to %.3f)", dt_rev$wald_beta, dt_rev$ci_lo, dt_rev$ci_hi)
dt_rev$pval_str  <- formatC(dt_rev$wald_pval, format = "e", digits = 2)

# Align: forward by outcome, reverse by exposure
outcome_order <- c("Dry_AMD", "Wet_AMD", "Any_AMD")
dt_fwd <- dt_fwd[match(outcome_order, dt_fwd$outcome), ]
dt_rev <- dt_rev[match(outcome_order, dt_rev$exposure), ]

# Build plot data frame with dual forest columns
plot_df <- data.frame(
  `Exposure`                   = "AD",
  `Outcome`                    = c("Dry AMD", "Wet AMD", "Any AMD"),
  `AD → AMD (95% CI)`          = dt_fwd$effect_str,
  `  Forest Plot (AD → AMD)  ` = paste(rep(" ", 20), collapse = ""),
  `AMD → AD (95% CI)`          = dt_rev$effect_str,
  `  Forest Plot (AMD → AD)  ` = paste(rep(" ", 20), collapse = ""),
  check.names = FALSE
)

# Theme (SMR style)
tm <- forest_theme(
  base_size = 10,
  core = list(bg_params = list(fill = c("white", "#E2EBF5"))),
  ci_pch    = 15,
  ci_alpha  = 1,
  refline_col = "black",
  refline_lwd = 1.2,
  xaxis_lwd   = 1.2
)

# Dual forest: ci_column = c(3, 5) for two independent forest panels
p <- forest(
  plot_df,
  est      = list(dt_fwd$wald_beta, dt_rev$wald_beta),
  lower    = list(dt_fwd$ci_lo,     dt_rev$ci_lo),
  upper    = list(dt_fwd$ci_hi,     dt_rev$ci_hi),
  sizes    = 0.5,
  ci_column = c(4, 6),
  ref_line  = c(0, 0),
  xlim      = list(c(-0.30, 0.05), c(-7, 2)),
  theme     = tm
)

# Widen both forest columns to prevent tick label overlap
p <- edit_plot(p, col = 4, width = unit(9, "cm"))
p <- edit_plot(p, col = 6, width = unit(9, "cm"))

# 三线表 borders
p <- add_border(p, part = "header", row = 1, where = "top",    gp = gpar(lwd = 2, col = "black"))
p <- add_border(p, part = "header", row = 1, where = "bottom", gp = gpar(lwd = 2, col = "black"))
p <- add_border(p, part = "body",   row = 3, where = "bottom", gp = gpar(lwd = 2, col = "black"))

# Color-code CI lines: forward = blue, reverse = red
for (i in 1:3) {
  p <- edit_plot(p, row = i, col = 4, which = "ci",
                 gp = gpar(col = "#2166AC", fill = "#2166AC"))
  p <- edit_plot(p, row = i, col = 6, which = "ci",
                 gp = gpar(col = "#B2182B", fill = "#B2182B"))
}

# Output
pdf("../figures/Fig4d_Wald_Ratio_Forest.pdf", width = 18, height = 1.8)
plot(p)
dev.off()

cat("Fig 3d saved\n")
