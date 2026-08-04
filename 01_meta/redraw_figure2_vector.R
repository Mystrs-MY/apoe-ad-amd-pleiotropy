###############################################################################
# Figure 2 vector redraw for Article 1
# Source logic: 01_meta/generate_figures_tables.R, primary REML + HKSJ meta-analysis
# Output: native vector PDF/SVG plus R-rendered preview and QA
###############################################################################

suppressPackageStartupMessages({
  library(meta)
  library(ggplot2)
  library(svglite)
  library(ragg)
  library(pdftools)
})

script_arg <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)][1])
script_dir <- dirname(normalizePath(script_file, winslash = "/", mustWork = TRUE))
root_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)

out_dir <- file.path(root_dir, "figures_submission")
main_dir <- file.path(out_dir, "main")
preview_dir <- file.path(out_dir, "previews")
source_dir <- file.path(out_dir, "source_data")
dir.create(main_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(preview_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

figure_pdf <- file.path(main_dir, "Figure_2.pdf")
figure_svg <- file.path(main_dir, "Figure_2.svg")
figure_tiff <- file.path(main_dir, "Figure_2.tiff")
preview_png <- file.path(preview_dir, "Figure_2_preview.png")
source_csv <- file.path(source_dir, "Figure_2_source_data.csv")
qa_csv <- file.path(out_dir, "Figure_2_vector_QA.csv")

# Tight page geometry preserves the original physical content size while
# removing the unused outer canvas.
width_mm <- 163
height_mm <- 108
width_in <- width_mm / 25.4
height_in <- height_mm / 25.4

theme_set(
  theme_void(base_size = 6.9, base_family = "Arial") +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(1.2, 1.2, 1.2, 1.2, unit = "mm")
    )
)

d <- data.frame(
  study    = c("Choi 2020", "Hwang PH 2021", "Keenan 2014", "Klaver 1999",
               "Lee CS 2019", "Shang 2021", "Tsai 2015", "Wen 2021",
               "Son 2025", "Lishinsky-Fischer 2025"),
  country  = c("South Korea", "United States", "United Kingdom", "Netherlands",
               "United States", "United Kingdom", "China (Taiwan)",
               "China (Taiwan)", "South Korea", "Multi-national"),
  effect_metric = c("aHR", "HR", "Rate Ratio", "RR", "HR", "HR", "HR", "aHR", "HR", "HR"),
  hr       = c(1.48, 1.87, 0.86, 1.50, 1.20, 1.08, 1.44, 1.23, 1.09, 1.02),
  ci_lower = c(1.25, 1.13, 0.67, 0.60, 1.02, 0.80, 1.26, 1.04, 1.05, 0.89),
  ci_upper = c(1.74, 3.09, 1.08, 3.50, 1.40, 1.44, 1.64, 1.46, 1.13, 1.16),
  n_amd    = c(2213, 668, 65894, 113, 1036, 3671, 4993, 10578, 21384, 56035),
  n_ctrl   = c(306127, 2707, 7700963, 811, 2284, 87602, 24965, 10578, 1662319, 56035),
  nos      = c(9, 9, 5, 8, 9, 6, 7, 8, 8, 8),
  stringsAsFactors = FALSE
)

d$logHR <- log(d$hr)
d$SE <- (log(d$ci_upper) - log(d$ci_lower)) / (2 * 1.96)

meta_all <- metagen(
  TE = logHR,
  seTE = SE,
  studlab = study,
  data = d,
  sm = "RR",
  method.tau = "REML",
  method.random.ci = "HK",
  common = FALSE,
  random = TRUE,
  prediction = TRUE
)

d$weight_random <- meta_all$w.random
d$weight_pct <- 100 * d$weight_random / sum(d$weight_random)
d$effect_label <- sprintf("%.2f (%.2f-%.2f)", d$hr, d$ci_lower, d$ci_upper)
d$weight_label <- sprintf("%.1f", d$weight_pct)

plot_d <- d[order(d$hr), ]
plot_d$y <- seq_len(nrow(plot_d))
plot_d$y <- max(plot_d$y) - plot_d$y + 1

x_min_hr <- 0.5
x_max_hr <- 4.0
forest_left <- 0.600
forest_right <- 0.795

map_hr <- function(x) {
  forest_left +
    (log(x) - log(x_min_hr)) / (log(x_max_hr) - log(x_min_hr)) *
    (forest_right - forest_left)
}

plot_d$x_hr <- map_hr(plot_d$hr)
plot_d$x_lo <- map_hr(pmax(plot_d$ci_lower, x_min_hr))
plot_d$x_hi <- map_hr(pmin(plot_d$ci_upper, x_max_hr))
plot_d$point_size <- 1.8 + 2.1 * sqrt(plot_d$weight_pct / max(plot_d$weight_pct))

summary_row <- data.frame(
  label = "Random-effects model",
  hr = exp(meta_all$TE.random),
  ci_lower = exp(meta_all$lower.random),
  ci_upper = exp(meta_all$upper.random),
  pred_lower = exp(meta_all$lower.predict),
  pred_upper = exp(meta_all$upper.predict),
  p = meta_all$pval.random,
  I2 = 100 * meta_all$I2,
  tau2 = meta_all$tau2,
  Q = meta_all$Q,
  p_Q = meta_all$pval.Q,
  stringsAsFactors = FALSE
)

diamond_y <- 0.25
diamond_h <- 0.34
diamond <- data.frame(
  x = map_hr(c(summary_row$ci_lower, summary_row$hr,
               summary_row$ci_upper, summary_row$hr)),
  y = c(diamond_y, diamond_y + diamond_h, diamond_y, diamond_y - diamond_h)
)

tick_values <- c(0.5, 1, 2, 4)
ticks <- data.frame(
  hr = tick_values,
  x = map_hr(tick_values),
  label = as.character(tick_values)
)

fmt_int <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)

header_y <- max(plot_d$y) + 1.15
axis_y <- -0.95
note_y1 <- -1.98
note_y2 <- -2.38
ylim <- c(-2.72, header_y + 0.45)

row_guides <- data.frame(y = c(plot_d$y - 0.5, diamond_y - 0.55))

stats_line <- sprintf(
  "REML + Hartung-Knapp random effects: pooled relative effect %.2f (95%% CI %.2f-%.2f), P = %.3f",
  summary_row$hr, summary_row$ci_lower, summary_row$ci_upper, summary_row$p
)
heterogeneity_line <- sprintf(
  "Heterogeneity: I² = %.1f%%, tau² = %.3f, Q = %.2f (P < 0.001); prediction interval %.2f-%.2f",
  summary_row$I2, summary_row$tau2, summary_row$Q,
  summary_row$pred_lower, summary_row$pred_upper
)

p <- ggplot() +
  geom_segment(
    data = row_guides,
    aes(x = 0.010, xend = 0.985, y = y, yend = y),
    linewidth = 0.18,
    colour = "#E6E6E6"
  ) +
  geom_segment(
    aes(x = 0.010, xend = 0.985, y = header_y - 0.38, yend = header_y - 0.38),
    linewidth = 0.35,
    colour = "#6B6B6B"
  ) +
  annotate("text", x = 0.012, y = header_y, label = "Study",
           hjust = 0, fontface = "bold", size = 2.12, family = "Arial") +
  annotate("text", x = 0.207, y = header_y, label = "Country",
           hjust = 0, fontface = "bold", size = 2.12, family = "Arial") +
  annotate("text", x = 0.386, y = header_y, label = "AMD (N)",
           hjust = 1, fontface = "bold", size = 2.12, family = "Arial") +
  annotate("text", x = 0.500, y = header_y, label = "Control (N)",
           hjust = 1, fontface = "bold", size = 2.12, family = "Arial") +
  annotate("text", x = 0.558, y = header_y, label = "NOS",
           hjust = 0.5, fontface = "bold", size = 2.12, family = "Arial") +
  annotate("text", x = (forest_left + forest_right) / 2, y = header_y,
           label = "Relative effect", hjust = 0.5,
           fontface = "bold", size = 2.12, family = "Arial") +
  annotate("text", x = 0.875, y = header_y, label = "Estimate (95% CI)",
           hjust = 0.5, fontface = "bold", size = 2.00, family = "Arial") +
  annotate("text", x = 0.976, y = header_y, label = "Weight",
           hjust = 0.5, fontface = "bold", size = 2.00, family = "Arial") +
  geom_text(
    data = plot_d,
    aes(x = 0.012, y = y, label = study),
    hjust = 0,
    size = 2.04,
    family = "Arial",
    colour = "#202020"
  ) +
  geom_text(
    data = plot_d,
    aes(x = 0.207, y = y, label = country),
    hjust = 0,
    size = 1.98,
    family = "Arial",
    colour = "#333333"
  ) +
  geom_text(
    data = plot_d,
    aes(x = 0.386, y = y, label = fmt_int(n_amd)),
    hjust = 1,
    size = 1.98,
    family = "Arial",
    colour = "#333333"
  ) +
  geom_text(
    data = plot_d,
    aes(x = 0.500, y = y, label = fmt_int(n_ctrl)),
    hjust = 1,
    size = 1.98,
    family = "Arial",
    colour = "#333333"
  ) +
  geom_text(
    data = plot_d,
    aes(x = 0.558, y = y, label = nos),
    hjust = 0.5,
    size = 1.98,
    family = "Arial",
    colour = "#333333"
  ) +
  geom_segment(
    aes(x = map_hr(1), xend = map_hr(1), y = -0.15, yend = max(plot_d$y) + 0.45),
    linewidth = 0.24,
    colour = "#707070",
    linetype = "dashed"
  ) +
  geom_segment(
    data = plot_d,
    aes(x = x_lo, xend = x_hi, y = y, yend = y),
    linewidth = 0.42,
    colour = "#2E75B6"
  ) +
  geom_point(
    data = plot_d,
    aes(x = x_hr, y = y, size = point_size),
    shape = 15,
    fill = "#2E75B6",
    colour = "#1F4E79",
    stroke = 0.2
  ) +
  scale_size_identity() +
  geom_polygon(
    data = diamond,
    aes(x = x, y = y),
    fill = "#2E75B6",
    colour = "#1F4E79",
    linewidth = 0.30
  ) +
  geom_segment(
    aes(x = map_hr(summary_row$pred_lower), xend = map_hr(summary_row$pred_upper),
        y = -0.55, yend = -0.55),
    linewidth = 0.38,
    colour = "#B23A3A"
  ) +
  geom_segment(
    aes(x = map_hr(summary_row$pred_lower), xend = map_hr(summary_row$pred_lower),
        y = -0.70, yend = -0.40),
    linewidth = 0.32,
    colour = "#B23A3A"
  ) +
  geom_segment(
    aes(x = map_hr(summary_row$pred_upper), xend = map_hr(summary_row$pred_upper),
        y = -0.70, yend = -0.40),
    linewidth = 0.32,
    colour = "#B23A3A"
  ) +
  annotate("text", x = 0.012, y = diamond_y, label = summary_row$label,
           hjust = 0, size = 2.12, family = "Arial", fontface = "bold") +
  annotate("text", x = 0.846, y = diamond_y,
           label = sprintf("%.2f (%.2f-%.2f)",
                           summary_row$hr, summary_row$ci_lower,
                           summary_row$ci_upper),
           hjust = 0, size = 2.12, family = "Arial", fontface = "bold") +
  annotate("text", x = 0.012, y = -0.55, label = "Prediction interval",
           hjust = 0, size = 1.90, family = "Arial", colour = "#6B2424") +
  annotate("text", x = 0.846, y = -0.55,
           label = sprintf("%.2f-%.2f", summary_row$pred_lower, summary_row$pred_upper),
           hjust = 0, size = 1.90, family = "Arial", colour = "#6B2424") +
  geom_text(
    data = plot_d,
    aes(x = 0.846, y = y, label = effect_label),
    hjust = 0,
    size = 1.98,
    family = "Arial",
    colour = "#202020"
  ) +
  geom_text(
    data = plot_d,
    aes(x = 0.978, y = y, label = weight_label),
    hjust = 1,
    size = 1.98,
    family = "Arial",
    colour = "#202020"
  ) +
  geom_segment(
    aes(x = forest_left, xend = forest_right, y = axis_y, yend = axis_y),
    linewidth = 0.30,
    colour = "black"
  ) +
  geom_segment(
    data = ticks,
    aes(x = x, xend = x, y = axis_y, yend = axis_y + 0.13),
    linewidth = 0.30,
    colour = "black"
  ) +
  geom_text(
    data = ticks,
    aes(x = x, y = axis_y - 0.26, label = label),
    size = 1.90,
    family = "Arial",
    colour = "black"
  ) +
  annotate("text", x = (forest_left + forest_right) / 2, y = axis_y - 0.62,
           label = "Relative effect (log scale)",
           hjust = 0.5, size = 1.94, family = "Arial") +
  annotate("text", x = 0.012, y = note_y1, label = stats_line,
           hjust = 0, size = 1.76, family = "Arial", colour = "#333333") +
  annotate("text", x = 0.012, y = note_y2, label = heterogeneity_line,
           hjust = 0, size = 1.76, family = "Arial", colour = "#333333") +
  coord_cartesian(xlim = c(0, 1), ylim = ylim, clip = "off", expand = FALSE) +
  theme(legend.position = "none")

svglite::svglite(figure_svg, width = width_in, height = height_in)
print(p)
dev.off()

grDevices::cairo_pdf(figure_pdf, width = width_in, height = height_in,
                     family = "Arial", bg = "white")
print(p)
dev.off()

ragg::agg_tiff(figure_tiff, width = width_in, height = height_in,
               units = "in", res = 600, background = "white",
               compression = "lzw")
print(p)
dev.off()

suppressWarnings(pdftools::pdf_convert(
  figure_pdf,
  format = "png",
  dpi = 200,
  pages = 1,
  filenames = preview_png,
  verbose = FALSE
))

dir.create(file.path(root_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
file.copy(figure_pdf, file.path(root_dir, "figures", "Figure_2.pdf"), overwrite = TRUE)

plot_d_out <- plot_d[, c("study", "country", "effect_metric", "hr", "ci_lower", "ci_upper",
                          "n_amd", "n_ctrl", "nos", "weight_pct")]
names(plot_d_out)[names(plot_d_out) == "hr"] <- "relative_effect"
write.csv(plot_d_out, source_csv, row.names = FALSE)

page <- pdftools::pdf_pagesize(figure_pdf)[1, ]
fonts <- pdftools::pdf_fonts(figure_pdf)
text_data <- pdftools::pdf_data(figure_pdf)[[1]]
qa <- data.frame(
  output_pdf = normalizePath(figure_pdf, winslash = "/", mustWork = FALSE),
  output_svg = normalizePath(figure_svg, winslash = "/", mustWork = FALSE),
  preview_png = normalizePath(preview_png, winslash = "/", mustWork = FALSE),
  width_mm = round(as.numeric(page[["width"]]) * 25.4 / 72, 1),
  height_mm = round(as.numeric(page[["height"]]) * 25.4 / 72, 1),
  text_objects = nrow(text_data),
  embedded_fonts = paste(unique(fonts$name), collapse = "; "),
  backend = "R",
  vector_pdf = TRUE,
  source_model = "meta::metagen(REML + Hartung-Knapp, random effects)",
  stringsAsFactors = FALSE
)
write.csv(qa, qa_csv, row.names = FALSE)

cat("Figure 2 vector redraw complete.\n")
cat("PDF: ", normalizePath(figure_pdf, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("SVG: ", normalizePath(figure_svg, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("Preview: ", normalizePath(preview_png, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("QA: ", normalizePath(qa_csv, winslash = "/", mustWork = FALSE), "\n", sep = "")
print(qa)
