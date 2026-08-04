#!/usr/bin/env Rscript

rm(list = ls())

required <- c("meta", "metafor", "pdftools")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R package(s): ", paste(missing, collapse = ", "))

suppressPackageStartupMessages({
  library(meta)
  library(metafor)
  library(pdftools)
})

script_arg <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)][1])
CODE_DIR <- dirname(normalizePath(script_file, winslash = "/", mustWork = TRUE))
SUBMISSION_DIR <- normalizePath(file.path(CODE_DIR, ".."), winslash = "/", mustWork = TRUE)
PROJECT_ROOT <- normalizePath(file.path(SUBMISSION_DIR, ".."), winslash = "/", mustWork = TRUE)
FIG_DIR <- file.path(SUBMISSION_DIR, "supplementary_figures")
PREVIEW_DIR <- file.path(SUBMISSION_DIR, "previews")
SOURCE_DIR <- file.path(SUBMISSION_DIR, "source_data")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PREVIEW_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SOURCE_DIR, recursive = TRUE, showWarnings = FALSE)

d <- data.frame(
  study = c(
    "Choi 2020", "Hwang PH 2021", "Keenan 2014", "Klaver 1999",
    "Lee CS 2019", "Shang 2021", "Tsai 2015", "Wen 2021",
    "Son 2025", "Lishinsky-Fischer 2025"
  ),
  country = c(
    "South Korea", "United States", "United Kingdom", "Netherlands",
    "United States", "United Kingdom", "China (Taiwan)", "China (Taiwan)",
    "South Korea", "Multi-national"
  ),
  effect_metric = c("aHR", "HR", "Rate Ratio", "RR", "HR", "HR", "HR", "aHR", "HR", "HR"),
  relative_effect = c(1.48, 1.87, 0.86, 1.50, 1.20, 1.08, 1.44, 1.23, 1.09, 1.02),
  ci_lower = c(1.25, 1.13, 0.67, 0.60, 1.02, 0.80, 1.26, 1.04, 1.05, 0.89),
  ci_upper = c(1.74, 3.09, 1.08, 3.50, 1.40, 1.44, 1.64, 1.46, 1.13, 1.16),
  n_amd = c(2213, 668, 65894, 113, 1036, 3671, 4993, 10578, 21384, 56035),
  n_ctrl = c(306127, 2707, 7700963, 811, 2284, 87602, 24965, 10578, 1662319, 56035),
  region = c(
    "Asia", "North America", "Europe", "Europe", "North America",
    "Europe", "Asia", "Asia", "Asia", "Multi-national"
  ),
  setting = c(
    "Community", "Community", "Hospital/treated", "Community", "Community",
    "Community", "Community", "Community", "Community", "Community"
  ),
  stringsAsFactors = FALSE
)
d$log_effect <- log(d$relative_effect)
d$SE <- (log(d$ci_upper) - log(d$ci_lower)) / (2 * 1.96)

fit_meta <- function(dat) {
  metagen(
    TE = log_effect, seTE = SE, studlab = study, data = dat,
    sm = "RR", method.tau = "REML", method.random.ci = "HK",
    common = FALSE, random = TRUE, prediction = TRUE
  )
}

meta_all <- fit_meta(d)
sg_setting <- update(meta_all, subgroup = setting)
sg_region <- update(meta_all, subgroup = region)
no_keenan <- d[d$study != "Keenan 2014", , drop = FALSE]
meta_no_keenan <- fit_meta(no_keenan)
rma_all <- rma(yi = d$log_effect, sei = d$SE, method = "REML", test = "knha")
egger <- regtest(rma_all, model = "lm")

common_forest_args <- list(
  rightcols = c("effect", "ci"),
  rightlabs = c("Estimate", "95% CI"),
  smlab = "Relative effect",
  text.random = "Random-effects model (HK)",
  col.diamond = "#2E75B6",
  col.diamond.lines = "#1F4E79",
  col.square = "#2E75B6",
  digits = 2,
  overall = TRUE,
  label.left = "Lower relative effect",
  label.right = "Higher relative effect"
)

draw_setting <- function() {
  do.call(forest, c(list(
    x = sg_setting,
    leftcols = c("studlab", "country", "n_amd"),
    leftlabs = c("Study", "Country", "AMD (N)"),
    sortvar = d$relative_effect,
    hetstat = TRUE, print.I2 = TRUE, test.subgroup = TRUE
  ), common_forest_args))
}

draw_region <- function() {
  do.call(forest, c(list(
    x = sg_region,
    leftcols = c("studlab", "country", "n_amd"),
    leftlabs = c("Study", "Country", "AMD (N)"),
    sortvar = d$relative_effect,
    hetstat = TRUE, print.I2 = TRUE, test.subgroup = TRUE
  ), common_forest_args))
}

draw_funnel <- function() {
  funnel(
    rma_all, main = "", xlab = "Log relative effect", yaxis = "sei",
    level = c(0.90, 0.95, 0.99), shade = c("white", "gray90", "gray75"),
    refline = 0, legend = TRUE
  )
  mtext(
    sprintf("Egger test: intercept = %.3f, P = %.3f", egger$est, egger$pval),
    side = 1, line = 3.5, cex = 0.75, col = "gray35"
  )
}

draw_no_keenan <- function() {
  do.call(forest, c(list(
    x = meta_no_keenan,
    leftcols = c("studlab", "country", "n_amd", "n_ctrl"),
    leftlabs = c("Study", "Country", "AMD (N)", "Control (N)"),
    sortvar = no_keenan$relative_effect,
    hetstat = TRUE, print.I2 = TRUE, print.tau2 = TRUE
  ), common_forest_args))
}

specs <- list(
  list(id = "FigS4_Subgroup_Setting", width = 10.5, height = 5.5, draw = draw_setting),
  list(id = "FigS5_Subgroup_Region", width = 8.5, height = 6.0, draw = draw_region),
  list(id = "FigS6_Funnel", width = 8.0, height = 6.2, draw = draw_funnel),
  list(id = "FigS9_Excluding_Keenan", width = 9.5, height = 3.9, draw = draw_no_keenan)
)

qa <- list()
for (spec in specs) {
  pdf_file <- file.path(FIG_DIR, paste0(spec$id, ".pdf"))
  preview_file <- file.path(PREVIEW_DIR, paste0(spec$id, ".png"))
  grDevices::cairo_pdf(pdf_file, width = spec$width, height = spec$height, family = "Arial", bg = "white")
  spec$draw()
  dev.off()
  suppressWarnings(pdf_convert(
    pdf_file, format = "png", dpi = 180, pages = 1,
    filenames = preview_file, verbose = FALSE
  ))
  page <- pdf_pagesize(pdf_file)[1, ]
  qa[[length(qa) + 1L]] <- data.frame(
    figure = spec$id,
    width_mm = round(as.numeric(page[["width"]]) * 25.4 / 72, 1),
    height_mm = round(as.numeric(page[["height"]]) * 25.4 / 72, 1),
    pdf_bytes = file.info(pdf_file)$size,
    preview_bytes = file.info(preview_file)$size,
    stringsAsFactors = FALSE
  )
}

write.table(
  d, file.path(SOURCE_DIR, "Meta_supplementary_relative_effect_source_data.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
write.table(
  do.call(rbind, qa), file.path(CODE_DIR, "FigS4_S5_S6_S9_export_QA.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
capture.output(sessionInfo(), file = file.path(CODE_DIR, "FigS4_S5_S6_S9_R_sessionInfo.txt"))
cat("Figs. S4, S5, S6, and S9 exported with relative-effect terminology.\n")
