#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(meta)
  library(ragg)
  library(pdftools)
})

script_arg <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
code_dir <- normalizePath(dirname(script_arg), winslash = "/", mustWork = TRUE)
submission <- normalizePath(file.path(code_dir, ".."), winslash = "/", mustWork = TRUE)
project <- normalizePath(file.path(submission, ".."), winslash = "/", mustWork = TRUE)
figure_dir <- file.path(submission, "supplementary_figures")
preview_dir <- file.path(submission, "previews")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(preview_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(project, "figures"), recursive = TRUE, showWarnings = FALSE)

figure_pdf <- file.path(figure_dir, "FigS8_LeaveOneOut.pdf")
preview_png <- file.path(preview_dir, "FigS8_LeaveOneOut.png")
source_tsv <- file.path(code_dir, "FigS8_LeaveOneOut_source_data.tsv")
qa_tsv <- file.path(code_dir, "FigS8_LeaveOneOut_export_QA.tsv")
session_file <- file.path(code_dir, "FigS8_LeaveOneOut_R_sessionInfo.txt")

d <- data.frame(
  study = c(
    "Choi 2020", "Hwang PH 2021", "Keenan 2014", "Klaver 1999",
    "Lee CS 2019", "Shang 2021", "Tsai 2015", "Wen 2021",
    "Son 2025", "Lishinsky-Fischer 2025"
  ),
  hr = c(1.48, 1.87, 0.86, 1.50, 1.20, 1.08, 1.44, 1.23, 1.09, 1.02),
  ci_lower = c(1.25, 1.13, 0.67, 0.60, 1.02, 0.80, 1.26, 1.04, 1.05, 0.89),
  ci_upper = c(1.74, 3.09, 1.08, 3.50, 1.40, 1.44, 1.64, 1.46, 1.13, 1.16),
  stringsAsFactors = FALSE
)
d$logHR <- log(d$hr)
d$SE <- (log(d$ci_upper) - log(d$ci_lower)) / (2 * 1.96)

fit_meta <- function(dat) {
  metagen(
    TE = logHR,
    seTE = SE,
    studlab = study,
    data = dat,
    sm = "RR",
    method.tau = "REML",
    method.random.ci = "HK",
    common = FALSE,
    random = TRUE,
    prediction = TRUE
  )
}

meta_all <- fit_meta(d)
loo <- metainf(meta_all, pooled = "random")

loo_source <- do.call(rbind, lapply(seq_len(nrow(d)), function(i) {
  fit <- fit_meta(d[-i, , drop = FALSE])
  data.frame(
    omitted_study = d$study[i],
    relative_effect = exp(fit$TE.random),
    ci_lower = exp(fit$lower.random),
    ci_upper = exp(fit$upper.random),
    p_value = fit$pval.random,
    I2_percent = fit$I2,
    tau2 = fit$tau2,
    stringsAsFactors = FALSE
  )
}))
write.table(loo_source, source_tsv, sep = "\t", row.names = FALSE, quote = FALSE)

draw_figure <- function() {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(y = grid::unit(0.5, "npc") + grid::unit(4, "mm")))
  forest(
    loo,
    leftcols = "studlab",
    leftlabs = "Study",
    rightcols = c("effect", "ci"),
    rightlabs = c("Estimate", "95% CI"),
    smlab = "Relative effect",
    text.random = "Random effects model (HK)",
    col.diamond = "#2E75B6",
    col.diamond.lines = "black",
    col.bg = "#BDBDBD",
    col.border = "#BDBDBD",
    col.inside = "black",
    hetstat = FALSE,
    print.I2 = FALSE,
    prediction = FALSE,
    new = FALSE,
    fontsize = 9.5,
    squaresize = 0.75,
    label.left = "Lower relative effect",
    label.right = "Higher relative effect"
  )
  grid::popViewport()
}

# Tight dimensions retain a narrow clipping-safe perimeter around the table.
width_mm <- 148
height_mm <- 90
w <- width_mm / 25.4
h <- height_mm / 25.4

grDevices::cairo_pdf(figure_pdf, width = w, height = h, family = "Arial", bg = "white")
draw_figure()
dev.off()

ragg::agg_png(preview_png, width = w, height = h, units = "in", res = 200,
              background = "white")
draw_figure()
dev.off()

file.copy(figure_pdf, file.path(project, "figures", "FigS8_LeaveOneOut.pdf"), overwrite = TRUE)

page <- pdftools::pdf_pagesize(figure_pdf)[1, ]
fonts <- pdftools::pdf_fonts(figure_pdf)
text_data <- pdftools::pdf_data(figure_pdf)[[1]]
qa <- data.frame(
  check = c("width_mm", "height_mm", "text_objects", "embedded_fonts", "pdf_bytes", "preview_bytes"),
  value = c(
    round(as.numeric(page[["width"]]) * 25.4 / 72, 1),
    round(as.numeric(page[["height"]]) * 25.4 / 72, 1),
    nrow(text_data),
    paste(unique(fonts$name), collapse = "; "),
    file.info(figure_pdf)$size,
    file.info(preview_png)$size
  ),
  stringsAsFactors = FALSE
)
write.table(qa, qa_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
capture.output(sessionInfo(), file = session_file)
cat("Fig. S8 leave-one-out forest exported with tight page bounds.\n")
