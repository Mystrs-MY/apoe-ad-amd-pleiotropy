#!/usr/bin/env Rscript

rm(list = ls())

required <- c("data.table", "readxl", "meta")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R package(s): ", paste(missing, collapse = ", "))

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(meta)
})

script_arg <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)][1])
SCRIPT_DIR <- dirname(normalizePath(script_file, winslash = "/", mustWork = TRUE))
ROOT <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), winslash = "/", mustWork = TRUE)
WORKBOOK <- file.path(ROOT, "01_meta", "纳入的文章", "AMD_AD_Final_10Studies.xlsx")
TABLE_DIR <- file.path(ROOT, "tables_submission", "supplementary_tables")
RESULT_DIR <- file.path(ROOT, "01_meta", "results")
dir.create(RESULT_DIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(WORKBOOK)) stop("Missing extraction workbook: ", WORKBOOK)

raw <- as.data.table(read_excel(WORKBOOK, sheet = "2_Effect_Estimates"))
raw <- raw[!is.na(Study) & !is.na(`HR/RR`)][1:10]
setnames(raw, c("HR/RR", "CI Lower", "CI Upper", "Effect Metric"),
         c("estimate", "ci_lower", "ci_upper", "effect_metric"))
raw[, `:=`(
  study_label = paste(Study, as.integer(Year)),
  estimate = as.numeric(estimate),
  ci_lower = as.numeric(ci_lower),
  ci_upper = as.numeric(ci_upper),
  log_effect = log(as.numeric(estimate)),
  SE = (log(as.numeric(ci_upper)) - log(as.numeric(ci_lower))) / (2 * 1.96)
)]

if (nrow(raw) != 10L || any(!is.finite(raw$SE))) {
  stop("Effect-metric extraction did not yield 10 complete study estimates.")
}

hr_only <- raw[effect_metric %chin% c("HR", "aHR")]
if (nrow(hr_only) != 8L) stop("Expected eight HR/aHR studies, found ", nrow(hr_only))

fit <- metagen(
  TE = log_effect,
  seTE = SE,
  studlab = study_label,
  data = hr_only,
  sm = "RR",
  method.tau = "REML",
  method.random.ci = "HK",
  common = FALSE,
  random = TRUE,
  prediction = TRUE
)

result <- data.table(
  analysis = "HR_only_sensitivity",
  included_effect_metrics = "HR;aHR",
  excluded_studies = paste(raw[!effect_metric %chin% c("HR", "aHR")]$study_label, collapse = ";"),
  k = fit$k,
  pooled_relative_effect = exp(fit$TE.random),
  CI_lower = exp(fit$lower.random),
  CI_upper = exp(fit$upper.random),
  P_value = fit$pval.random,
  I2_percent = fit$I2 * 100,
  tau2 = fit$tau2,
  prediction_interval_lower = exp(fit$lower.predict),
  prediction_interval_upper = exp(fit$upper.predict),
  model = "REML random effects with Hartung-Knapp confidence interval"
)
fwrite(result, file.path(RESULT_DIR, "HR_only_sensitivity.tsv"), sep = "\t", na = "NA")

audit <- raw[, .(
  Study = study_label, Year, outcome = Outcome,
  effect_metric, estimate, ci_lower, ci_upper,
  included_in_primary = TRUE,
  included_in_HR_only = effect_metric %chin% c("HR", "aHR")
)]
fwrite(audit, file.path(RESULT_DIR, "effect_metric_audit.tsv"), sep = "\t", na = "NA")

table_s12_file <- file.path(TABLE_DIR, "TableS12_Study_Characteristics.csv")
table_s12 <- fread(table_s12_file)
metric_cols <- intersect(c("Effect metric", "Effect metric.x", "Effect metric.y"), names(table_s12))
if (length(metric_cols)) table_s12[, (metric_cols) := NULL]
metric_map <- raw[, .(Study = study_label, `Effect metric` = effect_metric)]
table_s12 <- merge(table_s12, metric_map, by = "Study", all.x = TRUE, sort = FALSE)
if (anyNA(table_s12$`Effect metric`)) stop("Failed to map an effect metric into Table S12.")
if ("AMD→AD (95% CI)" %in% names(table_s12)) {
  setnames(table_s12, "AMD→AD (95% CI)", "Relative effect (95% CI)")
}
preferred_order <- c(
  "Study", "Country", "Design", "AMD (N)", "Control (N)", "Follow-up (yr)",
  "AMD Ascertainment", "AD Ascertainment", "Effect metric",
  "Relative effect (95% CI)", "Weight (%)", "NOS"
)
setcolorder(table_s12, preferred_order)
fwrite(table_s12, table_s12_file, bom = TRUE)

table_s14_file <- file.path(TABLE_DIR, "TableS14_Sensitivity_Analyses.csv")
table_s14 <- fread(table_s14_file)
if ("HR (95% CI)" %in% names(table_s14)) {
  setnames(table_s14, "HR (95% CI)", "Relative effect (95% CI)")
}
hr_row <- data.table(
  Analysis = "HR/aHR-only (excluding Keenan rate ratio and Klaver RR)",
  k = as.character(result$k),
  `Relative effect (95% CI)` = sprintf(
    "%.2f (%.2f–%.2f)", result$pooled_relative_effect, result$CI_lower, result$CI_upper
  ),
  `I² (%)` = sprintf("%.1f", result$I2_percent),
  `P value` = sprintf("%.3f", result$P_value)
)
table_s14 <- rbind(table_s14[Analysis != hr_row$Analysis], hr_row, fill = TRUE)
fwrite(table_s14, table_s14_file, bom = TRUE)

table_s13_file <- file.path(TABLE_DIR, "TableS13_Subgroup_Analyses.csv")
table_s13 <- fread(table_s13_file)
if ("HR (95% CI)" %in% names(table_s13)) {
  setnames(table_s13, "HR (95% CI)", "Relative effect (95% CI)")
}
fwrite(table_s13, table_s13_file, bom = TRUE)

table_s15_file <- file.path(TABLE_DIR, "TableS15_GRADE_SoF.csv")
table_s15 <- fread(table_s15_file)
table_s15[Item == "Pooled HR (95% CI)", Item := "Pooled relative effect (95% CI)"]
table_s15[Item == "Large Effect", Assessment :=
            "Not applicable: pooled relative effect 1.20 is modest (<2.0)"]
fwrite(table_s15, table_s15_file, bom = TRUE)

writeLines(capture.output(sessionInfo()), file.path(RESULT_DIR, "effect_metric_sensitivity_sessionInfo.txt"))
print(result)
