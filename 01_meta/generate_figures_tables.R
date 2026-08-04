###############################################################################
# Generate ALL Figures & Tables for AMD→AD Meta-Analysis Manuscript
# Output: figures/ (PDF) and tables/ (CSV)
###############################################################################

library(meta)
library(metafor)

dir.create("figures", showWarnings = FALSE)
dir.create("tables", showWarnings = FALSE)

# ============================
# 0. DATA
# ============================
d <- data.frame(
  study    = c("Choi 2020","Hwang PH 2021","Keenan 2014","Klaver 1999",
               "Lee CS 2019","Shang 2021","Tsai 2015","Wen 2021",
               "Son 2025","Lishinsky-Fischer 2025"),
  country  = c("South Korea","United States","United Kingdom","Netherlands",
               "United States","United Kingdom","China (Taiwan)","China (Taiwan)",
               "South Korea","Multi-national"),
  hr       = c(1.48,1.87,0.86,1.50,1.20,1.08,1.44,1.23,1.09,1.02),
  ci_lower = c(1.25,1.13,0.67,0.60,1.02,0.80,1.26,1.04,1.05,0.89),
  ci_upper = c(1.74,3.09,1.08,3.50,1.40,1.44,1.64,1.46,1.13,1.16),
  n_amd    = c(2213,668,65894,113,1036,3671,4993,10578,21384,56035),
  n_ctrl   = c(306127,2707,7700963,811,2284,87602,24965,10578,1662319,56035),
  region   = c("Asia","North America","Europe","Europe","North America",
               "Europe","Asia","Asia","Asia","Multi-national"),
  design   = c("Retrospective","Prospective","Record linkage","Prospective",
               "Prospective","Prospective","Retrospective","Retrospective",
               "Retrospective","Retrospective"),
  setting  = c("Community","Community","Hospital/treated","Community","Community",
               "Community","Community","Community","Community","Community"),
  nos      = c(9,9,5,8,9,6,7,8,8,8),
  fupy     = c("8.0","4.1","NR (min >1)","2.1","~8.0","11.0","4.4","5.7","9.7","5–10"),
  amd_def  = c("ICD-10; ≥2 visits/hosp","ICD-9+Medicare+Hx","Hospital admission (ICD)",
               "Fundus photography; grading","ICD-9 from EMR","Self-report+hosp records",
               "ICD-9-CM; ≥2 visits","ICD-9-CM; ophthalmologist dx","ICD-10 H35.3; claims",
               "ICD-10 from EMR"),
  ad_def   = c("ICD-10 F00/G30+anti-dementia Rx","NINCDS+DSM; panel+MRI",
               "Hospital admission (ICD)","NINCDS; 3-step screen","NINCDS; consensus panel",
               "ICD; hosp+death+self-report","ICD-9-CM 331.0/290.x; ≥2 dx",
               "ICD-9-CM 331","ICD-10 F00/G30+anti-dementia Rx","ICD-10 G30+medication records"),
  stringsAsFactors = FALSE
)
d$logHR <- log(d$hr)
d$SE    <- (log(d$ci_upper) - log(d$ci_lower)) / (2 * 1.96)
d$quality <- ifelse(d$nos >= 7, "NOS ≥7 (moderate-high)", "NOS <7 (lower)")
d$weight <- 1 / d$SE^2
d$weight_pct <- round(d$weight / sum(d$weight) * 100, 1)
d$effect_str <- sprintf("%.2f (%.2f–%.2f)", d$hr, d$ci_lower, d$ci_upper)

# ============================
# 1. PRIMARY META
# ============================
meta_all <- metagen(TE=logHR, seTE=SE, studlab=study, data=d,
                    sm="HR", method.tau="REML", method.random.ci="HK",
                    common=FALSE, random=TRUE, prediction=TRUE)

meta_noKeenan <- metagen(TE=logHR, seTE=SE, studlab=study,
                          data=subset(d, study!="Keenan 2014"),
                          sm="HR", method.tau="REML", method.random.ci="HK",
                          common=FALSE, random=TRUE)

# PRISMA screening counts are retained in machine-readable form, while the
# bespoke Word-template layout workflow is intentionally outside this package.
cat("PRISMA layout code is not distributed; verified counts are retained in 01_meta/PRISMA/Article1_PRISMA_content.json\n")

# ============================
# FIGURE 2 — Primary Forest Plot
# ============================
cat("Generating Fig2_Forest_Primary.pdf...\n")
pdf("../figures/Fig2_Forest_Primary.pdf", width=10.3, height=4.8)
forest(meta_all,
       leftcols=c("studlab","country","n_amd","n_ctrl","nos"),
       leftlabs=c("Study","Country","AMD (N)","Control (N)","NOS"),
       rightcols=c("effect","ci"),
       rightlabs=c("HR","95% CI"),
       smlab="Hazard Ratio (IV, Random, 95% CI)",
       sortvar=d$hr,
       col.diamond="#2E75B6", col.diamond.lines="#1F4E79",
       col.predict="#C00000", col.square="#2E75B6",
       digits=2, overall=TRUE, hetstat=TRUE,
       print.I2=TRUE, print.tau2=TRUE, print.Q=TRUE, print.pval.Q=TRUE,
       plotwidth="8cm", spacing=1.3,
       fs.lr=9, ff.lr="bold",
       label.left="Favours AMD\nnot associated with AD",
       label.right="Favours AMD\nassociated with AD")
dev.off()
cat("  -> figures/Fig2_Forest_Primary.pdf\n")

# ============================
# FIGURE 3 — Subgroup: AMD Treatment Setting
# ============================
cat("Generating Fig3_Subgroup_Setting.pdf...\n")
sg_setting <- update(meta_all, subgroup=setting)
pdf("../figures/Fig3_Subgroup_Setting.pdf", width=8.1, height=4.5)
forest(sg_setting,
       leftcols=c("studlab","country","n_amd"),
       leftlabs=c("Study","Country","AMD (N)"),
       rightcols=c("effect","ci"), rightlabs=c("HR","95% CI"),
       col.diamond="#2E75B6", sortvar=d$hr,
       overall=TRUE, hetstat=TRUE, print.I2=TRUE, test.subgroup=TRUE)
dev.off()
cat("  -> figures/Fig3_Subgroup_Setting.pdf\n")

# ============================
# FIGURE 4 — Subgroup: Geographic Region
# ============================
cat("Generating Fig4_Subgroup_Region.pdf...\n")
sg_region <- update(meta_all, subgroup=region)
pdf("../figures/Fig4_Subgroup_Region.pdf", width=8.5, height=6)
forest(sg_region,
       leftcols=c("studlab","country","n_amd"),
       leftlabs=c("Study","Country","AMD (N)"),
       rightcols=c("effect","ci"), rightlabs=c("HR","95% CI"),
       col.diamond="#2E75B6", sortvar=d$hr,
       overall=TRUE, hetstat=TRUE, print.I2=TRUE, test.subgroup=TRUE)
dev.off()
cat("  -> figures/Fig4_Subgroup_Region.pdf\n")

# ============================
# FIGURE 5 — Funnel Plot
# ============================
cat("Generating Fig5_Funnel.pdf...\n")
rma_all <- rma(yi=d$logHR, sei=d$SE, method="REML", test="knha")
pdf("../figures/Fig5_Funnel.pdf", width=10, height=8)
funnel(rma_all,
       main="Funnel Plot: AMD and Alzheimer's Disease (10 studies)",
       xlab="Log Hazard Ratio", yaxis="sei",
       level=c(0.90,0.95,0.99), shade=c("white","gray90","gray75"),
       refline=0, legend=TRUE)
egger <- regtest(rma_all, model="lm")
mtext(sprintf("Egger test: intercept=%.3f, P=%.3f", egger$est, egger$pval),
      side=1, line=3.5, cex=0.7, col="gray40")
dev.off()
cat("  -> figures/Fig5_Funnel.pdf\n")

# ============================
# FIGURE S7 — Baujat Plot
# ============================
cat("Generating FigS7_Baujat.pdf...\n")
pdf("../figures/FigS7_Baujat.pdf", width=9, height=7)
baujat(rma_all,
       main="Baujat Plot: Influence on Overall Result",
       xlab="Contribution to Heterogeneity (Cochran's Q)",
       ylab="Influence on Pooled Result")
dev.off()
cat("  -> figures/FigS7_Baujat.pdf\n")

# ============================
# FIGURE S1 — Leave-One-Out
# ============================
cat("Generating FigS1_LeaveOneOut.pdf...\n")
loo <- metainf(meta_all, pooled="random")
pdf("../figures/FigS1_LeaveOneOut.pdf", width=8, height=6)
forest(loo,
       leftcols="studlab", rightcols=c("effect","ci"),
       rightlabs=c("HR","95% CI"),
       col.diamond="#2E75B6",
       hetstat=TRUE, print.I2=TRUE,
       label.left="Favours AMD not associated with AD",
       label.right="Favours AMD associated with AD")
dev.off()
cat("  -> figures/FigS1_LeaveOneOut.pdf\n")

# ============================
# FIGURE S2 — Sensitivity Excluding Keenan 2014
# ============================
cat("Generating FigS2_Excluding_Keenan.pdf...\n")
pdf("../figures/FigS2_Excluding_Keenan.pdf", width=9.5, height=3.9)
forest(meta_noKeenan,
       leftcols=c("studlab","country","n_amd","n_ctrl"),
       leftlabs=c("Study","Country","AMD (N)","Control (N)"),
       rightcols=c("effect","ci"), rightlabs=c("HR","95% CI"),
       smlab="Sensitivity (9 studies)",
       sortvar=subset(d, study!="Keenan 2014")$hr,
       col.diamond="#2E75B6", col.square="#2E75B6",
       digits=2, overall=TRUE, hetstat=TRUE,
       print.I2=TRUE, print.tau2=TRUE,
       colgap.left="3mm", spacing=1.3)
dev.off()
cat("  -> figures/FigS2_Excluding_Keenan.pdf\n")

# ============================
# TABLES
# ============================

# --- Table 1: Study Characteristics ---
cat("\nGenerating tables...\n")
t1 <- data.frame(
  Study = d$study,
  Country = d$country,
  Design = d$design,
  `AMD (N)` = d$n_amd,
  `Control (N)` = d$n_ctrl,
  `Follow-up (yr)` = d$fupy,
  `AMD Ascertainment` = d$amd_def,
  `AD Ascertainment` = d$ad_def,
  `AMD→AD (95% CI)` = d$effect_str,
  `Weight (%)` = d$weight_pct,
  NOS = d$nos,
  check.names = FALSE
)
write.csv(t1, "tables/Table1_Study_Characteristics.csv", row.names=FALSE)

# --- Table 2: Subgroup Analyses ---
subgroup_results <- list(
  c("Overall (primary)", 10, "1.20 (1.04–1.38)", "77.5", "0.018", "—"),
  c("AMD Treatment Setting", "", "", "", "", "0.007"),
  c("  Community-based", 9, "1.23 (1.09–1.40)", "77.2", "0.005", ""),
  c("  Hospital/treated", 1, "0.86 (0.68–1.09)", "—", "0.210", ""),
  c("Geographic Region", "", "", "", "", "0.043"),
  c("  Asia", 4, "1.29 (1.02–1.63)", "89.2", "0.036", ""),
  c("  North America", 2, "1.40 (0.10–20.58)", "63.2", "0.514", ""),
  c("  Europe", 3, "0.97 (0.61–1.56)", "17.2", "0.894", ""),
  c("  Multi-national", 1, "1.02 (0.89–1.16)", "—", "0.792", ""),
  c("Study Design", "", "", "", "", "0.029"),
  c("  Retrospective", 5, "1.23 (1.00–1.51)", "86.7", "0.048", ""),
  c("  Prospective", 4, "1.22 (0.96–1.55)", "18.5", "0.079", ""),
  c("  Record linkage", 1, "0.86 (0.68–1.09)", "—", "0.210", ""),
  c("NOS Quality", "", "", "", "", "0.030"),
  c("  NOS ≥7 (moderate-high)", 8, "1.25 (1.08–1.44)", "80.0", "0.008", ""),
  c("  NOS <7 (lower)", 2, "0.95 (0.23–3.97)", "28.1", "0.922", "")
)
t2 <- as.data.frame(do.call(rbind, subgroup_results))
colnames(t2) <- c("Subgroup","k","HR (95% CI)","I² (%)","P (within)","P (between)")
write.csv(t2, "tables/Table2_Subgroup_Analyses.csv", row.names=FALSE)

# --- Table 3: Sensitivity Analyses ---
sens_results <- list(
  c("Primary (all 10 studies)", 10, "1.20 (1.04–1.38)", "77.5", "0.018"),
  c("− Keenan 2014", 9, "1.23 (1.09–1.40)", "77.2", "0.005"),
  c("− Shang 2021", 9, "1.21 (1.03–1.42)", "79.9", "0.024"),
  c("− Lishinsky-Fischer 2025", 9, "1.23 (1.05–1.43)", "78.8", "0.016"),
  c("NOS ≥7 only", 8, "1.25 (1.08–1.44)", "80.0", "0.008"),
  c("Prospective only", 4, "1.22 (0.96–1.55)", "18.5", "0.079"),
  c("Conservative (−Keenan,−Shang,−Lishinsky)", 7, "1.29 (1.12–1.49)", "81.6", "0.004"),
  c("Lishinsky-Fischer with US network (post hoc)", 10, "1.17 (0.98–1.39)", "87.7", "0.083")
)
t3 <- as.data.frame(do.call(rbind, sens_results))
colnames(t3) <- c("Analysis","k","HR (95% CI)","I² (%)","P value")
write.csv(t3, "tables/Table3_Sensitivity_Analyses.csv", row.names=FALSE)

# --- Table 4: GRADE SoF ---
t4 <- data.frame(
  Item = c("Outcome","Studies (Participants)","Study Design","Risk of Bias",
           "Inconsistency","Indirectness","Imprecision","Publication Bias",
           "Dose-Response","Large Effect","Plausible Confounding",
           "Pooled HR (95% CI)","Starting Certainty","Final Certainty"),
  Assessment = c(
    "Incident Alzheimer's disease",
    "10 cohort studies (10,658,212 participants; AMD: 166,574; non-AMD: 10,491,638)",
    "Cohort studies (4 prospective, 5 retrospective, 1 record linkage)",
    "Serious (−1): ICD-based diagnoses in most studies without standardized exams; missing key confounders in 2 studies; 1 study NOS=5",
    "Serious (−1): I²=77.5% (considerable); prediction interval 0.82–1.74 crosses null; effect sizes 0.86–1.87; no single factor fully explains heterogeneity",
    "Not serious: PICO well-matched across all 10 studies",
    "Not serious: 95% CI 1.04–1.38 excludes null; >10M participants; OIS met",
    "Undetected: Egger P=0.215; symmetric funnel plot; comprehensive search strategy",
    "Not upgraded: No formal within-study dose-response test; between-stratum CIs overlap",
    "Not applicable: HR 1.20 is modest (<2.0)",
    "Not applicable: direction uncertain",
    "1.20 (1.04–1.38), P=0.018, I²=77.5%",
    "LOW (observational studies)",
    "VERY LOW ⓉⓋⓋⓋ"
  ),
  check.names=FALSE
)
write.csv(t4, "tables/Table4_GRADE_SoF.csv", row.names=FALSE)

# --- Supplementary Tables ---
# S6: Excluded Studies
# The updated full-text exclusion table is manually curated from the 88 full-text
# candidates and contains one row per excluded report. Keep it intact when this
# batch script is rerun, rather than replacing it with an outdated summary table.
excluded_table <- "tables/TableS6_Excluded_Studies.csv"
if (file.exists(excluded_table)) {
  tS6 <- read.csv(excluded_table, check.names=FALSE)
  cat(sprintf("  -> %s retained (%d excluded reports)\n", excluded_table, nrow(tS6)))
} else {
  warning("Updated TableS6_Excluded_Studies.csv not found; no deprecated summary table was written.")
}

# S7: GRADE Evidence Profile
tS7 <- data.frame(
  Domain = c("Starting certainty","Risk of Bias","Inconsistency","Indirectness",
             "Imprecision","Publication Bias","Large Effect","Dose-Response",
             "Plausible Confounding","Final certainty"),
  Judgment = c("LOW","Serious (−1)","Serious (−1)","Not serious","Not serious",
               "Undetected","Not applicable","Not upgraded","Not applicable","VERY LOW ⓉⓋⓋⓋ"),
  Rationale = c(
    "Observational studies",
    "ICD-based diagnoses without standardized exams; missing confounders (smoking, education, APOE) in 2 studies; 1 study NOS=5",
    "I²=77.5% (considerable); prediction interval 0.82–1.74 crosses null; no single factor fully explains heterogeneity",
    "PICO well-matched across all studies",
    "95% CI excludes null; >10M participants; OIS met",
    "Egger P=0.215; funnel plot symmetric; pre-registered protocol",
    "HR 1.20 is modest (<2.0)",
    "No formal within-study dose-response test; between-stratum CIs overlap (e.g., Tsai 2015: dry vs wet AMD)",
    "Direction of residual confounding uncertain",
    "LOW → downgraded 2 levels → VERY LOW"
  ),
  check.names=FALSE
)
write.csv(tS7, "tables/TableS7_GRADE_Evidence_Profile.csv", row.names=FALSE)

# ============================
# DONE
# ============================
cat("\n========== ALL FILES GENERATED ==========\n")
cat("Figures (figures/):\n")
for(f in c("Fig2_Forest_Primary.pdf",
           "Fig3_Subgroup_Setting.pdf","Fig4_Subgroup_Region.pdf",
           "Fig5_Funnel.pdf","FigS7_Baujat.pdf",
           "FigS1_LeaveOneOut.pdf","FigS2_Excluding_Keenan.pdf")) {
  cat(sprintf("  %s\n", f))
}
cat("\nTables (tables/):\n")
for(f in c("Table1_Study_Characteristics.csv","Table2_Subgroup_Analyses.csv",
           "Table3_Sensitivity_Analyses.csv","Table4_GRADE_SoF.csv",
           "TableS6_Excluded_Studies.csv","TableS7_GRADE_Evidence_Profile.csv")) {
  cat(sprintf("  %s\n", f))
}
cat("\nDone.\n")
