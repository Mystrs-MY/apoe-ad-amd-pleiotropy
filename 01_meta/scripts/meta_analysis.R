###############################################################################
# REVISED Meta-Analysis: AMD → AD (All 10 Studies in Primary)
# PROSPERO: CRD420261339269
# Addressing peer review concerns re: Keenan 2014, heterogeneity, NOS, GRADE
###############################################################################

library(meta)
library(metafor)

# ============================
# 1. DATA (10 studies, revised NOS)
# ============================
d <- data.frame(
  study    = c("Choi 2020", "Hwang PH 2021", "Keenan 2014", "Klaver 1999",
               "Lee CS 2019", "Shang 2021", "Tsai 2015", "Wen 2021",
               "Son 2025", "Lishinsky-Fischer 2025"),
  country  = c("South Korea","United States","United Kingdom","Netherlands",
               "United States","United Kingdom","China (Taiwan)","China (Taiwan)",
               "South Korea","Multi-national"),
  hr       = c(1.48, 1.87, 0.86, 1.50, 1.20, 1.08, 1.44, 1.23, 1.09, 1.02),
  ci_lower = c(1.25, 1.13, 0.67, 0.60, 1.02, 0.80, 1.26, 1.04, 1.05, 0.89),
  ci_upper = c(1.74, 3.09, 1.08, 3.50, 1.40, 1.44, 1.64, 1.46, 1.13, 1.16),
  n_amd    = c(2213, 668, 65894, 113, 1036, 3671, 4993, 10578, 21384, 56035),
  n_ctrl   = c(306127, 2707, 7700963, 811, 2284, 87602, 24965, 10578, 1662319, 56035),
  region   = c("Asia","North America","Europe","Europe","North America",
               "Europe","Asia","Asia","Asia","Multi-national"),
  design   = c("Retro.","Prosp.","Record link.","Prosp.","Prosp.",
               "Prosp.","Retro.","Retro.","Retro.","Retro."),
  setting  = c("Community","Community","Hospital/treated","Community","Community",
               "Community","Community","Community","Community","Community"),
  # Revised NOS (peer-review corrections):
  # Klaver 1999: O2=0 (2.1yr <5yr) -> NOS=8 (S4+C2+O2=4+2+2=8)
  # Shang 2021: S1=0 (healthy volunteer), S3=0 (self-report) -> NOS=6
  # Tsai 2015: C2=0 (missing smoking, education, APOE) -> NOS=7
  nos      = c(9, 9, 5, 8, 9, 6, 7, 8, 8, 8),
  stringsAsFactors = FALSE
)
d$logHR <- log(d$hr)
d$SE    <- (log(d$ci_upper) - log(d$ci_lower)) / (2 * 1.96)
d$quality <- ifelse(d$nos >= 7, "NOS >= 7", "NOS < 7")

cat(strrep("=",80),"\n")
cat("REVISED META-ANALYSIS: AMD and Risk of AD (10 studies)\n")
cat("PROSPERO: CRD420261339269\n")
cat("Changes: Keenan 2014 included in primary; NOS revised; Keenan exclusion as sensitivity\n")
cat(strrep("=",80),"\n\n")

# ============================
# 2. PRIMARY (10 studies)
# ============================
meta10 <- metagen(TE = logHR, seTE = SE, studlab = study, data = d,
                  sm = "HR", method.tau = "REML", method.random.ci = "HK",
                  common = FALSE, random = TRUE, prediction = TRUE)

cat("=== PRIMARY (k=10) ===\n")
cat(sprintf("HR=%.2f (%.2f-%.2f) P=%.4f I2=%.1f%% Tau2=%.4f\n",
            exp(meta10$TE.random), exp(meta10$lower.random),
            exp(meta10$upper.random), meta10$pval.random,
            meta10$I2*100, meta10$tau2))
cat(sprintf("Prediction interval: %.2f-%.2f\n\n",
            exp(meta10$lower.predict), exp(meta10$upper.predict)))

# Forest
pdf("Fig1_Primary_Revised.pdf", width=14, height=8)
forest(meta10,
       leftcols=c("studlab","country","n_amd","n_ctrl","nos"),
       leftlabs=c("Study","Country","AMD (N)","Control (N)","NOS"),
       rightcols=c("effect","ci"), rightlabs=c("HR","95% CI"),
       smlab="Hazard Ratio (IV, Random, 95% CI)",
       sortvar=d$hr, col.diamond="#2E75B6", col.diamond.lines="#1F4E79",
       col.predict="#C00000", col.square="#2E75B6",
       digits=2, overall=TRUE, hetstat=TRUE,
       print.I2=TRUE, print.tau2=TRUE, print.Q=TRUE, print.pval.Q=TRUE,
       calcwidth.prediction=TRUE, fs.lr=9)
dev.off()
cat("Saved: Fig1_Primary_Revised.pdf\n")

# ============================
# 3. SUBGROUPS
# ============================
helper_sg <- function(var, title) {
  sg <- update(meta10, subgroup = d[[var]], title = title)
  cat(sprintf("\n--- Subgroup: %s ---\n", var))
  cat(sprintf("Between-group P = %.4f\n", sg$pval.Q.b.random))
  for (i in seq_along(sg$TE.random.w)) {
    cat(sprintf("  %s: k=%d HR=%.2f (%.2f-%.2f) I2=%.1f%%\n",
                sg$subgroup.levels[i], sg$k.w[i],
                exp(sg$TE.random.w[i]), exp(sg$lower.random.w[i]),
                exp(sg$upper.random.w[i]), sg$I2.w[i]*100))
  }
  sg
}

sg_region  <- helper_sg("region",  "AMD-AD by Region")
sg_design  <- helper_sg("design",  "AMD-AD by Design")
sg_setting <- helper_sg("setting", "AMD-AD by Treatment Setting")
sg_quality <- helper_sg("quality", "AMD-AD by NOS Quality")

# Forest: Setting (key new subgroup)
pdf("Fig2_Setting_Revised.pdf", width=14, height=7.5)
forest(sg_setting,
       leftcols=c("studlab","country","n_amd"),
       leftlabs=c("Study","Country","AMD (N)"),
       rightcols=c("effect","ci"), rightlabs=c("HR","95% CI"),
       col.diamond="#2E75B6", sortvar=d$hr,
       overall=TRUE, hetstat=TRUE, print.I2=TRUE, test.subgroup=TRUE)
dev.off()
cat("Saved: Fig2_Setting_Revised.pdf\n")

pdf("Fig3_Region_Revised.pdf", width=14, height=8)
forest(sg_region, leftcols=c("studlab","country","n_amd"),
       leftlabs=c("Study","Country","AMD (N)"),
       rightcols=c("effect","ci"), rightlabs=c("HR","95% CI"),
       col.diamond="#2E75B6", sortvar=d$hr,
       overall=TRUE, hetstat=TRUE, print.I2=TRUE, test.subgroup=TRUE)
dev.off()
cat("Saved: Fig3_Region_Revised.pdf\n")

# ============================
# 4. META-REGRESSION (exploratory)
# ============================
cat("\n=== EXPLORATORY META-REGRESSION ===\n")
rma_base <- rma(yi=logHR, sei=SE, data=d, method="REML", test="knha")

for (mod in c("region","design","setting","nos")) {
  mr <- rma(yi=logHR, sei=SE, mods=as.formula(paste("~", mod)),
            data=d, method="REML")
  R2 <- (rma_base$tau2 - mr$tau2) / rma_base$tau2 * 100
  cat(sprintf("Moderator: %-12s  QM=%.2f  P=%.4f  R2=%.1f%%\n",
              mod, mr$QM, mr$QMp, R2))
}

# Follow-up time (where available)
d$fupy <- c(8.0, 4.1, NA, 2.1, 8.0, 11.0, 4.4, 5.7, 9.7, 5.0)
if (sum(!is.na(d$fupy)) >= 6) {
  mr_fu <- rma(yi=logHR, sei=SE, mods=~fupy, data=d, method="REML")
  R2_fu <- (rma_base$tau2 - mr_fu$tau2) / rma_base$tau2 * 100
  cat(sprintf("Moderator: follow-up    QM=%.2f  P=%.4f  R2=%.1f%%\n",
              mr_fu$QM, mr_fu$QMp, R2_fu))
  # Dose-response test: AMD duration (Lee 2019 internal)
  cat("\nLee 2019 internal dose-response: established(>5yr) HR 1.50 vs recent(<=5yr) HR 1.20\n")
  cat("Son 2025 internal dose-response: AMD+VD HR 1.27 vs AMD-VD HR 1.09\n")
}

# ============================
# 5. SENSITIVITY
# ============================
cat("\n=== SENSITIVITY ANALYSES ===\n")
cat(sprintf("%-40s %2s  %-22s  %6s  %6s\n","Analysis","k","HR (95% CI)","I^2(%)","P"))
cat(strrep("-",90),"\n")

S <- function(lab, m) {
  cat(sprintf("%-40s %2d  %.2f (%.2f-%.2f)  %5.1f  %6.4f\n",
              lab, m$k, exp(m$TE.random), exp(m$lower.random),
              exp(m$upper.random), m$I2*100, m$pval.random))
}

S("PRIMARY (all 10)", meta10)

# S1: Exclude Keenan
m1 <- metagen(TE=logHR, seTE=SE, studlab=study,
              data=subset(d, study!="Keenan 2014"),
              sm="HR", method.tau="REML", method.random.ci="HK",
              common=FALSE, random=TRUE)
S("- Keenan 2014 (anti-VEGF cohort)", m1)

# S2: Exclude Shang
m2 <- metagen(TE=logHR, seTE=SE, studlab=study,
              data=subset(d, study!="Shang 2021"),
              sm="HR", method.tau="REML", method.random.ci="HK",
              common=FALSE, random=TRUE)
S("- Shang 2021 (self-report)", m2)

# S3: Exclude Lishinsky-Fischer
m3 <- metagen(TE=logHR, seTE=SE, studlab=study,
              data=subset(d, study!="Lishinsky-Fischer 2025"),
              sm="HR", method.tau="REML", method.random.ci="HK",
              common=FALSE, random=TRUE)
S("- Lishinsky-Fischer 2025", m3)

# S4: High Quality (NOS>=7, excluding Keenan NOS=5)
m4 <- metagen(TE=logHR, seTE=SE, studlab=study,
              data=subset(d, nos>=7),
              sm="HR", method.tau="REML", method.random.ci="HK",
              common=FALSE, random=TRUE)
S("High Quality (NOS>=7)", m4)

# S5: Prospective only
m5 <- metagen(TE=logHR, seTE=SE, studlab=study,
              data=subset(d, design=="Prosp."),
              sm="HR", method.tau="REML", method.random.ci="HK",
              common=FALSE, random=TRUE)
S("Prospective only", m5)

# S6: Conservative (no Keenan, no Shang, no Lishinsky)
m6 <- metagen(TE=logHR, seTE=SE, studlab=study,
              data=subset(d, !study %in% c("Keenan 2014","Shang 2021","Lishinsky-Fischer 2025")),
              sm="HR", method.tau="REML", method.random.ci="HK",
              common=FALSE, random=TRUE)
S("Conservative (-Keenan,-Shang,-Lishinsky)", m6)

# S7: US network of Lishinsky-Fischer 2025 (HR 0.79 instead of 1.02)
d_us <- d
d_us$hr[d_us$study=="Lishinsky-Fischer 2025"] <- 0.79
d_us$ci_lower[d_us$study=="Lishinsky-Fischer 2025"] <- 0.71
d_us$ci_upper[d_us$study=="Lishinsky-Fischer 2025"] <- 0.89
d_us$logHR <- log(d_us$hr)
d_us$SE <- (log(d_us$ci_upper)-log(d_us$ci_lower))/(2*1.96)
m7 <- metagen(TE=logHR, seTE=SE, studlab=study, data=d_us,
              sm="HR", method.tau="REML", method.random.ci="HK",
              common=FALSE, random=TRUE)
S("Lishinsky-Fischer using US network (HR 0.79)", m7)

# Leave-one-out
cat("\n--- Leave-One-Out ---\n")
loo <- metainf(meta10, pooled="random")
pdf("Fig4_LOO_Revised.pdf", width=11, height=7.5)
forest(loo, leftcols="studlab", rightcols=c("effect","ci"),
       rightlabs=c("HR","95% CI"), col.diamond="#2E75B6",
       hetstat=TRUE, print.I2=TRUE)
dev.off()
cat("LOO range: ")
cat(sprintf("HR %.2f-%.2f (all P<0.05? check output)\n",
            exp(min(loo$lower.random, na.rm=TRUE)),
            exp(max(loo$upper.random, na.rm=TRUE))))

# ============================
# 6. PUBLICATION BIAS
# ============================
cat("\n=== PUBLICATION BIAS ===\n")
rma10 <- rma(yi=d$logHR, sei=d$SE, method="REML", test="knha")
pdf("Fig5_Funnel_Revised.pdf", width=10, height=8)
funnel(rma10, main="Funnel Plot: AMD and AD (10 studies)",
       xlab="Log Hazard Ratio", yaxis="sei",
       level=c(0.90,0.95,0.99), shade=c("white","gray90","gray75"),
       refline=0, legend=TRUE)
dev.off()

egger <- regtest(rma10, model="lm")
cat(sprintf("Egger: intercept=%.4f t=%.2f P=%.4f\n", egger$est, egger$zval, egger$pval))

pdf("FigS7_Baujat_Revised.pdf", width=9, height=7)
baujat(rma10, main="Baujat Plot")
dev.off()

# ============================
# 7. FINAL SUMMARY
# ============================
cat("\n",strrep("=",80),"\n")
cat("FINAL SUMMARY FOR REVISED MANUSCRIPT\n")
cat(strrep("=",80),"\n\n")
cat("Primary (10 studies, including Keenan 2014):\n")
cat(sprintf("  HR=%.2f (%.2f-%.2f) P=%.4f I2=%.1f%% Tau2=%.4f\n",
            exp(meta10$TE.random), exp(meta10$lower.random),
            exp(meta10$upper.random), meta10$pval.random,
            meta10$I2*100, meta10$tau2))
cat(sprintf("  Prediction interval: %.2f-%.2f\n",
            exp(meta10$lower.predict), exp(meta10$upper.predict)))
cat(sprintf("  Egger P=%.4f | GRADE: VERY LOW\n\n", egger$pval))

cat("Revised NOS: Klaver=7, Shang=6, Keenan=5, others unchanged.\n")
cat("Keenan exclusion sensitivity: HR=%.2f (%.2f-%.2f)\n",
    exp(m1$TE.random), exp(m1$lower.random), exp(m1$upper.random))

save.image("Meta_Revised_Workspace.RData")
cat("\nDone. Workspace saved.\n")
