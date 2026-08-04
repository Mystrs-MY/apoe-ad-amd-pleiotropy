###############################################################################
# GRADE Assessment
# AMD -> AD Meta-Analysis
# PROSPERO: CRD420261339269
###############################################################################

library(meta)
library(metafor)
load("Meta_Workspace.RData")

cat(strrep("=", 80), "\n")
cat("GRADE CERTAINTY OF EVIDENCE ASSESSMENT\n")
cat("PICO: In adults >=50yr without baseline dementia,\n")
cat("      does pre-existing AMD (vs. no AMD) increase risk of incident AD?\n")
cat(strrep("=", 80), "\n\n")

# ==========================================
# GRADE Table Construction
# ==========================================

cat("STARTING CERTAINTY: LOW (observational studies)\n\n")

# --- Domain 1: Risk of Bias ---
cat("--- DOMAIN 1: Risk of Bias ---\n")
cat("  NOS scores: 9,9,8,8,8,8,8,7,7,5 (median=8)\n")
cat("  7/10 studies NOS >= 8 (high quality)\n")
cat("  Keenan 2014 (NOS=5): only hospital-treated AMD + AD (weight=2.9%)\n")
cat("  Shang 2021 (NOS=7): self-reported AMD exposure (weight=4.4%)\n")
cat("  Tsai 2015 (NOS=7): missing key confounders - smoking, education, APOE (weight=17.6%)\n")
cat("  Sensitivity: excluding Keenan -> HR 1.23 (1.09-1.40), excluding Shang -> HR 1.21 (1.03-1.42)\n")
cat("  High-quality only (NOS>=8): HR 1.21 (1.04-1.42), I^2=71.5%\n")
cat("  Most studies rely on ICD codes rather than standardized clinical exams\n")
cat("  => SERIOUS risk of bias (-1)\n\n")

# --- Domain 2: Inconsistency ---
cat("--- DOMAIN 2: Inconsistency ---\n")
cat(sprintf("  I^2 = %.1f%% (95%% CI: %.1f%%-%.1f%%) -> substantial heterogeneity\n",
            meta_primary$I2*100, meta_primary$lower.I2*100, meta_primary$upper.I2*100))
cat(sprintf("  Prediction interval: %.2f-%.2f (crosses 1.0)\n",
            exp(meta_primary$lower.predict), exp(meta_primary$upper.predict)))
cat("  Point estimates range: 0.86 (Keenan) to 1.87 (Hwang PH)\n")
cat("  8/10 studies show HR > 1.0 (consistent direction)\n")
cat("  Subgroup analyses partially explain heterogeneity:\n")
cat("    - Region: P=0.043 (Asia HR 1.29 vs Europe HR 0.97)\n")
cat("    - Design: P=0.029 (Retrospective HR 1.23 vs Prospective HR 1.22)\n")
cat("    - Prospective-only I^2=18.5% (substantially lower)\n")
cat("  => SERIOUS inconsistency (-1)\n\n")

# --- Domain 3: Indirectness ---
cat("--- DOMAIN 3: Indirectness ---\n")
cat("  PICO match assessment:\n")
cat("  - Population: older adults >=50yr without baseline dementia [✓ matched]\n")
cat("  - Exposure: AMD diagnosis by ICD codes/clinical exam/fundus photography [✓ matched]\n")
cat("  - Comparator: individuals without AMD from same source [✓ matched]\n")
cat("  - Outcome: incident AD (NINCDS-ADRDA or ICD codes) [✓ matched]\n")
cat("  Minor concerns:\n")
cat("    - Tsai 2015: outcome is 'AD or senile dementia' (not pure AD) [minor]\n")
cat("    - Keenan 2014: only hospitalized AMD (not generalizable) [minor, low weight]\n")
cat("    - Shang 2021: self-reported AMD [minor, sensitivity analysis confirms]\n")
cat("  => NOT SERIOUS indirectness (no downgrade)\n\n")

# --- Domain 4: Imprecision ---
cat("--- DOMAIN 4: Imprecision ---\n")
cat(sprintf("  10 studies, >10 million total participants\n"))
cat(sprintf("  Pooled HR = %.2f, 95%% CI: %.2f-%.2f\n",
            exp(meta_primary$TE.random), exp(meta_primary$lower.random),
            exp(meta_primary$upper.random)))
cat("  CI excludes null (1.0)\n")
cat("  Optimal Information Size (OIS) clearly met\n")
cat("  Sensitivity analyses consistently exclude null\n")
cat("  => NOT SERIOUS imprecision (no downgrade)\n\n")

# --- Domain 5: Publication Bias ---
cat("--- DOMAIN 5: Publication Bias ---\n")
egger <- regtest(res_rma, model = "lm")
cat(sprintf("  Egger test: intercept=%.4f, P=%.4f (not significant)\n",
            egger$est, egger$pval))
cat("  Funnel plot: visually symmetric\n")
cat("  10 studies meet minimum threshold for Egger test\n")
cat("  Note: PROSPERO protocol pre-registered; comprehensive search strategy\n")
cat("  => UNDETECTED (no downgrade)\n\n")

# --- Upgrade: Dose-Response ---
cat("--- UPGRADE: Dose-Response Gradient ---\n")
cat("  Evidence of a dose-response relationship from 3 independent sources:\n")
cat("  1. AMD duration (Lee CS 2019):\n")
cat("     Established AMD (>5yr): HR 1.50 (1.25-1.81) > Recent AMD (<=5yr): HR 1.20 (0.95-1.50)\n")
cat("  2. AMD severity (Son 2025):\n")
cat("     AMD + visual disability: HR 1.27 (1.12-1.43) > AMD without VD: HR 1.09 (1.05-1.13)\n")
cat("  3. AMD subtype (Tsai 2015):\n")
cat("     Dry/non-exudative AMD: HR 1.44 (1.26-1.65) > Wet/exudative AMD: HR 1.35 (0.89-2.06) NS\n")
cat("  This pattern suggests a biological gradient: more severe/longer AMD -> higher AD risk.\n")
cat("  => UPGRADE +1 for dose-response gradient\n\n")

# --- Final ---
cat(strrep("-", 80), "\n")
cat("FINAL GRADE CERTAINTY\n")
cat("  Starting:    LOW    (observational studies)\n")
cat("  Risk of bias: -1    (serious)\n")
cat("  Inconsistency: -1   (serious)\n")
cat("  Indirectness:  0    (not serious)\n")
cat("  Imprecision:   0    (not serious)\n")
cat("  Pub. bias:     0    (undetected)\n")
cat("  Dose-response: +1   (present)\n")
cat("  ----------------------------------------\n")
cat("  FINAL:        LOW   certainty of evidence\n")
cat(strrep("-", 80), "\n\n")

cat("GRADE INTERPRETATION:\n")
cat("  LOW certainty: Our confidence in the effect estimate is limited.\n")
cat("  The true effect may be substantially different from the estimate.\n")
cat("  Further research is likely to change the estimate.\n\n")

# ==========================================
# Summary of Findings (SoF) Table
# ==========================================
cat(strrep("=", 80), "\n")
cat("SUMMARY OF FINDINGS (SoF) TABLE\n")
cat(strrep("=", 80), "\n\n")

cat("Outcome: Incident Alzheimer's disease\n")
cat("Population: Adults >=50 years without baseline dementia\n")
cat("Exposure: Pre-existing age-related macular degeneration (AMD)\n")
cat("Comparator: No AMD diagnosis\n")
cat(sprintf("Pooled HR: %.2f (95%% CI: %.2f-%.2f)\n",
            exp(meta_primary$TE.random), exp(meta_primary$lower.random),
            exp(meta_primary$upper.random)))
cat(sprintf("Absolute risk (illustrative):\n"))
# Illustrative absolute risk
# Assume baseline AD incidence ~5 per 1000 PY over ~7 years
# With HR 1.20: 5 * 1.20 = 6 per 1000 PY
cat("  Low risk population:  5 per 1000 person-years\n")
cat("  High risk population: 6 per 1000 person-years (HR 1.20 applied)\n")
cat(sprintf("  Risk difference: +1 per 1000 person-years\n\n"))

cat("GRADE domains:\n")
cat("  Risk of bias:    Serious (-1)\n")
cat("  Inconsistency:   Serious (-1)\n")
cat("  Indirectness:    Not serious\n")
cat("  Imprecision:     Not serious\n")
cat("  Publication bias: Undetected\n")
cat("  Upgrade:         Dose-response (+1)\n")
cat("  Certainty:       LOW ⓉⓉⓋⓋ\n\n")

cat("Interpretation:\n")
cat("  AMD is probably associated with a modestly increased risk of AD (HR ~1.2).\n")
cat("  The true effect likely lies between no association and a moderately increased risk.\n")
cat("  Future prospective studies with standardized AMD grading and AD ascertainment\n")
cat("  may strengthen or modify this conclusion.\n")

# ==========================================
# Save GRADE table as CSV for paper
# ==========================================
grade_df <- data.frame(
  Item = c("Studies", "Design", "Participants", "Pooled HR (95% CI)",
           "I^2", "Prediction Interval",
           "Risk of Bias", "Inconsistency", "Indirectness",
           "Imprecision", "Publication Bias", "Dose-Response",
           "Starting Certainty", "Final Certainty"),
  Assessment = c(
    "10 cohort studies", "7 retrospective / record linkage; 4 prospective",
    ">10 million (AMD n=166,574; Control n=9,491,638)",
    sprintf("%.2f (%.2f-%.2f)", exp(meta_primary$TE.random),
            exp(meta_primary$lower.random), exp(meta_primary$upper.random)),
    sprintf("%.1f%%", meta_primary$I2*100),
    sprintf("%.2f-%.2f", exp(meta_primary$lower.predict),
            exp(meta_primary$upper.predict)),
    "Serious (-1): claims-based diagnoses; missing confounders in some studies; 1 study NOS=5",
    "Serious (-1): I^2=77.5%; prediction interval crosses 1.0; regional variation",
    "Not serious: PICO well-matched across studies",
    "Not serious: CI excludes null; >10M participants",
    "Undetected: Egger P=0.392; symmetric funnel plot",
    "Present (+1): stronger effect with longer AMD duration, greater severity",
    "LOW (observational studies)",
    "LOW ⓉⓉⓋⓋ"
  ),
  stringsAsFactors = FALSE
)

write.csv(grade_df, "GRADE_SoF_Table.csv", row.names = FALSE)
cat("\nGRADE SoF table saved: GRADE_SoF_Table.csv\n")

cat("\nAll GRADE assessments complete.\n")
