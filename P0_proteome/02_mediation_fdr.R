# Quick mediation FDR integration
library(data.table)

alpha <- fread("../04_protein_mr/C6_rs429358_effects.csv")
beta <- fread("../04_protein_mr/C6_ukbppp_mr_all.csv")
beta <- beta[method == "Inverse variance weighted", ]

# Rename for merge
setnames(alpha, "Gene", "exposure")

# Merge
med <- merge(alpha[, .(exposure, Category, BETA, SE)],
             beta[, .(exposure, outcome, b, se, pval)], by = "exposure", all.x = TRUE)

med <- med[!is.na(b), ]

# Mediation = alpha * beta
med[, med_eff := BETA * b]
med[, med_se := sqrt(BETA^2 * se^2 + b^2 * SE^2)]
med[, med_z := med_eff / med_se]
med[, med_p := 2 * pnorm(-abs(med_z))]

# FDR per outcome
med[, med_fdr := p.adjust(med_p, method = "BH"), by = outcome]

cat("\n=== FDR-Significant Mediation (FDR<0.05) ===\n")
sig <- med[med_fdr < 0.05, ]
if (nrow(sig) > 0) {
  print(sig[, .(exposure, outcome, Category, BETA, b, med_eff, med_p, med_fdr)])
} else {
  cat("  No paths pass FDR<0.05 threshold\n")
}

cat("\n=== Top 5 Nominal Mediation Paths ===\n")
print(med[order(med_p)][1:min(5, nrow(med)),
      .(exposure, outcome, Category, BETA, b, med_eff, med_p, med_fdr)])

cat(sprintf("\nTotal paths: %d, Nominal P<0.05: %d, FDR<0.05: %d\n",
  nrow(med), sum(med$med_p < 0.05), sum(med$med_fdr < 0.05, na.rm = TRUE)))

fwrite(med, "./results/mediation_with_fdr.csv")
cat("[Done] Saved to ./results/mediation_with_fdr.csv\n")
