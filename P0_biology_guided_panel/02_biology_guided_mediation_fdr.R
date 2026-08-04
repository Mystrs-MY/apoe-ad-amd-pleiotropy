# FDR sensitivity for two-step mediation in the biology-guided panel.
library(data.table)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this file with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = TRUE)
module_dir <- dirname(script_path)
project_root <- normalizePath(file.path(module_dir, ".."), winslash = "/", mustWork = TRUE)
result_dir <- file.path(module_dir, "results")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

alpha <- fread(file.path(project_root, "04_protein_mr", "C6_rs429358_effects.csv"))[Gene != "APOE"]
beta <- fread(file.path(project_root, "04_protein_mr", "C6_ukbppp_mr_all.csv"))
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

output <- file.path(result_dir, "biology_guided_mediation_with_fdr.csv")
fwrite(med, output)
cat("[Done] Saved to ", output, "\n", sep = "")
