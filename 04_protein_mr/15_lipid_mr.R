# ============================================================
# I5_c3_lipid_rerun.R — Simplified lipid MR re-run
# 5 lipids (GLGC 2021) → AD / 3 AMD
# Uses direct PLINK clump (bypasses old mr_helpers.R)
# ============================================================

library(data.table)
library(TwoSampleMR)
library(dplyr)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = TRUE)
SCRIPT_DIR <- dirname(script_path)
PROJECT_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), winslash = "/", mustWork = TRUE)
RESOURCE_ROOT <- Sys.getenv("A1_RESOURCE_ROOT", unset = file.path(PROJECT_ROOT, "data", "external"))
PLINK_BIN <- Sys.getenv("PLINK_BIN", unset = plinkbinr::get_plink_exe())
LD_REF <- Sys.getenv("A1_LD_PREFIX", unset = file.path(RESOURCE_ROOT, "EUR", "EUR"))
TMP_DIR   <- tempdir()
OUT_DIR <- file.path(SCRIPT_DIR, "results")
dir.create(OUT_DIR, showWarnings = FALSE)

LIPID_DIR <- file.path(RESOURCE_ROOT, "lipid_glgc2021")
GWAS_DIR <- file.path(RESOURCE_ROOT, "GWAS")

# ---- Lipids ----
lipid_files <- list(
  LDL_C   = file.path(LIPID_DIR, "LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz"),
  HDL_C   = file.path(LIPID_DIR, "HDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz"),
  logTG   = file.path(LIPID_DIR, "logTG_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz"),
  TC      = file.path(LIPID_DIR, "TC_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz"),
  nonHDL_C= file.path(LIPID_DIR, "nonHDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz")
)

# ---- Outcomes ----
outcome_files <- list(
  AD      = file.path(GWAS_DIR, "AD_Wightman_cleaned_hg19.tsv.gz"),
  Dry_AMD = file.path(GWAS_DIR, "AMD_Dry_R12_cleaned_hg19.tsv.gz"),
  Wet_AMD = file.path(GWAS_DIR, "AMD_Wet_R12_cleaned_hg19.tsv.gz"),
  Any_AMD = file.path(GWAS_DIR, "AMD_H7_R12_cleaned_hg19.tsv.gz")
)

outcomes <- list()
for (nm in names(outcome_files)) {
  dat <- fread(outcome_files[[nm]], select = c("SNP","CHR","BP","A1","A2","FREQ","BETA","SE","P","N"))
  setnames(dat, c("CHR","BP","A1","A2","FREQ","BETA","SE","P","N"),
           c("chr.outcome","pos.outcome","effect_allele.outcome","other_allele.outcome",
             "eaf.outcome","beta.outcome","se.outcome","pval.outcome","samplesize.outcome"))
  dat$outcome <- nm; dat$id.outcome <- nm
  outcomes[[nm]] <- dat
}

# ---- PLINK clump helper ----
plink_clump <- function(dat, p_thresh, tmp_tag) {
  clump_file <- file.path(TMP_DIR, paste0(tmp_tag, "_clump.txt"))
  fwrite(data.frame(SNP = dat$SNP, P = dat$pval.exposure), clump_file, sep="\t", quote=FALSE)
  out_prefix <- file.path(TMP_DIR, tmp_tag)
  cmd <- sprintf('"%s" --bfile "%s" --clump "%s" --clump-p1 %g --clump-r2 0.001 --clump-kb 10000 --out "%s"',
                 PLINK_BIN, LD_REF, clump_file, p_thresh, out_prefix)
  system(cmd, ignore.stdout=TRUE, ignore.stderr=TRUE)
  f <- paste0(out_prefix, ".clumped")
  if (!file.exists(f)) return(NULL)
  clumped <- fread(f, select="SNP", skip="SNP")
  dat[SNP %in% clumped$SNP, ]
}

# ---- Process each lipid ----
all_results <- list()

for (lname in names(lipid_files)) {
  cat(sprintf("\n=== %s ===\n", lname))

  # Read lipid GWAS
  gwas <- fread(lipid_files[[lname]], select = c("rsID","CHROM","POS_b37","ALT","REF",
                                                   "POOLED_ALT_AF","EFFECT_SIZE","SE","pvalue","N"))
  setnames(gwas, c("rsID","CHROM","POS_b37","ALT","REF","POOLED_ALT_AF","EFFECT_SIZE","SE","pvalue","N"),
           c("SNP","chr.exposure","pos.exposure","effect_allele.exposure","other_allele.exposure",
             "eaf.exposure","beta.exposure","se.exposure","pval.exposure","samplesize.exposure"))
  gwas$exposure <- lname; gwas$id.exposure <- lname

  cat(sprintf("  SNPs: %s\n", format(nrow(gwas), big.mark=",")))

  # Clump
  p_thresh <- ifelse(any(gwas$pval.exposure < 5e-8, na.rm=TRUE), 5e-8, 1e-5)
  gwas_clumped <- plink_clump(gwas[gwas$pval.exposure < p_thresh, ], p_thresh, lname)
  if (is.null(gwas_clumped) || nrow(gwas_clumped) < 3) { cat("  SKIP\n"); next }

  # APOE exclusion
  gwas_clumped <- gwas_clumped[!(gwas_clumped$chr.exposure == 19 &
                                 gwas_clumped$pos.exposure >= 44000000 &
                                 gwas_clumped$pos.exposure <= 46500000), ]
  if (nrow(gwas_clumped) < 3) { cat("  SKIP after APOE\n"); next }

  cat(sprintf("  IVs: %d\n", nrow(gwas_clumped)))

  # MR vs 4 outcomes
  for (oname in names(outcomes)) {
    harm <- harmonise_data(gwas_clumped, outcomes[[oname]], action = 2)
    if (nrow(harm) < 3) next
    mr_res <- mr(harm, method_list = c("mr_ivw","mr_weighted_median","mr_egger_regression"))
    mr_res$exposure <- lname; mr_res$outcome <- oname
    all_results[[length(all_results) + 1]] <- mr_res

    ivw <- mr_res[mr_res$method == "Inverse variance weighted", ]
    if (nrow(ivw) > 0)
      cat(sprintf("  %-8s: b=%.4f se=%.4f P=%.2e\n", oname, ivw$b, ivw$se, ivw$pval))
  }
}

# ---- Save ----
all_mr <- rbindlist(all_results)
fwrite(all_mr, file.path(OUT_DIR, "C3_lipid_mr_rerun.csv"))
cat(sprintf("\nSaved: results/C3_lipid_mr_rerun.csv (%d rows)\n", nrow(all_mr)))

# Summary
ivw_all <- all_mr[method == "Inverse variance weighted", ]
cat("\n=== Lipid MR Summary (IVW) ===\n")
for (lname in names(lipid_files)) {
  for (oname in names(outcomes)) {
    row <- ivw_all[exposure == lname & outcome == oname, ]
    if (nrow(row) == 0) next
    sig <- if(row$pval < 0.05) "*" else if(row$pval < 0.01) "**" else ""
    cat(sprintf("%-8s → %-8s: b=%.4f P=%.2e %s\n", lname, oname, row$b, row$pval, sig))
  }
}
