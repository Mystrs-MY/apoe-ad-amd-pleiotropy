# ============================================================
# C6_load_ukbppp.R
# UKB-PPP protein mediators: lazy RSID lookup -> PLINK clump -> MR
# Lazy approach: only load chromosome RSID maps for needed chromosomes
# ============================================================

library(data.table)
library(TwoSampleMR)
library(dplyr)

set.seed(20240603)

# ---- Config ----
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = TRUE)
SCRIPT_DIR <- dirname(script_path)
PROJECT_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), winslash = "/", mustWork = TRUE)
RESOURCE_ROOT <- Sys.getenv("A1_RESOURCE_ROOT", unset = file.path(PROJECT_ROOT, "data", "external"))
GWAS_DIR <- file.path(RESOURCE_ROOT, "GWAS")
MERGED_DIR <- file.path(RESOURCE_ROOT, "ukbppp_proteins", "ukbppp_merged")
RSID_DIR <- file.path(RESOURCE_ROOT, "ukbppp_proteins", "rsid_maps")
RESULTS_DIR <- file.path(SCRIPT_DIR, "results")
PLINK_BIN <- Sys.getenv("PLINK_BIN", unset = plinkbinr::get_plink_exe())
LD_REF <- Sys.getenv("A1_LD_PREFIX", unset = file.path(RESOURCE_ROOT, "EUR", "EUR"))
TMP_DIR     <- tempdir()

APOE_CHR <- 19; APOE_START <- 44000000; APOE_END <- 46500000
dir.create(RESULTS_DIR, showWarnings = FALSE)

# ---- Helpers ----
load_outcome_gwas <- function(path, trait_name) {
  dat <- fread(path, select=c("SNP","CHR","BP","A1","A2","FREQ","BETA","SE","P","N"))
  setnames(dat, c("CHR","BP","A1","A2","FREQ","BETA","SE","P","N"),
           c("chr.outcome","pos.outcome","effect_allele.outcome","other_allele.outcome",
             "eaf.outcome","beta.outcome","se.outcome","pval.outcome","samplesize.outcome"))
  dat$outcome <- trait_name; dat$id.outcome <- trait_name
  return(dat)
}

exclude_apoe <- function(dat) {
  dat[!(dat$chr.exposure == APOE_CHR &
        dat$pos.exposure >= APOE_START &
        dat$pos.exposure <= APOE_END), ]
}

# PLINK clump via direct system call
plink_clump <- function(dat, p_thresh, plink_bin, bfile, tmp_tag) {
  clump_file <- file.path(TMP_DIR, paste0(tmp_tag, "_clump.txt"))
  fwrite(data.frame(SNP=dat$SNP, P=dat$pval.exposure), clump_file, sep="\t", quote=FALSE)
  out_prefix <- file.path(TMP_DIR, tmp_tag)

  cmd <- sprintf('"%s" --bfile "%s" --clump "%s" --clump-p1 %g --clump-r2 0.001 --clump-kb 10000 --out "%s"',
                 plink_bin, bfile, clump_file, p_thresh, out_prefix)
  system(cmd, ignore.stdout=TRUE, ignore.stderr=TRUE)

  clumped_file <- paste0(out_prefix, ".clumped")
  if (!file.exists(clumped_file)) return(NULL)

  clumped <- fread(clumped_file, select="SNP", skip="SNP")
  if (nrow(clumped) == 0) return(NULL)

  dat[dat$SNP %in% clumped$SNP, ]
}

# ---- Load outcomes (no RSID pre-load!) ----
cat("\n=== Loading 4 GWAS outcomes ===\n")
outcome_files <- list(
  AD      = paste0(GWAS_DIR, "AD_Wightman_cleaned_hg19.tsv.gz"),
  Dry_AMD = paste0(GWAS_DIR, "AMD_Dry_R12_cleaned_hg19.tsv.gz"),
  Wet_AMD = paste0(GWAS_DIR, "AMD_Wet_R12_cleaned_hg19.tsv.gz"),
  Any_AMD = paste0(GWAS_DIR, "AMD_H7_R12_cleaned_hg19.tsv.gz")
)

outcomes <- list()
for (name in names(outcome_files)) {
  cat(sprintf("  Loading %s...\n", name))
  outcomes[[name]] <- load_outcome_gwas(outcome_files[[name]], name)
  cat(sprintf("    %s SNPs\n", format(nrow(outcomes[[name]]), big.mark=",")))
}

# ---- Categories ----
protein_categories <- list(
  Complement   = c("C1QA","C2","C3","C5","CFB","CFD","CFH","CFI","CFP","SERPING1"),
  Inflammatory = c("CSF1","IFNG","IL10","IL18","IL1B","IL6","TGFB1","TNF"),
  Chemokine    = c("CCL2","CCL5","CX3CL1","CXCL10","CXCL12"),
  Immune       = c("CD14","CD40","TLR4","TREM2"),
  Lipid        = c("APOE","LPA","PON1")
)
gene_cat <- setNames(rep(names(protein_categories), lengths(protein_categories)),
                     unlist(protein_categories))

# ---- Process each protein ----
merged_files <- list.files(MERGED_DIR, pattern="_merged\\.txt\\.gz$", full.names=TRUE)
cat(sprintf("\n=== Processing %d proteins ===\n", length(merged_files)))

all_exposures <- list()
rs429358_report <- data.frame()
failed_genes <- c()

# Pre-index RSID mapping file list by chromosome
rsid_files <- list.files(RSID_DIR, pattern="olink_rsid_map.*_chr.*patched_v2\\.tsv\\.gz$", full.names=TRUE)
rsid_by_chr <- list()
for (cf in rsid_files) {
  bn <- basename(cf)
  # Extract chr number from filename: ..._chr10_... or ..._chr1_... or ..._chrXY_...
  chr_str <- gsub(".*_chr([0-9XY]+)_.*", "\\1", bn)
  rsid_by_chr[[chr_str]] <- cf
}
cat(sprintf("Indexed %d chromosome RSID maps\n", length(rsid_by_chr)))

for (f in merged_files) {
  gene <- toupper(gsub("_merged\\.txt\\.gz$", "", basename(f)))
  cat(sprintf("[%s] ", gene))

  # Read GWAS
  gwas <- tryCatch(fread(f, header=TRUE, showProgress=FALSE, fill=TRUE), error=function(e) NULL)
  if (is.null(gwas) || nrow(gwas) < 1000) { cat("READ ERR\n"); failed_genes <- c(failed_genes, gene); next }

  # rs429358 check
  rs_row <- gwas[grepl("45411941", ID), ]
  if (nrow(rs_row) > 0) {
    r <- rs_row[1,]; pv <- 10^(-r$LOG10P)
    rs429358_report <- rbind(rs429358_report, data.frame(
      Gene=gene, Category=ifelse(gene %in% names(gene_cat), gene_cat[[gene]], "Unknown"),
      BETA=r$BETA, SE=r$SE, P=pv, A1=r$ALLELE1, A0=r$ALLELE0,
      A1FREQ=r$A1FREQ, N=r$N, stringsAsFactors=FALSE))
  } else { cat("NO rs429358\n"); failed_genes <- c(failed_genes, gene); next }

  # Pre-filter to P<1e-5 (16M → ~5K rows)
  gwas$Pval <- 10^(-gwas$LOG10P)
  gwas <- gwas[gwas$Pval < 1e-5, ]
  if (nrow(gwas) < 2) { cat("no P<1e-5\n"); failed_genes <- c(failed_genes, gene); next }

  # Parse chr from ID, identify needed chromosomes
  parts <- strsplit(gwas$ID, ":")
  gwas$chr <- as.integer(sapply(parts, `[`, 1))
  gwas$pos <- as.integer(sapply(parts, `[`, 2))
  needed_chrs <- unique(as.character(gwas$chr))

  # Lazy load only the RSID maps we need
  gwas$rsid <- NA_character_
  for (chr_str in needed_chrs) {
    map_file <- rsid_by_chr[[chr_str]]
    if (is.null(map_file)) next  # XY or other

    # Load only this chromosome's mapping
    chr_map <- fread(map_file, select=c("ID","rsid"), showProgress=FALSE)
    # Match and assign rsIDs
    idx <- which(gwas$chr == as.integer(chr_str))
    m <- match(gwas$ID[idx], chr_map$ID)
    gwas$rsid[idx] <- ifelse(is.na(m), gwas$ID[idx], chr_map$rsid[m])
  }
  gwas$rsid[is.na(gwas$rsid)] <- gwas$ID[is.na(gwas$rsid)]

  # Format for TwoSampleMR (use rsID as SNP)
  gwas$SNP <- gwas$rsid
  gwas <- as.data.frame(gwas)
  cat(sprintf("%d SNPs -> ", nrow(gwas)), flush=TRUE)

  exposure <- format_data(gwas, type="exposure",
    snp_col="SNP", beta_col="BETA", se_col="SE",
    effect_allele_col="ALLELE1", other_allele_col="ALLELE0",
    eaf_col="A1FREQ", pval_col="Pval",
    samplesize_col="N", chr_col="chr", pos_col="pos")

  # PLINK clump
  p_thresh <- ifelse(any(exposure$pval.exposure < 5e-8), 5e-8, 1e-5)
  exp_clumped <- plink_clump(exposure[exposure$pval.exposure < p_thresh, ],
                              p_thresh, PLINK_BIN, LD_REF, gene)

  if (is.null(exp_clumped) || nrow(exp_clumped) < 2) {
    cat("SKIP (<2 IVs)\n"); failed_genes <- c(failed_genes, gene); next
  }

  exp_clean <- exclude_apoe(exp_clumped)
  if (nrow(exp_clean) < 2) { cat("SKIP (all APOE)\n"); failed_genes <- c(failed_genes, gene); next }

  cat(sprintf("%d IVs MR ", nrow(exp_clean)))
  all_exposures[[gene]] <- exp_clean

  # MR vs 4 outcomes
  for (outcome_name in names(outcomes)) {
    harm <- tryCatch(harmonise_data(exp_clean, outcomes[[outcome_name]], action=2), error=function(e) NULL)
    if (is.null(harm) || nrow(harm) < 2) next
    mr_res <- tryCatch(mr(harm, method_list=c("mr_ivw","mr_weighted_median","mr_egger_regression")),
                       error=function(e) NULL)
    if (is.null(mr_res)) next
    mr_res$exposure <- gene; mr_res$outcome <- outcome_name
    mr_res$category <- ifelse(gene %in% names(gene_cat), gene_cat[[gene]], "Unknown")
    fwrite(mr_res, file.path(RESULTS_DIR, sprintf("mr_%s_%s.csv", gene, outcome_name)))
  }
  cat("done\n")
}

# ---- Save ----
if (length(all_exposures) > 0) {
  saveRDS(all_exposures, file.path(RESULTS_DIR, "ukbppp_exposures.rds"))
  mr_files <- list.files(RESULTS_DIR, pattern="^mr_.*\\.csv$", full.names=TRUE)
  if (length(mr_files) > 0) {
    mr_all <- rbindlist(lapply(mr_files, fread), fill=TRUE)
    fwrite(mr_all, file.path(RESULTS_DIR, "C6_ukbppp_mr_all.csv"))
    cat(sprintf("\nCombined MR: %d rows\n", nrow(mr_all)))
  }
}
fwrite(rs429358_report, file.path(RESULTS_DIR, "C6_rs429358_effects.csv"))

if (length(failed_genes) > 0)
  cat(sprintf("Failed: %s\n", paste(failed_genes, collapse=", ")))
cat(sprintf("C6 done: %d processed, %d with IVs\n",
            length(merged_files) - length(failed_genes), length(all_exposures)))
