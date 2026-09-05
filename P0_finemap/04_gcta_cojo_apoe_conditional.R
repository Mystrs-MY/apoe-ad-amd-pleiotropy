# Reproducible APOE-region approximate conditional analysis with GCTA-COJO.
# This is a reference-panel-limited sensitivity analysis because the available
# 1000 Genomes EUR LD panel contains 503 individuals.

rm(list = ls())

required <- c("data.table")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R package(s): ", paste(missing, collapse = ", "))

suppressPackageStartupMessages(library(data.table))

script_arg <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)][1])
SCRIPT_DIR <- dirname(normalizePath(script_file, winslash = "/", mustWork = TRUE))
PROJECT_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), winslash = "/", mustWork = TRUE)
RESOURCE_ROOT <- Sys.getenv("A1_RESOURCE_ROOT", unset = file.path(PROJECT_ROOT, "data", "external"))
GWAS_DIR <- file.path(RESOURCE_ROOT, "GWAS")
LD_PREFIX <- Sys.getenv("A1_LD_PREFIX", unset = file.path(RESOURCE_ROOT, "EUR", "EUR"))
GCTA_BIN <- Sys.getenv("GCTA_BIN", unset = "gcta64")
RUN_TAG <- Sys.getenv("A1_COJO_RUN_TAG", unset = "correctedN_20260903")
FORCE_RERUN <- tolower(Sys.getenv("A1_COJO_FORCE", unset = "false")) %in%
  c("1", "true", "yes")
OUT_DIR <- file.path(
  PROJECT_ROOT, "P0_finemap", "results",
  paste0("gcta_cojo_apoe_", RUN_TAG)
)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

gcta_resolved <- if (file.exists(GCTA_BIN)) GCTA_BIN else Sys.which(GCTA_BIN)
if (!nzchar(gcta_resolved)) stop("GCTA binary not found. Set GCTA_BIN or add gcta64 to PATH.")
GCTA_BIN <- gcta_resolved
required_files <- paste0(LD_PREFIX, c(".bed", ".bim", ".fam"))
if (any(!file.exists(required_files))) {
  stop("Missing required file(s):\n", paste(required_files[!file.exists(required_files)], collapse = "\n"))
}

APOE_CHR <- 19L
APOE_START <- 44000000L
APOE_END <- 46500000L
REFERENCE_N <- length(readLines(paste0(LD_PREFIX, ".fam")))
if (REFERENCE_N != 503L) stop("Unexpected LD reference sample size: ", REFERENCE_N)

traits <- data.table(
  trait_id = c("AD", "Dry_AMD", "Wet_AMD", "Any_AMD"),
  trait_label = c("AD", "Dry AMD", "Wet AMD", "Any AMD"),
  file_name = c(
    "AD_Wightman_cleaned_hg19.tsv.gz",
    "AMD_Dry_R12_cleaned_hg19.tsv.gz",
    "AMD_Wet_R12_cleaned_hg19.tsv.gz",
    "AMD_H7_R12_cleaned_hg19.tsv.gz"
  )
)

scenarios <- list(
  rs429358 = "rs429358",
  rs429358_rs7412 = c("rs429358", "rs7412")
)

bim <- fread(
  paste0(LD_PREFIX, ".bim"),
  col.names = c("CHR_ref", "SNP", "CM_ref", "BP_ref", "A1_ref", "A2_ref")
)
bim <- bim[CHR_ref == APOE_CHR & BP_ref >= APOE_START & BP_ref <= APOE_END]
if (anyDuplicated(bim$SNP)) stop("Duplicated SNP identifiers in APOE LD reference subset.")

harmonize_to_reference <- function(gwas, trait_id) {
  merged <- merge(gwas, bim, by = "SNP", all = FALSE)
  merged <- merged[CHR == CHR_ref & BP == BP_ref]
  merged[, allele_relation := fifelse(
    A1 == A1_ref & A2 == A2_ref, "same",
    fifelse(A1 == A2_ref & A2 == A1_ref, "swapped", "mismatch")
  )]

  mismatch <- merged[allele_relation == "mismatch"]
  if (nrow(mismatch)) {
    fwrite(mismatch, file.path(OUT_DIR, paste0(trait_id, "_allele_mismatch.tsv")), sep = "\t", na = "NA")
  }
  merged <- merged[allele_relation != "mismatch"]
  if (!nrow(merged)) stop("No reference-matched APOE variants for ", trait_id)

  # GCTA accepts summary alleles as supplied. FREQ is the frequency of A1.
  merged[, P := pmax(as.numeric(P), 1e-300)]
  merged[, N := as.numeric(N)]
  merged <- merged[is.finite(BETA) & is.finite(SE) & SE > 0 & is.finite(P) &
                     is.finite(FREQ) & FREQ > 0 & FREQ < 1 & is.finite(N) & N > 0]
  merged <- unique(merged, by = "SNP")

  missing_cond <- setdiff(c("rs429358", "rs7412"), merged$SNP)
  if (length(missing_cond)) {
    stop(trait_id, " is missing required condition variant(s): ", paste(missing_cond, collapse = ", "))
  }
  merged[]
}

run_gcta <- function(ma_file, condition_file, prefix) {
  cma_file <- paste0(prefix, ".cma.cojo")
  if (!FORCE_RERUN && file.exists(cma_file) && file.info(cma_file)$size > 0) {
    message("Reusing completed GCTA output: ", cma_file)
    return(cma_file)
  }
  args <- c(
    "--bfile", normalizePath(LD_PREFIX, winslash = "/", mustWork = FALSE),
    "--chr", as.character(APOE_CHR),
    "--maf", "0.01",
    "--cojo-file", normalizePath(ma_file, winslash = "/", mustWork = TRUE),
    "--cojo-cond", normalizePath(condition_file, winslash = "/", mustWork = TRUE),
    "--thread-num", "4",
    "--out", normalizePath(prefix, winslash = "/", mustWork = FALSE)
  )
  stdout_file <- paste0(prefix, ".stdout.log")
  stderr_file <- paste0(prefix, ".stderr.log")
  status <- system2(GCTA_BIN, args = args, stdout = stdout_file, stderr = stderr_file)
  if (!identical(status, 0L)) {
    stop("GCTA-COJO failed (status ", status, "). Inspect ", stdout_file, " and ", stderr_file)
  }
  if (!file.exists(cma_file)) stop("Expected GCTA output is missing: ", cma_file)
  cma_file
}

all_results <- list()
all_summary <- list()
command_manifest <- list()

for (i in seq_len(nrow(traits))) {
  trait <- traits[i]
  gwas_path <- file.path(GWAS_DIR, trait$file_name)
  if (!file.exists(gwas_path)) stop("Missing GWAS: ", gwas_path)

  gwas <- fread(gwas_path, select = c("SNP", "CHR", "BP", "A1", "A2", "FREQ", "BETA", "SE", "P", "N"))
  gwas <- gwas[CHR == APOE_CHR & BP >= APOE_START & BP <= APOE_END]
  gwas <- harmonize_to_reference(gwas, trait$trait_id)

  ma <- gwas[, .(SNP, A1, A2, freq = FREQ, b = BETA, se = SE, p = P, N)]
  ma_file <- file.path(OUT_DIR, paste0(trait$trait_id, "_APOE_region.ma"))
  fwrite(ma, ma_file, sep = "\t", quote = FALSE, na = "NA")

  harmonization_audit <- gwas[, .(
    trait = trait$trait_label,
    SNP, CHR, BP, A1, A2, FREQ, BETA, SE, P, N,
    A1_ref, A2_ref, allele_relation
  )]
  fwrite(
    harmonization_audit,
    file.path(OUT_DIR, paste0(trait$trait_id, "_APOE_harmonization_audit.tsv")),
    sep = "\t", na = "NA"
  )

  for (scenario_name in names(scenarios)) {
    condition_snps <- scenarios[[scenario_name]]
    condition_file <- file.path(OUT_DIR, paste0("condition_", scenario_name, ".txt"))
    writeLines(condition_snps, condition_file)

    prefix <- file.path(OUT_DIR, paste0(trait$trait_id, "_cond_", scenario_name))
    cma_file <- run_gcta(ma_file, condition_file, prefix)
    cma <- fread(cma_file)

    required_cols <- c("Chr", "SNP", "bp", "refA", "freq", "b", "se", "p", "bC", "bC_se", "pC")
    if (!all(required_cols %in% names(cma))) {
      stop("Unexpected GCTA .cma.cojo schema: ", paste(names(cma), collapse = ", "))
    }
    cma[, `:=`(
      trait_id = trait$trait_id,
      trait = trait$trait_label,
      scenario = scenario_name,
      conditioned_snps = paste(condition_snps, collapse = "+"),
      ld_reference = "1000_Genomes_EUR",
      ld_reference_n = REFERENCE_N,
      reference_panel_limited_sensitivity = TRUE,
      sample_size_field = "N_EFFECTIVE",
      gcta_version = "1.95.1",
      analysis_run_tag = RUN_TAG
    )]
    all_results[[length(all_results) + 1L]] <- cma

    unconditioned <- copy(cma)[order(p)][1]
    residual <- copy(cma)[!SNP %in% condition_snps][order(pC)][1]
    key_rs10414043 <- cma[SNP == "rs10414043"][1]
    key_rs7412 <- cma[SNP == "rs7412"][1]
    all_summary[[length(all_summary) + 1L]] <- data.table(
      trait_id = trait$trait_id,
      trait = trait$trait_label,
      scenario = scenario_name,
      conditioned_snps = paste(condition_snps, collapse = "+"),
      n_variants = nrow(cma),
      top_unconditioned_snp = unconditioned$SNP,
      top_unconditioned_p = unconditioned$p,
      top_conditional_snp = residual$SNP,
      top_conditional_p = residual$pC,
      rs10414043_beta_conditional = key_rs10414043$bC,
      rs10414043_se_conditional = key_rs10414043$bC_se,
      rs10414043_p_conditional = key_rs10414043$pC,
      rs7412_beta_conditional = key_rs7412$bC,
      rs7412_se_conditional = key_rs7412$bC_se,
      rs7412_p_conditional = key_rs7412$pC,
      ld_reference = "1000_Genomes_EUR",
      ld_reference_n = REFERENCE_N,
      reference_panel_limited_sensitivity = TRUE
    )
    command_manifest[[length(command_manifest) + 1L]] <- data.table(
      trait = trait$trait_label,
      scenario = scenario_name,
      gwas_resource = file.path("data/external/GWAS", trait$file_name),
      ma_file = file.path("P0_finemap/results", basename(OUT_DIR), basename(ma_file)),
      condition_file = file.path("P0_finemap/results", basename(OUT_DIR), basename(condition_file)),
      output_prefix = file.path("P0_finemap/results", basename(OUT_DIR), basename(prefix))
    )
  }
}

results <- rbindlist(all_results, fill = TRUE, use.names = TRUE)
summary <- rbindlist(all_summary, fill = TRUE, use.names = TRUE)
manifest <- rbindlist(command_manifest, fill = TRUE, use.names = TRUE)

fwrite(results, file.path(OUT_DIR, "APOE_COJO_conditional_source.tsv"), sep = "\t", na = "NA")
fwrite(summary, file.path(OUT_DIR, "APOE_COJO_conditional_summary.tsv"), sep = "\t", na = "NA")
fwrite(manifest, file.path(OUT_DIR, "APOE_COJO_run_manifest.tsv"), sep = "\t", na = "NA")
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
writeLines(
  c(
    paste0("analysis_run_tag=", RUN_TAG),
    paste0("force_rerun=", FORCE_RERUN),
    "gwas_directory=data/external/GWAS",
    "ld_reference=data/external/EUR/EUR",
    paste0("ld_reference_n=", REFERENCE_N)
    ,"sample_size_field=N (explicit copy of N_EFFECTIVE)"
    ,"AD_N_definition=variant-specific_N_effective supplied by GCST013196"
    ,"AMD_N_definition=4/(1/N_cases+1/N_controls) from final FinnGen R12 counts"
  ),
  file.path(OUT_DIR, "run_metadata.txt")
)

cat("GCTA-COJO APOE conditional sensitivity completed.\n")
print(summary)
