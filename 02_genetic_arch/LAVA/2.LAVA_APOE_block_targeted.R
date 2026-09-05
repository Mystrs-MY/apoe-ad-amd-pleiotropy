#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(LAVA)
  library(data.table)
})

script_arg <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)][1])
script_dir <- dirname(normalizePath(script_file, winslash = "/", mustWork = TRUE))
run_tag <- Sys.getenv("A1_LAVA_RUN_TAG", unset = "correctedN_20260903")
input_dir <- Sys.getenv("A1_LAVA_INPUT_DIR", unset = file.path(script_dir, paste0("inputs_", run_tag)))
ref_prefix <- Sys.getenv("A1_LAVA_REF_PREFIX", unset = "")
output_dir <- Sys.getenv("A1_LAVA_OUTPUT_DIR", unset = file.path(script_dir, paste0("results_", run_tag)))
locus_file <- Sys.getenv("A1_LAVA_LOCUS_FILE", unset = file.path(script_dir, "blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile"))
target_locus_id <- as.integer(Sys.getenv("A1_LAVA_TARGET_LOCUS", unset = "2351"))
seed <- as.integer(Sys.getenv("A1_LAVA_SEED", unset = "20260903"))

if (!nzchar(ref_prefix)) stop("A1_LAVA_REF_PREFIX is not set.")
if (!file.exists(locus_file)) stop("Missing LAVA locus file: ", locus_file)
if (packageVersion("LAVA") < "0.1.5") stop("LAVA >= 0.1.5 is required.")
input_info_file <- file.path(input_dir, "LAVA_input_info.tsv")
if (!file.exists(input_info_file)) stop("Missing LAVA input info: ", input_info_file)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(seed)

loci <- read.loci(locus_file)
target <- loci[loci$LOC == target_locus_id, , drop = FALSE]
if (nrow(target) != 1L) stop("Expected one target locus, found ", nrow(target))
input <- process.input(input.info.file = input_info_file, sample.overlap.file = NULL, ref.prefix = ref_prefix,
                       phenos = c("AD", "Dry_AMD", "Wet_AMD", "Any_AMD"))
locus <- process.locus(target, input)
if (is.null(locus)) stop("Target locus could not be processed.")
biv <- run.bivar(locus, target = "AD")
if (is.null(biv) || nrow(biv) == 0L) stop("No bivariate result returned for target locus.")

result <- cbind(data.table(locus = locus$id, CHR = locus$chr, START = locus$start, STOP = locus$stop), as.data.table(biv))
fwrite(result, file.path(output_dir, paste0("LAVA_APOE_block_targeted_", run_tag, ".csv")), quote = FALSE)
fwrite(data.table(run_tag = run_tag, lava_version = as.character(packageVersion("LAVA")), reference_prefix = ref_prefix,
                  target_locus = target_locus_id, input_info = normalizePath(input_info_file, winslash = "/", mustWork = TRUE)),
       file.path(output_dir, "LAVA_APOE_targeted_manifest.tsv"), sep = "\t", quote = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "LAVA_APOE_targeted_sessionInfo.txt"))
print(result)
