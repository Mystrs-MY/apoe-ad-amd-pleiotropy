#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(LAVA)
  library(data.table)
  library(foreach)
  library(doParallel)
  library(doSNOW)
})

script_arg <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)][1])
script_dir <- dirname(normalizePath(script_file, winslash = "/", mustWork = TRUE))
run_tag <- Sys.getenv("A1_LAVA_RUN_TAG", unset = "correctedN_20260903")
input_dir <- Sys.getenv("A1_LAVA_INPUT_DIR", unset = file.path(script_dir, paste0("inputs_", run_tag)))
ref_prefix <- Sys.getenv("A1_LAVA_REF_PREFIX", unset = "")
output_dir <- Sys.getenv("A1_LAVA_OUTPUT_DIR", unset = file.path(script_dir, paste0("results_", run_tag)))
locus_file <- Sys.getenv("A1_LAVA_LOCUS_FILE", unset = file.path(script_dir, "blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile"))
n_cores <- max(1L, as.integer(Sys.getenv("A1_LAVA_CORES", unset = "4")))
force <- identical(tolower(Sys.getenv("A1_LAVA_FORCE", unset = "false")), "true")
seed <- as.integer(Sys.getenv("A1_LAVA_SEED", unset = "20260903"))

if (!nzchar(ref_prefix)) stop("A1_LAVA_REF_PREFIX is not set.")
if (!file.exists(locus_file)) stop("Missing locus file: ", locus_file)
if (packageVersion("LAVA") < "0.1.5") stop("LAVA >= 0.1.5 is required for the UKB binary LD reference.")
input_info_file <- file.path(input_dir, "LAVA_input_info.tsv")
if (!file.exists(input_info_file)) stop("Missing LAVA input info: ", input_info_file)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
checkpoint_dir <- file.path(output_dir, "checkpoints")
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(seed)
loci <- read.loci(locus_file)
input <- process.input(input.info.file = input_info_file, sample.overlap.file = NULL, ref.prefix = ref_prefix,
                       phenos = c("AD", "Dry_AMD", "Wet_AMD", "Any_AMD"))
checkpoint_path <- function(i) file.path(checkpoint_dir, sprintf("locus_%04d.rds", i))
pending <- seq_len(nrow(loci))
if (!force) pending <- pending[!file.exists(vapply(pending, checkpoint_path, character(1)))]

message("LAVA version: ", packageVersion("LAVA"))
message("Reference prefix: ", ref_prefix)
message("Loci total: ", nrow(loci), "; pending: ", length(pending))

if (length(pending)) {
  n_cores <- min(n_cores, length(pending))
  cl <- makeCluster(n_cores)
  on.exit(try(stopCluster(cl), silent = TRUE), add = TRUE)
  registerDoSNOW(cl)
  pb <- txtProgressBar(max = length(pending), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  foreach(i = pending, .packages = c("LAVA", "data.table"), .options.snow = opts,
          .export = c("loci", "input", "checkpoint_path")) %dopar% {
    result <- tryCatch({
      locus <- process.locus(loci[i, ], input)
      if (is.null(locus)) {
        list(status = "not_estimable", index = i, result = NULL, message = NA_character_)
      } else {
        biv <- run.bivar(locus, target = "AD")
        if (is.null(biv) || !nrow(biv)) {
          list(status = "not_estimable", index = i, result = NULL, message = NA_character_)
        } else {
          loc_base <- data.table(locus = locus$id, CHR = locus$chr, START = locus$start, STOP = locus$stop)
          list(status = "ok", index = i, result = cbind(loc_base, as.data.table(biv)), message = NA_character_)
        }
      }
    }, error = function(e) {
      list(status = "error", index = i, result = NULL, message = conditionMessage(e))
    })
    saveRDS(result, checkpoint_path(i))
    result$status
  }
  close(pb)
  stopCluster(cl)
}

checkpoint_files <- vapply(seq_len(nrow(loci)), checkpoint_path, character(1))
if (any(!file.exists(checkpoint_files))) stop("Missing LAVA checkpoint files.")
objects <- lapply(checkpoint_files, readRDS)
status <- rbindlist(lapply(objects, function(x) data.table(locus_index = x$index, status = x$status, message = x$message)), fill = TRUE)
fwrite(status, file.path(output_dir, "LAVA_locus_status.tsv"), sep = "\t", quote = FALSE)
all_results <- rbindlist(lapply(objects, function(x) x$result), fill = TRUE)
if (!nrow(all_results)) stop("No bivariate LAVA results were returned.")

extract_pair <- function(other) {
  ans <- all_results[(phen1 == "AD" & phen2 == other) | (phen1 == other & phen2 == "AD")]
  if ("p" %in% names(ans)) setorder(ans, p)
  ans
}
outputs <- list(DryAMD = extract_pair("Dry_AMD"), WetAMD = extract_pair("Wet_AMD"), AnyAMD = extract_pair("Any_AMD"))
for (nm in names(outputs)) fwrite(outputs[[nm]], file.path(output_dir, paste0("LAVA_FullScan_AD_vs_", nm, "_", run_tag, ".csv")), quote = FALSE)

manifest <- data.table(run_tag = run_tag, lava_version = as.character(packageVersion("LAVA")), reference_prefix = ref_prefix,
                       input_info = normalizePath(input_info_file, winslash = "/", mustWork = TRUE),
                       locus_file = normalizePath(locus_file, winslash = "/", mustWork = TRUE), seed = seed, cores = n_cores,
                       loci_total = nrow(loci), loci_ok = sum(status$status == "ok"),
                       loci_not_estimable = sum(status$status == "not_estimable"), loci_error = sum(status$status == "error"))
fwrite(manifest, file.path(output_dir, "LAVA_run_manifest.tsv"), sep = "\t", quote = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
if (manifest$loci_error > 0L) message("LAVA completed with recorded locus-level errors; inspect LAVA_locus_status.tsv.")
message("LAVA full scan completed: ", output_dir)
