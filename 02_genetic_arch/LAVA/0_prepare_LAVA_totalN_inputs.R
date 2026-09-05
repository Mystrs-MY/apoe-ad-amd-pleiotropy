#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

script_arg <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)][1])
script_dir <- dirname(normalizePath(script_file, winslash = "/", mustWork = TRUE))
project_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)
resource_root <- Sys.getenv("A1_RESOURCE_ROOT", unset = file.path(project_root, "data", "external"))
run_tag <- Sys.getenv("A1_LAVA_RUN_TAG", unset = "correctedN_20260903")
output_dir <- Sys.getenv(
  "A1_LAVA_INPUT_DIR",
  unset = file.path(script_dir, paste0("inputs_", run_tag))
)
force <- identical(tolower(Sys.getenv("A1_LAVA_FORCE", unset = "false")), "true")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

traits <- data.table(
  phenotype = c("AD", "Dry_AMD", "Wet_AMD", "Any_AMD"),
  source_file = file.path(
    resource_root, "GWAS",
    c(
      "AD_Wightman_cleaned_hg19.tsv.gz",
      "AMD_Dry_R12_cleaned_hg19.tsv.gz",
      "AMD_Wet_R12_cleaned_hg19.tsv.gz",
      "AMD_H7_R12_cleaned_hg19.tsv.gz"
    )
  ),
  output_file = file.path(
    output_dir,
    c("AD.tsv.gz", "Dry_AMD.tsv.gz", "Wet_AMD.tsv.gz", "Any_AMD.tsv.gz")
  ),
  cases = c(39918L, 8570L, 6699L, 12495L),
  controls = c(358140L, 329258L, 331070L, 461686L)
)

required <- c("SNP", "A1", "A2", "BETA", "P", "N_TOTAL")
manifest_rows <- vector("list", nrow(traits))

for (i in seq_len(nrow(traits))) {
  x <- traits[i]
  if (!file.exists(x$source_file)) stop("Missing cleaned GWAS: ", x$source_file)
  reused <- !force && file.exists(x$output_file) &&
    file.info(x$output_file)$mtime >= file.info(x$source_file)$mtime

  if (!reused) {
    message("Preparing LAVA total-N input: ", x$phenotype)
    dat <- fread(x$source_file, select = required, showProgress = TRUE)
    missing_cols <- setdiff(required, names(dat))
    if (length(missing_cols)) {
      stop("Missing required columns in ", x$source_file, ": ", paste(missing_cols, collapse = ", "))
    }
    if (anyNA(dat$N_TOTAL) || any(dat$N_TOTAL <= 0)) stop("Invalid N_TOTAL values in ", x$source_file)
    out <- dat[, .(SNP, A1, A2, BETA, P, N = N_TOTAL)]
    fwrite(out, x$output_file, sep = "\t", quote = FALSE, compress = "gzip")
    rm(dat, out)
    invisible(gc())
  } else {
    message("Reusing current LAVA input: ", x$output_file)
  }

  qc <- fread(x$output_file, select = c("SNP", "N"), showProgress = FALSE)
  manifest_rows[[i]] <- data.table(
    run_tag = run_tag,
    phenotype = x$phenotype,
    cases = x$cases,
    controls = x$controls,
    source_file = normalizePath(x$source_file, winslash = "/", mustWork = TRUE),
    output_file = normalizePath(x$output_file, winslash = "/", mustWork = TRUE),
    rows = nrow(qc),
    n_min = min(qc$N),
    n_median = median(qc$N),
    n_max = max(qc$N),
    md5 = unname(tools::md5sum(x$output_file))
  )
  rm(qc)
  invisible(gc())
}

manifest <- rbindlist(manifest_rows)
fwrite(manifest, file.path(output_dir, "LAVA_input_manifest.tsv"), sep = "\t", quote = FALSE)
fwrite(
  manifest[, .(phenotype, cases, controls, filename = output_file)],
  file.path(output_dir, "LAVA_input_info.tsv"), sep = "\t", quote = FALSE
)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
message("LAVA inputs prepared in: ", output_dir)
