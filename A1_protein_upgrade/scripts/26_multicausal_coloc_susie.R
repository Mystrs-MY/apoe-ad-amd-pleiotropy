#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
  library(susieR)
})

set.seed(20260712)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
region_dir <- file.path(root, "data_processed", "multicausal_coloc_regions")
out_dir <- file.path(root, "tables")
log_dir <- file.path(root, "logs", "multicausal_coloc")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

project_root <- normalizePath(file.path(root, ".."), winslash = "/", mustWork = TRUE)
resource_root <- Sys.getenv("A1_RESOURCE_ROOT", unset = file.path(project_root, "data", "external"))
plink <- Sys.getenv("PLINK_BIN", unset = "plink")
plink <- if (file.exists(plink)) plink else Sys.which(plink)
if (!nzchar(plink)) stop("PLINK binary not found. Set PLINK_BIN or add plink to PATH.")
bfile <- Sys.getenv("A1_LD_PREFIX", unset = file.path(resource_root, "EUR", "EUR"))
bim <- fread(paste0(bfile, ".bim"), header = FALSE,
             col.names = c("CHR", "ref_snp", "CM", "BP", "ref_A1", "ref_A2"))

genes <- data.table(
  gene = c("TREM2", "IL6R", "CFB", "CFH", "CFP", "PCSK9"),
  chr = c(6L, 1L, 6L, 1L, 23L, 1L),
  start = c(40126786L, 153377669L, 30913724L, 195621008L, 46413002L, 54505221L),
  end = c(42126786L, 155377669L, 32913724L, 197621008L, 48413002L, 56505221L),
  high_complexity = c(FALSE, FALSE, TRUE, TRUE, FALSE, FALSE)
)
outcomes <- data.table(
  outcome = c("AD", "Dry", "Wet", "Any"),
  prevalence = c(0.22, 0.037, 0.028, 0.036)
)
prespecified_pairs <- data.table(
  gene = c("TREM2", rep(c("IL6R", "CFB", "CFH", "PCSK9"), each = 3)),
  outcome = c("AD", rep(c("Dry", "Wet", "Any"), times = 4))
)
max_susie_snps <- 400L

complement <- function(x) chartr("ACGT", "TGCA", x)

orient_to_reference <- function(dt, effect, other, beta, freq, prefix) {
  direct <- dt[[effect]] == dt$ref_A1 & dt[[other]] == dt$ref_A2
  swapped <- dt[[effect]] == dt$ref_A2 & dt[[other]] == dt$ref_A1
  comp_direct <- complement(dt[[effect]]) == dt$ref_A1 & complement(dt[[other]]) == dt$ref_A2
  comp_swapped <- complement(dt[[effect]]) == dt$ref_A2 & complement(dt[[other]]) == dt$ref_A1
  sign <- fifelse(direct | comp_direct, 1, fifelse(swapped | comp_swapped, -1, NA_real_))
  dt[, (paste0(prefix, "_beta_ref")) := get(beta) * sign]
  dt[, (paste0(prefix, "_freq_ref")) := fifelse(sign == 1, get(freq), 1 - get(freq))]
  dt[, (paste0(prefix, "_orientation")) := fifelse(
    direct, "direct", fifelse(swapped, "swapped", fifelse(comp_direct, "complement", fifelse(comp_swapped, "complement_swapped", NA_character_)))
  )]
  dt[!is.na(sign)]
}

read_pqtl <- function(gene, ref) {
  x <- fread(file.path(region_dir, paste0(gene, "_pqtl_region.tsv.gz")))
  x <- x[TEST == "ADD" & is.finite(BETA) & is.finite(SE) & SE > 0]
  x[, parsed_BP := as.integer(tstrsplit(ID, ":", fixed = TRUE)[[2]])]
  x <- x[, .(BP = parsed_BP, pqtl_snp = ID, pqtl_A1 = ALLELE1,
             pqtl_A2 = ALLELE0, pqtl_freq = A1FREQ, pqtl_beta = BETA,
             pqtl_se = SE, pqtl_N = N, pqtl_INFO = INFO)]
  x <- merge(x, ref, by = "BP", allow.cartesian = TRUE)
  x <- orient_to_reference(x, "pqtl_A1", "pqtl_A2", "pqtl_beta", "pqtl_freq", "pqtl")
  x[pqtl_INFO >= 0.8 & pqtl_freq_ref >= 0.01 & pqtl_freq_ref <= 0.99]
}

read_gwas <- function(gene, outcome, ref) {
  path <- file.path(region_dir, paste0(gene, "_", outcome, "_gwas_region.tsv.gz"))
  if (!file.exists(path) || file.info(path)$size < 50) return(NULL)
  x <- fread(path)
  if (!nrow(x)) return(NULL)
  x <- x[is.finite(BETA) & is.finite(SE) & SE > 0]
  required_n <- c("N_TOTAL", "CASE_FRACTION")
  if (!all(required_n %in% names(x))) {
    stop("GWAS regional extract lacks explicit total-N/case-fraction fields: ", path)
  }
  x <- x[, .(BP = as.integer(BP), gwas_snp = SNP, gwas_A1 = A1,
             gwas_A2 = A2, gwas_freq = FREQ, gwas_beta = BETA,
             gwas_se = SE, gwas_N_total = N_TOTAL,
             gwas_case_fraction = CASE_FRACTION)]
  x <- merge(x, ref, by = "BP", allow.cartesian = TRUE)
  x <- orient_to_reference(x, "gwas_A1", "gwas_A2", "gwas_beta", "gwas_freq", "gwas")
  x[gwas_freq_ref >= 0.01 & gwas_freq_ref <= 0.99]
}

run_ld <- function(gene, snps) {
  prefix <- file.path(log_dir, paste0(gene, "_common"))
  list_path <- paste0(prefix, ".extract")
  fwrite(data.table(snp = unique(snps)), list_path, col.names = FALSE)
  args <- c("--bfile", bfile, "--extract", list_path, "--maf", "0.01", "--geno", "0.05",
            "--write-snplist", "--r", "square", "gz", "--threads", "4", "--memory", "12000",
            "--out", prefix)
  status <- system2(plink, args = args, stdout = paste0(prefix, ".stdout.txt"), stderr = paste0(prefix, ".stderr.txt"))
  if (status != 0) stop("PLINK LD failed for ", gene)
  order <- fread(paste0(prefix, ".snplist"), header = FALSE)[[1]]
  R <- as.matrix(fread(cmd = paste("gzip -cd", shQuote(paste0(prefix, ".ld.gz"))), header = FALSE))
  if (nrow(R) != length(order) || ncol(R) != length(order)) stop("LD dimension mismatch for ", gene)
  diag(R) <- 1
  rownames(R) <- order
  colnames(R) <- order
  list(order = order, R = R)
}

status_rows <- list()
signal_rows <- list()

for (i in seq_len(nrow(genes))) {
  g <- genes[i]
  ref <- bim[CHR == g$chr & BP >= g$start & BP <= g$end]
  if (!nrow(ref)) {
    status_rows[[length(status_rows) + 1L]] <- data.table(gene = g$gene, outcome = outcomes$outcome,
      status = "no_LD_reference_region", n_common = 0L, high_complexity = g$high_complexity)
    next
  }
  pqtl <- read_pqtl(g$gene, ref)
  gwas_list <- setNames(lapply(outcomes$outcome, function(o) read_gwas(g$gene, o, ref)), outcomes$outcome)
  available <- names(gwas_list)[vapply(gwas_list, function(x) !is.null(x) && nrow(x) >= 50, logical(1))]
  if (!length(available)) {
    status_rows[[length(status_rows) + 1L]] <- data.table(gene = g$gene, outcome = outcomes$outcome,
      status = "not_estimable_no_outcome_region", n_common = 0L, high_complexity = g$high_complexity)
    next
  }

  common <- Reduce(intersect, c(list(unique(pqtl$ref_snp)), lapply(gwas_list[available], function(x) unique(x$ref_snp))))
  common <- unique(common)
  if (length(common) < 50) {
    status_rows[[length(status_rows) + 1L]] <- data.table(gene = g$gene, outcome = available,
      status = "not_estimable_fewer_than_50_common_SNPs", n_common = length(common), high_complexity = g$high_complexity)
    next
  }

  selection_method <- "all_strictly_harmonized_common_SNPs"
  if (length(common) > max_susie_snps) {
    candidates <- pqtl[ref_snp %in% common,
      .(ref_snp, BP, min_p = 2 * pnorm(-abs(pqtl_beta_ref / pqtl_se)))]
    candidates <- candidates[!duplicated(ref_snp)]
    for (outcome_name in available) {
      tmp <- gwas_list[[outcome_name]][ref_snp %in% common,
        .(ref_snp, p_tmp = 2 * pnorm(-abs(gwas_beta_ref / gwas_se)))]
      tmp <- tmp[!duplicated(ref_snp)]
      candidates <- merge(candidates, tmp, by = "ref_snp", all.x = TRUE)
      candidates[, min_p := pmin(min_p, p_tmp, na.rm = TRUE)][, p_tmp := NULL]
    }
    top_n <- max_susie_snps %/% 2L
    top <- candidates[order(min_p)][seq_len(top_n), ref_snp]
    remaining <- candidates[!ref_snp %in% top][order(BP)]
    grid_index <- unique(round(seq(1, nrow(remaining), length.out = max_susie_snps - length(top))))
    common <- unique(c(top, remaining[grid_index, ref_snp]))
    selection_method <- "top_signal_plus_genomic_grid_sensitivity"
  }

  ld <- tryCatch(run_ld(g$gene, common), error = identity)
  if (inherits(ld, "error")) {
    status_rows[[length(status_rows) + 1L]] <- data.table(gene = g$gene, outcome = available,
      status = paste0("LD_error: ", conditionMessage(ld)), n_common = length(common), high_complexity = g$high_complexity)
    next
  }
  common <- ld$order
  pdat <- pqtl[match(common, ref_snp)]
  if (anyNA(pdat$ref_snp)) stop("pQTL ordering failure for ", g$gene)
  d1 <- list(beta = pdat$pqtl_beta_ref, varbeta = pdat$pqtl_se^2, snp = common,
             position = pdat$BP, type = "quant", N = median(pdat$pqtl_N, na.rm = TRUE),
             MAF = pdat$pqtl_freq_ref, LD = ld$R)
  s1 <- tryCatch(runsusie(d1, L = 3, maxit = 50, coverage = 0.95), error = identity)

  for (outcome in outcomes$outcome) {
    if (!any(prespecified_pairs$gene == g$gene & prespecified_pairs$outcome == outcome)) {
      status_rows[[length(status_rows) + 1L]] <- data.table(gene = g$gene, outcome = outcome,
        status = "not_prespecified_multicausal_pair", n_common = length(common),
        SNP_selection = selection_method, high_complexity = g$high_complexity)
      next
    }
    odat <- gwas_list[[outcome]]
    if (is.null(odat) || nrow(odat) < 50) {
      status_rows[[length(status_rows) + 1L]] <- data.table(gene = g$gene, outcome = outcome,
        status = "not_estimable_no_outcome_region", n_common = 0L,
        SNP_selection = selection_method, high_complexity = g$high_complexity)
      next
    }
    odat <- odat[match(common, ref_snp)]
    if (anyNA(odat$ref_snp)) {
      status_rows[[length(status_rows) + 1L]] <- data.table(gene = g$gene, outcome = outcome,
        status = "not_estimable_common_set_mismatch", n_common = sum(!is.na(odat$ref_snp)),
        SNP_selection = selection_method, high_complexity = g$high_complexity)
      next
    }
    prev <- outcomes$prevalence[outcomes$outcome == outcome]
    observed_case_fraction <- median(odat$gwas_case_fraction, na.rm = TRUE)
    if (!is.finite(observed_case_fraction) || observed_case_fraction <= 0 || observed_case_fraction >= 1) {
      stop("Invalid case fraction for ", g$gene, " / ", outcome)
    }
    d2 <- list(beta = odat$gwas_beta_ref, varbeta = odat$gwas_se^2, snp = common,
               position = odat$BP, type = "cc", s = observed_case_fraction,
               N = median(odat$gwas_N_total, na.rm = TRUE),
               MAF = odat$gwas_freq_ref, LD = ld$R)
    s2 <- tryCatch(runsusie(d2, L = 3, maxit = 50, coverage = 0.95), error = identity)
    if (inherits(s1, "error") || inherits(s2, "error")) {
      msg <- paste(c(if (inherits(s1, "error")) paste0("pQTL: ", conditionMessage(s1)),
                     if (inherits(s2, "error")) paste0("GWAS: ", conditionMessage(s2))), collapse = " | ")
      status_rows[[length(status_rows) + 1L]] <- data.table(gene = g$gene, outcome = outcome,
        status = paste0("SuSiE_error: ", msg), n_common = length(common),
        SNP_selection = selection_method, high_complexity = g$high_complexity)
      next
    }
    coloc_result <- tryCatch(coloc.susie(s1, s2), error = identity)
    if (inherits(coloc_result, "error")) {
      status_rows[[length(status_rows) + 1L]] <- data.table(gene = g$gene, outcome = outcome,
        status = paste0("coloc_susie_error: ", conditionMessage(coloc_result)), n_common = length(common),
        SNP_selection = selection_method, high_complexity = g$high_complexity)
      next
    }
    summary <- as.data.table(coloc_result$summary)
    summary[, `:=`(gene = g$gene, outcome = outcome, n_common = length(common),
                   gwas_N_definition = "per-variant total analysed N",
                   gwas_case_fraction = observed_case_fraction,
                   SNP_selection = selection_method, high_complexity = g$high_complexity)]
    signal_rows[[length(signal_rows) + 1L]] <- summary
    status_rows[[length(status_rows) + 1L]] <- data.table(gene = g$gene, outcome = outcome,
      status = "estimated", n_common = length(common), n_pqtl_signals = length(s1$sets$cs),
      n_gwas_signals = length(s2$sets$cs), n_signal_pairs = nrow(summary),
      SNP_selection = selection_method, high_complexity = g$high_complexity)
  }
}

status <- rbindlist(status_rows, fill = TRUE)
signals <- rbindlist(signal_rows, fill = TRUE)
fwrite(status, file.path(out_dir, "multicausal_coloc_status.tsv"), sep = "\t")
fwrite(signals, file.path(out_dir, "multicausal_coloc_signal_pairs.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo.txt"))
stopifnot(nrow(status) == nrow(genes) * nrow(outcomes))
message("Multi-causal coloc status rows: ", nrow(status), "; signal-pair rows: ", nrow(signals))
