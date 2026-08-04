###############################################################################
# parse_mixer_json.R — 从 MiXeR JSON 自动生成 CSV
# 用法: 把此脚本放在 JSON 文件同目录, Rscript parse_mixer_json.R
# 输出: MiXeR_univariate.csv, MiXeR_bivariate.csv
# Public release: downstream multi-panel Venn assembly is intentionally excluded.
###############################################################################

library(jsonlite)
# Run from MiXeR/ directory (or set working directory accordingly)
# setwd("...02_genetic_arch/MiXeR")
OUT_TAB <- "../../tables"
dir.create(OUT_TAB, showWarnings=FALSE, recursive=TRUE)

# ---- 1. 扫描目录, 区分 univariate / bivariate ----------------------------
files <- list.files(pattern = "\\.json$")

univ_files <- c()
biv_files  <- c()

for (f in files) {
  js <- fromJSON(f)

  # MiXeR JSON 里: univariate → analysis=="univariate"
  #                bivariate  → analysis=="bivariate"
  analysis_type <- js$analysis
  if (is.null(analysis_type)) next

  if (analysis_type == "univariate") {
    univ_files <- c(univ_files, f)
  } else if (analysis_type == "bivariate") {
    biv_files <- c(biv_files, f)
  }
}

cat(sprintf("Found %d univariate, %d bivariate JSON files\n",
            length(univ_files), length(biv_files)))

# ---- 2. 解析 univariate --------------------------------------------------
univ_rows <- list()
for (f in univ_files) {
  js <- fromJSON(f)
  p  <- js$params                     # {pi, sig2_beta, sig2_zero}

  # 提取 AIC/BIC (可能在 optimize 列表的最后一个元素中)
  opt <- js$optimize
  if (is.data.frame(opt) || is.list(opt)) {
    if (is.list(opt) && !is.data.frame(opt) && length(opt) > 0) {
      last <- opt[[length(opt)]]
      aic  <- last$AIC
      bic  <- last$BIC
    } else {
      aic <- opt$AIC
      bic <- opt$BIC
    }
  } else {
    aic <- NA; bic <- NA
  }

  trait_name <- gsub("\\.json$", "", f)
  univ_rows[[f]] <- data.frame(
    Trait      = trait_name,
    pi         = p$pi,
    sig2_beta  = p$sig2_beta,
    sig2_zero  = p$sig2_zero,
    AIC        = if (is.null(aic)) NA else aic,
    BIC        = if (is.null(bic)) NA else bic,
    stringsAsFactors = FALSE
  )
}
univ_df <- do.call(rbind, univ_rows)
rownames(univ_df) <- NULL

# ---- 3. 解析 bivariate ---------------------------------------------------
biv_rows <- list()
for (f in biv_files) {
  js <- fromJSON(f)
  p  <- js$params
  # p$pi      = [pi1, pi2, pi12]
  # p$sig2_beta = [s1, s2]
  # p$sig2_zero = [z1, z2]
  # p$rho_beta, p$rho_zero

  pi_vec <- p$pi
  pi1  <- pi_vec[1]; pi2 <- pi_vec[2]; pi12 <- pi_vec[3]

  overlap_pct <- pi12 / min(pi1, pi2) * 100

  # 尝试从文件名推断 trait 名称
  fname <- gsub("\\.json$", "", f)
  # 文件名格式通常是 "AD_vs_AMD_Dry" 或 "Trait1_vs_Trait2"
  # 但也可能包含 univariate 的前缀, 需要从 params 推断
  # 简单策略: 使用文件名中的 vs 拆分
  parts <- strsplit(fname, "_vs_")[[1]]
  t1 <- if (length(parts) >= 2) parts[1] else "Trait1"
  t2 <- if (length(parts) >= 2) parts[2] else "Trait2"

  # 如果 t1 或 t2 是 univariate JSON 中的 trait 名, 尝试匹配
  if (exists("univ_df") && nrow(univ_df) > 0) {
    for (j in seq_len(nrow(univ_df))) {
      tn <- univ_df$Trait[j]
      if (grepl(tn, t1, fixed=TRUE) || grepl(t1, tn, fixed=TRUE)) t1 <- tn
      if (grepl(tn, t2, fixed=TRUE) || grepl(t2, tn, fixed=TRUE)) t2 <- tn
    }
  }

  opt <- js$optimize
  if (is.list(opt) && !is.data.frame(opt) && length(opt) > 0) {
    last <- opt[[length(opt)]]
    aic  <- last$AIC
    bic  <- last$BIC
  } else if (is.list(opt)) {
    aic <- opt$AIC; bic <- opt$BIC
  } else {
    aic <- NA; bic <- NA
  }

  comparison <- paste(t1, "vs.", t2)

  biv_rows[[f]] <- data.frame(
    Comparison   = comparison,
    Trait1       = t1,
    Trait2       = t2,
    pi1          = pi1,
    pi2          = pi2,
    pi12         = pi12,
    Overlap_pct  = round(overlap_pct, 1),
    rho_beta     = p$rho_beta,
    rho_zero     = p$rho_zero,
    AIC          = if (is.null(aic)) NA else aic,
    BIC          = if (is.null(bic)) NA else bic,
    stringsAsFactors = FALSE
  )
}
biv_df <- do.call(rbind, biv_rows)
rownames(biv_df) <- NULL

# ---- 4. 写 CSV -----------------------------------------------------------
if (nrow(univ_df) > 0) {
  write.csv(univ_df, file.path(OUT_TAB, "TableS3a_MiXeR_Univariate.csv"), row.names = FALSE, quote = FALSE)
  cat("TableS3a_MiXeR_Univariate.csv (", nrow(univ_df), "rows)\n", sep="")
}

if (nrow(biv_df) > 0) {
  write.csv(biv_df, file.path(OUT_TAB, "TableS3b_MiXeR_Bivariate.csv"), row.names = FALSE, quote = FALSE)
  cat("TableS3b_MiXeR_Bivariate.csv (", nrow(biv_df), "rows)\n", sep="")
}

# ---- 5. 打印预览 ----------------------------------------------------------
if (nrow(univ_df) > 0) {
  cat("\n--- Univariate ---\n")
  print(univ_df[, c("Trait","pi","sig2_beta","sig2_zero")], row.names=FALSE)
}
if (nrow(biv_df) > 0) {
  cat("\n--- Bivariate ---\n")
  print(biv_df[, c("Comparison","Overlap_pct","rho_beta","rho_zero")], row.names=FALSE)
}

cat("\nDone. MiXeR machine-readable summaries are ready.\n")
