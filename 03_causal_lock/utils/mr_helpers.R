# ============================================================
# NetMR Article 1 — MR 通用工具函数
# ============================================================

# APOE 区域 (hg19/GRCh37)
APOE_CHR <- 19
APOE_START <- 44000000
APOE_END <- 46500000

# ---------- 数据加载 ----------

#' 加载本地 GWAS 并格式化为 TwoSampleMR 暴露格式
#' @param path GWAS 文件路径 (tsv.gz)
#' @param trait_name 性状名称
#' @param p_threshold P 值过滤阈值 (默认 1, 不过滤)
#' @return TwoSampleMR 格式的 data.table
load_local_gwas <- function(path, trait_name, p_threshold = 1) {
  cat(sprintf("[LOAD] %s <- %s\n", trait_name, basename(path)))
  dat <- data.table::fread(path)

  # 标准化列名
  data.table::setnames(dat,
    old = c("CHR", "BP", "A1", "A2", "FREQ", "BETA", "SE", "P", "N"),
    new = c("chr.exposure", "pos.exposure", "effect_allele.exposure",
            "other_allele.exposure", "eaf.exposure", "beta.exposure",
            "se.exposure", "pval.exposure", "samplesize.exposure"),
    skip_absent = TRUE
  )

  dat$exposure <- trait_name
  dat$id.exposure <- trait_name
  dat$SNP <- dat$SNP  # 保持 rsID

  # P 值过滤
  if (p_threshold < 1) {
    dat <- dat[dat$pval.exposure < p_threshold, ]
    cat(sprintf("  [Filter] P < %.0e: %d SNPs\n", p_threshold, nrow(dat)))
  }

  cat(sprintf("  [Done] %d SNPs loaded\n", nrow(dat)))
  return(dat)
}

#' 加载本地 GWAS 并格式化为 TwoSampleMR 结局格式
#' @param gwas data.table (已 load_local_gwas)
#' @param outcome_name 结局名称
#' @param snp_list SNP 列表
#' @return TwoSampleMR 结局格式的 data.table
prepare_outcome <- function(gwas, outcome_name, snp_list) {
  out <- gwas[gwas$SNP %in% snp_list, ]
  colnames(out) <- gsub("exposure", "outcome", colnames(out))
  out$outcome <- outcome_name
  out$id.outcome <- outcome_name
  return(out)
}

#' 验证 rs429358 等位基因方向
#' @param gwases 命名 list of GWAS data.tables
verify_rs429358 <- function(gwases) {
  cat("\n=== rs429358 Allele Verification ===\n")
  for (nm in names(gwases)) {
    row <- gwases[[nm]][gwases[[nm]]$SNP == "rs429358", ]
    if (nrow(row) == 0) {
      cat(sprintf("[WARN] %s: rs429358 NOT FOUND!\n", nm))
    } else {
      cat(sprintf("  %-10s: A1=%s, A2=%s, EAF=%.3f, BETA=%+.4f, SE=%.4f, P=%.2e\n",
        nm, row$effect_allele.exposure, row$other_allele.exposure,
        row$eaf.exposure, row$beta.exposure, row$se.exposure, row$pval.exposure))
    }
  }
  cat("  Expected AD:     A1=T, A2=C, BETA=-1.1275\n")
  cat("  Expected AMD(s): A1=C, A2=T, BETA≈-0.21\n")
  cat("  -> AD and AMD have FLIPPED effect alleles!\n")
  cat("  -> harmonise_data(action=2) should auto-handle this.\n\n")
}

# ---------- Clumping ----------

#' LD Clumping
#' @param gwas TwoSampleMR 格式 data.table
#' @param p_threshold P 阈值
#' @param plink_bin PLINK binary path
#' @param bfile LD reference path
#' @return clumped data.table
ld_clump_local <- function(gwas, p_threshold = 5e-8,
                            plink_bin = plinkbinr::get_plink_exe(),
                            bfile = "../../Resource/EUR/EUR") {

  sig <- gwas[gwas$pval.exposure < p_threshold, ]
  if (nrow(sig) == 0) {
    warning("No significant SNPs after P filter")
    return(NULL)
  }

  cat(sprintf("  Clumping: P < %.0e, %d SNPs -> ", p_threshold, nrow(sig)))

  clumped <- ieugwasr::ld_clump(
    dat = dplyr::tibble(rsid = sig$SNP, pval = sig$pval.exposure, id = sig$id.exposure),
    clump_kb = 10000,
    clump_r2 = 0.001,
    clump_p = p_threshold,
    plink_bin = plink_bin,
    bfile = bfile,
    pop = "EUR"
  )

  keep <- paste(gwas$SNP, gwas$id.exposure) %in% paste(clumped$rsid, clumped$id)
  result <- gwas[keep, ]
  cat(sprintf("%d IVs\n", nrow(result)))
  return(result)
}

# ---------- APOE Exclusion ----------

#' 排除 APOE 区域
#' @param dat harmonised 数据或 exposure 数据
#' @return 排除后的 data.table
exclude_apoe <- function(dat) {
  # 始终使用 exposure 侧的 chr/pos 来判断 APOE 区域
  # harmonised 数据同时有 chr.exposure 和 chr.outcome
  # 排除依据: 工具变量所在的染色体位置 (即 exposure 侧)
  chr_col <- "chr.exposure"
  pos_col <- "pos.exposure"

  if (!chr_col %in% names(dat) || !pos_col %in% names(dat)) {
    warning("chr.exposure or pos.exposure not found — cannot exclude APOE region")
    return(dat)
  }

  n_before <- nrow(dat)
  keep <- !(dat[[chr_col]] == APOE_CHR &
            dat[[pos_col]] >= APOE_START &
            dat[[pos_col]] <= APOE_END)
  result <- dat[keep, ]
  n_removed <- n_before - nrow(result)
  if (n_removed > 0) {
    cat(sprintf("  [APOE Exclusion] chr19:44-46.5Mb: %d -> %d IVs (removed %d)\n",
        n_before, nrow(result), n_removed))
  }
  return(result)
}

# ---------- MR 执行 ----------

#' 执行标准 MR 分析流程
#' @param exp_dat 暴露数据 (clumped)
#' @param out_gwas 结局 GWAS 完整 data.table
#' @param exp_label 暴露标签
#' @param out_label 结局标签
#' @param exclude_apoe 是否排除 APOE 区域
#' @param n_min_iv 最小 IV 数量
#' @return list(dat, mr_results, hetero, pleio, presso, steiger, n_iv)
run_standard_mr <- function(exp_dat, out_gwas, exp_label, out_label,
                             exclude_apoe = FALSE, n_min_iv = 3L) {

  cat(sprintf("\n--- %s -> %s ---\n", exp_label, out_label))

  # 提取结局 SNP
  snps <- exp_dat$SNP
  out_dat <- out_gwas[out_gwas$SNP %in% snps, ]

  if (nrow(out_dat) == 0) {
    cat("  [SKIP] No overlapping SNPs in outcome\n")
    return(NULL)
  }

  # 格式化为 outcome 格式
  out_formatted <- out_dat
  colnames(out_formatted) <- gsub("exposure", "outcome", colnames(out_formatted))
  out_formatted$outcome <- out_label
  out_formatted$id.outcome <- out_label

  # Harmonise (action=2 尝试正向链推断)
  dat <- TwoSampleMR::harmonise_data(
    exposure_dat = exp_dat,
    outcome_dat = out_formatted,
    action = 2
  )

  # ⚠️ 验证 mr_keep
  dat <- dat[dat$mr_keep, ]

  if (nrow(dat) < n_min_iv) {
    cat(sprintf("  [SKIP] Only %d IVs after harmonise (need >= %d)\n", nrow(dat), n_min_iv))
    return(NULL)
  }

  # APOE 排除 (在 harmonise 之后)
  if (exclude_apoe) {
    dat <- exclude_apoe(dat)
    if (nrow(dat) < n_min_iv) {
      cat(sprintf("  [SKIP] Only %d IVs after APOE exclusion\n", nrow(dat)))
      return(NULL)
    }
  }

  # F 统计量
  dat$F_stat <- (dat$beta.exposure / dat$se.exposure)^2
  mean_f <- mean(dat$F_stat, na.rm = TRUE)
  cat(sprintf("  IVs=%d, Mean F=%.1f\n", nrow(dat), mean_f))

  # MR 分析 (4 种方法)
  mr_res <- TwoSampleMR::mr(dat, method_list = c(
    "mr_ivw", "mr_egger_regression",
    "mr_weighted_median", "mr_weighted_mode"
  ))

  # 异质性检验
  hetero <- TwoSampleMR::mr_heterogeneity(dat)

  # 多效性检验
  pleio <- TwoSampleMR::mr_pleiotropy_test(dat)

  # Steiger 方向性
  steiger <- TwoSampleMR::directionality_test(dat)

  # MR-PRESSO
  presso <- NULL
  if (nrow(dat) >= 4) {
    tryCatch({
      presso <- MRPRESSO::mr_presso(
        BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
        SdOutcome = "se.outcome", SdExposure = "se.exposure",
        OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
        data = dat, NbDistribution = 1000, SignifThreshold = 0.05
      )
    }, error = function(e) {
      cat("  [WARN] MR-PRESSO failed:", e$message, "\n")
    })
  }

  # 提取关键统计量
  ivw <- mr_res[mr_res$method == "Inverse variance weighted", ]
  egger <- mr_res[mr_res$method == "MR Egger", ]

  summary <- data.table::data.table(
    exposure = exp_label,
    outcome = out_label,
    apoe_excluded = exclude_apoe,
    n_iv = nrow(dat),
    mean_f = mean_f,
    ivw_beta = ivw$b,
    ivw_se = ivw$se,
    ivw_pval = ivw$pval,
    egger_beta = if (nrow(egger) > 0) egger$b else NA_real_,
    egger_intercept = if (nrow(pleio) > 0) pleio$egger_intercept[1] else NA_real_,
    egger_intercept_pval = if (nrow(pleio) > 0) pleio$pval[1] else NA_real_,
    cochran_q = if (nrow(hetero) > 0) hetero$Q[hetero$method == "Inverse variance weighted"] else NA_real_,
    cochran_q_pval = if (nrow(hetero) > 0) hetero$Q_pval[hetero$method == "Inverse variance weighted"] else NA_real_,
    steiger_pval = if (nrow(steiger) > 0) steiger$pval[1] else NA_real_,
    steiger_correct_direction = if (nrow(steiger) > 0) steiger$correct_direction[1] else NA
  )

  # 添加 MR-PRESSO 信息
  if (!is.null(presso)) {
    summary$presso_global_pval <- presso$`MR-PRESSO results`$`Global Test`$Pvalue
    summary$presso_n_outliers <- ifelse(
      is.null(presso$`MR-PRESSO results`$`Outlier Test`),
      0, nrow(presso$`MR-PRESSO results`$`Outlier Test`)
    )
  }

  return(list(
    summary = summary,
    dat = dat,
    mr_res = mr_res,
    hetero = hetero,
    pleio = pleio,
    presso = presso,
    steiger = steiger
  ))
}

#' 执行 Wald Ratio (单 SNP)
#' @param exp_gwas 暴露 GWAS
#' @param out_gwas 结局 GWAS
#' @param snp SNP rsID
#' @param exp_label 暴露标签
#' @param out_label 结局标签
run_wald_ratio <- function(exp_gwas, out_gwas, snp = "rs429358",
                            exp_label = "", out_label = "") {

  exp_snp <- exp_gwas[exp_gwas$SNP == snp, ]
  out_snp <- out_gwas[out_gwas$SNP == snp, ]

  if (nrow(exp_snp) == 0 || nrow(out_snp) == 0) {
    cat(sprintf("[SKIP] %s not found in one or both GWAS\n", snp))
    return(NULL)
  }

  exp_snp$id.exposure <- exp_label
  exp_snp$exposure <- exp_label

  out_formatted <- out_snp
  colnames(out_formatted) <- gsub("exposure", "outcome", colnames(out_formatted))
  out_formatted$outcome <- out_label
  out_formatted$id.outcome <- out_label

  dat <- TwoSampleMR::harmonise_data(
    exposure_dat = exp_snp,
    outcome_dat = out_formatted,
    action = 2
  )

  stopifnot(dat$mr_keep == TRUE)

  mr_res <- TwoSampleMR::mr(dat, method_list = c("mr_wald_ratio"))

  # 手动验算
  wald_manual <- dat$beta.outcome / dat$beta.exposure

  return(list(
    summary = data.table::data.table(
      exposure = exp_label,
      outcome = out_label,
      snp = snp,
      beta_exposure = dat$beta.exposure,
      se_exposure = dat$se.exposure,
      beta_outcome = dat$beta.outcome,
      se_outcome = dat$se.outcome,
      wald_beta = mr_res$b,
      wald_se = mr_res$se,
      wald_pval = mr_res$pval,
      wald_manual = wald_manual,
      exposure_a1 = dat$effect_allele.exposure,
      outcome_a1 = dat$effect_allele.outcome,
      alleles_flipped = dat$effect_allele.exposure != dat$effect_allele.outcome,
      mr_keep = dat$mr_keep
    ),
    dat = dat,
    mr_res = mr_res
  ))
}

#' 汇总多个 MR 结果为 data.table
#' @param mr_list 命名 list of MR results
summarise_mr_list <- function(mr_list) {
  summaries <- lapply(mr_list, function(x) x$summary)
  result <- data.table::rbindlist(summaries, fill = TRUE)
  return(result)
}

cat("[mr_helpers.R] Utility functions loaded.\n")
cat("  load_local_gwas() | prepare_outcome() | verify_rs429358()\n")
cat("  ld_clump_local() | exclude_apoe()\n")
cat("  run_standard_mr() | run_wald_ratio() | summarise_mr_list()\n")
cat("  APOE region: chr19:", APOE_START, "-", APOE_END, "\n")
