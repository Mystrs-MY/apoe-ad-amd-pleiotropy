# =====================================================================
# 脚本名称: 1.LAVA_Ultimate_Upgrade_Fixed.R
# 核心功能: 全基因组局部遗传相关性扫描 (Local Genetic Correlation)
# 终极升级: 对接 Wightman/R12 终极纯净数据 + 顶刊 Any_AMD 命名规范
# =====================================================================

library(LAVA)
library(data.table)
library(foreach)
library(doParallel)
library(doSNOW)

# ================= 1. 环境与参数配置 =================
message("==================================================")
message("📡 启动 LAVA 全基因组局部遗传重叠雷达...")

# 确保路径指向你存放数据的真实文件夹
ref_prefix <- "../EUR/EUR"
locus_file <- "blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile"

# ================= 2. 读取区块雷达 =================
message("--- 第一步：加载全基因组 LD 独立区块 ---")
if(!file.exists(locus_file)) stop("❌ 找不到区块文件 (.locfile)！")
loci <- read.loci(locus_file)

# ================= 3. 创建配置表 (严格吻合 Wightman & FinnGen R12) =================
message("--- 第二步：配置多性状输入参数 (精准校正真实样本量) ---")

# FinnGen R12 总人数约 412,181 (推算 Control 约 40w)
input_info_df <- data.frame(
  phenotype = c("AD", "Dry_AMD", "Wet_AMD", "Any_AMD"),
  cases     = c(21982, 5272, 7529, 8913),
  controls  = c(41944, 406909, 404652, 403268),
  filename  = c("../GWAS/AD_Wightman_cleaned_hg19.tsv.gz",
                "../GWAS/AMD_Dry_R12_cleaned_hg19.tsv.gz",
                "../GWAS/AMD_Wet_R12_cleaned_hg19.tsv.gz",
                "../GWAS/AMD_H7_R12_cleaned_hg19.tsv.gz")
)

fwrite(input_info_df, "LAVA_input_info_R12.txt", sep="\t", quote=FALSE)

# ================= 4. 初始化 LAVA 数据对象 =================
message("--- 第三步：初始化 LAVA 对象 (底层自动映射并提取交叉位点) ---")
# 这里会自动与 1000G 欧洲人参考面板进行 LD 校准
input <- process.input(input.info.file = "LAVA_input_info_R12.txt",
                       sample.overlap.file = NULL,
                       ref.prefix = ref_prefix,
                       phenos = c("AD", "Dry_AMD", "Wet_AMD", "Any_AMD"))

# ================= 5. 多核全基因组并行扫描 =================
message("\n--- 第四步：启动多核并发扫描 (全力抓取全基因组热点) ---")

cl <- makeCluster(8) # 根据你的电脑性能可调为 6 或 10
registerDoSNOW(cl)

# 优美的进度条设置
pb <- txtProgressBar(max = nrow(loci), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# 并行轮询计算局部 rg
results <- foreach(i = 1:nrow(loci), .packages = c("LAVA", "data.table"), .options.snow = opts) %dopar% {

  locus <- process.locus(loci[i, ], input)
  res_list <- list(Dry = NULL, Wet = NULL, Any = NULL)

  if (!is.null(locus)) {
    # 尝试运行二元因果分析，tryCatch 防止没有信号的“死区块”阻断进程
    tryCatch({
      biv_out <- run.bivar(locus)
      if (!is.null(biv_out)) {

        # 追加该区块的绝对物理坐标
        loc_base <- data.frame(locus = locus$id, CHR = locus$chr, START = locus$start, STOP = locus$stop)
        biv_all <- cbind(loc_base, biv_out)

        # 提取 AD vs Dry
        dry_row <- biv_all[(biv_all$phen1 == "AD" & biv_all$phen2 == "Dry_AMD") |
                             (biv_all$phen1 == "Dry_AMD" & biv_all$phen2 == "AD"), ]
        if(nrow(dry_row) > 0) res_list$Dry <- dry_row

        # 提取 AD vs Wet
        wet_row <- biv_all[(biv_all$phen1 == "AD" & biv_all$phen2 == "Wet_AMD") |
                             (biv_all$phen1 == "Wet_AMD" & biv_all$phen2 == "AD"), ]
        if(nrow(wet_row) > 0) res_list$Wet <- wet_row

        # 提取 AD vs Any (原 H7)
        any_row <- biv_all[(biv_all$phen1 == "AD" & biv_all$phen2 == "Any_AMD") |
                             (biv_all$phen1 == "Any_AMD" & biv_all$phen2 == "AD"), ]
        if(nrow(any_row) > 0) res_list$Any <- any_row
      }
    }, error = function(e) {})
  }

  return(res_list)
}

close(pb)
stopCluster(cl)

# ================= 6. 合并导出并智能检索 APOE 目标 =================
message("\n--- 第五步：汇总落盘并定向捕捉 Chr19 的终极密码 ---")

df_dry <- rbindlist(lapply(results, function(x) x$Dry))
df_wet <- rbindlist(lapply(results, function(x) x$Wet))
df_any <- rbindlist(lapply(results, function(x) x$Any))

# 按 P 值排序，确保最强信号在顶端
df_dry <- df_dry[order(p)]
df_wet <- df_wet[order(p)]
df_any <- df_any[order(p)]

# 保存终极结果
fwrite(df_dry, "LAVA_FullScan_AD_vs_DryAMD_Final.csv")
fwrite(df_wet, "LAVA_FullScan_AD_vs_WetAMD_Final.csv")
fwrite(df_any, "LAVA_FullScan_AD_vs_AnyAMD_Final.csv")

# 自动检索 19 号染色体上的最强信号 (APOE 定位器)
find_apoe_block <- function(res_df, target_name) {
  top_block <- res_df[CHR == 19][order(p)][1, ]
  message(sprintf(">> [%s] 最强区块 Locus %d (P = %.2e) | 局部效应 rho = %.3f | 物理坐标: %d - %d",
                  target_name, top_block$locus, top_block$p, top_block$rho, top_block$START, top_block$STOP))
}

message("\n🚨🚨🚨 锁定以下坐标！这将是你论文中选定 APOE 区域的绝对理由：🚨🚨🚨")
find_apoe_block(df_dry, "Dry AMD")
find_apoe_block(df_wet, "Wet AMD")
find_apoe_block(df_any, "Any AMD")

message("\n🎉 LAVA 扫描全部收官！看看是不是所有的 rho 都是负数(完美拮抗)，且 P 值爆表！")