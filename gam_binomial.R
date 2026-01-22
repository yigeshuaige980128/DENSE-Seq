#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(mgcv)
  library(ggplot2)
  library(foreach)
  library(doParallel)
  library(tools)
})

# ============================================================
# 用法：
#   Rscript fit_pool_ctrl_corrected_af_manhattan_extrema.R combo_paths.txt [out_dir] [n_cores]
#
# combo_paths.txt: 每行一个输入文件路径（允许行内带多列，默认只取第一列）
#
# 每个 combo 输出：
#   1) POOL_corrected_by_CTRL.AF_fit.12chr.pdf
#   2) POOL_vs_CTRL.manhattan_-log10p.12chr.pdf
#   3) Chr*/fitted_values.ChrX.txt.gz
#   4) fit_info.txt
#   5) extrema_top5_max_min.perChr.correctedAF.txt   (每样本每Chr前5极大+前5极小)
#
# 全局输出（out_root）：
#   ALL_samples.extrema_top5_max_min.perChr.correctedAF.txt
# ============================================================

# --------------------------
# CLI 参数
# --------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop(paste0(
    "用法：Rscript fit_pool_ctrl_corrected_af_manhattan_extrema.R combo_paths.txt [out_dir] [n_cores]\n",
    "  combo_paths.txt: 每行一个输入文件路径\n",
    "  out_dir: 输出目录（默认 af_bam_output）\n",
    "  n_cores: 并行处理文件数（默认 4）\n"
  ))
}
paths_list_file <- args[1]
out_root <- ifelse(length(args) >= 2, args[2], "af_bam_output")
n_cores  <- ifelse(length(args) >= 3, as.integer(args[3]), 4)

if (!dir.exists(out_root)) dir.create(out_root, recursive = TRUE)

# --------------------------
# 读取路径清单（防呆：只取第一列 + 跳过不存在）
# --------------------------
combo_paths <- readLines(paths_list_file, warn = FALSE)
combo_paths <- trimws(combo_paths)
combo_paths <- combo_paths[nzchar(combo_paths)]
combo_paths <- sub("[[:space:]\t].*$", "", combo_paths)   # 只取第一列

ok <- file.exists(combo_paths)
if (!all(ok)) {
  warning("以下路径不存在，将被跳过：\n", paste(combo_paths[!ok], collapse = "\n"))
}
combo_paths <- combo_paths[ok]
if (length(combo_paths) == 0) stop("路径清单里没有任何有效文件路径：", paths_list_file)

# --------------------------
# 染色体长度（bp）
# --------------------------
chr_len_bp <- c(
  Chr1=43270923, Chr2=35937250, Chr3=36413819, Chr4=35502694,
  Chr5=29958434, Chr6=31248787, Chr7=29697621, Chr8=28443022,
  Chr9=23012720, Chr10=23207287, Chr11=29021106, Chr12=27531856
)
keep_chr <- names(chr_len_bp)

# --------------------------
# 参数（可调）
# --------------------------
pool_success_default <- "POOL_T65.np7"
pool_failure_default <- "POOL_GLA4.np7"
ctrl_success_default <- "CTRL_T65.np7"
ctrl_failure_default <- "CTRL_GLA4.np7"

min_depth <- 15
target_bw_mb <- 0.5
k_cap <- 200
grid_n <- 4000

# 两张图统一 12chr 3x4 版式
pdf_width  <- 16
pdf_height <- 12
facet_ncol <- 3
pdf_dpi    <- 300

# 曼哈顿阈值 & 抽样
p_thresh <- 1e-6
max_points_manhattan <- 2e6

# 图1拟合线颜色（深蓝）
fit_line_color <- "darkblue"

# extrema：每Chr取前5 max/min
extrema_topN <- 5

# --------------------------
# 工具函数
# --------------------------
parse_combo_pool <- function(fpath) {
  base <- basename(fpath)
  combo <- sub("\\.[^.]+$", "", base)
  pool  <- sub("\\.minus.*$", "", combo, perl = TRUE)
  list(combo = combo, pool = pool)
}

choose_k <- function(chr, chr_len_bp, target_bw_mb = 0.5, k_cap = 200) {
  L <- chr_len_bp[[chr]]
  if (is.na(L) || is.null(L)) return(80L)
  k0 <- ceiling((L / 1e6) / target_bw_mb) + 10
  k0 <- max(40L, min(as.integer(k0), as.integer(k_cap)))
  k0
}

.sign05 <- function(x, eps = 1e-12) {
  ifelse(x > 0.5 + eps,  1L, ifelse(x < 0.5 - eps, -1L, 0L))
}

# 0.5基线方向校正（你的规则）
correct_pool_by_ctrl <- function(pool, ctrl) {
  dirP <- .sign05(pool)
  dirC <- .sign05(ctrl)
  d    <- abs(ctrl - 0.5)
  corr <- pool + ifelse(dirP == 0L, 0,
                        ifelse(dirP == dirC, -dirP * d, dirP * d))
  pmin(1, pmax(0, corr))
}

# 从曲线 y(x) 提取局部极值（max/min），各取 topN
find_extrema_topN <- function(x, y, topN = 5, eps = 1e-10) {
  if (length(x) < 5 || length(y) < 5) {
    return(data.frame(type=character(), pos_mb=numeric(), value=numeric()))
  }
  ord <- order(x); x <- x[ord]; y <- y[ord]

  dy <- diff(y)
  s  <- ifelse(abs(dy) < eps, 0L, sign(dy))

  # forward-fill 0（平台）
  for (i in seq_along(s)) {
    if (s[i] == 0L && i > 1) s[i] <- s[i-1]
  }
  # backward-fill 头部仍为0
  if (length(s) > 0 && s[1] == 0L) {
    j <- which(s != 0L)
    if (length(j) > 0) s[1:j[1]] <- s[j[1]]
  }

  max_idx <- which(s[-length(s)] > 0 & s[-1] < 0) + 1
  min_idx <- which(s[-length(s)] < 0 & s[-1] > 0) + 1

  ext <- rbind(
    data.frame(type="max", pos_mb=x[max_idx], value=y[max_idx]),
    data.frame(type="min", pos_mb=x[min_idx], value=y[min_idx])
  )
  if (nrow(ext) == 0) return(ext)

  ext_max <- ext[ext$type=="max", , drop=FALSE]
  ext_min <- ext[ext$type=="min", , drop=FALSE]

  if (nrow(ext_max) > 0) ext_max <- ext_max[order(-ext_max$value), ][seq_len(min(topN, nrow(ext_max))), , drop=FALSE]
  if (nrow(ext_min) > 0) ext_min <- ext_min[order( ext_min$value), ][seq_len(min(topN, nrow(ext_min))), , drop=FALSE]

  rbind(ext_max, ext_min)
}

# --------------------------
# 图1：12chr合一，仅画“校正后拟合线”（无均值标注），线为深蓝
# --------------------------
plot_corrected_af_12chr <- function(pred_corr, combo, out_pdf,
                                    facet_ncol = 3,
                                    width = 16, height = 12, dpi = 300,
                                    line_color = "darkblue") {
  pred_corr$Chr <- factor(pred_corr$Chr, levels = paste0("Chr", 1:12))

  p <- ggplot(pred_corr, aes(x = pos_mb, y = af_corr_fit)) +
    geom_hline(yintercept = 0.5, linetype = "11", linewidth = 0.4) +
    geom_line(linewidth = 0.7, color = line_color) +
    facet_wrap(~ Chr, ncol = facet_ncol, scales = "free_x") +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      title = paste0(combo, " | corrected POOL AF (CTRL-based, baseline=0.5)"),
      x = "Genomic Position (Mb)",
      y = "Corrected AF (0–1)"
    ) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "white"),
      strip.text = element_text(face = "bold", size = 10),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
    )

  ggsave(out_pdf, p, width = width, height = height, units = "in",
         dpi = dpi, limitsize = FALSE)
}

# --------------------------
# 图2：曼哈顿（12chr facet 3x4，每Chr一个小框）
# two-proportion z test, y = -log10(p)
# --------------------------
plot_manhattan_pool_vs_ctrl_12chr <- function(dt, combo, out_pdf,
                                             p_thresh = 1e-6, max_points = 2e6,
                                             facet_ncol = 3,
                                             width = 16, height = 12, dpi = 300) {
  keep_chr <- paste0("Chr", 1:12)
  dt <- dt[Chr %in% keep_chr]
  dt[, Chr := factor(Chr, levels = keep_chr)]

  # two-proportion z-test（向量化）
  n1 <- dt$pool_depth
  n2 <- dt$ctrl_depth
  p1 <- dt$pool_s / n1
  p2 <- dt$ctrl_s / n2
  p  <- (dt$pool_s + dt$ctrl_s) / (n1 + n2)

  se <- sqrt(p * (1 - p) * (1/n1 + 1/n2))
  se[!is.finite(se) | se <= 0] <- NA_real_

  z <- (p1 - p2) / se
  pval  <- 2 * pnorm(-abs(z))
  mlogp <- -log10(pval + 1e-300)

  dt[, `:=`(pval = pval, mlogp = mlogp, x_mb = POS / 1e6)]

  # 点太多就抽样
  if (nrow(dt) > max_points) {
    set.seed(1)
    dt <- dt[sample.int(nrow(dt), max_points)]
  }

  thr_y <- -log10(p_thresh)

  pplot <- ggplot(dt, aes(x = x_mb, y = mlogp)) +
    geom_point(size = 0.25, alpha = 0.65, color = "grey35", na.rm = TRUE) +
    geom_hline(yintercept = thr_y, linetype = "11", linewidth = 0.45) +
    facet_wrap(~ Chr, ncol = facet_ncol, scales = "free_x") +
    labs(
      title = paste0(combo, " | Manhattan: POOL vs CTRL (two-proportion z test)"),
      x = "Genomic Position (Mb)",
      y = expression(-log[10](p))
    ) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "white"),
      strip.text = element_text(face = "bold", size = 10),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
    )

  ggsave(out_pdf, pplot, width = width, height = height, units = "in",
         dpi = dpi, limitsize = FALSE)
}

# ============================================================
# 单文件处理
# ============================================================
process_one_file <- function(fpath, out_root) {
  if (!file.exists(fpath)) {
    warning("文件不存在，跳过：", fpath)
    return(NULL)
  }

  meta <- parse_combo_pool(fpath)
  combo <- meta$combo
  pool  <- meta$pool

  sample_dir <- file.path(out_root, combo)
  if (!dir.exists(sample_dir)) dir.create(sample_dir, recursive = TRUE)

  dt <- tryCatch(
    fread(fpath, sep = "\t", header = TRUE, data.table = TRUE, showProgress = FALSE),
    error = function(e) NULL
  )
  if (is.null(dt) || ncol(dt) < 6) {
    warning("无法读取或列数不足：", fpath)
    return(NULL)
  }

  # 默认：第1列Chr，第2列POS
  chr_col <- names(dt)[1]
  pos_col <- names(dt)[2]

  nms <- names(dt)
  tail4 <- tail(nms, 4)

  ctrl_success_col <- if (ctrl_success_default %in% nms) ctrl_success_default else tail4[1]
  ctrl_failure_col <- if (ctrl_failure_default %in% nms) ctrl_failure_default else tail4[2]
  pool_success_col <- if (pool_success_default %in% nms) pool_success_default else tail4[3]
  pool_failure_col <- if (pool_failure_default %in% nms) pool_failure_default else tail4[4]

  # 转数值 + 清洗
  dt[, (pos_col) := suppressWarnings(as.numeric(get(pos_col)))]
  dt[, ctrl_s := suppressWarnings(as.numeric(get(ctrl_success_col)))]
  dt[, ctrl_f := suppressWarnings(as.numeric(get(ctrl_failure_col)))]
  dt[, pool_s := suppressWarnings(as.numeric(get(pool_success_col)))]
  dt[, pool_f := suppressWarnings(as.numeric(get(pool_failure_col)))]

  dt <- dt[is.finite(get(pos_col)) & is.finite(ctrl_s) & is.finite(ctrl_f) & is.finite(pool_s) & is.finite(pool_f)]

  dt[, ctrl_depth := ctrl_s + ctrl_f]
  dt[, pool_depth := pool_s + pool_f]
  dt <- dt[ctrl_depth >= min_depth & pool_depth >= min_depth]
  if (nrow(dt) < 1000) {
    warning("有效点太少（<1000）跳过：", combo)
    return(NULL)
  }

  dt[, CTRL_af_obs := ctrl_s / ctrl_depth]
  dt[, POOL_af_obs := pool_s / pool_depth]
  dt[, Chr := as.character(get(chr_col))]
  dt[, POS := as.integer(round(get(pos_col)))]
  dt[, pos_mb := POS / 1e6]

  dt <- dt[Chr %in% keep_chr]
  if (nrow(dt) < 1000) {
    warning("过滤Chr1-12后有效点太少：", combo)
    return(NULL)
  }
  dt[, Chr := factor(Chr, levels = keep_chr)]

  # 图1：收集“校正后的拟合线”
  line_corr_list <- list()

  # fit_info
  fit_info_list <- list()

  # extrema（每Chr前5 max/min）
  extrema_list <- list()

  for (cc in levels(dt$Chr)) {
    dtc <- dt[Chr == cc][order(pos_mb)]
    if (nrow(dtc) < 500) next

    chr_dir <- file.path(sample_dir, as.character(cc))
    if (!dir.exists(chr_dir)) dir.create(chr_dir, recursive = TRUE)

    k <- choose_k(as.character(cc), chr_len_bp, target_bw_mb, k_cap)

    m_pool <- tryCatch(
      bam(cbind(pool_s, pool_f) ~ s(pos_mb, bs = "ps", k = k),
          family = binomial(), data = dtc,
          method = "fREML", discrete = TRUE,
          nthreads = 1, chunk.size = 100000),
      error = function(e) NULL
    )

    m_ctrl <- tryCatch(
      bam(cbind(ctrl_s, ctrl_f) ~ s(pos_mb, bs = "ps", k = k),
          family = binomial(), data = dtc,
          method = "fREML", discrete = TRUE,
          nthreads = 1, chunk.size = 100000),
      error = function(e) NULL
    )

    if (is.null(m_pool) || is.null(m_ctrl)) {
      warning("拟合失败：", combo, " ", cc)
      next
    }

    # SNP点 fitted（用于写文件）
    pool_fit <- as.numeric(predict(m_pool, type = "response"))
    ctrl_fit <- as.numeric(predict(m_ctrl, type = "response"))
    corr_fit <- correct_pool_by_ctrl(pool_fit, ctrl_fit)

    dtc[, POOL_af_fit := pool_fit]
    dtc[, CTRL_af_fit := ctrl_fit]
    dtc[, POOL_af_corr_fit := corr_fit]

    # 网格预测（用于图1 + extrema）
    xg <- seq(min(dtc$pos_mb), max(dtc$pos_mb), length.out = grid_n)
    pg_pool <- as.numeric(predict(m_pool, newdata = data.frame(pos_mb = xg), type = "response"))
    pg_ctrl <- as.numeric(predict(m_ctrl, newdata = data.frame(pos_mb = xg), type = "response"))
    pg_corr <- correct_pool_by_ctrl(pg_pool, pg_ctrl)

    # 图1拟合线收集
    line_corr_list[[as.character(cc)]] <- data.frame(
      Chr = as.character(cc),
      pos_mb = xg,
      af_corr_fit = pg_corr
    )

    # extrema：在 pg_corr 上找局部极值，取前5 max/min
    ext_cc <- find_extrema_topN(xg, pg_corr, topN = extrema_topN)
    if (nrow(ext_cc) > 0) {
      ext_cc$Combo <- combo
      ext_cc$Pool  <- pool
      ext_cc$Chr   <- as.character(cc)
      ext_cc$pos_bp <- round(ext_cc$pos_mb * 1e6)
      ext_cc$metric <- "corrected_AF_fit"
      ext_cc <- ext_cc[, c("Combo","Pool","Chr","type","pos_mb","pos_bp","value","metric")]
      extrema_list[[as.character(cc)]] <- ext_cc
    }

    # fit_info
    smp <- summary(m_pool)
    smc <- summary(m_ctrl)
    fit_info_list[[as.character(cc)]] <- data.frame(
      Chr = as.character(cc),
      N_points = nrow(dtc),
      k = k,
      target_bw_mb = target_bw_mb,
      AIC_pool = AIC(m_pool),
      AIC_ctrl = AIC(m_ctrl),
      edf_pool = sum(smp$edf),
      edf_ctrl = sum(smc$edf),
      dev_expl_pool = smp$dev.expl,
      dev_expl_ctrl = smc$dev.expl,
      Method = "bam_binomial_ps (pool & ctrl) + corrected around 0.5",
      stringsAsFactors = FALSE
    )

    # 写每条Chr点级输出
    fitted_path <- file.path(chr_dir, sprintf("fitted_values.%s.txt.gz", as.character(cc)))
    con <- gzfile(fitted_path, "wt")
    writeLines(paste(
      c("Chr","POS",
        "POOL_success","POOL_failure","POOL_depth","POOL_af_obs","POOL_af_fit",
        "CTRL_success","CTRL_failure","CTRL_depth","CTRL_af_obs","CTRL_af_fit",
        "POOL_af_corr_fit"),
      collapse = "\t"
    ), con = con)

    write.table(
      data.frame(
        Chr = as.character(dtc$Chr),
        POS = dtc$POS,
        POOL_success = dtc$pool_s,
        POOL_failure = dtc$pool_f,
        POOL_depth   = dtc$pool_depth,
        POOL_af_obs  = dtc$POOL_af_obs,
        POOL_af_fit  = dtc$POOL_af_fit,
        CTRL_success = dtc$ctrl_s,
        CTRL_failure = dtc$ctrl_f,
        CTRL_depth   = dtc$ctrl_depth,
        CTRL_af_obs  = dtc$CTRL_af_obs,
        CTRL_af_fit  = dtc$CTRL_af_fit,
        POOL_af_corr_fit = dtc$POOL_af_corr_fit
      ),
      file = con, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
    )
    close(con)
  }

  # 写 fit_info
  fit_info <- if (length(fit_info_list) > 0) do.call(rbind, fit_info_list) else NULL
  if (!is.null(fit_info) && nrow(fit_info) > 0) {
    write.table(fit_info, file.path(sample_dir, "fit_info.txt"),
                sep = "\t", row.names = FALSE, quote = FALSE)
  } else {
    file.create(file.path(sample_dir, "fit_info.txt"))
  }

  # 图1：12chr合一（仅拟合线，深蓝；无均值标注）
  if (length(line_corr_list) > 0) {
    pred_corr <- do.call(rbind, line_corr_list)
    out_pdf1 <- file.path(sample_dir, "POOL_corrected_by_CTRL.AF_fit.12chr.pdf")
    plot_corrected_af_12chr(
      pred_corr = pred_corr,
      combo = combo,
      out_pdf = out_pdf1,
      facet_ncol = facet_ncol,
      width = pdf_width,
      height = pdf_height,
      dpi = pdf_dpi,
      line_color = fit_line_color
    )
  } else {
    warning("没有可用的拟合线，跳过图1：", combo)
  }

  # 图2：曼哈顿（12chr facet 3x4）
  out_pdf2 <- file.path(sample_dir, "POOL_vs_CTRL.manhattan_-log10p.12chr.pdf")
  plot_manhattan_pool_vs_ctrl_12chr(
    dt = dt[, .(Chr=as.character(Chr), POS=POS,
                pool_s=pool_s, pool_depth=pool_depth,
                ctrl_s=ctrl_s, ctrl_depth=ctrl_depth)],
    combo = combo,
    out_pdf = out_pdf2,
    p_thresh = p_thresh,
    max_points = max_points_manhattan,
    facet_ncol = facet_ncol,
    width = pdf_width,
    height = pdf_height,
    dpi = pdf_dpi
  )

  # 每样本 extrema 输出 + 返回给全局汇总
  extrema_df <- if (length(extrema_list) > 0) rbindlist(extrema_list, use.names = TRUE, fill = TRUE) else NULL
  out_ext <- file.path(sample_dir, "extrema_top5_max_min.perChr.correctedAF.txt")
  if (!is.null(extrema_df) && nrow(extrema_df) > 0) {
    fwrite(extrema_df, out_ext, sep = "\t")
  } else {
    file.create(out_ext)
  }

  return(extrema_df)
}

# ============================================================
# 并行：按文件并行
# ============================================================
cl <- makeCluster(n_cores)
registerDoParallel(cl)

res_list <- foreach(f = combo_paths,
                    .packages = c("data.table","mgcv","ggplot2","foreach","doParallel","tools"),
                    .errorhandling = "pass") %dopar% {
  tryCatch(process_one_file(f, out_root), error = function(e) NULL)
}

stopCluster(cl)

# ============================================================
# 全局 extrema 汇总
# ============================================================
all_ext <- rbindlist(res_list, use.names = TRUE, fill = TRUE)
out_all <- file.path(out_root, "ALL_samples.extrema_top5_max_min.perChr.correctedAF.txt")
if (!is.null(all_ext) && nrow(all_ext) > 0) {
  fwrite(all_ext, out_all, sep = "\t")
} else {
  file.create(out_all)
}

cat("完成：输出根目录 => ", normalizePath(out_root), "\n",
    "每个 combo 子目录：\n",
    "  - POOL_corrected_by_CTRL.AF_fit.12chr.pdf\n",
    "  - POOL_vs_CTRL.manhattan_-log10p.12chr.pdf\n",
    "  - extrema_top5_max_min.perChr.correctedAF.txt\n",
    "  - Chr1-12 子目录下：fitted_values.<Chr>.txt.gz\n",
    "  - fit_info.txt\n",
    "全局 extrema 汇总：\n",
    "  - ", out_all, "\n", sep = "")
