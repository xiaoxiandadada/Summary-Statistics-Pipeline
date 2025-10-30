#!/usr/bin/env Rscript

# Pipeline 主分析流程
# GhostKnockoff 生成与 LAVA-Knock 局部遗传相关分析

# 加载所需 R 包
pipeline_packages <- c(
  "corpcor", "dplyr", "readr",
  "SKAT", "SPAtest", "CompQuadForm",
  "GhostKnockoff", "irlba", "matrixsampling",
  "LAVAKnock", "snpStats",
  "parallel"
)
missing_packages <- pipeline_packages[
  !vapply(pipeline_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(sprintf("缺少依赖: %s。请先运行 install_packages.R", paste(missing_packages, collapse = ", ")))
}

suppressPackageStartupMessages(invisible(lapply(pipeline_packages, library, character.only = TRUE)))

# 源配置脚本
source("configure_parameters.R")

# Python 加速（直接从 accelerator.py 载入，可选）
py_env <- new.env(parent = emptyenv())
py_ok <- FALSE
if (requireNamespace("reticulate", quietly = TRUE)) {
  try({ reticulate::source_python("accelerator.py", envir = py_env); py_ok <- TRUE }, silent = TRUE)
}
if (!isTRUE(py_ok)) {
  stop('Python 加速模块加载失败，请确认已安装 reticulate 并能导入 accelerator.py')
}

py_fast_correlation <- function(X, verbose = FALSE) {
  reticulate::py_to_r(py_env$fast_correlation_matrix(X))
}

py_ld_pruning <- function(corr, threshold = 0.75, verbose = FALSE) {
  idx0 <- reticulate::py_to_r(py_env$ld_pruning(corr, threshold))
  as.integer(idx0) + 1L
}

py_qc_genotype_with_mask <- function(genotype_matrix, maf_threshold = 0.0, mac_threshold = 25.0) {
  reticulate::py_to_r(py_env$qc_genotype_with_mask(genotype_matrix, maf_threshold = maf_threshold, mac_threshold = mac_threshold))
}

py_auto_prepare_inputs <- function(zscore, panel_default, trait_arg = NULL, verbose = FALSE) {
  reticulate::py_to_r(py_env$auto_prepare_inputs(zscore, panel_default, trait_arg, verbose))
}

py_load_variant_info <- function(path) {
  reticulate::py_to_r(py_env$load_variant_info(path))
}

py_load_gwas_table <- function(path) {
  reticulate::py_to_r(py_env$load_gwas_table(path))
}

py_load_multi_gwas_table <- function(path, zcols) {
  reticulate::py_to_r(py_env$load_multi_gwas_table(path, zcols))
}

py_load_genotype_csv <- function(path, fmt = 'samples_by_snps') {
  reticulate::py_to_r(py_env$load_genotype_csv(path, fmt = fmt))
}

py_load_gene_catalog <- function(coord, base_dir = '.') {
  reticulate::py_to_r(py_env$load_gene_catalog(coord, base_dir))
}

py_load_genotype_plink <- function(prefix, snp_ids = NULL, chrpos_ids = NULL) {
  reticulate::py_to_r(py_env$load_genotype_plink(prefix, snp_ids, chrpos_ids))
}

py_align_genotype_to_gwas <- function(geno_matrix, colnames_geno, chrpos_attr, rsid_attr, gwas_df, prefer_rsid = TRUE) {
  reticulate::py_to_r(py_env$align_genotype_to_gwas(
    geno_matrix,
    colnames_geno,
    chrpos_attr,
    rsid_attr,
    gwas_df,
    prefer_rsid
  ))
}

load_gene_catalog <- function(coord, base_dir = '.') {
  genes <- py_load_gene_catalog(coord, base_dir)
  df <- as.data.frame(genes, stringsAsFactors = FALSE)
  col_lower <- tolower(colnames(df))
  rename_map <- c(
    id = "id",
    chr = "chr",
    start = "start",
    end = "end",
    stop = "end",
    stop_bp = "end",
    pos_end = "end",
    position_end = "end",
    tss = "tss"
  )
  for (nm in names(rename_map)) {
    idx <- which(col_lower == nm)
    if (length(idx) == 1) {
      colnames(df)[idx] <- rename_map[[nm]]
      col_lower[idx] <- rename_map[[nm]]
    }
  }
  if (!"id" %in% colnames(df)) {
    df$id <- paste0("gene_", seq_len(nrow(df)))
  }
  if (!"end" %in% colnames(df)) {
    stop("gene catalog 缺少 end 列")
  }
  df
}

partition_ld_blocks <- function(chr_vec, pos_vec, coord, workers = NULL, fallback = TRUE) {
  args <- list(
    chromosomes = as.integer(chr_vec),
    positions = as.numeric(pos_vec),
    coord_version = coord,
    fallback = fallback
  )
  if (!is.null(workers)) {
    args$workers <- as.integer(workers)
  }
  py_res <- do.call(py_env$partition_ld_blocks, args)
  reticulate::py_to_r(py_res)
}

annotate_nearest_gene <- function(chrs, positions, gene_catalog) {
  if (is.null(gene_catalog) || nrow(gene_catalog) == 0) {
    return(rep(NA_character_, length(positions)))
  }
  py_gene <- reticulate::r_to_py(gene_catalog)
  reticulate::py_to_r(py_env$annotate_nearest_gene(as.integer(chrs), as.numeric(positions), py_gene))
}

load_genotype_plink <- function(prefix, snp_ids = NULL, chrpos_ids = NULL) {
  res <- py_load_genotype_plink(prefix, snp_ids, chrpos_ids)
  geno_matrix <- res$matrix
  if (!is.null(res$colnames)) colnames(geno_matrix) <- res$colnames
  if (!is.null(res$chrpos)) attr(geno_matrix, 'chrpos') <- res$chrpos
  if (!is.null(res$rsid)) attr(geno_matrix, 'rsid') <- res$rsid
  map_df <- if (!is.null(res$map)) as.data.frame(res$map, stringsAsFactors = FALSE) else NULL
  list(matrix = geno_matrix, map = map_df)
}

align_genotype_to_gwas <- function(geno_matrix, gwas_df, prefer_rsid = TRUE) {
  res <- py_align_genotype_to_gwas(
    geno_matrix,
    colnames(geno_matrix),
    attr(geno_matrix, 'chrpos'),
    attr(geno_matrix, 'rsid'),
    gwas_df,
    prefer_rsid
  )
  geno_selected <- res$matrix
  if (!is.null(res$variant_ids)) {
    colnames(geno_selected) <- res$variant_ids
    attr(geno_selected, 'chrpos') <- res$variant_ids
  }
  if (!is.null(res$rsid)) {
    attr(geno_selected, 'rsid') <- res$rsid
  }
  list(
    matrix = geno_selected,
    gwas_indices = as.integer(res$gwas_indices),
    matched = as.integer(res$matched),
    variant_ids = res$variant_ids,
    positions = as.numeric(res$positions),
    rsid = res$rsid
  )
}

parse_trait_option <- function(x) {
  if (is.null(x) || is.na(x) || x == '') return(NULL)
  trimws(unlist(strsplit(x, ',')))
}

infer_gene_catalog_key <- function(ld_coord_arg, default = "GRCh37") {
  if (is.null(ld_coord_arg) || !nzchar(ld_coord_arg)) return(default)
  upper_val <- toupper(ld_coord_arg)
  if (upper_val %in% c("GRCH37", "GRCH38")) return(upper_val)
  file_name <- basename(ld_coord_arg)
  if (grepl("38|HG38|GRCH38", file_name, ignore.case = TRUE)) return("GRCh38")
  if (grepl("37|HG19|GRCH37", file_name, ignore.case = TRUE)) return("GRCh37")
  default
}

load_genotype_rds <- function(path) {
  obj <- readRDS(path)
  if (!is.matrix(obj)) stop('RDS文件未包含matrix对象')
  storage.mode(obj) <- 'numeric'
  obj
}

# 窗口管理函数 — 构建窗口区间

#' 创建分析窗口
create_analysis_windows <- function(variant_info, window_size = 100000, window_step = 50000) {
  if (nrow(variant_info) == 0) {
    return(data.frame())
  }

  chr <- unique(variant_info$chr)[1]

  # 检查pos_bp列是否存在
  if ("pos_bp" %in% colnames(variant_info)) {
    positions <- variant_info$pos_bp
  } else if ("pos" %in% colnames(variant_info)) {
    positions <- variant_info$pos
  } else if ("POS" %in% colnames(variant_info)) {
    positions <- variant_info$POS
  } else {
    warning("未找到位置信息列 (pos_bp, pos, POS)")
    return(data.frame())
  }

  # 移除NA值
  valid_indices <- !is.na(positions)
  if (sum(valid_indices) == 0) {
    warning("所有位置信息都是NA")
    return(data.frame())
  }

  positions <- positions[valid_indices]
  min_pos <- min(positions)
  max_pos <- max(positions)
  windows <- data.frame()
  window_id <- 1

  start_pos <- min_pos
  while (start_pos < max_pos) {
    end_pos <- start_pos + window_size - 1

    # 检查窗口内是否有变异
    variants_in_window <- which(positions >= start_pos & positions <= end_pos)

    if (length(variants_in_window) >= 5) {  # 至少5个变异
      windows <- rbind(windows, data.frame(
        window_id = window_id,
        chr = chr,
        start = start_pos,
        end = end_pos,
        n_variants = length(variants_in_window),
        variant_indices = I(list(variants_in_window))
      ))
      window_id <- window_id + 1
    }
    start_pos <- start_pos + window_step
  }

  return(windows)
}

# GhostKnockoff 相关函数 — knockoff 分数生成

#' 使用GhostKnockoff生成knockoff Z-scores（可接入真实基因型/LD）
#' @param genotype_matrix 数值矩阵，维度: 样本×SNP；列名建议为 chr:pos 或 rsid（并与zscore对齐）
#' @param zscore_matrix 数值矩阵，维度: SNP×表型；行顺序需与基因型列顺序一致
#' @param ld_threshold 数值，LD聚类阈值（0-1），用于去冗余
#' @param n_samples 整数，GWAS样本量
#' @param n_knockoffs 整数，knockoff副本数
#' @param variant_positions 可选，数值向量，对应每个SNP的基因组位置（优先用于排序）
#' @param verbose 逻辑，是否打印日志
generate_ghostknockoff_scores <- function(genotype_matrix, zscore_matrix,
                                          ld_threshold = 0.75, n_samples = 20000,
                                          n_knockoffs = 5, variant_positions = NULL,
                                          verbose = FALSE) {

  if (verbose) cat("开始生成GhostKnockoff分数...\n")

  # 数据预处理与 QC（由 Python 侧完成）
  orig_colnames <- colnames(genotype_matrix)
  qc <- py_qc_genotype_with_mask(genotype_matrix)
  keep_index <- as.integer(qc$indices) + 1L
  genotype_matrix <- qc$matrix
  if (!is.null(orig_colnames)) {
    colnames(genotype_matrix) <- orig_colnames[keep_index]
  }
  zscore_matrix <- zscore_matrix[keep_index, , drop = FALSE]
  if (!is.null(variant_positions)) {
    variant_positions <- as.numeric(variant_positions)
    variant_positions <- variant_positions[keep_index]
  }

  if (ncol(genotype_matrix) < 5) {
    stop("过滤后变异数量过少")
  }

  # 按位置排序（优先使用提供的位置向量）
  if (!is.null(variant_positions)) {
    pos <- as.numeric(variant_positions)
  } else {
    pos <- suppressWarnings(as.numeric(gsub("^.*:", "", colnames(genotype_matrix))))
  }
  order_idx <- order(pos, na.last = TRUE)
  genotype_matrix <- genotype_matrix[, order_idx, drop = FALSE]
  zscore_matrix <- zscore_matrix[order_idx, , drop = FALSE]
  keep_index <- keep_index[order_idx]
  variant_positions <- pos[order_idx]

  # LD过滤（Python实现）
  if (ncol(genotype_matrix) > 1) {
    corr_matrix <- py_fast_correlation(genotype_matrix, verbose)
    keep_cols <- py_ld_pruning(corr_matrix, threshold = ld_threshold, verbose = verbose)
    genotype_matrix <- genotype_matrix[, keep_cols, drop = FALSE]
    zscore_matrix <- zscore_matrix[keep_cols, , drop = FALSE]
    keep_index <- keep_index[keep_cols]
  }

  if (verbose) cat(paste("LD过滤后变异数量:", ncol(genotype_matrix), "\n"))

  # 计算收缩LD矩阵
corr_shrink <- corpcor::cor.shrink(genotype_matrix, verbose = FALSE)

  # 生成knockoff
  tryCatch({
    set.seed(12345)
    fit_prelim <- GhostKnockoff::GhostKnockoff.prelim(corr_shrink, M = n_knockoffs,
                                      method = 'asdp',
                                      max.size = ncol(genotype_matrix))

    # 为每个表型生成knockoff
    knockoff_results <- list()

    for (pheno in seq_len(ncol(zscore_matrix))) {
      gk_fit <- GhostKnockoff::GhostKnockoff.fit(as.matrix(zscore_matrix[, pheno]),
                                  n.study = n_samples,
                                  fit.prelim = fit_prelim,
                                  gamma = 1,
                                  weight.study = NULL)

      # 组合原始和knockoff Z-scores
      combined_scores <- cbind(zscore_matrix[, pheno], gk_fit$GK.Zscore_k)
      colnames(combined_scores) <- c("org", paste0("knock", 1:n_knockoffs))

      knockoff_results[[paste0("pheno", pheno)]] <- combined_scores
    }

    knockoff_results$index <- keep_index
    if (verbose) cat("GhostKnockoff分数生成完成\n")
    return(knockoff_results)

  }, error = function(e) {
    warning(paste("GhostKnockoff生成失败:", e$message))
    return(NULL)
  })
}

# LAVA-Knock 分析函数 — 单/双变量窗口统计

#' 单变量遗传性检验
run_univariate_analysis <- function(genotype_matrix, zscore_matrix,
                                   chr, locus_start, locus_end,
                                   n_samples = 20000, prune_thresh = 99,
                                   verbose = FALSE) {

  if (verbose) cat("运行单变量遗传性分析...\n")

  tryCatch({
    result <- LAVAKnock::LAVAKnock_univariate(
      Zscore = zscore_matrix,
      G_locus = genotype_matrix,
      chr = chr,
      locus_start = locus_start,
      locus_end = locus_end,
      n = n_samples,
      prune.thresh = prune_thresh
    )

    if (verbose) cat("单变量分析完成\n")
    return(result)

  }, error = function(e) {
    warning(paste("单变量分析失败:", e$message))
    return(NULL)
  })
}

# GhostKnockoff 变量选择（单表型）— knockoff + FDR 控制

#' 基于GhostKnockoff执行变量选择（单表型）
#' @param genotype_matrix 矩阵，参考基因型，维度: 样本×SNP
#' @param z_vector 数值向量，长度为SNP数，对应单表型原始Z-score
#' @param n_samples 整数，GWAS样本量
#' @param n_knockoffs 整数，knockoff副本数量
#' @param fdr 数值，目标FDR
#' @param variant_ids 可选，字符向量，长度=SNP数，用作输出的ID（默认使用列名或索引）
#' @param verbose 逻辑值，是否打印详细日志
#' @return 列表，包含p值矩阵、LAVAKnock选择结果，以及按变异汇总的结果表
run_ghostknockoff_selection <- function(genotype_matrix, z_vector,
                                        n_samples = 20000,
                                        n_knockoffs = 5,
                                        fdr = 0.1,
                                        variant_ids = NULL,
                                        verbose = FALSE) {
  if (is.null(genotype_matrix) || is.null(z_vector)) {
    stop("genotype_matrix 或 z_vector 为空")
  }
  if (is.matrix(z_vector) || is.data.frame(z_vector)) {
    z_vector <- as.numeric(z_vector[, 1])
  }
  # 统一SNP维度（使用列-变异的约定）
  n_snps <- ncol(genotype_matrix)
  if (length(z_vector) != n_snps) {
    # 如果行是变异（某些用户可能提供SNP×样本），尝试转置
    if (nrow(genotype_matrix) == length(z_vector)) {
      genotype_matrix <- t(genotype_matrix)
      n_snps <- ncol(genotype_matrix)
    } else {
      stop("z_vector长度与变异数量不一致")
    }
  }

  zscore_matrix <- cbind(z_vector)
  colnames(zscore_matrix) <- c("z1")

  # 生成GhostKnockoff Z-scores
  gk <- generate_ghostknockoff_scores(
    genotype_matrix = genotype_matrix,
    zscore_matrix = zscore_matrix,
    n_samples = n_samples,
    n_knockoffs = n_knockoffs,
    verbose = verbose
  )
  if (is.null(gk) || is.null(gk$pheno1)) {
    warning("GhostKnockoff分数生成失败")
    return(NULL)
  }

  scores <- as.data.frame(gk$pheno1)
  # 计算双侧p值
  p0 <- 2 * pnorm(-abs(scores$org))
  pko <- sapply(1:n_knockoffs, function(k) {
    2 * pnorm(-abs(scores[[paste0("knock", k)]]))
  })
  colnames(pko) <- paste0("pval.knockoff", 1:n_knockoffs)
  pko <- as.data.frame(pko)

  # 变异ID
  if (is.null(variant_ids)) {
    if (!is.null(colnames(genotype_matrix))) {
      variant_ids <- colnames(genotype_matrix)
    } else {
      variant_ids <- as.character(seq_len(ncol(genotype_matrix)))
    }
  }
  if (!is.null(gk$index)) {
    valid_idx <- as.integer(gk$index)
    valid_idx <- valid_idx[valid_idx >= 1 & valid_idx <= length(variant_ids)]
    if (length(valid_idx) == nrow(scores)) {
      variant_ids <- variant_ids[valid_idx]
    } else if (length(variant_ids) >= nrow(scores)) {
      variant_ids <- variant_ids[seq_len(nrow(scores))]
    }
  } else if (length(variant_ids) >= nrow(scores)) {
    variant_ids <- variant_ids[seq_len(nrow(scores))]
  }

  # 使用LAVAKnock进行FDR控制的选择
  res <- LAVAKnock::LAVAKnock(
    M = n_knockoffs,
    p0 = p0,
    p_ko = as.matrix(pko),
    fdr = fdr,
    window_id = variant_ids,
    Rej.Bound = max(20000, length(p0))
  )

  threshold_val <- if (!is.null(res$W.threshold)) as.numeric(res$W.threshold) else NA_real_

  # 汇总为结果表
  out_tbl <- data.frame(
    id = variant_ids,
    pval.orginal = p0,
    pko,
    W = res$W,
    Qvalue = res$Qvalue,
    selected = res$W >= res$W.threshold,
    W_threshold = threshold_val,
    stringsAsFactors = FALSE
  )

  if (verbose) {
    cat(paste0("变量选择完成: 选中 ", sum(out_tbl$selected), "/", nrow(out_tbl), " 变异\n"))
  }

  list(
    pvals = list(p0 = p0, pko = pko),
    lavaknock = res,
    table = out_tbl,
    threshold = threshold_val
  )
}

plot_manhattan_w <- function(df, out_path, title = 'Manhattan plot (W statistic)', threshold = NULL) {
  df <- df[!is.na(df$chr) & !is.na(df$pos) & !is.na(df$W), , drop = FALSE]
  if (nrow(df) == 0) return(invisible())
  df <- df[order(df$chr, df$pos), , drop = FALSE]

  df$CHR <- as.integer(df$chr)
  df$BP <- as.numeric(df$pos)
  df$SNP <- ifelse(!is.na(df$id) & nzchar(df$id), df$id, paste0(df$CHR, ':', df$BP))
  if ('W_threshold' %in% colnames(df) && (is.null(threshold) || length(threshold) == 0)) {
    thr_candidates <- unique(stats::na.omit(df$W_threshold))
    if (length(thr_candidates) > 0) {
      threshold <- max(thr_candidates)
    }
  }

  colors <- c('#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b')
  highlight_color <- '#e41a1c'

  chr_levels <- sort(unique(df$CHR))
  df$color <- colors[(match(df$CHR, chr_levels) - 1) %% length(colors) + 1]
  df$cum_pos <- rep(NA_real_, nrow(df))
  tick_pos <- numeric(length(chr_levels))
  cumulative <- 0
  for (i in seq_along(chr_levels)) {
    chr_val <- chr_levels[i]
    idx <- which(df$CHR == chr_val)
    if (length(idx) == 0) next
    df$cum_pos[idx] <- df$BP[idx] + cumulative
    tick_pos[i] <- cumulative + stats::median(df$BP[idx])
    cumulative <- cumulative + max(df$BP[idx])
  }

  y_vals <- df$W
  if (!is.null(threshold) && length(threshold) > 0 && is.finite(threshold)) {
    y_vals <- c(y_vals, threshold)
  }
  y_lim <- range(y_vals, finite = TRUE)
  y_pad <- diff(y_lim)
  if (!is.finite(y_pad) || y_pad == 0) y_pad <- abs(y_lim[1]) * 0.1 + 1
  y_lim[1] <- y_lim[1] - 0.05 * y_pad
  y_lim[2] <- y_lim[2] + 0.05 * y_pad

  png(out_path, width = 1600, height = 600)
  old_par <- par(no.readonly = TRUE)
  on.exit({par(old_par); dev.off()}, add = TRUE)
  par(mar = c(5, 5, 3, 1))

  plot(
    df$cum_pos, df$W,
    col = df$color, pch = 16, cex = 1.1,
    xaxt = 'n', xlab = 'Chromosome', ylab = 'W statistic',
    main = title, cex.axis = 1.2, cex.lab = 1.4, cex.main = 1.6,
    las = 1, font.main = 2, font.lab = 2, ylim = y_lim
  )
  axis(1, at = tick_pos, labels = chr_levels, cex.axis = 1.2)
  axis(2, las = 1, cex.axis = 1.2)
  box()

  highlight_idx <- which(df$selected %in% TRUE)
  if (length(highlight_idx) > 0) {
    points(df$cum_pos[highlight_idx], df$W[highlight_idx],
           col = highlight_color, pch = 16, cex = 1.3)
  }

  thr_line <- NULL
  if (!is.null(threshold) && length(threshold) > 0 && is.finite(threshold)) {
    thr_line <- as.numeric(threshold[1])
    abline(h = thr_line, col = highlight_color, lwd = 2, lty = 2)
  }

  legend_entries <- character()
  legend_cols <- character()
  legend_pch <- numeric()
  legend_lty <- numeric()
  legend_lwd <- numeric()

  if (length(highlight_idx) > 0) {
    legend_entries <- c(legend_entries, 'Selected variants')
    legend_cols <- c(legend_cols, highlight_color)
    legend_pch <- c(legend_pch, 16)
    legend_lty <- c(legend_lty, NA)
    legend_lwd <- c(legend_lwd, NA)
  }

  if (!is.null(thr_line)) {
    legend_entries <- c(legend_entries, sprintf('Threshold (W = %.3f)', thr_line))
    legend_cols <- c(legend_cols, highlight_color)
    legend_pch <- c(legend_pch, NA)
    legend_lty <- c(legend_lty, 2)
    legend_lwd <- c(legend_lwd, 2)
  }

  if (length(legend_entries) > 0) {
    legend('topright',
           legend = legend_entries,
           col = legend_cols,
           pch = legend_pch,
           lty = legend_lty,
           lwd = legend_lwd,
           bty = 'n',
           cex = 1.0)
  }

  invisible()
}

# 核心处理函数（供外部调用）— 染色体/窗口批处理

#' 处理单个染色体的单个文件
#' @param chr 染色体编号
#' @param info_file 变异信息文件路径
#' @param gwas_data GWAS汇总统计数据
#' @param genotype_data 基因型数据
#' @param params 参数列表
#' @param verbose 是否详细输出
#' @return 分析结果列表
process_single_locus <- function(chr, info_file, gwas_data, genotype_data, params, verbose = FALSE) {

  if (verbose) cat(paste("处理文件:", basename(info_file), "\n"))

  # 加载变异信息
  variant_info <- as.data.frame(py_load_variant_info(info_file), stringsAsFactors = FALSE)
  if (is.null(variant_info) || nrow(variant_info) == 0) {
    return(NULL)
  }

  results <- list(
    univariate_results = data.frame(),
    bivariate_results = data.frame(),
    knockoff_results = NULL
  )

  # 准备ID与对齐
  if (!all(c("chr","pos_bp") %in% colnames(gwas_data))) {
    stop("gwas_data 缺少 chr/pos_bp 列")
  }
  gid <- paste0(gwas_data$chr, ":", gwas_data$pos_bp)
  if (is.null(colnames(genotype_data))) stop("genotype_data 缺少列名用于对齐")
  col_ids <- colnames(genotype_data)
  # 对齐Z分数到基因型顺序
  idx <- match(col_ids, gid)
  keep <- which(!is.na(idx))
  if (length(keep) < 2) {
    warning("基因型与GWAS变异对齐失败或匹配过少")
    return(NULL)
  }
  genotype_data <- genotype_data[, keep, drop = FALSE]
  z_mat <- as.matrix(gwas_data[idx[keep], c("zscore_pheno1", "zscore_pheno2")])
  rownames(z_mat) <- colnames(genotype_data)

  # 单变量分析
  if (params$run_mode %in% c("full", "lavaknock_only")) {
    locus_start <- min(variant_info$pos_bp)
    locus_end <- max(variant_info$pos_bp)

    univ_result <- run_univariate_analysis(
      genotype_data, z_mat, chr, locus_start, locus_end,
      params$sample_size, params$prune_threshold, verbose
    )

    if (!is.null(univ_result)) {
      results$univariate_results <- univ_result
    }
  }

  # 生成knockoff并进行双变量分析
  if (params$run_mode %in% c("full", "ghostknockoff_only")) {

    # 创建分析窗口
    windows <- create_analysis_windows(variant_info, params$window_size, params$window_step)

    if (nrow(windows) > 0) {

      # 生成knockoff（记录选择后的索引以便窗口映射）
      knockoff_scores <- generate_ghostknockoff_scores(
        genotype_matrix = genotype_data,
        zscore_matrix = z_mat,
        ld_threshold = params$ld_threshold,
        n_samples = params$sample_size,
        n_knockoffs = params$knockoffs,
        verbose = verbose
      )

      if (!is.null(knockoff_scores)) {

        # 对每个窗口进行双变量分析
        window_results <- data.frame()

        # 为variant_info构建至基因型列的索引映射
        v_ids <- paste0(variant_info$chr, ":", variant_info$pos_bp)
        map_idx <- match(v_ids, colnames(genotype_data))
        for (i in seq_len(nrow(windows))) {
          window <- windows[i, ]
          variant_indices <- window$variant_indices[[1]]
          geno_indices <- map_idx[variant_indices]
          geno_indices <- geno_indices[!is.na(geno_indices)]
          if (length(geno_indices) < 2) next
          window_genotype <- genotype_data[, geno_indices, drop = FALSE]

          # 将基因型列索引映射到GhostKnockoff选择后的索引空间
          if (is.null(knockoff_scores$index)) {
            stop("GhostKnockoff返回缺少index映射，无法对齐窗口Z-scores")
          }
          ko_idx <- match(geno_indices, knockoff_scores$index)
          ko_idx <- ko_idx[!is.na(ko_idx)]
          if (length(ko_idx) < 2) next
          # 提取窗口的knockoff分数（使用映射后的行索引）
          zscore_pheno1_window <- as.data.frame(knockoff_scores$pheno1[ko_idx, ])
          zscore_pheno2_window <- as.data.frame(knockoff_scores$pheno2[ko_idx, ])

          # 双变量分析
          bivar_result <- run_bivariate_analysis(
            window_genotype, zscore_pheno1_window, zscore_pheno2_window,
            chr, window$start, window$end, params$sample_size,
            params$prune_threshold, verbose
          )

          if (!is.null(bivar_result)) {
            window_results <- rbind(window_results, bivar_result)
          }
        }

        results$bivariate_results <- window_results

        # 对该文件进行knockoff过滤
        if (nrow(window_results) > 0) {
          results$knockoff_results <- run_knockoff_filter(
            window_results, params$fdr, params$knockoffs,
            verbose = verbose
          )
        }
      }
    }
  }

  return(results)
}

#' 双变量局部遗传相关性分析
run_bivariate_analysis <- function(genotype_matrix, zscore_pheno1, zscore_pheno2,
                                  chr, window_start, window_end,
                                  n_samples = 20000, prune_thresh = 99,
                                  verbose = FALSE) {

  if (verbose) cat("运行双变量局部遗传相关性分析...\n")

  # 数据验证
  if (is.null(genotype_matrix) || nrow(genotype_matrix) == 0 || ncol(genotype_matrix) == 0) {
    warning("基因型矩阵为空或维度为0")
    return(NULL)
  }

  if (is.null(zscore_pheno1) || nrow(zscore_pheno1) == 0) {
    warning("表型1 Z-score数据为空")
    return(NULL)
  }

  if (is.null(zscore_pheno2) || nrow(zscore_pheno2) == 0) {
    warning("表型2 Z-score数据为空")
    return(NULL)
  }

  # 确保zscore数据是data.frame格式
  if (!is.data.frame(zscore_pheno1)) {
    zscore_pheno1 <- as.data.frame(zscore_pheno1)
  }
  if (!is.data.frame(zscore_pheno2)) {
    zscore_pheno2 <- as.data.frame(zscore_pheno2)
  }

  # 强化数据检查和预处理
  # 1. 保证基因型矩阵为数值型，且无全零/无变异行
  genotype_matrix <- as.matrix(genotype_matrix)
  if (!is.numeric(genotype_matrix)) {
    genotype_matrix <- apply(genotype_matrix, c(1,2), as.numeric)
  }
  # 剔除无变异SNP
  snp_var <- apply(genotype_matrix, 1, var, na.rm=TRUE)
  valid_snps <- which(snp_var > 0 & !is.na(snp_var))
  if (length(valid_snps) < 2) {
    warning("有效SNP数量不足")
    return(NULL)
  }
  genotype_matrix <- genotype_matrix[valid_snps,,drop=FALSE]
  zscore_pheno1 <- zscore_pheno1[valid_snps,,drop=FALSE]
  zscore_pheno2 <- zscore_pheno2[valid_snps,,drop=FALSE]

  # 2. 保证Z-score数据为data.frame，行数与基因型矩阵一致，列名标准
  zscore_pheno1 <- as.data.frame(zscore_pheno1)
  zscore_pheno2 <- as.data.frame(zscore_pheno2)
  # 自动修正列名
  knockoff_names <- paste0("knock", seq_len(ncol(zscore_pheno1)-1))
  colnames(zscore_pheno1)[1] <- "org"
  if (length(knockoff_names) > 0) colnames(zscore_pheno1)[2:ncol(zscore_pheno1)] <- knockoff_names
  colnames(zscore_pheno2)[1] <- "org"
  if (length(knockoff_names) > 0) colnames(zscore_pheno2)[2:ncol(zscore_pheno2)] <- knockoff_names

  # 3. 检查维度完全一致
  if (nrow(genotype_matrix) != nrow(zscore_pheno1) || nrow(genotype_matrix) != nrow(zscore_pheno2)) {
    warning(paste("最终维度不匹配: 基因型", nrow(genotype_matrix), "表型1", nrow(zscore_pheno1), "表型2", nrow(zscore_pheno2)))
    return(NULL)
  }
  if (any(is.na(genotype_matrix)) || any(is.na(zscore_pheno1)) || any(is.na(zscore_pheno2))) {
    warning("输入数据存在NA值")
    return(NULL)
  }

  # 检查是否有足够的knockoff列
  required_cols <- c("org", paste0("knock", 1:5))  # 假设5个knockoffs

  if (!all(required_cols %in% colnames(zscore_pheno1))) {
    if (verbose) cat("表型1缺少knockoff列，尝试基本格式...\n")
    # 如果没有标准的knockoff列，检查是否有数字列
    if (ncol(zscore_pheno1) < 2) {
      warning("表型1 Z-score数据列数不足")
      return(NULL)
    }
  }

  if (!all(required_cols %in% colnames(zscore_pheno2))) {
    if (verbose) cat("表型2缺少knockoff列，尝试基本格式...\n")
    if (ncol(zscore_pheno2) < 2) {
      warning("表型2 Z-score数据列数不足")
      return(NULL)
    }
  }

  tryCatch({
    # 确保数据格式正确
    # LAVAKnock期望Z-score数据有特定的列格式
    if (!"org" %in% colnames(zscore_pheno1)) {
      # 如果没有org列，将第一列作为原始分数
      if (ncol(zscore_pheno1) >= 1) {
        colnames(zscore_pheno1)[1] <- "org"
      }
    }

    if (!"org" %in% colnames(zscore_pheno2)) {
      if (ncol(zscore_pheno2) >= 1) {
        colnames(zscore_pheno2)[1] <- "org"
      }
    }

    # 确保基因型矩阵是数值型
    if (!is.numeric(genotype_matrix)) {
      genotype_matrix <- as.matrix(as.numeric(genotype_matrix))
      dim(genotype_matrix) <- c(nrow(genotype_matrix), ncol(genotype_matrix))
    }

    # 检查基因型矩阵是否有变异
    genotype_var <- apply(genotype_matrix, 1, var, na.rm = TRUE)
    if (any(is.na(genotype_var)) || any(genotype_var == 0)) {
      if (verbose) cat("移除无变异的SNPs...\n")
      valid_snps <- !is.na(genotype_var) & genotype_var > 0
      genotype_matrix <- genotype_matrix[valid_snps, , drop = FALSE]
      zscore_pheno1 <- zscore_pheno1[valid_snps, , drop = FALSE]
      zscore_pheno2 <- zscore_pheno2[valid_snps, , drop = FALSE]

      if (nrow(genotype_matrix) < 2) {
        warning("有效SNPs数量不足")
        return(NULL)
      }
    }

    if (verbose) {
      cat("调试信息:\n")
      cat("  最终基因型矩阵维度:", paste(dim(genotype_matrix), collapse = "x"), "\n")
      cat("  最终表型1维度:", paste(dim(zscore_pheno1), collapse = "x"), "\n")
      cat("  最终表型2维度:", paste(dim(zscore_pheno2), collapse = "x"), "\n")
      cat("  表型1列名:", paste(colnames(zscore_pheno1), collapse = ", "), "\n")
      cat("  表型2列名:", paste(colnames(zscore_pheno2), collapse = ", "), "\n")
    }

    # 只尝试LAVAKnock_bivariate，失败时直接返回NULL
    result <- tryCatch({
      LAVAKnock::LAVAKnock_bivariate(
        Zscore_pheno1_window = zscore_pheno1,
        Zscore_pheno2_window = zscore_pheno2,
        G_window = genotype_matrix,
        chr = chr,
        window_start = window_start,
        window_end = window_end,
        n = n_samples,
        prune.thresh = prune_thresh
      )
    }, error = function(e1) {
      warning(paste("LAVAKnock_bivariate失败:", e1$message))
      if (verbose) {
        cat("详细错误调试信息:\n")
        cat("  错误消息:", e1$message, "\n")
        cat("  基因型矩阵维度:", paste(dim(genotype_matrix), collapse = "x"), "\n")
        cat("  基因型矩阵类型:", class(genotype_matrix), "\n")
        cat("  表型1维度:", paste(dim(zscore_pheno1), collapse = "x"), "\n")
        cat("  表型2维度:", paste(dim(zscore_pheno2), collapse = "x"), "\n")
        cat("  表型1列名:", paste(colnames(zscore_pheno1), collapse = ", "), "\n")
        cat("  表型2列名:", paste(colnames(zscore_pheno2), collapse = ", "), "\n")
        cat("  基因型矩阵前几行前几列:\n")
        print(genotype_matrix[seq_len(min(3, nrow(genotype_matrix))), seq_len(min(3, ncol(genotype_matrix)))])
      }
      return(NULL)
    })
    if (verbose && !is.null(result)) cat("双变量分析完成\n")
    return(result)

  }, error = function(e) {
    warning(paste("双变量分析失败:", e$message))
    if (verbose) {
      cat("详细错误调试信息:\n")
      cat("  错误消息:", e$message, "\n")
      cat("  基因型矩阵维度:", paste(dim(genotype_matrix), collapse = "x"), "\n")
      cat("  基因型矩阵类型:", class(genotype_matrix), "\n")
      cat("  表型1维度:", paste(dim(zscore_pheno1), collapse = "x"), "\n")
      cat("  表型2维度:", paste(dim(zscore_pheno2), collapse = "x"), "\n")
      cat("  表型1列名:", paste(colnames(zscore_pheno1), collapse = ", "), "\n")
      cat("  表型2列名:", paste(colnames(zscore_pheno2), collapse = ", "), "\n")
      cat("  基因型矩阵前几行前几列:\n")
      print(genotype_matrix[1:min(3, nrow(genotype_matrix)), 1:min(3, ncol(genotype_matrix))])
    }
    return(NULL)
  })
}

#' Knockoff 多重检验校正
run_knockoff_filter <- function(bivariate_results, fdr = 0.1, n_knockoffs = 5,
                               rej_bound = 20000, verbose = FALSE) {

  if (verbose) cat("运行knockoff多重检验校正...\n")

  if (nrow(bivariate_results) == 0) {
    warning("没有双变量分析结果")
    return(NULL)
  }

  tryCatch({
    # 提取原始p值
    p0 <- bivariate_results$pval.orginal

    # 提取knockoff p值
    p_ko <- bivariate_results[, paste0("pval.knockoff", 1:n_knockoffs), drop = FALSE]

    # 创建窗口ID
    window_id <- paste(bivariate_results$chr,
                      bivariate_results$window.start,
                      bivariate_results$window.end, sep = ":")

    result <- LAVAKnock::LAVAKnock(
      M = n_knockoffs,
      p0 = p0,
      p_ko = as.matrix(p_ko),
      fdr = fdr,
      window_id = window_id,
      Rej.Bound = rej_bound
    )

    if (verbose) {
      cat(paste("检测到", length(result$window_sign), "个显著窗口\n"))
      cat(paste("W阈值:", round(result$W.threshold, 4), "\n"))
    }

    return(result)

  }, error = function(e) {
    warning(paste("Knockoff过滤失败:", e$message))
    return(NULL)
  })
}

execute_pipeline <- function(opts) {
  auto_temp_paths <- character()
  auto_trait_names <- NULL

  if (!is.null(opts$panel) && is.null(opts$ref_plink)) {
    opts$ref_plink <- opts$panel
  }
  panel_default <- if (!is.null(opts$panel)) opts$panel else if (!is.null(opts$ref_plink)) opts$ref_plink else file.path('g1000_eur', 'g1000_eur')

  trait_option_vec <- NULL
  if (!is.null(opts$zscore_traits)) {
    trait_option_vec <- parse_trait_option(opts$zscore_traits)
  }

  if (is.null(opts$info) && !is.null(opts$zscore)) {
    trait_arg <- if (is.null(trait_option_vec)) NULL else as.list(trait_option_vec)
    auto_res <- py_auto_prepare_inputs(opts$zscore, panel_default, trait_arg, opts$verbose)
    info_df <- as.data.frame(auto_res$variant_info, stringsAsFactors = FALSE)
    safe_base <- gsub('[^A-Za-z0-9_]+', '_', tools::file_path_sans_ext(basename(opts$zscore)))
    info_path <- tempfile(paste0('auto_info_', safe_base, '_'), fileext = '.csv')
    readr::write_csv(info_df, info_path)

    gwas_paths <- character()
    trait_names <- if (is.null(auto_res$traits)) character() else as.character(auto_res$traits)
    tables <- auto_res$gwas_tables
    if (!is.null(tables)) {
      if (length(trait_names) < length(tables)) {
        auto_names <- paste0('trait', seq_along(tables))
        trait_names <- head(c(trait_names, auto_names), length(tables))
      }
      for (i in seq_along(tables)) {
        gwas_df <- as.data.frame(tables[[i]], stringsAsFactors = FALSE)
        gwas_df <- gwas_df[!is.na(gwas_df$CHR) & !is.na(gwas_df$POS) & !is.na(gwas_df$Z), , drop = FALSE]
        gwas_path <- tempfile(paste0('auto_gwas_', safe_base, '_', trait_names[i], '_'), fileext = '.tsv')
        readr::write_tsv(gwas_df, gwas_path)
        gwas_paths <- c(gwas_paths, gwas_path)
      }
    }
    if (opts$verbose) {
      cat(sprintf('自动生成Info与GWAS：%d 个变异\n', nrow(info_df)))
    }
    opts$info <- info_path
    if (is.null(opts$gwas1) && length(gwas_paths) >= 1) opts$gwas1 <- gwas_paths[1]
    if (is.null(opts$gwas2) && length(gwas_paths) >= 2) opts$gwas2 <- gwas_paths[2]
    if (is.null(opts$ref_plink)) opts$ref_plink <- panel_default
    auto_temp_paths <- c(auto_temp_paths, info_path, gwas_paths)
    auto_trait_names <- trait_names
  }

  if (is.null(opts$mode)) {
    stop('必须提供 --mode')
  }

  if (is.null(opts$info)) {
    stop('必须提供 --info 或同时提供 --zscore 与 --panel/--ref_plink')
  }

  if (!is.null(auto_trait_names)) {
    if (length(auto_trait_names) >= 1 && (is.null(opts$pheno1_name) || opts$pheno1_name == 'trait1')) {
      opts$pheno1_name <- auto_trait_names[1]
    }
    if (length(auto_trait_names) >= 2 && (is.null(opts$pheno2_name) || opts$pheno2_name == 'trait2')) {
      opts$pheno2_name <- auto_trait_names[2]
    }
  }

  dir.create(opts$outdir, showWarnings = FALSE, recursive = TRUE)

  if (is.null(opts$ld_coord) || !nzchar(opts$ld_coord)) {
    stop('必须提供 --ld_coord (GRCh37, GRCh38 或自定义LD block文件路径)')
  }

  gene_catalog <- if (!is.null(opts$gene_catalog) && nzchar(opts$gene_catalog)) {
    load_gene_catalog(opts$gene_catalog)
  } else {
    gene_key <- infer_gene_catalog_key(opts$ld_coord, default = "GRCh37")
    tryCatch(
      load_gene_catalog(gene_key),
      error = function(e) {
        warning(sprintf("gene catalog 加载失败 (%s)，尝试使用默认 GRCh37: %s", gene_key, e$message))
        load_gene_catalog("GRCh37")
      }
    )
  }

  info <- as.data.frame(py_load_variant_info(opts$info), stringsAsFactors = FALSE)
  if (is.null(info) || nrow(info) == 0) stop('变异信息文件读取失败或为空')

  g1 <- NULL
  g2 <- NULL
  if (!is.null(opts$gwas1)) {
    gwas1_df <- as.data.frame(py_load_gwas_table(opts$gwas1), stringsAsFactors = FALSE)
    gwas1_df$CHR <- as.integer(gwas1_df$CHR)
    gwas1_df$POS <- as.numeric(gwas1_df$POS)
    gwas1_df$Z <- as.numeric(gwas1_df$Z)
    gwas1_df <- gwas1_df[!is.na(gwas1_df$CHR) & !is.na(gwas1_df$POS) & !is.na(gwas1_df$Z), , drop = FALSE]
    g1 <- merge(info, gwas1_df, by.x = c('chr', 'pos_bp'), by.y = c('CHR', 'POS'))
    if (nrow(g1) == 0) stop('变异信息与GWAS在CHR+POS上没有交集')
    if (!is.null(opts$gwas2)) {
      gwas2_df <- as.data.frame(py_load_gwas_table(opts$gwas2), stringsAsFactors = FALSE)
      gwas2_df$CHR <- as.integer(gwas2_df$CHR)
      gwas2_df$POS <- as.numeric(gwas2_df$POS)
      gwas2_df$Z <- as.numeric(gwas2_df$Z)
      gwas2_df <- gwas2_df[!is.na(gwas2_df$CHR) & !is.na(gwas2_df$POS) & !is.na(gwas2_df$Z), , drop = FALSE]
      g2 <- merge(info, gwas2_df, by.x = c('chr', 'pos_bp'), by.y = c('CHR', 'POS'))
      if (nrow(g2) == 0) stop('变异信息与GWAS在CHR+POS上没有交集 (gwas2)')
    }
  } else if (!is.null(opts$multi_gwas) && !is.null(opts$zcols)) {
    mg <- py_load_multi_gwas_table(opts$multi_gwas, opts$zcols)
    m <- as.data.frame(mg$data, stringsAsFactors = FALSE)
    m$CHR <- as.integer(m$CHR)
    m$POS <- as.numeric(m$POS)
    m <- m[!is.na(m$CHR) & !is.na(m$POS), , drop = FALSE]
    merged <- merge(info, m, by.x = c('chr', 'pos_bp'), by.y = c('CHR', 'POS'))
    if (nrow(merged) == 0) stop('变异信息与multi_gwas无法匹配')
    zc <- as.character(mg$zcols)
    if (length(zc) < 1) stop('未在 multi_gwas 中找到指定的Z列')
    g1 <- merged
    g1$Z <- as.numeric(merged[[zc[1]]])
    if (length(zc) >= 2) {
      g2 <- merged
      g2$Z <- as.numeric(merged[[zc[2]]])
    }
  } else {
    stop('未提供 GWAS 输入。请使用 --gwas1/--gwas2 或 --multi_gwas 配合 --zcols')
  }

  if (!is.null(g2)) {
    shared_idx <- !is.na(g1$Z) & !is.na(g2$Z)
    g1 <- g1[shared_idx, ]
    g2 <- g2[shared_idx, ]
  } else {
    g1 <- g1[!is.na(g1$Z), ]
  }
  if (nrow(g1) == 0) {
    stop('未能获得有效的GWAS记录，请检查输入文件是否与变异信息匹配')
  }

  gwas_data <- data.frame(
    rsid = g1$rsid,
    chr = g1$chr,
    pos_bp = g1$pos_bp,
    zscore_pheno1 = g1$Z,
    zscore_pheno2 = if (!is.null(g2)) g2$Z else rep(NA_real_, nrow(g1)),
    stringsAsFactors = FALSE
  )
  attr(gwas_data, 'pheno1_name') <- opts$pheno1_name
  attr(gwas_data, 'pheno2_name') <- if (!is.null(g2)) opts$pheno2_name else NA_character_

  target_snp_ids <- unique(stats::na.omit(gwas_data$rsid))
  chrpos_ids <- paste0(gwas_data$chr, ':', gwas_data$pos_bp)
  using_reference_panel <- FALSE

  geno_raw <- NULL

  if (!is.null(opts$geno_rds)) {
    if (!file.exists(opts$geno_rds)) stop('geno_rds文件不存在')
    # RDS 文件仍由 R 端读取
    geno_raw <- load_genotype_rds(opts$geno_rds)
  } else if (!is.null(opts$geno_csv)) {
    if (!file.exists(opts$geno_csv)) stop('geno_csv文件不存在')
    geno_raw <- py_load_genotype_csv(opts$geno_csv, fmt = opts$geno_format)
  } else if (!is.null(opts$geno_plink)) {
    geno_raw <- load_genotype_plink(opts$geno_plink, snp_ids = target_snp_ids, chrpos_ids = chrpos_ids)$matrix
  } else {
    prefix <- if (!is.null(opts$ref_plink)) opts$ref_plink else file.path('g1000_eur', 'g1000_eur')
    using_reference_panel <- TRUE
    message('未提供真实基因型，自动使用参考面板: ', prefix)
    geno_raw <- load_genotype_plink(prefix, snp_ids = target_snp_ids, chrpos_ids = chrpos_ids)$matrix
  }

  geno_source_type <- if (using_reference_panel) 'reference' else 'real'

  pheno1_attr <- attr(gwas_data, 'pheno1_name')
  pheno2_attr <- attr(gwas_data, 'pheno2_name')

  ali <- align_genotype_to_gwas(geno_raw, gwas_data, prefer_rsid = TRUE)
  geno_matrix <- ali$matrix
  match_idx <- as.integer(ali$gwas_indices) + 1L
  if (length(match_idx) < 2 || any(is.na(match_idx))) {
    stop('基因型与GWAS变异匹配数量不足')
  }
  gwas_selected <- gwas_data[match_idx, , drop = FALSE]
  variant_ids <- colnames(geno_matrix)
  if (is.null(variant_ids)) {
    variant_ids <- paste0(gwas_selected$chr, ':', gwas_selected$pos_bp)
    colnames(geno_matrix) <- variant_ids
  }
  variant_positions <- as.numeric(ali$positions)
  if (length(variant_positions) != length(variant_ids) || any(is.na(variant_positions))) {
    variant_positions <- gwas_selected$pos_bp
  }
  z_vec <- as.numeric(gwas_selected$zscore_pheno1)
  rsid_lookup <- ali$rsid
  if (is.null(rsid_lookup) || length(rsid_lookup) != length(variant_ids)) {
    rsid_lookup <- gwas_selected$rsid
  }
  rsid_lookup <- as.character(rsid_lookup)
  attr(gwas_selected, 'pheno1_name') <- pheno1_attr
  attr(gwas_selected, 'pheno2_name') <- pheno2_attr
  gwas_data <- gwas_selected

  if (opts$mode == 'select') {
    if (length(variant_ids) < 5) {
      stop('匹配到的变异数量不足以执行变量选择')
    }

    selection_table <- NULL
    chunk_summary <- NULL
    chunk_tables_data <- list()
    chunk_records <- list()
    ld_reference_used <- NA_character_

    variant_chr_vec <- gwas_selected$chr
    variant_pos_vec <- variant_positions

    chunk_plan <- partition_ld_blocks(
      variant_chr_vec,
      variant_pos_vec,
      opts$ld_coord,
      workers = opts$threads
    )
    chunk_list <- c(
      if (!is.null(chunk_plan$chunks)) chunk_plan$chunks else list(),
      if (!is.null(chunk_plan$fallback)) chunk_plan$fallback else list()
    )

    if (length(chunk_list) > 0) {
      # n_cores <- min(detectCores() - 1, 8)  # 自动选择核心数，最多 8
      # if (opts$verbose) {
        cat(sprintf("Using %d cores for parallel processing...\n", opts$threads))
      # }

      # 定义单个 chunk 的处理函数
      process_chunk <- function(chunk, chunk_id) {
        cat(sprintf("[Core %d] Processing chunk %d on PID %d\n",
              as.integer(Sys.getpid()), chunk_id, Sys.getpid()))

        idx0 <- as.integer(chunk$indices)
        if (is.null(idx0) || length(idx0) < 5) return(NULL)
        idx <- idx0 + 1L
        sub_geno <- geno_matrix[, idx, drop = FALSE]
        sub_z <- z_vec[idx]
        sub_ids <- variant_ids[idx]
        # if (opts$verbose) {
          cat(sprintf('Chunk %d: %d SNPs (chr %s: %.0f-%.0f)\n',
                      chunk_id, length(idx), as.character(chunk$chr),
                      as.numeric(chunk$start), as.numeric(chunk$end)))
        # }
        res <- run_ghostknockoff_selection(
          genotype_matrix = sub_geno,
          z_vector = sub_z,
          n_samples = opts$n,
          n_knockoffs = opts$knockoffs,
          fdr = opts$fdr,
          variant_ids = sub_ids,
          verbose = opts$verbose
        )
        if (is.null(res)) return(NULL)

        tbl <- res$table
        chunk_chr <- as.integer(chunk$chr)
        chunk_start <- as.numeric(chunk$start)
        chunk_end <- as.numeric(chunk$end)
        chunk_source <- if (!is.null(chunk$source)) as.character(chunk$source) else 'ld_block'
        chunk_size_bp <- chunk_end - chunk_start + 1
        tbl$chunk_id <- chunk_id
        tbl$chunk_chr <- chunk_chr
        tbl$chunk_start <- chunk_start
        tbl$chunk_end <- chunk_end
        tbl$chunk_size_bp <- chunk_size_bp
        tbl$partition_source <- chunk_source
        chunk_record <- data.frame(
          chunk_id = chunk_id,
          chr = chunk_chr,
          start = chunk_start,
          end = chunk_end,
          n_snps = length(idx),
          selected = sum(tbl$selected),
          source = chunk_source,
          chunk_size_bp = chunk_size_bp,
          stringsAsFactors = FALSE
        )
        list(table = tbl, record = chunk_record)
      }

      # 并行处理所有 chunk
      results <- mclapply(seq_along(chunk_list), function(i) {
        process_chunk(chunk_list[[i]], i)
      }, mc.cores = opts$threads)

      # 去除空结果
      results <- Filter(Negate(is.null), results)

      quit()

      if (length(results) > 0) {
        chunk_tables_data <- lapply(results, `[[`, "table")
        chunk_records <- lapply(results, `[[`, "record")

        selection_table <- do.call(rbind, chunk_tables_data)
        chunk_summary <- do.call(rbind, chunk_records)
        ld_reference_used <- opts$ld_coord
      }
    }

    if (is.null(selection_table)) {
      sel <- run_ghostknockoff_selection(
        genotype_matrix = geno_matrix,
        z_vector = z_vec,
        n_samples = opts$n,
        n_knockoffs = opts$knockoffs,
        fdr = opts$fdr,
        variant_ids = variant_ids,
        verbose = opts$verbose
      )
      if (is.null(sel)) stop('变量选择失败')
      selection_table <- sel$table
      chunk_summary <- NULL
      ld_reference_used <- NA_character_
    }

    if (!is.null(rsid_lookup)) {
      rsid_map <- rsid_lookup
      names(rsid_map) <- variant_ids
      selection_table$rsid <- rsid_map[selection_table$id]
    }
    parts <- strsplit(selection_table$id, ':', fixed = TRUE)
    selection_table$chr <- suppressWarnings(as.integer(vapply(parts, function(x) if (length(x) >= 1) x[1] else NA_character_, character(1), USE.NAMES = FALSE)))
    pos_raw <- vapply(parts, function(x) if (length(x) >= 2) x[2] else NA_character_, character(1), USE.NAMES = FALSE)
    selection_table$pos <- suppressWarnings(as.integer(gsub('[^0-9]', '', pos_raw)))
    selection_table$geno_source <- geno_source_type

    if (!is.null(gene_catalog) && nrow(gene_catalog) > 0) {
      selection_table$nearest_gene <- annotate_nearest_gene(selection_table$chr, selection_table$pos, gene_catalog)
    } else {
      selection_table$nearest_gene <- NA_character_
    }

    plot_cols <- intersect(c('id', 'chr', 'pos', 'W', 'selected', 'W_threshold'), colnames(selection_table))
    plot_df <- selection_table[, plot_cols, drop = FALSE]
    drop_cols <- intersect(c('chunk_id', 'chunk_chr', 'chunk_start', 'chunk_end', 'chunk_size_bp',
                             'partition_source', 'rsid', 'chr', 'pos', 'geno_source'), colnames(selection_table))
    final_table <- selection_table
    if (length(drop_cols) > 0) {
      final_table <- final_table[, setdiff(colnames(final_table), drop_cols), drop = FALSE]
    }

    out_dir <- file.path(opts$outdir, 'selection')
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    old_files <- list.files(out_dir, full.names = TRUE)
    if (length(old_files) > 0) unlink(old_files, recursive = TRUE, force = TRUE)

    prefix_base <- if (!is.null(opts$zscore)) {
      gsub('[^A-Za-z0-9_]+', '_', tools::file_path_sans_ext(basename(opts$zscore)))
    } else {
      tools::file_path_sans_ext(basename(opts$info))
    }
    prefix <- paste0(prefix_base, '_', gsub('[^A-Za-z0-9_]+', '_', attr(gwas_data, 'pheno1_name')))

    out_csv <- file.path(out_dir, paste0(prefix, '_selection.csv'))
    readr::write_csv(final_table, out_csv)

    summary_df <- data.frame(
      phenotype = attr(gwas_data, 'pheno1_name'),
      total_snps = nrow(selection_table),
      selected = sum(selection_table$selected),
      total_chunks = if (!is.null(chunk_summary)) nrow(chunk_summary) else NA_integer_,
      ld_reference = if (!is.null(chunk_summary) && nrow(chunk_summary) > 0) ld_reference_used else NA_character_,
      stringsAsFactors = FALSE
    )
    summary_csv <- file.path(out_dir, paste0(prefix, '_summary.csv'))
    readr::write_csv(summary_df, summary_csv)

    if (!is.null(chunk_summary) && nrow(chunk_summary) > 0) {
      if (!is.na(ld_reference_used)) {
        chunk_summary$ld_reference <- ld_reference_used
      }
      chunk_csv <- file.path(out_dir, paste0(prefix, '_chunk_summary.csv'))
      readr::write_csv(chunk_summary, chunk_csv)

      chunk_dir <- file.path(out_dir, paste0(prefix, '_chunks'))
      if (dir.exists(chunk_dir)) unlink(chunk_dir, recursive = TRUE, force = TRUE)
      dir.create(chunk_dir, recursive = TRUE)
      if (length(chunk_tables_data) > 0) {
        for (tbl in chunk_tables_data) {
          chunk_file <- file.path(chunk_dir, sprintf('%s_chunk_%02d.csv', prefix, tbl$chunk_id[1]))
          readr::write_csv(tbl, chunk_file)
        }
      }
    }

    manhattan_path <- file.path(out_dir, paste0(prefix, '_manhattan.png'))
    threshold_val <- NULL
    if ('W_threshold' %in% colnames(plot_df)) {
      thr_candidates <- unique(stats::na.omit(plot_df$W_threshold))
      if (length(thr_candidates) > 0) threshold_val <- max(thr_candidates)
    }
    plot_manhattan_w(plot_df, manhattan_path, title = paste0(prefix_base, ' (', attr(gwas_data, 'pheno1_name'), ')'), threshold = threshold_val)

    if (opts$verbose) cat(sprintf('已保存结果: %s\n', out_csv))

  } else if (opts$mode == 'correl') {
    params <- list(
      data_dir = dirname(opts$info),
      output_dir = opts$outdir,
      sample_size = opts$n,
      window_size = opts$win,
      window_step = opts$step,
      knockoffs = opts$knockoffs,
      fdr = opts$fdr,
      ld_threshold = 0.75,
      prune_threshold = 99,
      run_mode = 'ghostknockoff_only'
    )
    chr <- unique(gwas_data$chr)[1]
    res <- process_single_locus(
      chr = chr,
      info_file = opts$info,
      gwas_data = gwas_data,
      genotype_data = geno_matrix,
      params = params,
      verbose = opts$verbose
    )
    if (is.null(res)) stop('相关性分析失败')

    out_dir <- file.path(opts$outdir, 'correlation')
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    prefix <- paste0(
      tools::file_path_sans_ext(basename(opts$info)),
      '_', gsub('[^A-Za-z0-9_]+', '_', attr(gwas_data, 'pheno1_name')),
      '__', gsub('[^A-Za-z0-9_]+', '_', attr(gwas_data, 'pheno2_name'))
    )
    if (is.null(res$bivariate_results) || nrow(res$bivariate_results) == 0) {
      empty_biv <- data.frame(
        chr = integer(),
        window.start = integer(),
        window.end = integer(),
        pval.orginal = numeric(),
        stringsAsFactors = FALSE
      )
      readr::write_csv(empty_biv, file.path(out_dir, paste0(prefix, '_bivariate.csv')))
    } else {
      readr::write_csv(res$bivariate_results, file.path(out_dir, paste0(prefix, '_bivariate.csv')))
    }
    if (!is.null(res$knockoff_results) && !is.null(res$knockoff_results$window_sign) && length(res$knockoff_results$window_sign) > 0) {
      ko <- res$knockoff_results
      sig <- data.frame(
        window = ko$window_sign,
        W = ko$W[match(ko$window_sign, names(ko$W))],
        Qvalue = ko$Qvalue[match(ko$window_sign, names(ko$Qvalue))],
        stringsAsFactors = FALSE
      )
      readr::write_csv(sig, file.path(out_dir, paste0(prefix, '_significant_windows.csv')))
    } else {
      empty_sig <- data.frame(window = character(), W = numeric(), Qvalue = numeric(), stringsAsFactors = FALSE)
      readr::write_csv(empty_sig, file.path(out_dir, paste0(prefix, '_significant_windows.csv')))
    }

    if (!is.null(res$bivariate_results) && nrow(res$bivariate_results) > 0 && !is.null(res$knockoff_results)) {
      ko <- res$knockoff_results
      window_ids <- paste(res$bivariate_results$chr,
                          res$bivariate_results$window.start,
                          res$bivariate_results$window.end,
                          sep = ":")
      W_vec <- ko$W
      if (is.null(names(W_vec))) {
        names(W_vec) <- window_ids[seq_along(W_vec)]
      }
      w_match <- W_vec[match(window_ids, names(W_vec))]
      plot_df <- data.frame(
        id = window_ids,
        chr = res$bivariate_results$chr,
        pos = floor((res$bivariate_results$window.start + res$bivariate_results$window.end) / 2),
        W = as.numeric(w_match),
        selected = window_ids %in% ko$window_sign,
        W_threshold = if (!is.null(ko$W.threshold)) as.numeric(ko$W.threshold) else NA_real_,
        stringsAsFactors = FALSE
      )
      plot_df <- plot_df[!is.na(plot_df$chr) & !is.na(plot_df$pos) & !is.na(plot_df$W), , drop = FALSE]
      if (nrow(plot_df) > 0) {
        thr_val <- if (!is.null(ko$W.threshold)) as.numeric(ko$W.threshold) else NULL
        manhattan_path <- file.path(out_dir, paste0(prefix, '_manhattan.png'))
        plot_manhattan_w(
          plot_df,
          manhattan_path,
          title = paste0(prefix, ' (W statistic)'),
          threshold = thr_val
        )
      }
    }

    if (opts$verbose) cat('相关性分析完成并已保存结果\n')
  } else {
    stop('未知的mode。请使用 select 或 correl')
  }

  list(auto_temp_paths = auto_temp_paths)
}
