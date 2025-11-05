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

if (file.exists("plotting.R")) {
  source("plotting.R")
}

# Python 加速（直接从 accelerator.py 载入）
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

py_build_info_from_gwas <- function(path) {
  reticulate::py_to_r(py_env$build_info_from_gwas(path))
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

py_align_dual_traits <- function(primary_df, secondary_df = NULL) {
  sec <- if (is.null(secondary_df)) NULL else reticulate::r_to_py(secondary_df)
  res <- reticulate::py_to_r(py_env$align_dual_traits(reticulate::r_to_py(primary_df), sec))
  as.data.frame(res, stringsAsFactors = FALSE)
}

py_align_positions <- function(variant_ids, variant_info, fallback_positions = NULL) {
  info_py <- reticulate::r_to_py(as.data.frame(variant_info, stringsAsFactors = FALSE))
  reticulate::py_to_r(py_env$align_positions_with_info(variant_ids, info_py, fallback_positions))
}

merge_gwas_with_info <- function(info_df, gwas_df) {
  if (is.null(gwas_df)) return(NULL)
  df <- as.data.frame(gwas_df, stringsAsFactors = FALSE)
  required <- c("CHR", "POS", "Z")
  if (!all(required %in% names(df))) {
    stop("GWAS 文件缺少 CHR/POS/Z 列")
  }
  df$CHR <- as.integer(df$CHR)
  df$POS <- as.numeric(df$POS)
  df$Z <- as.numeric(df$Z)
  df <- df[!is.na(df$CHR) & !is.na(df$POS) & !is.na(df$Z), , drop = FALSE]
  if (nrow(df) == 0) return(NULL)

  merged <- merge(info_df, df, by.x = c("chr", "pos_bp"), by.y = c("CHR", "POS"), all = FALSE)
  if (!"rsid" %in% names(merged) || all(is.na(merged$rsid))) {
    merged$rsid <- paste0(merged$chr, ":", merged$pos_bp)
  }
  merged
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

align_positions_with_info <- function(variant_ids, variant_info, fallback_positions = NULL) {
  if (is.null(variant_info) || nrow(variant_info) == 0) {
    return(if (is.null(fallback_positions)) rep(NA_real_, length(variant_ids)) else fallback_positions)
  }

  key <- NULL
  if ("rsid" %in% names(variant_info)) {
    key <- as.character(variant_info$rsid)
  } else {
    chr_col <- if ("chr" %in% names(variant_info)) variant_info$chr else variant_info[[which(grepl("^chr$", names(variant_info), ignore.case = TRUE))[1]]]
    pos_col <- NULL
    if ("pos_bp" %in% names(variant_info)) {
      pos_col <- variant_info$pos_bp
    } else {
      pos_idx <- which(grepl("pos", names(variant_info), ignore.case = TRUE))
      if (length(pos_idx) > 0) pos_col <- variant_info[[pos_idx[1]]]
    }
    if (!is.null(chr_col) && !is.null(pos_col)) {
      key <- paste0(chr_col, ":", pos_col)
    }
  }

  pos_values <- NULL
  if ("pos_bp" %in% names(variant_info)) {
    pos_values <- as.numeric(variant_info$pos_bp)
  } else {
    pos_idx <- which(grepl("pos", names(variant_info), ignore.case = TRUE))
    if (length(pos_idx) > 0) pos_values <- as.numeric(variant_info[[pos_idx[1]]])
  }

  aligned <- rep(NA_real_, length(variant_ids))
  if (!is.null(key) && !is.null(pos_values)) {
    match_idx <- match(variant_ids, key)
    take <- !is.na(match_idx)
    if (any(take)) {
      aligned[take] <- pos_values[match_idx[take]]
    }
  }

  if (!is.null(fallback_positions)) {
    aligned[is.na(aligned)] <- fallback_positions[is.na(aligned)]
  }

  aligned
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

  # 移除NA值，同时保留原始行索引
  valid_indices <- which(!is.na(positions))
  if (length(valid_indices) == 0) {
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
    idx_in_window <- which(positions >= start_pos & positions <= end_pos)

    if (length(idx_in_window) >= 5) {  # 至少5个变异
      original_rows <- valid_indices[idx_in_window]
      windows <- rbind(windows, data.frame(
        window_id = window_id,
        chr = chr,
        start = start_pos,
        end = end_pos,
        n_variants = length(original_rows),
        variant_indices = I(list(original_rows))
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

  G <- as.matrix(genotype_matrix)
  if (nrow(G) < 2 || ncol(G) < 2) {
    if (verbose) message("跳过单变量分析: 有效样本或变异数量不足")
    return(NULL)
  }

  Z <- as.matrix(zscore_matrix)
  if (nrow(Z) != ncol(G)) {
    stop("单变量分析失败: Zscore行数与基因型变异数量不一致")
  }

  effective_n <- n_samples
  geno_n <- nrow(G)
  if (!is.na(geno_n) && geno_n > 1 && geno_n != n_samples) {
    if (verbose) {
      message(sprintf(
        "  注意: 提供的样本量(n=%s)与基因型行数(n=%s)不一致，单变量分析将使用后者以匹配矩阵维度",
        n_samples, geno_n
      ))
    }
    effective_n <- geno_n
  }

  tryCatch({
    result <- LAVAKnock::LAVAKnock_univariate(
      Zscore = Z,
      G_locus = G,
      chr = chr,
      locus_start = locus_start,
      locus_end = locus_end,
      n = effective_n,
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

  knock_cols <- paste0("knock", seq_len(n_knockoffs))
  T0 <- abs(scores$org)
  Tk <- abs(as.matrix(scores[, knock_cols, drop = FALSE]))

  filter_res <- GhostKnockoff::GhostKnockoff.filter(T0, Tk)
  q_vals <- as.numeric(filter_res$q)
  tau_vals <- as.numeric(filter_res$tau)
  kappa_vals <- as.integer(filter_res$kappa)

  selected_flag <- !is.na(q_vals) & q_vals <= fdr
  threshold_val <- if (any(selected_flag, na.rm = TRUE)) {
    min(tau_vals[selected_flag], na.rm = TRUE)
  } else {
    NA_real_
  }

  leading_source <- ifelse(
    is.na(kappa_vals) | kappa_vals == 0,
    "original",
    paste0("knock", kappa_vals)
  )

  out_tbl <- data.frame(
    id = variant_ids,
    pval.orginal = p0,
    pko,
    W = tau_vals,
    Qvalue = q_vals,
    selected = selected_flag,
    W_threshold = threshold_val,
    leading_source = leading_source,
    stringsAsFactors = FALSE
  )

  if (verbose) {
    cat(paste0("变量选择完成: 选中 ", sum(out_tbl$selected), "/", nrow(out_tbl), " 变异\n"))
  }

  list(
    pvals = list(p0 = p0, pko = pko),
    filter = filter_res,
    table = out_tbl,
    threshold = threshold_val
  )
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

  lava_sample_n <- nrow(genotype_data)
  if (!is.null(params$sample_size) && !is.na(params$sample_size) &&
      params$sample_size != lava_sample_n) {
    if (verbose) {
      message(sprintf(
        "样本量提醒: GWAS 提供 n=%s, 基因型行数 n=%s；LAVAKnock 将使用后者以避免矩阵维度不一致 (GhostKnockoff 仍使用 n=%s)",
        params$sample_size, lava_sample_n, params$sample_size
      ))
    }
  }

  v_ids <- paste0(variant_info$chr, ":", variant_info$pos_bp)
  matched_rows <- match(colnames(genotype_data), v_ids)
  matched_rows <- matched_rows[!is.na(matched_rows)]
  if (length(matched_rows) != ncol(genotype_data)) {
    stop('变异信息与基因型列无法一一匹配')
  }
  variant_info <- variant_info[matched_rows, , drop = FALSE]
  variant_pos_all <- if ("pos_bp" %in% colnames(variant_info)) {
    variant_info$pos_bp
  } else if ("pos" %in% colnames(variant_info)) {
    variant_info$pos
  } else {
    variant_info[[which(grepl('pos', names(variant_info), ignore.case = TRUE))[1]]]
  }

  # 单变量分析
  if (params$run_mode %in% c("full", "lavaknock_only")) {
    locus_start <- min(variant_pos_all, na.rm = TRUE)
    locus_end <- max(variant_pos_all, na.rm = TRUE)

    univ_result <- run_univariate_analysis(
      genotype_data, z_mat, chr, locus_start, locus_end,
      lava_sample_n, params$prune_threshold, verbose
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

      window_variants <- lapply(windows$variant_indices, function(x) {
        if (is.list(x)) as.integer(x[[1]]) else as.integer(x)
      })

      window_results <- data.frame()
      result_list <- vector("list", length = nrow(windows))

      for (w_idx in seq_len(nrow(windows))) {
        window_cols <- window_variants[[w_idx]]
        window_cols <- window_cols[window_cols %in% seq_len(ncol(genotype_data))]
        if (length(window_cols) < 5) next

        sub_geno <- genotype_data[, window_cols, drop = FALSE]
        sub_z <- z_mat[window_cols, , drop = FALSE]

        knockoff_scores <- generate_ghostknockoff_scores(
          genotype_matrix = sub_geno,
          zscore_matrix = sub_z,
          ld_threshold = params$ld_threshold,
          n_samples = params$sample_size,
          n_knockoffs = params$knockoffs,
          variant_positions = variant_pos_all[window_cols],
          verbose = verbose
        )

        if (is.null(knockoff_scores) || is.null(knockoff_scores$index)) next

        ko_idx <- as.integer(knockoff_scores$index)
        ko_idx <- ko_idx[ko_idx >= 1 & ko_idx <= ncol(sub_geno)]
        if (length(ko_idx) < 2) next

        window_genotype <- sub_geno[, ko_idx, drop = FALSE]
        zscore_pheno1_window <- knockoff_scores$pheno1
        zscore_pheno2_window <- knockoff_scores$pheno2

        window <- windows[w_idx, ]

        bivar_result <- run_bivariate_analysis(
          window_genotype,
          zscore_pheno1_window,
          zscore_pheno2_window,
          chr = window$chr,
          window_start = window$start,
          window_end = window$end,
          n_samples = lava_sample_n,
          prune_thresh = params$prune_threshold,
          verbose = verbose
        )

        if (!is.null(bivar_result)) {
          result_list[[w_idx]] <- bivar_result
        }
      }

      result_list <- Filter(Negate(is.null), result_list)

      if (length(result_list) > 0) {
        window_results <- do.call(rbind, result_list)
        results$knockoff_results <- run_knockoff_filter(
          window_results, params$fdr, params$knockoffs,
          verbose = verbose
        )

        if (!is.null(results$knockoff_results)) {
          ko_res <- results$knockoff_results
          window_ids <- paste(window_results$chr,
                              window_results$window.start,
                              window_results$window.end,
                              sep = ":")

          W_vec <- ko_res$W
          if (!is.null(W_vec)) {
            if (is.null(names(W_vec))) {
              names(W_vec) <- window_ids[seq_along(W_vec)]
            }
            window_results$W <- as.numeric(W_vec[match(window_ids, names(W_vec))])
          }

          Q_vec <- ko_res$Qvalue
          if (!is.null(Q_vec)) {
            if (is.null(names(Q_vec))) {
              names(Q_vec) <- window_ids[seq_along(Q_vec)]
            }
            window_results$Qvalue <- as.numeric(Q_vec[match(window_ids, names(Q_vec))])
          }

          thr_val <- if (!is.null(ko_res$W.threshold)) as.numeric(ko_res$W.threshold) else NA_real_
          window_results$W_threshold <- thr_val
          window_results$selected <- ifelse(!is.na(window_results$W) & !is.na(thr_val),
                                            window_results$W >= thr_val, FALSE)
        }

        if (!'W' %in% names(window_results)) window_results$W <- NA_real_
        if (!'Qvalue' %in% names(window_results)) window_results$Qvalue <- NA_real_
        if (!'W_threshold' %in% names(window_results)) window_results$W_threshold <- NA_real_
        if (!'selected' %in% names(window_results)) window_results$selected <- FALSE

        results$bivariate_results <- window_results
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

  G <- as.matrix(genotype_matrix)
  if (nrow(G) < 2 || ncol(G) < 2) return(NULL)

  Z1 <- as.data.frame(zscore_pheno1)
  Z2 <- as.data.frame(zscore_pheno2)
  if (nrow(Z1) != ncol(G) || nrow(Z2) != ncol(G)) return(NULL)

  colnames(Z1)[1] <- "org"
  colnames(Z2)[1] <- "org"
  if (ncol(Z1) < 2 || ncol(Z2) < 2) return(NULL)
  knock_cols <- paste0("knock", seq_len(ncol(Z1) - 1))
  colnames(Z1)[-1] <- knock_cols
  colnames(Z2)[-1] <- knock_cols

  keep_cols <- which(apply(G, 2, var, na.rm = TRUE) > 0)
  if (length(keep_cols) < 2) return(NULL)
  G <- G[, keep_cols, drop = FALSE]
  Z1 <- Z1[keep_cols, , drop = FALSE]
  Z2 <- Z2[keep_cols, , drop = FALSE]

  if (ncol(G) < 2 || qr(G)$rank < ncol(G)) {
    if (verbose) message("跳过窗口: 基因型矩阵秩不足导致Sigma非正定")
    return(NULL)
  }

  effective_n <- n_samples
  geno_n <- nrow(G)
  if (!is.na(geno_n) && geno_n > 1 && geno_n != n_samples) {
    if (verbose) {
      message(sprintf(
        "  注意: 提供的样本量(n=%s)与基因型行数(n=%s)不一致，双变量分析将使用后者以匹配矩阵维度",
        n_samples, geno_n
      ))
    }
    effective_n <- geno_n
  }

  res <- tryCatch(
    LAVAKnock::LAVAKnock_bivariate(
      Zscore_pheno1_window = Z1,
      Zscore_pheno2_window = Z2,
      G_window = G,
      chr = chr,
      window_start = window_start,
      window_end = window_end,
      n = effective_n,
      prune.thresh = prune_thresh
    ),
    error = function(e) {
      if (verbose) message("LAVAKnock_bivariate失败: ", e$message)
      NULL
    }
  )

  res
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
  should_plot <- exists("plot_manhattan")

  if (is.null(opts$info) && !is.null(opts$gwas1)) {
    info_df <- as.data.frame(py_build_info_from_gwas(opts$gwas1), stringsAsFactors = FALSE)
    if (is.null(info_df) || nrow(info_df) == 0) {
      stop("自动生成Info失败: gwas1无有效CHR/POS")
    }
    info_path <- tempfile(paste0("auto_info_", gsub("[^A-Za-z0-9_]+", "_", tools::file_path_sans_ext(basename(opts$gwas1))), "_"), fileext = ".csv")
    readr::write_csv(info_df, info_path)
    opts$info <- info_path
    auto_temp_paths <- c(auto_temp_paths, info_path)
  }

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
    opts$info_label <- safe_base
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

  if (is.null(opts$info_label) || !nzchar(opts$info_label)) {
    opts$info_label <- gsub('[^A-Za-z0-9_]+', '_', tools::file_path_sans_ext(basename(opts$info)))
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

  gwas_primary <- NULL
  gwas_secondary <- NULL

  if (!is.null(opts$gwas1)) {
    gwas_primary <- merge_gwas_with_info(info, as.data.frame(py_load_gwas_table(opts$gwas1), stringsAsFactors = FALSE))
    if (is.null(gwas_primary) || nrow(gwas_primary) == 0) stop('变异信息与GWAS1在CHR+POS上没有交集')
    if (!is.null(opts$gwas2)) {
      gwas_secondary <- merge_gwas_with_info(info, as.data.frame(py_load_gwas_table(opts$gwas2), stringsAsFactors = FALSE))
      if (is.null(gwas_secondary) || nrow(gwas_secondary) == 0) stop('变异信息与GWAS2在CHR+POS上没有交集')
    }
  } else if (!is.null(opts$multi_gwas) && !is.null(opts$zcols)) {
    mg <- py_load_multi_gwas_table(opts$multi_gwas, opts$zcols)
    merged <- merge(info, as.data.frame(mg$data, stringsAsFactors = FALSE),
                    by.x = c('chr', 'pos_bp'), by.y = c('CHR', 'POS'),
                    all = FALSE, sort = FALSE)
    if (nrow(merged) == 0) stop('变异信息与multi_gwas无法匹配')
    zc <- trimws(as.character(mg$zcols))
    if (length(zc) < 1) stop('未在 multi_gwas 中找到指定的Z列')
    gwas_primary <- merged[, c(names(info), zc[1]), drop = FALSE]
    names(gwas_primary)[ncol(gwas_primary)] <- "Z"
    if (length(zc) >= 2) {
      gwas_secondary <- merged[, c(names(info), zc[2]), drop = FALSE]
      names(gwas_secondary)[ncol(gwas_secondary)] <- "Z"
    }
  } else {
    stop('未提供 GWAS 输入。请使用 --gwas1/--gwas2 或 --multi_gwas 配合 --zcols')
  }

  aligned_gwas <- py_align_dual_traits(gwas_primary, gwas_secondary)
  if (nrow(aligned_gwas) == 0) {
    stop('未能获得有效的GWAS记录，请检查输入文件是否与变异信息匹配')
  }

  gwas_data <- data.frame(
    rsid = aligned_gwas$rsid,
    chr = aligned_gwas$chr,
    pos_bp = aligned_gwas$pos_bp,
    zscore_pheno1 = aligned_gwas$zscore_pheno1,
    zscore_pheno2 = aligned_gwas$zscore_pheno2,
    stringsAsFactors = FALSE
  )
  attr(gwas_data, 'pheno1_name') <- opts$pheno1_name
  has_trait2 <- any(!is.na(gwas_data$zscore_pheno2))
  attr(gwas_data, 'pheno2_name') <- if (has_trait2) opts$pheno2_name else NA_character_

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
  variant_positions <- py_align_positions(
    variant_ids,
    info,
    fallback_positions = as.numeric(ali$positions)
  )
  if (length(variant_positions) != length(variant_ids) || all(is.na(variant_positions))) {
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

    prefix_base <- opts$info_label
    prefix <- paste0(prefix_base, '_', gsub('[^A-Za-z0-9_]+', '_', attr(gwas_data, 'pheno1_name')))

    if (should_plot && exists("plot_manhattan")) {
      plot_cols <- intersect(c('id', 'chr', 'pos', 'W', 'selected', 'W_threshold'), colnames(selection_table))
      plot_df <- selection_table[, plot_cols, drop = FALSE]
      if (all(c('chr', 'pos', 'W') %in% names(plot_df))) {
        thr_val <- NULL
        if ('W_threshold' %in% names(plot_df)) {
          thr_candidates <- unique(stats::na.omit(plot_df$W_threshold))
          if (length(thr_candidates) > 0) thr_val <- max(thr_candidates)
        }
        manhattan_path <- file.path(out_dir, paste0(prefix, '_manhattan.png'))
        try(plot_manhattan(
          plot_df,
          manhattan_path,
          title = paste0(prefix_base, ' (', attr(gwas_data, 'pheno1_name'), ')'),
          threshold = thr_val
        ), silent = TRUE)
      }
    }

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

    if (opts$verbose) cat(sprintf('已保存结果: %s\n', out_csv))

  } else if (opts$mode == 'correl') {
    params <- list(
      data_dir = dirname(opts$info),
      output_dir = opts$outdir,
      sample_size = opts$n,
      window_size = 5000000L,
      window_step = 5000000L,
      knockoffs = opts$knockoffs,
      fdr = opts$fdr,
      ld_threshold = 0.75,
      prune_threshold = 99,
      run_mode = 'ghostknockoff_only',
      ld_coord = opts$ld_coord,
      threads = opts$threads
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
      opts$info_label,
      '_', gsub('[^A-Za-z0-9_]+', '_', attr(gwas_data, 'pheno1_name')),
      '__', gsub('[^A-Za-z0-9_]+', '_', attr(gwas_data, 'pheno2_name'))
    )
    n_kn <- if (!is.null(opts$knockoffs) && opts$knockoffs > 0) opts$knockoffs else 5
    desired_kn_cols <- paste0('pval.knockoff', seq_len(n_kn))
    desired_cols <- c('chr', 'window.start', 'window.end', 'pval.orginal', desired_kn_cols, 'W', 'Q')

    if (is.null(res$bivariate_results) || nrow(res$bivariate_results) == 0) {
      empty_biv <- as.data.frame(setNames(replicate(length(desired_cols), numeric(), simplify = FALSE), desired_cols))
      readr::write_csv(empty_biv, file.path(out_dir, paste0(prefix, '_bivariate.csv')))
    } else {
      biv_export <- res$bivariate_results
      if (!'W' %in% colnames(biv_export)) biv_export$W <- NA_real_
      if (!'Qvalue' %in% colnames(biv_export)) biv_export$Qvalue <- NA_real_
      for (col in desired_kn_cols) {
        if (!col %in% colnames(biv_export)) {
          biv_export[[col]] <- NA_real_
        }
      }
      biv_export$Q <- as.numeric(biv_export$Qvalue)
      biv_export$Qvalue <- NULL
      missing_cols <- setdiff(desired_cols, colnames(biv_export))
      for (col in missing_cols) {
        biv_export[[col]] <- if (col %in% c('chr')) NA_integer_ else NA_real_
      }
      biv_export$chr <- as.integer(biv_export$chr)
      biv_export$window.start <- as.numeric(biv_export$window.start)
      biv_export$window.end <- as.numeric(biv_export$window.end)
      biv_export$pval.orginal <- as.numeric(biv_export$pval.orginal)
      for (col in desired_kn_cols) {
        biv_export[[col]] <- as.numeric(biv_export[[col]])
      }
      biv_export$W <- as.numeric(biv_export$W)
      biv_export$Q <- as.numeric(biv_export$Q)
      biv_export <- biv_export[, desired_cols, drop = FALSE]
      readr::write_csv(biv_export, file.path(out_dir, paste0(prefix, '_bivariate.csv')))
    }

    if (!is.null(res$knockoff_results) && !is.null(res$knockoff_results$window_sign) && length(res$knockoff_results$window_sign) > 0) {
      ko <- res$knockoff_results
      win_split <- strsplit(ko$window_sign, ':', fixed = TRUE)
      chr_vec <- suppressWarnings(as.integer(vapply(win_split, function(x) x[1], character(1))))
      start_vec <- suppressWarnings(as.numeric(vapply(win_split, function(x) ifelse(length(x) >= 2, x[2], NA_character_), character(1))))
      end_vec <- suppressWarnings(as.numeric(vapply(win_split, function(x) ifelse(length(x) >= 3, x[3], NA_character_), character(1))))
      W_vals <- ko$W[match(ko$window_sign, names(ko$W))]
      Q_vals <- ko$Qvalue[match(ko$window_sign, names(ko$Qvalue))]
      sig <- data.frame(
        chr = chr_vec,
        window.start = start_vec,
        window.end = end_vec,
        W = as.numeric(W_vals),
        Q = as.numeric(Q_vals),
        W_threshold = if (!is.null(ko$W.threshold)) as.numeric(ko$W.threshold) else NA_real_,
        stringsAsFactors = FALSE
      )
      sig <- sig[!is.na(sig$chr) & !is.na(sig$window.start) & !is.na(sig$window.end), , drop = FALSE]
      if (nrow(sig) > 0) sig <- unique(sig)
      readr::write_csv(sig, file.path(out_dir, paste0(prefix, '_significant_windows.csv')))
    } else {
      empty_sig <- data.frame(chr = integer(), window.start = numeric(), window.end = numeric(), W = numeric(), Q = numeric(), W_threshold = numeric(), stringsAsFactors = FALSE)
      readr::write_csv(empty_sig, file.path(out_dir, paste0(prefix, '_significant_windows.csv')))
    }

    if (should_plot && exists("plot_manhattan") && !is.null(res$bivariate_results) && nrow(res$bivariate_results) > 0) {
      plot_df <- res$bivariate_results
      if (!'pos' %in% colnames(plot_df) || all(is.na(plot_df$pos))) {
        plot_df$pos <- floor((as.numeric(plot_df$window.start) + as.numeric(plot_df$window.end)) / 2)
      }
      needed_cols <- intersect(c('id', 'chr', 'pos', 'W', 'selected', 'W_threshold'), colnames(plot_df))
      plot_df <- plot_df[, needed_cols, drop = FALSE]
      if (all(c('chr', 'pos', 'W') %in% names(plot_df))) {
        thr_val <- NULL
        if ('W_threshold' %in% names(plot_df)) {
          thr_candidates <- unique(stats::na.omit(plot_df$W_threshold))
          if (length(thr_candidates) > 0) thr_val <- max(thr_candidates)
        }
        manhattan_path <- file.path(out_dir, paste0(prefix, '_manhattan.png'))
        try(plot_manhattan(
          plot_df,
          manhattan_path,
          title = paste0(prefix, ' (W statistic)'),
          threshold = thr_val
        ), silent = TRUE)
      }
    }

    if (opts$verbose) cat('相关性分析完成并已保存结果\n')
  } else {
    stop('未知的mode。请使用 select 或 correl')
  }

  list(auto_temp_paths = auto_temp_paths)
}
