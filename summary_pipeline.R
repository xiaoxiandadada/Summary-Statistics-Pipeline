#!/usr/bin/env Rscript

# =============================================================================
# Pipeline - 主分析流程
# =============================================================================
# 包括GhostKnockoff生成和LAVA-Knock局部遗传相关性分析

suppressPackageStartupMessages({
  library(Matrix)
  library(MASS)
  library(corpcor)
  library(parallel)
  library(doParallel)
  library(foreach)
  library(data.table)
  library(dplyr)
  library(readr)
  library(stringr)
})

# 尝试加载必需的包
required_packages <- c("SKAT", "SPAtest", "CompQuadForm", "GhostKnockoff",
                      "irlba", "matrixsampling", "LAVAKnock", "snpStats")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste("包", pkg, "未安装。请先运行 install_packages.R"))
  }
}

library(SKAT)
library(SPAtest)
library(CompQuadForm)
library(GhostKnockoff)
library(irlba)
library(matrixsampling)
library(LAVAKnock)
library(snpStats)

# 源配置脚本
source("configure_parameters.R")

# Python 加速（直接从 accelerator.py 载入，可选）
py_env <- new.env(parent = emptyenv())
py_ok <- FALSE
if (requireNamespace("reticulate", quietly = TRUE)) {
  try({ reticulate::source_python("accelerator.py", envir = py_env); py_ok <- TRUE }, silent = TRUE)
}

py_fast_correlation <- function(X, verbose = FALSE) {
  if (isTRUE(py_ok) && !is.null(py_env$fast_correlation_matrix)) {
    return(py_env$fast_correlation_matrix(X))
  }
  if (verbose) message("Python accel unavailable; using R cor()")
  return(stats::cor(X))
}

py_ld_pruning <- function(corr, threshold = 0.75, verbose = FALSE) {
  if (isTRUE(py_ok) && !is.null(py_env$ld_pruning)) {
    idx0 <- py_env$ld_pruning(corr, threshold)
    # convert 0-based numpy indices to 1-based R indices
    return(as.integer(idx0) + 1L)
  }
  if (verbose) message("Python accel unavailable; using hclust pruning")
  Sigma.distance <- stats::as.dist(1 - abs(corr))
  fit <- stats::hclust(Sigma.distance, method = "single")
  clusters <- stats::cutree(fit, h = 1 - threshold)
  return(match(unique(clusters), clusters))
}

# =============================================================================
# 数据加载和预处理函数
# =============================================================================

standardize_variant_info <- function(info) {
  if (is.null(info) || nrow(info) == 0) {
    return(info)
  }
  df <- as.data.frame(info, stringsAsFactors = FALSE)
  colnames(df) <- tolower(colnames(df))

  chr_col <- intersect(c('chr', 'chrom', 'chromosome'), colnames(df))
  if (length(chr_col) == 0) {
    stop('变异信息缺少染色体列 (chr/chrom/chromosome)')
  }
  pos_col <- intersect(c('pos_bp', 'pos', 'position', 'bp', 'poshg38', 'poshg19'), colnames(df))
  if (length(pos_col) == 0) {
    stop('变异信息缺少位点列 (pos, position, pos_bp, poshg38, poshg19)')
  }

  raw_chr <- df[[chr_col[1]]]
  chr_str <- toupper(as.character(raw_chr))
  chr_str <- sub('^CHR', '', chr_str, ignore.case = TRUE)
  chr_num <- suppressWarnings(as.integer(chr_str))
  chr_num[is.na(chr_num) & chr_str == 'X'] <- 23L
  chr_num[is.na(chr_num) & chr_str == 'Y'] <- 24L
  chr_num[is.na(chr_num) & chr_str %in% c('M', 'MT', 'MITO')] <- 25L
  df$chr <- chr_num
  df$pos_bp <- suppressWarnings(as.numeric(df[[pos_col[1]]]))

  if (!'rsid' %in% colnames(df)) {
    df$rsid <- paste0('chr', df$chr, ':', df$pos_bp)
  }
  df$rsid <- as.character(df$rsid)
  df <- df[!is.na(df$chr) & !is.na(df$pos_bp), , drop = FALSE]
  df <- df[order(df$chr, df$pos_bp), , drop = FALSE]
  unique(df)
}

parse_trait_option <- function(x) {
  if (is.null(x) || is.na(x) || x == '') return(NULL)
  trimws(unlist(strsplit(x, ',')))
}

normalize_gwas_table <- function(gwas_df) {
  if (is.null(gwas_df) || nrow(gwas_df) == 0) {
    return(gwas_df)
  }
  df <- as.data.frame(gwas_df, stringsAsFactors = FALSE)
  colnames(df) <- tolower(colnames(df))
  if (!all(c('chr', 'pos') %in% colnames(df))) {
    id_col <- intersect(c('rsid', 'snp', 'id', 'variant_id'), colnames(df))
    if (length(id_col) > 0) {
      id_vals <- as.character(df[[id_col[1]]])
      parts <- strsplit(id_vals, ':', fixed = TRUE)
      chr_raw <- vapply(parts, function(x) if (length(x) >= 1) x[1] else NA_character_, character(1), USE.NAMES = FALSE)
      pos_raw <- vapply(parts, function(x) if (length(x) >= 2) x[2] else NA_character_, character(1), USE.NAMES = FALSE)
      chr_clean <- toupper(chr_raw)
      chr_clean <- sub('^CHR', '', chr_clean, ignore.case = TRUE)
      chr_num <- suppressWarnings(as.integer(chr_clean))
      chr_num[is.na(chr_num) & chr_clean == 'X'] <- 23L
      chr_num[is.na(chr_num) & chr_clean == 'Y'] <- 24L
      chr_num[is.na(chr_num) & chr_clean %in% c('M', 'MT', 'MITO')] <- 25L
      pos_num <- suppressWarnings(as.numeric(gsub('[^0-9]', '', pos_raw)))
      df$chr <- chr_num
      df$pos <- pos_num
    }
  }
  if (!'chr' %in% colnames(df)) {
    alt_chr <- intersect(c('chrom', 'chromosome'), colnames(df))
    if (length(alt_chr) > 0) df$chr <- suppressWarnings(as.integer(df[[alt_chr[1]]]))
  }
  if (!'pos' %in% colnames(df)) {
    alt_pos <- intersect(c('position', 'bp', 'pos_bp'), colnames(df))
    if (length(alt_pos) > 0) df$pos <- suppressWarnings(as.numeric(df[[alt_pos[1]]]))
  }
  df$chr <- suppressWarnings(as.integer(df$chr))
  df$pos <- suppressWarnings(as.numeric(df$pos))
  df <- df[!is.na(df$chr) & !is.na(df$pos), , drop = FALSE]
  colnames(df) <- toupper(colnames(df))
  df
}

detect_single_z_column <- function(df) {
  if ('Z' %in% colnames(df)) return('Z')
  numeric_cols <- names(df)[vapply(df, is.numeric, logical(1))]
  numeric_cols <- setdiff(numeric_cols, c('CHR', 'POS', 'BP', 'N', 'SE', 'P', 'BETA'))
  if (length(numeric_cols) == 1) {
    return(numeric_cols[1])
  }
  NULL
}

auto_prepare_inputs <- function(zscore_path, panel_prefix, trait_cols = NULL, verbose = FALSE) {
  if (!file.exists(zscore_path)) stop(paste('zscore文件不存在:', zscore_path))
  bim <- paste0(panel_prefix, '.bim')
  if (!file.exists(bim)) stop(paste('参考面板缺少 BIM 文件:', bim))

  if (verbose) cat('读取zscore文件...\n')
  z_dt <- data.table::fread(zscore_path, data.table = FALSE)
  cols <- colnames(z_dt)
  if (!all(c('CHR', 'POS') %in% cols)) {
    id_col <- intersect(c('rsid', 'RSID', 'SNP', 'ID', 'variant_id'), cols)
    if (length(id_col) == 0) stop('zscore文件缺少CHR/POS列，也无法从rsid推断')
    id_vals <- as.character(z_dt[[id_col[1]]])
    pieces <- strsplit(id_vals, ':', fixed = TRUE)
    chr_raw <- vapply(pieces, function(x) if (length(x) >= 1) x[1] else NA_character_, character(1), USE.NAMES = FALSE)
    pos_raw <- vapply(pieces, function(x) if (length(x) >= 2) x[2] else NA_character_, character(1), USE.NAMES = FALSE)
    chr_clean <- toupper(chr_raw)
    chr_clean <- sub('^CHR', '', chr_clean, ignore.case = TRUE)
    chr_num <- suppressWarnings(as.integer(chr_clean))
    chr_num[is.na(chr_num) & chr_clean == 'X'] <- 23L
    chr_num[is.na(chr_num) & chr_clean == 'Y'] <- 24L
    chr_num[is.na(chr_num) & chr_clean %in% c('M', 'MT', 'MITO')] <- 25L
    pos_num <- suppressWarnings(as.integer(gsub('[^0-9]', '', pos_raw)))
    z_dt$CHR <- chr_num
    z_dt$POS <- pos_num
  }
  z_dt$CHR <- suppressWarnings(as.integer(z_dt$CHR))
  z_dt$POS <- suppressWarnings(as.integer(z_dt$POS))
  if (!all(c('CHR', 'POS') %in% colnames(z_dt))) stop('zscore文件仍缺少CHR/POS列')

  numeric_cols <- names(z_dt)[vapply(z_dt, is.numeric, logical(1))]
  numeric_cols <- setdiff(numeric_cols, c('CHR', 'POS', 'BP', 'N', 'SE', 'P', 'BETA'))
  if (!is.null(trait_cols)) {
    trait_cols <- trimws(unlist(strsplit(trait_cols, ',')))
    if (!all(trait_cols %in% numeric_cols)) {
      missing_cols <- setdiff(trait_cols, numeric_cols)
      stop(paste('zscore文件缺少指定的trait列:', paste(missing_cols, collapse = ', ')))
    }
  } else {
    trait_cols <- numeric_cols
  }
  trait_cols <- trait_cols[seq_len(min(length(trait_cols), 2))]
  if (length(trait_cols) == 0) stop('无法确定zscore中的trait列')

  for (col in trait_cols) {
    z_dt[[col]] <- as.numeric(z_dt[[col]])
  }
  primary_col <- trait_cols[1]
  z_dt$Z <- as.numeric(z_dt[[primary_col]])
  z_dt <- z_dt[!is.na(z_dt$CHR) & !is.na(z_dt$POS) & !is.na(z_dt$Z), ]
  z_dt <- z_dt[!duplicated(z_dt[c('CHR', 'POS')]), ]
  if (nrow(z_dt) == 0) stop('zscore文件有效记录为0')

  if (verbose) cat('读取参考面板BIM并匹配...\n')
  map_dt <- data.table::fread(bim, col.names = c('CHR', 'RSID', 'CM', 'POS', 'A1', 'A2'), data.table = FALSE)
  map_dt$CHR <- suppressWarnings(as.integer(map_dt$CHR))
  merged <- merge(map_dt, z_dt, by = c('CHR', 'POS'))
  if (nrow(merged) == 0) stop('参考面板与zscore在CHR+POS上没有重叠')

  info_df <- merged[, c('RSID', 'CHR', 'POS')]
  colnames(info_df) <- c('rsid', 'chr', 'pos_bp')

  safe_base <- gsub('[^A-Za-z0-9_]+', '_', tools::file_path_sans_ext(basename(zscore_path)))
  info_path <- tempfile(paste0('auto_info_', safe_base, '_'), fileext = '.csv')
  data.table::fwrite(info_df, info_path)

  gwas_paths <- character()
  trait_names <- character()
  for (trait in trait_cols) {
    gwas_df <- merged[, c('CHR', 'POS', trait)]
    colnames(gwas_df) <- c('CHR', 'POS', 'Z')
    gwas_path <- tempfile(paste0('auto_gwas_', safe_base, '_', trait, '_'), fileext = '.tsv')
    data.table::fwrite(gwas_df, gwas_path, sep = '\t')
    gwas_paths <- c(gwas_paths, gwas_path)
    trait_names <- c(trait_names, trait)
  }
  if (verbose) cat(sprintf('自动生成Info与GWAS：%d 个变异\n', nrow(info_df)))
  list(info = info_path, gwas = list(paths = gwas_paths, traits = trait_names))
}

load_gwas_for_info <- function(gwas_file, variant_info) {
  if (is.null(gwas_file)) return(NULL)
  if (!file.exists(gwas_file)) stop(paste('GWAS文件不存在:', gwas_file))
  gwas <- read.table(gwas_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  gwas <- normalize_gwas_table(gwas)
  if (!all(c('CHR', 'POS', 'Z') %in% colnames(gwas))) {
    stop('GWAS文件需要包含 CHR, POS, Z 列或可推断的位置信息')
  }
  gwas <- gwas[!is.na(gwas$CHR) & !is.na(gwas$POS), , drop = FALSE]
  merged <- merge(variant_info, gwas, by.x = c('chr', 'pos_bp'), by.y = c('CHR', 'POS'))
  if (nrow(merged) == 0) {
    stop('变异信息与GWAS在CHR+POS上没有交集')
  }
  merged
}

load_multi_gwas_for_info <- function(multi_file, zcols, variant_info) {
  if (is.null(multi_file) || is.null(zcols)) return(NULL)
  if (!file.exists(multi_file)) stop(paste('multi_gwas 文件不存在:', multi_file))
  gwas <- read.table(multi_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  gwas <- normalize_gwas_table(gwas)
  if (!all(c('CHR', 'POS') %in% colnames(gwas))) {
    stop('multi_gwas 缺少CHR/POS信息，且无法从rsid推断')
  }
  gwas <- gwas[!is.na(gwas$CHR) & !is.na(gwas$POS), ]
  gwas <- gwas[!duplicated(gwas[c('CHR', 'POS')]), ]
  zlist <- unlist(strsplit(zcols, ','))
  zlist <- trimws(zlist)
  if (length(zlist) < 1) stop('--zcols 至少需要提供一个列名')
  if (!all(zlist %in% colnames(gwas))) {
    stop(paste('--zcols 列不存在:', paste(setdiff(zlist, colnames(gwas)), collapse = ', ')))
  }
  m <- merge(variant_info, gwas, by.x = c('chr', 'pos_bp'), by.y = c('CHR', 'POS'))
  if (nrow(m) == 0) stop('变异信息与multi_gwas无法匹配')
  list(data = m, zcols = zlist)
}

#' 读取遗传变异信息文件
load_variant_info <- function(info_file) {
  tryCatch({
    if (file.exists(info_file)) {
      info <- read_csv(info_file, show_col_types = FALSE)
      info <- standardize_variant_info(info)
      return(info)
    } else {
      warning(paste("文件不存在:", info_file))
      return(NULL)
    }
  }, error = function(e) {
    warning(paste("读取文件出错:", info_file, "-", e$message))
    return(NULL)
  })
}

#' 加载真实GWAS汇总统计数据
load_gwas_summary_stats <- function(gwas_file, variant_info = NULL) {
  tryCatch({
    if (!file.exists(gwas_file)) {
      stop(paste("GWAS文件不存在:", gwas_file))
    }
    gwas_data <- read.table(gwas_file, header = TRUE, stringsAsFactors = FALSE)
    required_cols <- c("CHR", "POS", "Z")
    if (!all(required_cols %in% colnames(gwas_data))) {
      stop(paste("GWAS文件缺少必需列:", paste(setdiff(required_cols, colnames(gwas_data)), collapse = ", ")))
    }
    gwas_data$CHR <- suppressWarnings(as.integer(gwas_data$CHR))
    gwas_data$POS <- suppressWarnings(as.numeric(gwas_data$POS))
    gwas_data <- gwas_data[!is.na(gwas_data$CHR) & !is.na(gwas_data$POS), , drop = FALSE]

    if (!is.null(variant_info)) {
      matched <- merge(variant_info, gwas_data, by.x = c("chr", "pos_bp"), by.y = c("CHR", "POS"))
      if (nrow(matched) == 0) {
        stop("变异信息与GWAS在CHR+POS上没有交集")
      }
      result <- matched[, c("rsid", "chr", "pos_bp"), drop = FALSE]
      result$zscore_pheno1 <- matched$Z
      if ("Z2" %in% colnames(matched)) {
        result$zscore_pheno2 <- matched$Z2
      } else if ("Z_pheno2" %in% colnames(matched)) {
        result$zscore_pheno2 <- matched$Z_pheno2
      } else if ("zscore_pheno2" %in% colnames(matched)) {
        result$zscore_pheno2 <- matched$zscore_pheno2
      } else {
        result$zscore_pheno2 <- rep(NA_real_, nrow(result))
      }
      return(result)
    }
    gwas_data <- gwas_data[order(gwas_data$CHR, gwas_data$POS), , drop = FALSE]
    return(gwas_data)
  }, error = function(e) {
    warning(paste("加载GWAS数据失败:", e$message))
    return(NULL)
  })
}

load_ld_blocks <- function(build, base_dir = '.') {
  if (is.null(build) || is.na(build) || build == '') {
    stop('必须指定LD block坐标版本 (GRCh37 或 GRCh38)')
  }
  key <- toupper(trimws(build))
  file_map <- c(
    GRCH37 = 'LAVA_s2500_m25_f1_w200.blocks',
    GRCH38 = 'deCODE_EUR_LD_blocks.bed'
  )
  if (!key %in% names(file_map)) {
    stop("未知的LD坐标版本: ", build, ". 仅支持 GRCh37 或 GRCh38")
  }
  candidate <- file_map[[key]]
  paths <- c(candidate, file.path(base_dir, candidate))
  path <- paths[file.exists(paths)][1]
  if (is.na(path)) {
    stop('无法找到LD blocks文件: ', candidate)
  }
  blocks <- data.table::fread(path, data.table = FALSE)
  if (nrow(blocks) == 0) {
    stop('LD blocks文件为空: ', path)
  }
  colnames(blocks) <- tolower(colnames(blocks))
  if (!'chr' %in% colnames(blocks)) {
    stop('LD blocks文件缺少chr列: ', path)
  }
  if (!('stop' %in% colnames(blocks) || 'end' %in% colnames(blocks))) {
    stop('LD blocks文件缺少stop/end列: ', path)
  }
  blocks$chr <- gsub('^chr', '', blocks$chr, ignore.case = TRUE)
  blocks$chr <- suppressWarnings(as.integer(blocks$chr))
  if ('stop' %in% colnames(blocks)) {
    blocks$end <- blocks$stop
  }
  blocks$start <- suppressWarnings(as.numeric(blocks$start))
  blocks$end <- suppressWarnings(as.numeric(blocks$end))
  blocks <- blocks[!is.na(blocks$chr) & !is.na(blocks$start) & !is.na(blocks$end), , drop = FALSE]
  blocks <- blocks[order(blocks$chr, blocks$start, blocks$end), , drop = FALSE]
  blocks$source_build <- key
  blocks
}

load_genotype_rds <- function(path) {
  obj <- readRDS(path)
  if (!is.matrix(obj)) stop('RDS未包含matrix对象')
  storage.mode(obj) <- 'numeric'
  obj
}

load_genotype_csv <- function(path, fmt = c('samples_by_snps', 'snps_by_samples')) {
  fmt <- match.arg(fmt)
  m <- suppressWarnings(as.data.frame(data.table::fread(path)))
  M <- as.matrix(m)
  storage.mode(M) <- 'numeric'
  if (fmt == 'snps_by_samples') {
    M <- t(M)
  }
  M
}

load_genotype_plink <- function(prefix, snp_ids = NULL, chrpos_ids = NULL) {
  if (!requireNamespace('snpStats', quietly = TRUE)) {
    stop("需要安装snpStats包以读取PLINK文件，请运行 install.packages('snpStats')")
  }
  bed <- paste0(prefix, '.bed')
  bim <- paste0(prefix, '.bim')
  fam <- paste0(prefix, '.fam')
  if (!file.exists(bed) || !file.exists(bim) || !file.exists(fam)) {
    stop(paste0('PLINK文件缺失，请检查前缀: ', prefix))
  }
  if (!is.null(snp_ids)) {
    snp_ids <- unique(stats::na.omit(snp_ids))
    snp_ids <- snp_ids[nchar(snp_ids) > 0]
  }
  if (!is.null(chrpos_ids)) {
    chrpos_ids <- unique(stats::na.omit(chrpos_ids))
    chrpos_ids <- chrpos_ids[nchar(chrpos_ids) > 0]
    if (length(chrpos_ids) > 0) {
      map_dt <- data.table::fread(bim, col.names = c('chr', 'rsid', 'cm', 'pos', 'a1', 'a2'), showProgress = FALSE)
      map_dt[, chrpos_key := paste0(chr, ':', pos)]
      matched_ids <- unique(map_dt[chrpos_key %in% chrpos_ids, rsid])
      snp_ids <- unique(c(snp_ids, matched_ids))
    }
  }
  if (!is.null(snp_ids) && length(snp_ids) == 0) {
    snp_ids <- NULL
  }
  plink <- snpStats::read.plink(bed, bim, fam, select.snps = snp_ids)
  geno <- plink$genotypes
  if (is.null(geno) || ncol(geno) == 0) {
    stop('PLINK文件未读取到任何SNP (可能过滤条件过严)')
  }
  G <- as(geno, 'numeric')
  storage.mode(G) <- 'numeric'
  if (anyNA(G)) {
    col_means <- colMeans(G, na.rm = TRUE)
    na_idx <- which(is.na(G), arr.ind = TRUE)
    if (length(na_idx) > 0) {
      G[na_idx] <- col_means[na_idx[, 2]]
      G[is.na(G)] <- 0
    }
  }
  map <- plink$map
  chrpos <- paste0(map$chromosome, ':', map$position)
  rsid <- map$snp.name
  colnames(G) <- ifelse(!is.na(rsid) & rsid != '', rsid, chrpos)
  attr(G, 'chrpos') <- chrpos
  attr(G, 'rsid') <- rsid
  list(matrix = G, map = map)
}

align_genotype_to_gwas <- function(G, gwas_df, prefer_rsid = TRUE) {
  if (is.null(colnames(G)) && is.null(attr(G, 'chrpos'))) {
    stop('基因型矩阵缺少列名或位置信息，无法对齐')
  }
  col_ids <- colnames(G)
  chrpos <- paste0(gwas_df$chr, ':', gwas_df$pos_bp)
  chrpos_attr <- attr(G, 'chrpos')
  rsid_attr <- attr(G, 'rsid')

  if (prefer_rsid && !is.null(gwas_df$rsid) && !is.null(col_ids)) {
    idx <- match(gwas_df$rsid, col_ids)
    matched <- which(!is.na(idx))
    if (length(matched) >= 2) {
      G2 <- G[, idx[matched], drop = FALSE]
      new_ids <- chrpos[matched]
      colnames(G2) <- new_ids
      if (!is.null(rsid_attr)) {
        attr(G2, 'rsid') <- rsid_attr[idx[matched]]
      }
      return(list(G = G2, matched = length(matched), ids = new_ids, pos = gwas_df$pos_bp[matched]))
    }
  }
  if (!is.null(col_ids)) {
    idx <- match(chrpos, col_ids)
    matched <- which(!is.na(idx))
    if (length(matched) >= 2) {
      G2 <- G[, idx[matched], drop = FALSE]
      new_ids <- chrpos[matched]
      colnames(G2) <- new_ids
      if (!is.null(rsid_attr)) {
        attr(G2, 'rsid') <- rsid_attr[idx[matched]]
      }
      return(list(G = G2, matched = length(matched), ids = new_ids, pos = gwas_df$pos_bp[matched]))
    }
  }
  if (!is.null(chrpos_attr)) {
    idx <- match(chrpos, chrpos_attr)
    matched <- which(!is.na(idx))
    if (length(matched) >= 2) {
      G2 <- G[, idx[matched], drop = FALSE]
      new_ids <- chrpos[matched]
      colnames(G2) <- new_ids
      if (!is.null(rsid_attr)) {
        attr(G2, 'rsid') <- rsid_attr[idx[matched]]
      }
      return(list(G = G2, matched = length(matched), ids = new_ids, pos = gwas_df$pos_bp[matched]))
    }
  }
  stop('基因型列名与GWAS变异无法对齐 (尝试rsid与chr:pos均失败)')
}

# =============================================================================
# 窗口管理函数
# =============================================================================

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
  valid_variant_info <- variant_info[valid_indices, ]
  
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

# =============================================================================
# GhostKnockoff 相关函数
# =============================================================================

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
                                        use_python_accel = FALSE,
                                        verbose = FALSE) {
  
  if (verbose) cat("开始生成GhostKnockoff分数...\n")
  
  # 数据预处理
  genotype_matrix[genotype_matrix < 0 | genotype_matrix > 2] <- 0
  
  # 计算MAF和MAC
  MAF <- colMeans(genotype_matrix) / 2
  MAC <- colSums(genotype_matrix)
  s <- colMeans(genotype_matrix^2) - colMeans(genotype_matrix)^2
  
  # 过滤变异
  SNP.index <- which(MAF > 0 & MAC >= 25 & s != 0 & !is.na(MAF))
  
  if (length(SNP.index) < 5) {
    warning("过滤后变异数量过少")
    return(NULL)
  }
  
  genotype_matrix <- genotype_matrix[, SNP.index, drop = FALSE]
  zscore_matrix <- zscore_matrix[SNP.index, , drop = FALSE]
  keep_index <- SNP.index
  
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
  
  # LD过滤（支持Python加速）
  if (ncol(genotype_matrix) > 1) {
    cor.X <- if (use_python_accel) py_fast_correlation(genotype_matrix, verbose) else stats::cor(genotype_matrix)
    keep_cols <- py_ld_pruning(cor.X, threshold = ld_threshold, verbose = verbose)
    genotype_matrix <- genotype_matrix[, keep_cols, drop = FALSE]
    zscore_matrix <- zscore_matrix[keep_cols, , drop = FALSE]
    keep_index <- keep_index[keep_cols]
  }
  
  if (verbose) cat(paste("LD过滤后变异数量:", ncol(genotype_matrix), "\n"))
  
  # 计算收缩LD矩阵
  cor.G <- cor.shrink(genotype_matrix, verbose = FALSE)
  
  # 生成knockoff
  tryCatch({
    set.seed(12345)
    fit.prelim <- GhostKnockoff.prelim(cor.G, M = n_knockoffs, 
                                      method = 'asdp', 
                                      max.size = ncol(genotype_matrix))
    
    # 为每个表型生成knockoff
    knockoff_results <- list()
    
    for (pheno in seq_len(ncol(zscore_matrix))) {
      GK.stat <- GhostKnockoff.fit(as.matrix(zscore_matrix[, pheno]), 
                                  n.study = n_samples,
                                  fit.prelim = fit.prelim, 
                                  gamma = 1, 
                                  weight.study = NULL)
      
      # 组合原始和knockoff Z-scores
      combined_scores <- cbind(zscore_matrix[, pheno], GK.stat$GK.Zscore_k)
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

# =============================================================================
# LAVA-Knock 分析函数
# =============================================================================

#' 单变量遗传性检验
run_univariate_analysis <- function(genotype_matrix, zscore_matrix, 
                                   chr, locus_start, locus_end, 
                                   n_samples = 20000, prune_thresh = 99,
                                   verbose = FALSE) {
  
  if (verbose) cat("运行单变量遗传性分析...\n")
  
  tryCatch({
    result <- LAVAKnock_univariate(
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

# =============================================================================
# GhostKnockoff 变量选择（单表型）
# =============================================================================

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
                                        use_python_accel = FALSE,
                                        verbose = FALSE) {
  if (is.null(genotype_matrix) || is.null(z_vector)) {
    warning("genotype_matrix 或 z_vector 为空")
    return(NULL)
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
    use_python_accel = isTRUE(use_python_accel),
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
  res <- LAVAKnock(
    M = n_knockoffs,
    p0 = p0,
    p_ko = as.matrix(pko),
    fdr = fdr,
    window_id = variant_ids,
    Rej.Bound = max(20000, length(p0))
  )

  # 汇总为结果表
  out_tbl <- data.frame(
    id = variant_ids,
    pval.orginal = p0,
    pko,
    W = res$W,
    Qvalue = res$Qvalue,
    selected = res$W >= res$W.threshold,
    stringsAsFactors = FALSE
  )

  if (verbose) {
    cat(paste0("变量选择完成: 选中 ", sum(out_tbl$selected), "/", nrow(out_tbl), " 变异\n"))
  }

  list(
    pvals = list(p0 = p0, pko = pko),
    lavaknock = res,
    table = out_tbl
  )
}

run_block_selection <- function(genotype_matrix, z_vector, variant_ids, variant_chr, variant_pos,
                                blocks_df, n_samples, n_knockoffs, fdr, use_python_accel, verbose) {
  if (ncol(genotype_matrix) != length(z_vector) || length(z_vector) != length(variant_ids)) {
    stop('变异信息与基因型矩阵长度不一致')
  }

  variant_chr <- suppressWarnings(as.integer(variant_chr))
  variant_pos <- suppressWarnings(as.numeric(variant_pos))
  valid <- !is.na(variant_chr) & !is.na(variant_pos)
  if (sum(valid) < 5) {
    warning('可用于LD blocks的变异数量不足')
    return(NULL)
  }
  if (!all(valid)) {
    invalid <- which(!valid)
    genotype_matrix <- genotype_matrix[, valid, drop = FALSE]
    z_vector <- z_vector[valid]
    variant_ids <- variant_ids[valid]
    variant_chr <- variant_chr[valid]
    variant_pos <- variant_pos[valid]
    warning(sprintf('已忽略 %d 个缺少chr/pos的变异', length(invalid)))
  }

  ord <- order(variant_chr, variant_pos, na.last = TRUE)
  genotype_matrix <- genotype_matrix[, ord, drop = FALSE]
  z_vector <- z_vector[ord]
  variant_ids <- variant_ids[ord]
  variant_chr <- variant_chr[ord]
  variant_pos <- variant_pos[ord]

  blocks_df <- blocks_df[blocks_df$chr %in% unique(variant_chr), , drop = FALSE]
  blocks_df <- blocks_df[order(blocks_df$chr, blocks_df$start, blocks_df$end), , drop = FALSE]

  combined_tbl <- data.frame()
  chunk_summary <- data.frame(chunk_id = integer(), chr = integer(), start = numeric(), end = numeric(),
                              n_snps = integer(), selected = integer(), source = character(),
                              stringsAsFactors = FALSE)
  chunk_tables <- list()
  assigned <- rep(FALSE, length(variant_ids))

  for (i in seq_len(nrow(blocks_df))) {
    block_chr <- blocks_df$chr[i]
    block_start <- blocks_df$start[i]
    block_end <- blocks_df$end[i]
    idx <- which(!assigned & variant_chr == block_chr & variant_pos >= block_start & variant_pos <= block_end)
    if (length(idx) == 0) next
    if (length(idx) < 5) {
      next
    }
    chunk_id <- length(chunk_tables) + 1L
    sub_G <- genotype_matrix[, idx, drop = FALSE]
    sub_z <- z_vector[idx]
    sub_ids <- variant_ids[idx]
    if (verbose) {
      cat(sprintf('Chunk %d: %d SNPs (chr %d: %d-%d)\n', chunk_id, length(idx), block_chr, floor(block_start), floor(block_end)))
      cat(sprintf('  调用 GhostKnockoff (M=%d)\n', n_knockoffs))
    }
    res <- run_ghostknockoff_selection(
      genotype_matrix = sub_G,
      z_vector = sub_z,
      n_samples = n_samples,
      n_knockoffs = n_knockoffs,
      fdr = fdr,
      variant_ids = sub_ids,
      use_python_accel = use_python_accel,
      verbose = verbose
    )
    if (is.null(res)) next
    assigned[idx] <- TRUE
    tbl <- res$table
    tbl$chunk_id <- chunk_id
    tbl$chunk_chr <- block_chr
    tbl$chunk_start <- block_start
    tbl$chunk_end <- block_end
    tbl$chunk_size_bp <- block_end - block_start + 1
    tbl$partition_source <- 'ld_block'
    combined_tbl <- rbind(combined_tbl, tbl)
    chunk_summary <- rbind(chunk_summary, data.frame(
      chunk_id = chunk_id,
      chr = block_chr,
      start = block_start,
      end = block_end,
      n_snps = length(idx),
      selected = sum(tbl$selected),
      source = 'ld_block',
      stringsAsFactors = FALSE
    ))
    chunk_tables[[length(chunk_tables) + 1L]] <- list(
      id = chunk_id,
      chr = block_chr,
      start = block_start,
      end = block_end,
      source = 'ld_block',
      data = tbl
    )
  }

  if (!all(assigned)) {
    warning('部分变异不在预定义LD blocks中，按染色体范围回退切分')
    leftover_idx <- which(!assigned)
    leftover_chr <- unique(variant_chr[leftover_idx])
    for (chr_val in leftover_chr) {
      idx_chr <- leftover_idx[variant_chr[leftover_idx] == chr_val]
      if (length(idx_chr) < 5) {
        assigned[idx_chr] <- TRUE
        next
      }
      block_start <- min(variant_pos[idx_chr], na.rm = TRUE)
      block_end <- max(variant_pos[idx_chr], na.rm = TRUE)
      chunk_id <- length(chunk_tables) + 1L
      sub_G <- genotype_matrix[, idx_chr, drop = FALSE]
      sub_z <- z_vector[idx_chr]
      sub_ids <- variant_ids[idx_chr]
      if (verbose) {
        cat(sprintf('Chunk %d: %d SNPs (chr %d: %d-%d) [fallback]\n', chunk_id, length(idx_chr), chr_val, floor(block_start), floor(block_end)))
        cat(sprintf('  调用 GhostKnockoff (M=%d)\n', n_knockoffs))
      }
      res <- run_ghostknockoff_selection(
        genotype_matrix = sub_G,
        z_vector = sub_z,
        n_samples = n_samples,
        n_knockoffs = n_knockoffs,
        fdr = fdr,
        variant_ids = sub_ids,
        use_python_accel = use_python_accel,
        verbose = verbose
      )
      if (is.null(res)) {
        assigned[idx_chr] <- TRUE
        next
      }
      assigned[idx_chr] <- TRUE
      tbl <- res$table
      tbl$chunk_id <- chunk_id
      tbl$chunk_chr <- chr_val
      tbl$chunk_start <- block_start
      tbl$chunk_end <- block_end
      tbl$chunk_size_bp <- block_end - block_start + 1
      tbl$partition_source <- 'fallback'
      combined_tbl <- rbind(combined_tbl, tbl)
      chunk_summary <- rbind(chunk_summary, data.frame(
        chunk_id = chunk_id,
        chr = chr_val,
        start = block_start,
        end = block_end,
        n_snps = length(idx_chr),
        selected = sum(tbl$selected),
        source = 'fallback',
        stringsAsFactors = FALSE
      ))
      chunk_tables[[length(chunk_tables) + 1L]] <- list(
        id = chunk_id,
        chr = chr_val,
        start = block_start,
        end = block_end,
        source = 'fallback',
        data = tbl
      )
    }
  }

  if (nrow(combined_tbl) == 0) {
    warning('LD 切块后未得到任何结果')
    return(NULL)
  }
  ord_idx <- match(combined_tbl$id, variant_ids)
  combined_tbl <- combined_tbl[order(combined_tbl$chunk_id, ord_idx), ]
  list(table = combined_tbl, chunks = chunk_summary, chunk_tables = chunk_tables)
}

execute_pipeline <- function(opts) {
  auto_temp_paths <- character()
  auto_trait_names <- NULL

  if (!is.null(opts$panel) && is.null(opts$ref_plink)) {
    opts$ref_plink <- opts$panel
  }
  panel_default <- if (!is.null(opts$panel)) opts$panel else if (!is.null(opts$ref_plink)) opts$ref_plink else file.path('g1000_eur', 'g1000_eur')

  trait_option_str <- NULL
  if (!is.null(opts$zscore_traits)) {
    parsed_traits <- parse_trait_option(opts$zscore_traits)
    if (!is.null(parsed_traits)) trait_option_str <- paste(parsed_traits, collapse = ',')
  }

  if (is.null(opts$info) && !is.null(opts$zscore)) {
    auto_inputs <- auto_prepare_inputs(opts$zscore, panel_default, trait_cols = trait_option_str, verbose = opts$verbose)
    opts$info <- auto_inputs$info
    if (is.null(opts$gwas1)) opts$gwas1 <- auto_inputs$gwas$paths[1]
    if (is.null(opts$gwas2) && length(auto_inputs$gwas$paths) >= 2) opts$gwas2 <- auto_inputs$gwas$paths[2]
    if (is.null(opts$ref_plink)) opts$ref_plink <- panel_default
    auto_temp_paths <- c(auto_temp_paths, auto_inputs$info, auto_inputs$gwas$paths)
    auto_trait_names <- auto_inputs$gwas$traits
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

  ld_blocks <- load_ld_blocks(opts$ld_coord)

  info <- load_variant_info(opts$info)
  if (is.null(info) || nrow(info) == 0) stop('变异信息文件读取失败或为空')

  g1 <- NULL
  g2 <- NULL
  if (!is.null(opts$gwas1)) {
    g1 <- load_gwas_for_info(opts$gwas1, info)
    if (!is.null(opts$gwas2)) {
      g2 <- load_gwas_for_info(opts$gwas2, info)
    }
  } else if (!is.null(opts$multi_gwas) && !is.null(opts$zcols)) {
    mg <- load_multi_gwas_for_info(opts$multi_gwas, opts$zcols, info)
    zc <- mg$zcols
    m <- mg$data
    g1 <- m
    g1$Z <- as.numeric(m[[zc[1]]])
    if (length(zc) >= 2) {
      g2 <- m
      g2$Z <- as.numeric(m[[zc[2]]])
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
  plink_map <- NULL
  using_reference_panel <- FALSE
  reference_prefix <- NULL

  if (!is.null(opts$geno_rds)) {
    if (!file.exists(opts$geno_rds)) stop('geno_rds文件不存在')
    G_raw <- load_genotype_rds(opts$geno_rds)
  } else if (!is.null(opts$geno_csv)) {
    if (!file.exists(opts$geno_csv)) stop('geno_csv文件不存在')
    G_raw <- load_genotype_csv(opts$geno_csv, fmt = opts$geno_format)
  } else if (!is.null(opts$geno_plink)) {
    plink_data <- load_genotype_plink(opts$geno_plink, snp_ids = target_snp_ids, chrpos_ids = chrpos_ids)
    G_raw <- plink_data$matrix
    plink_map <- plink_data$map
  } else {
    prefix <- if (!is.null(opts$ref_plink)) opts$ref_plink else file.path('g1000_eur', 'g1000_eur')
    using_reference_panel <- TRUE
    reference_prefix <- prefix
    message('未提供真实基因型，自动使用参考面板: ', prefix)
    plink_data <- load_genotype_plink(prefix, snp_ids = target_snp_ids, chrpos_ids = chrpos_ids)
    G_raw <- plink_data$matrix
    plink_map <- plink_data$map
  }

  geno_source_type <- if (using_reference_panel) 'reference' else 'real'

  ali <- align_genotype_to_gwas(G_raw, gwas_data, prefer_rsid = TRUE)
  G <- ali$G
  variant_positions <- ali$pos

  if (opts$mode == 'select') {
    variant_ids <- colnames(G)
    chrpos_all <- paste0(gwas_data$chr, ':', gwas_data$pos_bp)
    match_idx <- match(variant_ids, chrpos_all)
    if (any(is.na(match_idx))) {
      warning('部分变异无法在GWAS数据中找到匹配，将忽略这些变异')
      keep <- which(!is.na(match_idx))
      if (length(keep) < 5) stop('可匹配的变异数量不足')
      G <- G[, keep, drop = FALSE]
      variant_ids <- colnames(G)
      match_idx <- match(variant_ids, chrpos_all)
      if (!is.null(variant_positions)) {
        variant_positions <- variant_positions[keep]
      }
    }
    z_vec <- gwas_data$zscore_pheno1[match_idx]
    rsid_lookup <- gwas_data$rsid[match_idx]

    selection_table <- NULL
    chunk_summary <- NULL
    chunk_res <- NULL
    ld_reference_used <- NA_character_

    variant_chr_vec <- gwas_data$chr[match_idx]
    variant_pos_vec <- gwas_data$pos_bp[match_idx]
    if (!is.null(variant_positions) && length(variant_positions) == length(variant_ids)) {
      variant_pos_vec <- variant_positions
    }

    if (!is.null(ld_blocks) && nrow(ld_blocks) > 0) {
      chunk_res <- run_block_selection(
        genotype_matrix = G,
        z_vector = z_vec,
        variant_ids = variant_ids,
        variant_chr = variant_chr_vec,
        variant_pos = variant_pos_vec,
        blocks_df = ld_blocks,
        n_samples = opts$n,
        n_knockoffs = opts$knockoffs,
        fdr = opts$fdr,
        use_python_accel = opts$py_accel,
        verbose = opts$verbose
      )
      if (is.null(chunk_res)) {
        stop('变量选择失败 (LD block阶段)')
      }
      selection_table <- chunk_res$table
      chunk_summary <- chunk_res$chunks
      ld_reference_used <- opts$ld_coord
    }

    if (is.null(selection_table)) {
      sel <- run_ghostknockoff_selection(
        genotype_matrix = G,
        z_vector = z_vec,
        n_samples = opts$n,
        n_knockoffs = opts$knockoffs,
        fdr = opts$fdr,
        variant_ids = variant_ids,
        use_python_accel = opts$py_accel,
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
    readr::write_csv(selection_table, out_csv)

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
      chunk_summary$chunk_size_bp <- chunk_summary$end - chunk_summary$start + 1
      if (!is.na(ld_reference_used)) {
        chunk_summary$ld_reference <- ld_reference_used
      }
      chunk_csv <- file.path(out_dir, paste0(prefix, '_chunk_summary.csv'))
      readr::write_csv(chunk_summary, chunk_csv)

      chunk_dir <- file.path(out_dir, paste0(prefix, '_chunks'))
      if (dir.exists(chunk_dir)) unlink(chunk_dir, recursive = TRUE, force = TRUE)
      dir.create(chunk_dir, recursive = TRUE)
      chunk_tables <- if (!is.null(chunk_res)) chunk_res$chunk_tables else list()
      if (!is.null(chunk_tables) && length(chunk_tables) > 0) {
        for (chunk_info in chunk_tables) {
          chunk_tbl <- chunk_info$data
          chunk_file <- file.path(chunk_dir, sprintf('%s_chunk_%02d.csv', prefix, chunk_info$id))
          readr::write_csv(chunk_tbl, chunk_file)
        }
      }
    }

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
      run_mode = 'ghostknockoff_only',
      use_python_accel = opts$py_accel
    )
    chr <- unique(gwas_data$chr)[1]
    res <- process_single_locus(
      chr = chr,
      info_file = opts$info,
      gwas_data = gwas_data,
      genotype_data = G,
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
    if (opts$verbose) cat('相关性分析完成并已保存结果\n')
  } else {
    stop('未知的mode。请使用 select 或 correl')
  }

  list(auto_temp_paths = auto_temp_paths)
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
      LAVAKnock_bivariate(
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
    
    result <- LAVAKnock(
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

# =============================================================================
# 并行处理函数
# =============================================================================

#' 并行处理染色体
process_chromosome_parallel <- function(chr, params, verbose = FALSE) {
  
  if (verbose) cat(paste("处理染色体", chr, "...\n"))
  
  # 查找染色体数据文件
  chr_dir <- file.path(params$data_dir, paste0("chr", chr))
  
  if (!dir.exists(chr_dir)) {
    warning(paste("染色体", chr, "目录不存在:", chr_dir))
    return(NULL)
  }
  
  # 获取该染色体的所有info文件
  info_files <- list.files(chr_dir, pattern = "Info_start.*\\.csv$", full.names = TRUE)
  
  if (length(info_files) == 0) {
    warning(paste("染色体", chr, "没有找到数据文件"))
    return(NULL)
  }
  
  chr_results <- list(
    univariate_results = data.frame(),
    bivariate_results = data.frame(),
    knockoff_results = NULL
  )
  
  # 处理每个locus文件
  for (info_file in info_files) {
    if (verbose) cat(paste("  处理文件:", basename(info_file), "\n"))
    
    # 加载变异信息
    variant_info <- load_variant_info(info_file)
    if (is.null(variant_info) || nrow(variant_info) == 0) next
    
    # 模拟或加载GWAS数据
    if (params$use_example_data) {
      # 使用示例数据进行测试
      warning("使用示例数据模式，应该由测试脚本调用")
      return(NULL)
    } else {
      # 实际使用时应该由调用者提供数据加载函数
      warning("此函数应该由测试脚本调用，并传入真实数据")
      return(NULL)
    }
    
    # 准备Z-score矩阵
    zscore_matrix <- as.matrix(gwas_data[, c("zscore_pheno1", "zscore_pheno2")])
    rownames(zscore_matrix) <- gwas_data$rsid
    
    # 单变量分析
    if (params$run_mode %in% c("full", "lavaknock_only")) {
      locus_start <- min(variant_info$pos_bp)
      locus_end <- max(variant_info$pos_bp)
      
      univ_result <- run_univariate_analysis(
        genotype_data, zscore_matrix, chr, locus_start, locus_end,
        params$sample_size, params$prune_threshold, verbose
      )
      
      if (!is.null(univ_result)) {
        chr_results$univariate_results <- rbind(chr_results$univariate_results, univ_result)
      }
    }
    
    # 生成knockoff并进行双变量分析
    if (params$run_mode %in% c("full", "ghostknockoff_only")) {
      
      # 创建分析窗口
      windows <- create_analysis_windows(variant_info, params$window_size, params$window_step)
      
      if (nrow(windows) == 0) next
      
      # 生成knockoff
      knockoff_scores <- generate_ghostknockoff_scores(
        genotype_matrix = genotype_data,
        zscore_matrix = zscore_matrix,
        ld_threshold = params$ld_threshold,
        n_samples = params$sample_size,
        n_knockoffs = params$knockoffs,
        use_python_accel = isTRUE(params$use_python_accel),
        verbose = verbose
      )
      
      if (is.null(knockoff_scores)) next
      
      # 对每个窗口进行双变量分析
      for (i in 1:nrow(windows)) {
        window <- windows[i, ]
        variant_indices <- window$variant_indices[[1]]
        
        window_genotype <- genotype_data[, variant_indices, drop = FALSE]
        
        # 提取窗口的knockoff分数
        zscore_pheno1_window <- as.data.frame(knockoff_scores$pheno1[variant_indices, ])
        zscore_pheno2_window <- as.data.frame(knockoff_scores$pheno2[variant_indices, ])
        
        # 双变量分析
        bivar_result <- run_bivariate_analysis(
          window_genotype, zscore_pheno1_window, zscore_pheno2_window,
          chr, window$start, window$end, params$sample_size, 
          params$prune_threshold, verbose
        )
        
        if (!is.null(bivar_result)) {
          chr_results$bivariate_results <- rbind(chr_results$bivariate_results, bivar_result)
        }
      }
    }
  }
  
  # 对该染色体进行knockoff过滤
  if (nrow(chr_results$bivariate_results) > 0) {
    chr_results$knockoff_results <- run_knockoff_filter(
      chr_results$bivariate_results, params$fdr, params$knockoffs, 
      verbose = verbose
    )
  }
  
  if (verbose) cat(paste("染色体", chr, "处理完成\n"))
  return(chr_results)
}

# =============================================================================
# 核心处理函数 - 供外部调用
# =============================================================================

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
  variant_info <- load_variant_info(info_file)
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
        use_python_accel = isTRUE(params$use_python_accel),
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
