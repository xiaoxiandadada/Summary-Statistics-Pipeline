#!/usr/bin/env Rscript

cat('=== Summary Pipeline 功能回归测试 ===\n')
cat(sprintf('开始时间: %s\n\n', Sys.time()))

Sys.setenv(
  OMP_WAIT_POLICY = 'PASSIVE',
  KMP_DUPLICATE_LIB_OK = 'TRUE',
  OMP_NUM_THREADS = '1',
  KMP_INIT_AT_FORK = 'FALSE',
  KMP_AFFINITY = 'none',
  KMP_BLOCKTIME = '0'
)

source('summary_pipeline.R')

# 默认参数模板（模拟 optparse 输出结构）
default_opts <- function() {
  list(
    mode = 'select',
    info = NULL,
    multi_gwas = NULL,
    zcols = NULL,
    gwas1 = NULL,
    gwas2 = NULL,
    outdir = 'results',
    n = 20000L,
    knockoffs = 5L,
    fdr = 0.1,
    win = 100000L,
    step = 50000L,
    pheno1_name = 'trait1',
    pheno2_name = 'trait2',
    geno_rds = NULL,
    geno_csv = NULL,
    geno_plink = NULL,
    panel = NULL,
    ref_plink = NULL,
    zscore = NULL,
    zscore_traits = NULL,
    geno_format = 'samples_by_snps',
    py_accel = FALSE,
    ld_coord = 'GRCh37',
    verbose = TRUE
  )
}

run_case <- function(title, opts) {
  cat(sprintf('--- %s ---\n', title))
  on.exit(cat('\n'), add = TRUE)
  if (dir.exists(opts$outdir)) {
    unlink(opts$outdir, recursive = TRUE, force = TRUE)
  }
  res <- tryCatch(execute_pipeline(opts), error = identity)
  if (inherits(res, 'error')) {
    stop(sprintf('%s 失败: %s', title, res$message))
  }
  if (!is.null(res$auto_temp_paths) && length(res$auto_temp_paths) > 0) {
    suppressWarnings(file.remove(res$auto_temp_paths[file.exists(res$auto_temp_paths)]))
  }
  selection_dir <- file.path(opts$outdir, 'selection')
  if (!dir.exists(selection_dir)) {
    stop(sprintf('%s 未生成 selection 目录', title))
  }
  cat(sprintf('%s 完成，输出目录: %s\n', title, normalizePath(selection_dir, mustWork = FALSE)))
}

# 测试1：仅参考面板（模拟 --panel）
opts_panel <- default_opts()
opts_panel$mode <- 'select'
opts_panel$zscore <- 'test/zscore_demo.tsv'
opts_panel$panel <- 'test/demo_chr10'
opts_panel$outdir <- 'test/results_panel_test'
run_case('参考面板模式', opts_panel)

# 测试2：真实基因型（模拟 --geno_plink）
opts_geno <- default_opts()
opts_geno$mode <- 'select'
opts_geno$zscore <- 'test/zscore_demo.tsv'
opts_geno$geno_plink <- 'test/demo_chr10'
opts_geno$panel <- 'test/demo_chr10'  # 用同一个demo面板完成自动准备
opts_geno$outdir <- 'test/results_geno_test'
run_case('真实基因型模式', opts_geno)

cat('\n所有测试通过。\n')
