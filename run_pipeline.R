#!/usr/bin/env Rscript

# # 入口脚本：解析命令行，并将计算委托给 summary_pipeline.R
# Sys.setenv(
#   RETICULATE_USE_UV = "0",
#   RETICULATE_DISABLE_AUTOCONFIGURE = "1",
#   RETICULATE_PYTHON = file.path(Sys.getenv("CONDA_PREFIX", ""), "bin", "python")
# )
# suppressMessages(library(reticulate))
# if (nzchar(Sys.getenv("RETICULATE_PYTHON"))) {
#   use_python(Sys.getenv("RETICULATE_PYTHON"), required = TRUE)
# } else {
#   use_condaenv("R", required = TRUE)
# }

# py_config()

suppressPackageStartupMessages({
  if (!requireNamespace('optparse', quietly = TRUE)) {
    stop("需要安装optparse包: install.packages('optparse')")
  }
  library(optparse)
})

source('summary_pipeline.R')

opt_list <- list(
  make_option(c('-m', '--mode'), type = 'character', default = NULL,
              help = '运行模式: select(变量选择) 或 correl(局部遗传相关性)', metavar = 'MODE'),
  make_option(c('-i', '--info'), type = 'character', default = NULL,
              help = '变异信息CSV文件路径 (包含列rsid, chr, pos_bp)', metavar = 'FILE'),
  make_option(c('--multi_gwas'), type = 'character', default = NULL,
              help = '单个GWAS汇总统计文件，包含多表型Z列；配合 --zcols 使用', metavar = 'FILE'),
  make_option(c('--zcols'), type = 'character', default = NULL,
              help = '在 --multi_gwas 中选择的Z列名，逗号分隔。示例: Z_A,Z_B', metavar = 'STR'),
  make_option(c('--gwas1'), type = 'character', default = NULL,
              help = 'GWAS汇总统计文件1 (列: CHR, POS, Z)', metavar = 'FILE'),
  make_option(c('--gwas2'), type = 'character', default = NULL,
              help = 'GWAS汇总统计文件2 (可选; correl模式建议提供)', metavar = 'FILE'),
  make_option(c('-o', '--outdir'), type = 'character', default = 'results',
              help = '输出目录', metavar = 'DIR'),
  make_option(c('-n', '--n'), type = 'integer', default = 20000,
              help = 'GWAS样本量', metavar = 'INT'),
  make_option(c('-k', '--knockoffs'), type = 'integer', default = 5,
              help = 'GhostKnockoff副本数量', metavar = 'INT'),
  make_option(c('--fdr'), type = 'double', default = 0.1,
              help = 'FDR阈值', metavar = 'NUM'),
  make_option(c('--win'), type = 'integer', default = 100000,
              help = '窗口大小(bp) (correl模式)', metavar = 'INT'),
  make_option(c('--step'), type = 'integer', default = 50000,
              help = '窗口步长(bp) (correl模式)', metavar = 'INT'),
  make_option(c('--threads'), type = 'integer', default = NULL,
              help = 'Python 侧分块时的线程数 (默认自动)', metavar = 'INT'),
  make_option(c('--pheno1_name'), type = 'character', default = 'trait1',
              help = '表型1名称（用于日志/输出命名）', metavar = 'STR'),
  make_option(c('--pheno2_name'), type = 'character', default = 'trait2',
              help = '表型2名称（用于日志/输出命名）', metavar = 'STR'),
  make_option(c('--geno_rds'), type = 'character', default = NULL,
              help = '真实基因型矩阵RDS (矩阵: 样本×SNP; 列名为rsid或chr:pos)', metavar = 'FILE'),
  make_option(c('--geno_csv'), type = 'character', default = NULL,
              help = '真实基因型CSV (矩阵: 样本×SNP; 列名为rsid或chr:pos)', metavar = 'FILE'),
  make_option(c('--geno_plink'), type = 'character', default = NULL,
              help = '真实基因型PLINK前缀(无需扩展名)', metavar = 'PREFIX'),
  make_option(c('-P', '--panel'), type = 'character', default = NULL,
              help = '仅提供zscore时使用的参考面板PLINK前缀 (默认 g1000_eur/g1000_eur)', metavar = 'PREFIX'),
  make_option(c('--ref_plink'), type = 'character', default = NULL,
              help = '参考面板PLINK前缀(默认使用内置的 g1000_eur)', metavar = 'PREFIX'),
  make_option(c('--ld_coord'), type = 'character', default = NULL,
              help = '必填: 指定LD block坐标版本 (GRCh37 或 GRCh38)', metavar = 'STR'),
  make_option(c('--gene_catalog'), type = 'character', default = NULL,
              help = '自定义基因注释文件 (可选)', metavar = 'FILE'),
  make_option(c('-Z', '--zscore'), type = 'character', default = NULL,
              help = '单表型或多表型zscore文件 (列: CHR POS Z 或 rsid + 数值列)', metavar = 'FILE'),
  make_option(c('-T', '--zscore_traits'), type = 'character', default = NULL,
              help = '在 zscore 文件中选择的trait列名，逗号分隔 (最多取前两个)', metavar = 'STR'),
  make_option(c('--geno_format'), type = 'character', default = 'samples_by_snps',
              help = '基因型CSV格式: samples_by_snps 或 snps_by_samples', metavar = 'STR'),
  make_option(c('-v', '--verbose'), action = 'store_true', default = FALSE,
              help = '打印详细日志')
)

parser <- OptionParser(option_list = opt_list)
opts <- parse_args(parser)

if (is.null(opts$threads) || is.na(opts$threads) || opts$threads <= 0) {
  detected_threads <- tryCatch(parallel::detectCores(logical = TRUE), error = function(...) NA_integer_)
  if (is.na(detected_threads) || detected_threads < 1L) detected_threads <- 1L
  opts$threads <- detected_threads
}

execution <- execute_pipeline(opts)

if (!is.null(execution$auto_temp_paths) && length(execution$auto_temp_paths) > 0) {
  suppressWarnings(file.remove(execution$auto_temp_paths[file.exists(execution$auto_temp_paths)]))
}
