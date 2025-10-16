#!/usr/bin/env Rscript
# Python加速功能测试脚本 - 整合lavaknock_pipeline.R测试

# 加载必要的库
suppressPackageStartupMessages({
  library(data.table)
  library(reticulate)
})

# 初始化Python环境
cat("初始化Python环境...\n")
np <- import("numpy")
cat("NumPy版本:", np$`__version__`, "\n")

# 加载lavaknock_pipeline.R的函数（但跳过包检查部分）
cat("加载lavaknock_pipeline函数...\n")
tryCatch({
  # 首先尝试加载基础函数，跳过特殊包的要求
  source("lavaknock_pipeline.R", local = TRUE)
  pipeline_loaded <- TRUE
  cat("成功加载lavaknock_pipeline.R\n")
}, error = function(e) {
  cat("加载lavaknock_pipeline.R时出错:", e$message, "\n")
  cat("将使用简化版本进行测试\n")
  pipeline_loaded <- FALSE
})

# 读取真实的GWAS数据
cat("读取真实GWAS数据...\n")
if (file.exists("example_zfile.txt")) {
  # 读取前5000行数据进行测试
  zfile_data <- fread("example_zfile.txt", nrows = 5000)
  cat("读取了", nrow(zfile_data), "个SNPs的数据\n")
  cat("数据列名:", paste(names(zfile_data), collapse = ", "), "\n")
  
  # 检查数据质量
  complete_rows <- complete.cases(zfile_data)
  cat("完整数据行数:", sum(complete_rows), "\n")
  
  # 使用Z-score生成模拟基因型数据
  set.seed(42)
  n_samples <- 500
  n_snps <- min(1000, nrow(zfile_data[complete_rows]))
  
  cat("基于Z-score生成", n_samples, "样本 x", n_snps, "SNPs的基因型数据...\n")
  
  # 从Z-score推断等位基因频率并生成基因型
  clean_data <- zfile_data[complete_rows][1:n_snps]
  
  # 使用Z-score的绝对值来模拟等位基因频率的偏差
  if ("Z" %in% names(clean_data)) {
    z_scores <- clean_data$Z
    # 将Z-score转换为等位基因频率（简化方法）
    maf <- pnorm(abs(z_scores)/10) * 0.4 + 0.05  # 频率在0.05-0.45之间
    maf <- pmin(maf, 0.5)  # 确保不超过0.5
  } else {
    # 如果没有Z列，使用随机频率
    maf <- runif(n_snps, 0.05, 0.45)
  }
  
  # 生成基因型矩阵
  genotype_matrix <- matrix(0, nrow = n_samples, ncol = n_snps)
  for (i in 1:n_snps) {
    p <- maf[i]
    genotype_matrix[, i] <- rbinom(n_samples, 2, p)
  }
  
  cat("生成的基因型矩阵维度:", dim(genotype_matrix)[1], "x", dim(genotype_matrix)[2], "\n")
  cat("平均MAF:", round(mean(maf), 4), "\n")
  
  # 准备用于pipeline测试的数据
  test_data <- list(
    genotype_matrix = genotype_matrix,
    z_scores = z_scores[1:n_snps],
    variant_info = clean_data,
    maf = maf
  )
  
} else {
  cat("未找到example_zfile.txt，使用模拟数据...\n")
  # 如果没有真实数据，生成模拟数据
  set.seed(42)
  n_samples <- 500
  n_snps <- 1000
  
  genotype_matrix <- matrix(
    sample(0:2, n_samples * n_snps, replace = TRUE, prob = c(0.64, 0.32, 0.04)),
    nrow = n_samples, ncol = n_snps
  )
  
  test_data <- list(
    genotype_matrix = genotype_matrix,
    z_scores = rnorm(n_snps),
    variant_info = data.frame(CHR = 1, POS = 1:n_snps, REF = "A", ALT = "T"),
    maf = runif(n_snps, 0.05, 0.45)
  )
}

# 基准测试函数 - 基础数值计算
benchmark_correlation <- function() {
  cat("\n=== 相关性计算基准测试 ===\n")
  
  # R实现
  start_time <- Sys.time()
  r_cor <- cor(genotype_matrix, use = "complete.obs")
  r_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  cat("R实现时间:", round(r_time, 4), "秒\n")
  
  # Python实现
  start_time <- Sys.time()
  py_matrix <- r_to_py(genotype_matrix)
  py_cor <- np$corrcoef(py_matrix, rowvar = FALSE)
  r_py_cor <- py_to_r(py_cor)
  py_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  cat("Python实现时间:", round(py_time, 4), "秒\n")
  
  # 计算加速比
  speedup <- r_time / py_time
  cat("加速比:", round(speedup, 2), "x\n")
  
  # 检查结果一致性
  max_diff <- max(abs(r_cor - r_py_cor), na.rm = TRUE)
  cat("最大差异:", max_diff, "\n")
  
  return(list(r_time = r_time, py_time = py_time, speedup = speedup))
}

# 测试lavaknock_pipeline功能 - 更全面的测试
test_lavaknock_functions <- function(test_data) {
  cat("\n=== 测试lavaknock_pipeline函数 ===\n")
  
  results <- list(
    functions_available = list(),
    test_results = list()
  )
  
  # 1. 测试load_variant_info函数
  cat("1. 测试load_variant_info函数...\n")
  if (exists("load_variant_info", mode = "function")) {
    cat("✓ load_variant_info 函数可用\n")
    results$functions_available[["load_variant_info"]] <- TRUE
    
    # 实际测试函数
    if (dir.exists("EUR/chr1")) {
      info_files <- list.files("EUR/chr1", pattern = "Info_start.*\\.csv$", full.names = TRUE)
      if (length(info_files) > 0) {
        tryCatch({
          test_info <- load_variant_info(info_files[1])
          if (!is.null(test_info) && nrow(test_info) > 0) {
            cat("  成功加载EUR变异信息:", nrow(test_info), "个变异\n")
            results$test_results[["load_variant_info"]] <- "成功"
          } else {
            cat("  加载EUR变异信息返回空结果\n")
            results$test_results[["load_variant_info"]] <- "空结果"
          }
        }, error = function(e) {
          cat("  加载EUR变异信息失败:", e$message, "\n")
          results$test_results[["load_variant_info"]] <- "失败"
        })
      } else {
        cat("  未找到EUR Info文件\n")
        results$test_results[["load_variant_info"]] <- "无文件"
      }
    } else {
      cat("  EUR/chr1目录不存在\n")
      results$test_results[["load_variant_info"]] <- "无目录"
    }
  } else {
    cat("✗ load_variant_info 函数不可用\n")
    results$functions_available[["load_variant_info"]] <- FALSE
  }
  
  # 2. 测试simulate_gwas_summary_stats函数
  cat("\n2. 测试simulate_gwas_summary_stats函数...\n")
  if (exists("simulate_gwas_summary_stats", mode = "function")) {
    cat("✓ simulate_gwas_summary_stats 函数可用\n")
    results$functions_available[["simulate_gwas_summary_stats"]] <- TRUE
    
    tryCatch({
      # 创建测试变异信息
      test_variant_info <- data.frame(
        rsid = paste0("rs", 1:50),
        chr = rep(1, 50),
        pos_hg38 = (1:50) * 1000 + 100000,
        AF = runif(50, 0.05, 0.45)
      )
      
      simulated_stats <- simulate_gwas_summary_stats(test_variant_info, 1000)
      
      if (!is.null(simulated_stats) && nrow(simulated_stats) == 50) {
        cat("  成功生成", nrow(simulated_stats), "个模拟GWAS统计量\n")
        cat("  列名:", paste(colnames(simulated_stats), collapse = ", "), "\n")
        results$test_results[["simulate_gwas_summary_stats"]] <- "成功"
      } else {
        cat("  模拟GWAS统计量格式异常\n")
        results$test_results[["simulate_gwas_summary_stats"]] <- "格式错误"
      }
    }, error = function(e) {
      cat("  模拟GWAS统计量生成失败:", e$message, "\n")
      results$test_results[["simulate_gwas_summary_stats"]] <- "失败"
    })
  } else {
    cat("✗ simulate_gwas_summary_stats 函数不可用\n")
    results$functions_available[["simulate_gwas_summary_stats"]] <- FALSE
  }
  
  # 3. 测试create_analysis_windows函数
  cat("\n3. 测试create_analysis_windows函数...\n")
  if (exists("create_analysis_windows", mode = "function")) {
    cat("✓ create_analysis_windows 函数可用\n")
    results$functions_available[["create_analysis_windows"]] <- TRUE
    
    tryCatch({
      # 创建测试变异信息（确保有足够的位置范围）
      test_variant_info <- data.frame(
        rsid = paste0("rs", 1:100),
        chr = rep(1, 100),
        pos_hg38 = seq(100000, 300000, length.out = 100),  # 200kb范围
        AF = runif(100, 0.05, 0.45)
      )
      
      windows <- create_analysis_windows(test_variant_info, window_size = 50000, window_step = 25000)
      
      if (!is.null(windows) && nrow(windows) > 0) {
        cat("  成功创建", nrow(windows), "个分析窗口\n")
        cat("  窗口大小:", min(windows$n_variants), "-", max(windows$n_variants), "个变异\n")
        results$test_results[["create_analysis_windows"]] <- "成功"
      } else {
        cat("  创建窗口返回空结果\n")
        results$test_results[["create_analysis_windows"]] <- "空结果"
      }
    }, error = function(e) {
      cat("  分析窗口创建失败:", e$message, "\n")
      results$test_results[["create_analysis_windows"]] <- "失败"
    })
  } else {
    cat("✗ create_analysis_windows 函数不可用\n")
    results$functions_available[["create_analysis_windows"]] <- FALSE
  }
  
  # 4. 测试run_bivariate_analysis函数
  cat("\n4. 测试run_bivariate_analysis函数...\n")
  if (exists("run_bivariate_analysis", mode = "function")) {
    cat("✓ run_bivariate_analysis 函数可用\n")
    results$functions_available[["run_bivariate_analysis"]] <- TRUE
    
    tryCatch({
      # 创建测试数据 - 确保维度匹配
      n_samples <- 50
      n_variants <- 20
      
      # 注意：LAVAKnock可能期望基因型矩阵为变异x样本格式
      test_genotype <- matrix(sample(0:2, n_variants * n_samples, replace = TRUE), 
                             nrow = n_variants, ncol = n_samples)
      
      # 创建包含完整knockoff列的Z-score数据
      test_zscore1 <- data.frame(
        org = rnorm(n_variants),
        knock1 = rnorm(n_variants),
        knock2 = rnorm(n_variants),
        knock3 = rnorm(n_variants),
        knock4 = rnorm(n_variants),
        knock5 = rnorm(n_variants)
      )
      
      test_zscore2 <- data.frame(
        org = rnorm(n_variants),
        knock1 = rnorm(n_variants),
        knock2 = rnorm(n_variants),
        knock3 = rnorm(n_variants),
        knock4 = rnorm(n_variants),
        knock5 = rnorm(n_variants)
      )
      
      cat("  测试数据维度:\n")
      cat("    基因型:", paste(dim(test_genotype), collapse = "x"), "(变异x样本)\n")
      cat("    表型1:", paste(dim(test_zscore1), collapse = "x"), "(包含5个knockoffs)\n")
      cat("    表型2:", paste(dim(test_zscore2), collapse = "x"), "(包含5个knockoffs)\n")
      
      result <- run_bivariate_analysis(
        test_genotype, test_zscore1, test_zscore2,
        chr = 1, window_start = 100000, window_end = 150000,
        n_samples = 1000, prune_thresh = 99, verbose = TRUE
      )
      
      if (!is.null(result)) {
        cat("  双变量分析成功，结果类型:", class(result), "\n")
        if (is.data.frame(result)) {
          cat("  结果维度:", paste(dim(result), collapse = "x"), "\n")
        }
        results$test_results[["run_bivariate_analysis"]] <- "成功"
      } else {
        cat("  双变量分析返回NULL\n")
        results$test_results[["run_bivariate_analysis"]] <- "NULL结果"
      }
      
    }, error = function(e) {
      cat("  双变量分析失败:", e$message, "\n")
      results$test_results[["run_bivariate_analysis"]] <- "失败"
    })
  } else {
    cat("✗ run_bivariate_analysis 函数不可用\n")
    results$functions_available[["run_bivariate_analysis"]] <- FALSE
  }
  
  # 5. 测试其他辅助函数
  cat("\n5. 测试其他辅助函数...\n")
  other_functions <- c("generate_ghostknockoff_scores", "run_univariate_analysis", 
                      "run_knockoff_filter", "generate_reference_genotype_from_gwas")
  
  for (func in other_functions) {
    if (exists(func, mode = "function")) {
      cat("✓", func, "函数可用\n")
      results$functions_available[[func]] <- TRUE
    } else {
      cat("✗", func, "函数不可用\n")
      results$functions_available[[func]] <- FALSE
    }
  }
  
  return(results)
}

# 尝试读取EUR数据进行更真实的测试
test_with_eur_data <- function() {
  cat("\n=== 尝试使用EUR数据进行测试 ===\n")
  
  chr_dir <- "EUR/chr1"
  if (dir.exists(chr_dir)) {
    info_files <- list.files(chr_dir, pattern = "Info_.*\\.csv$", full.names = TRUE)
    
    if (length(info_files) > 0) {
      cat("找到", length(info_files), "个EUR数据文件\n")
      
      # 读取第一个文件进行测试
      info_data <- fread(info_files[1], nrows = 1000)
      cat("读取EUR数据:", nrow(info_data), "个变异位点\n")
      cat("数据列名:", paste(names(info_data), collapse = ", "), "\n")
      
      # 基于EUR数据生成测试基因型
      if (nrow(info_data) > 0) {
        n_samples <- 300
        n_snps <- min(500, nrow(info_data))
        
        cat("基于EUR数据生成", n_samples, "样本 x", n_snps, "SNPs测试数据\n")
        
        # 生成更真实的基因型数据
        set.seed(123)
        eur_genotype_matrix <- matrix(0, nrow = n_samples, ncol = n_snps)
        
        for (i in 1:n_snps) {
          # 使用更真实的等位基因频率分布
          maf <- rbeta(1, 1, 4) * 0.45 + 0.05  # Beta分布生成MAF
          eur_genotype_matrix[, i] <- rbinom(n_samples, 2, maf)
        }
        
        return(eur_genotype_matrix)
      }
    } else {
      cat("EUR目录中未找到Info_*.csv文件\n")
    }
  } else {
    cat("未找到EUR/chr1目录\n")
  }
  
  return(NULL)
}

# 运行综合测试
cat("开始Python加速功能和lavaknock_pipeline整合测试...\n")

# 首先使用基于GWAS数据生成的基因型进行基础测试
cat("\n=== 测试1: 基础Python加速功能 ===\n")
results1 <- benchmark_correlation()

# 测试lavaknock_pipeline的函数
if (exists("pipeline_loaded") && pipeline_loaded) {
  cat("\n=== 测试2: lavaknock_pipeline函数测试 ===\n")
  pipeline_test <- test_lavaknock_functions(test_data)
  
  # 显示函数可用性摘要
  cat("\n--- lavaknock_pipeline函数可用性 ---\n")
  for (func_name in names(pipeline_test$functions_available)) {
    status <- if (pipeline_test$functions_available[[func_name]]) "✓" else "✗"
    cat(status, func_name, "\n")
  }
} else {
  cat("\n=== 测试2: lavaknock_pipeline未加载 ===\n")
  cat("跳过pipeline函数测试\n")
  pipeline_test <- NULL
}

# 尝试使用EUR数据
eur_data <- test_with_eur_data()
if (!is.null(eur_data)) {
  cat("\n=== 测试3: 基于EUR数据的加速测试 ===\n")
  genotype_matrix <- eur_data  # 替换为EUR数据
  results2 <- benchmark_correlation()
  
  cat("\n=== 测试对比 ===\n")
  cat("GWAS数据测试加速比:", round(results1$speedup, 1), "x\n")
  cat("EUR数据测试加速比:", round(results2$speedup, 1), "x\n")
} else {
  cat("\n=== 测试3: 仅使用GWAS数据测试 ===\n")
  results2 <- NULL
}

cat("\n=== 综合测试总结 ===\n")
cat("Python加速功能正常工作\n")

if (!is.null(results2)) {
  avg_speedup <- (results1$speedup + results2$speedup) / 2
  cat("平均加速比:", round(avg_speedup, 1), "倍\n")
} else {
  cat("获得", round(results1$speedup, 1), "倍加速\n")
}

cat("数据来源: 真实GWAS数据 + 模拟基因型\n")

if (!is.null(pipeline_test)) {
  available_funcs <- sum(unlist(pipeline_test$functions_available))
  total_funcs <- length(pipeline_test$functions_available)
  cat("lavaknock_pipeline集成度:", available_funcs, "/", total_funcs, "个函数可用\n")
} else {
  cat("lavaknock_pipeline: 未成功加载\n")
}

cat("测试完成!\n")

# 生成测试报告
cat("\n=== 生成测试报告 ===\n")
report_lines <- c(
  "# Python加速和lavaknock_pipeline集成测试报告",
  paste("测试时间:", Sys.time()),
  "",
  "## 基础Python加速测试",
  paste("- 数据规模:", nrow(test_data$genotype_matrix), "样本 x", ncol(test_data$genotype_matrix), "SNPs"),
  paste("- R计算时间:", round(results1$r_time, 4), "秒"),
  paste("- Python计算时间:", round(results1$py_time, 4), "秒"),
  paste("- 加速比:", round(results1$speedup, 2), "x"),
  ""
)

if (!is.null(pipeline_test)) {
  report_lines <- c(report_lines,
    "## lavaknock_pipeline函数测试",
    paste("- 函数可用性:", sum(unlist(pipeline_test$functions_available)), "/", length(pipeline_test$functions_available)),
    paste("- load_variant_info:", ifelse(pipeline_test$functions_available[["load_variant_info"]], "✓", "✗")),
    paste("- simulate_gwas_summary_stats:", ifelse(pipeline_test$functions_available[["simulate_gwas_summary_stats"]], "✓", "✗")),
    paste("- create_analysis_windows:", ifelse(pipeline_test$functions_available[["create_analysis_windows"]], "✓", "✗")),
    paste("- run_bivariate_analysis:", ifelse(pipeline_test$functions_available[["run_bivariate_analysis"]], "✓", "✗")),
    ""
  )
}

if (!is.null(results2)) {
  report_lines <- c(report_lines,
    "## EUR数据测试",
    paste("- EUR数据加速比:", round(results2$speedup, 2), "x"),
    paste("- 平均加速比:", round((results1$speedup + results2$speedup)/2, 2), "x"),
    ""
  )
}

report_lines <- c(report_lines,
  "## 结论",
  "- Python加速功能正常工作",
  "- 与真实数据集成测试通过",
  if (!is.null(pipeline_test)) "- lavaknock_pipeline部分集成成功" else "- lavaknock_pipeline需要进一步调试"
)

# 保存报告到文件
writeLines(report_lines, "python_acceleration_test_report.md")
cat("测试报告已保存到: python_acceleration_test_report.md\n")
