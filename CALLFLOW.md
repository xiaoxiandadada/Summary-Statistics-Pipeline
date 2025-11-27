# 函数调用流程（详单）

> 列出主要函数的调用关系、功能、输入/输出，覆盖 select 与 correl 两种模式。涉及文件：`run_pipeline.R`、`summary_pipeline.R`、`accelerator.py`、`plotting.R`、`LAVAKnock`。

## 入口与通用准备
- `run_pipeline.R::parse_args()` → 解析 CLI。
- `execute_pipeline(opts)`（summary_pipeline.R）
  - 载入 R 包与 `plotting.R`（若存在）。
  - `source_python("accelerator.py")` 暴露 py_*。
  - GWAS/Info 读取：
    - `py_auto_prepare_inputs()`（仅给 zscore 时自动产出 info/gwas）
    - 或 `py_load_variant_info()` / `py_load_gwas_table()` / `py_load_multi_gwas_table()`.
  - `load_sample_overlap_matrix()`：可选 2×2 重叠矩阵（无则 NULL）。
  - 基因型加载：`load_genotype_rds()` / `py_load_genotype_csv()` / `load_genotype_plink()`（缺省用参考 panel）。
  - 对齐：`align_genotype_to_gwas()`（调用 `py_align_genotype_to_gwas`），输出等长的基因型矩阵、位点顺序。
  - 生成 `gwas_data`（列 rsid/chr/pos_bp/zscore_pheno*，附 pheno 名）。

## 选择模式（--mode select）
1) 分块：`partition_ld_blocks()`（Python 按 LD block 返回索引列表）。
2) 每块执行 `run_ghostknockoff_selection()`：
   - `generate_ghostknockoff_scores()`  
     输入：genotype (samples×SNP)、z_vector、n_samples、n_knockoffs、LD 阈值。  
     过程：`py_qc_genotype_with_mask` → `py_fast_correlation` → `py_ld_pruning` → `GhostKnockoff.prelim/fit`。  
     输出：原始+M 个 knockoff Z，保留索引。
   - 计算原始/knockoff p 值矩阵；调用 `GhostKnockoff::GhostKnockoff.filter`（包内函数）做 FDR，返回表格与阈值。
3) 汇总：拼接各块 selection 表 → `*_selection.csv`；若有 `plot_manhattan()`（plotting.R），绘制 `*_manhattan.png`。

## 相关模式（--mode correl）
1) 窗口生成：`create_analysis_windows()`（100 kb 窗、100 kb 步），输出窗口列表及变异索引。
2) 窗口循环：
   - 截取子基因型/子 Z（两表型）。
   - `generate_ghostknockoff_scores()`（同上，但输入为窗口子矩阵，输出两表型的原始+knockoff Z）。
   - `run_bivariate_analysis()`：
     - 去零方差/秩检查；若秩不足跳过。
     - 单变量过滤：`univ_p_thresh`（默认 0.1），任一表型 univariate p 超阈值则跳过窗口。
     - `lavaknock_bivariate_local()`  
       输入：窗口基因型、两表型 Z(含 knockoff)、n、prune_thresh、sample_overlap、univ_p_thresh。  
       过程：标准化基因型 → 对每组 Z 调用 `univ_bivariate_rg_overlap()`；内部沿用 LAVAKnock 的 `integral.p`/`bivariate.integral`/`bivar.cond.stats`/`conditional.norm` 计算 rg 与双变量 p，sigma/omega 用 nearPD/ridge 纠偏。  
       输出：每窗口的 rg.orginal/pval.orginal 及 knockoff 列。
3) 多重校正：`run_knockoff_filter()` → `LAVAKnock::LAVAKnock` 生成 W/Q/阈值，标记显著窗口。
4) 输出：`*_bivariate.csv`、`*_significant_windows.csv`，若有 `plot_manhattan()` 则自动绘图。

## Python 侧（accelerator.py）主要接口
- IO：`load_gwas_table` / `load_multi_gwas_table` / `load_variant_info` / `load_genotype_csv/plink`。
- 预处理：`qc_genotype_with_mask`、`fast_correlation_matrix`、`ld_pruning`。
- 对齐：`align_genotype_to_gwas`、`align_dual_traits`、`align_positions_with_info`。
- 自动构建：`auto_prepare_inputs`（zscore→info/gwas）。

## 关键参数（默认）
- correl 窗口：size=100,000；step=100,000。
- knockoffs=5；FDR=0.1；LD 阈值=0.75；prune_threshold=99。
- 样本重叠：若提供 `--sample_overlap`，用 2×2 重叠矩阵调整 sigma；否则视为独立。

## 输出位置
- select：`results/selection/*_selection.csv`（可选曼哈顿图）。
- correl：`results/correlation/*_bivariate.csv`、`*_significant_windows.csv`（可选曼哈顿图）。

## 未使用/重复说明
- 主要流程函数均被调用；冗余度低。未发现完全未用的核心函数，但若提供 A1/A2 等位基，可考虑改为直接使用 LAVA alignment 以减少自定义对齐逻辑。***
