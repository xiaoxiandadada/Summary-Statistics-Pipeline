# Pipeline 工作流概览

## 函数调用流程
```
run_pipeline.R
  └─ execute_pipeline(opts)
       ├─ py_auto_prepare_inputs()/merge_gwas_with_info()
       │     └─ 若 GWAS 无 Z 列，默认使用第一列数值数据
       │     └─ 缺少 Info 时通过 `build_info_from_gwas()` 自动生成
       ├─ align_dual_traits() → 输出齐整的 gwas_data
       ├─ load_genotype_plink()/load_genotype_rds()/load_genotype_csv()
       ├─ align_genotype_to_gwas()
       ├─ align_positions_with_info()
       ├─ partition_ld_blocks()  # Python 切块
       ├─ generate_ghostknockoff_scores()
       │     └─ accelerator.py: QC → corr → LD pruning → GhostKnockoff.prelim/fit
       ├─ run_ghostknockoff_selection()  # select 模式
       │     └─ GhostKnockoff.filter()
       └─ run_bivariate_analysis()/run_knockoff_filter()  # correl 模式
             ├─ LAVAKnock_bivariate()
             └─ LAVAKnock()
       （若启用 `--plot_manhattan` 且 plotting.R 存在，则最终调用 `plot_manhattan()` 输出 PNG）
```

## 关键脚本
- `summary_pipeline.R`：核心流程与工具函数
- `accelerator.py`：Python/NumPy/Numba 加速（相关矩阵、LD cut、PLINK 读取等）
- `run_pipeline.R`：命令行入口与参数解析
- `install_packages.R`：依赖安装（R 包 + Python numpy/numba）
- `scripts/`：Info/GWAS 构建与 demo 工具

## 输入输出（与 README 相同）
- **select** → `results/selection/*.csv`
- **correl** → `results/correlation/*.csv`

## 运行提示
- 一律指定 `--ld_coord`。
- 双表型需两份 GWAS 或 `--multi_gwas --zcols`。
- 默认使用全部 CPU，可通过 `--threads` 限制。
