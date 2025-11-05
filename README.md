# Summary Statistics Pipeline

用于处理 GWAS 汇总统计的端到端流程。入口脚本 `run_pipeline.R` 支持两种模式：
- `--mode select`：单表型变量选择（GhostKnockoff + LAVAKnock）
- `--mode correl`：双表型局部遗传相关（LAVAKnock_bivariate + FDR）

## 运行概览
1. 解析 Info CSV（`rsid, chr, pos_bp`）与 GWAS（`CHR, POS, Z`）。若提供两份 GWAS，则对齐共享的 `chr:pos`。
2. 加载基因型：优先使用真实数据（RDS / CSV / PLINK），否则回退参考面板 `--ref_plink`（默认 `g1000_eur/g1000_eur`）。
3. 调用 Python (`accelerator.py`) 完成 QC、相关矩阵、LD 去冗余与 knockoff Z 生成。
4. select 模式：使用 `GhostKnockoff.filter` 输出 τ / q / 领先来源，再整合最近基因信息。
5. correl 模式：按 LD block 切片运行 `LAVAKnock_bivariate`，随后 `LAVAKnock` 做窗口级 FDR。
6. 将结果写入 `results/selection/` 或 `results/correlation/` 目录。

## 快速开始
### 依赖安装
```bash
Rscript install_packages.R --python /path/to/python
```

### Demo：变量选择
```bash
RETICULATE_PYTHON=.venv/bin/python \
Rscript run_pipeline.R --mode select \
  --zscore data/demo/zscore_demo.tsv \
  --panel data/demo/demo_chr10 \
  --ld_coord GRCh37 \
  --fdr 0.1 --knockoffs 5 \
  -o results/select_demo -v
```

### Demo：局部遗传相关
```bash
RETICULATE_PYTHON=.venv/bin/python \
Rscript run_pipeline.R --mode correl \
  --info data/demo/demo_info.csv \
  --gwas1 data/demo/gwas_trait1.tsv \
  --gwas2 data/demo/gwas_trait2.tsv \
  --panel data/demo/demo_chr10 \
  --ld_coord GRCh37 \
  --fdr 0.1 --knockoffs 5 \
  -o results/correl_demo -v
```

## 输入要求
- Info CSV：`chr, pos_bp, rsid`（未提供时会依据 `--gwas1` 自动生成）
- GWAS：`CHR, POS, Z`（单表型或双表型；双表型时需两份文件或 `--multi_gwas --zcols`）
- 基因型：`--geno_rds`、`--geno_csv`（需指定 `--geno_format`）、`--geno_plink`
- `--ld_coord`：指定 LD block 坐标版本（`GRCh37` / `GRCh38` 或路径）

## 输出
- `results/selection/<info>_<pheno>_selection.csv`
- `results/selection/<prefix>_chunk_summary.csv`
- `results/correlation/<info>_<pheno1>__<pheno2>_bivariate.csv`
- `results/correlation/<info>_<pheno1>__<pheno2>_significant_windows.csv`

字段示例：`id, pval.orginal, pval.knockoff*, W, Qvalue, selected, leading_source, nearest_gene`。

自动会生成 `*_manhattan.png`（需确保根目录存在 `plotting.R`）。

## 目录速览
```
g1000_eur/             # 默认参考面板
test/                  # demo 数据与输出
summary_pipeline.R     # 核心 R 流程
accelerator.py         # Python 加速模块
run_pipeline.R         # CLI 入口
plotting.R             # 可选的曼哈顿图绘制工具
scripts/               # Info/GWAS 构建脚本
install_packages.R     # 依赖安装
results/               # 最新分析结果
```

## 注意事项
- 请设置 `RETICULATE_PYTHON` 指向包含 `numpy/numba` 的解释器。
- 双表型模式会对齐共享位点，不匹配的 SNP 将被丢弃。
- `--threads` 控制 Python 分块并行度，默认使用全部核心。
