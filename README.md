# Summary Statistics Pipeline

用于 GWAS 汇总统计（summary statistics）的端到端分析工具，入口脚本为 `run_pipeline.R`。

- **变量选择** (`--mode select`): GhostKnockoff + LAVAKnock FDR 控制，输出显著变异。
- **局部遗传相关** (`--mode correl`): LAVA-Knock 双表型窗口分析，输出窗口级局部相关与 knockoff 过滤结果。

## Pipeline 版本：Standard vs Fast

| 版本 | 触发方式 | 主要差异 |
| --- | --- | --- |
| Standard (默认) | 不加 `--py_accel` | 全部计算在 R 中完成，适合通用环境。 |
| Fast | 命令后追加 `--py_accel` | 使用 `accelerator.py` 中的 NumPy/Numba 加速相关系数和 LD 去冗余运算。 |

无法加载 Python 环境时，会自动回退到标准实现。

## 快速开始

以下命令可直接在仓库根目录执行。若只提供 `zscore.txt` 与参考面板 `g1000_eur`，pipeline 会自动生成所需的 Info/GWAS 文件，再依据用户指定的 LD 坐标（`--ld_coord`，需显式给出 `GRCh37` 或 `GRCh38`）进行切块。

### 1) 安装依赖

```bash
Rscript install_packages.R
```

### 2) 变量选择（参考面板 + 自动切块 + Python 加速）

```bash
Rscript run_pipeline.R --mode select   --zscore zscore.txt   --panel g1000_eur   --ld_coord GRCh37   --py_accel -v -o results
```

- 自动解析 `zscore.txt` 中的 `chr:pos:ref:alt` 或 `CHR/POS` 列，与参考面板按 `CHR+POS` 匹配。
- `--panel` 可省略（默认 `g1000_eur/g1000_eur`）。
- `--ld_coord` 必须显式指定 LD block 坐标（`GRCh37` 或 `GRCh38`）。

### 3) 变量选择（真实基因型）

```bash
Rscript run_pipeline.R --mode select   --zscore zscore.txt   --geno_plink genotype   --ld_coord GRCh37   -o results -v
```

### 5) Demo 小数据快速验证

仓库自带 `test/zscore_demo.tsv` 与 `test/demo_chr10`，可在本地快速跑通四种组合：

```bash
# 参考面板 + GRCh37 blocks
Rscript run_pipeline.R --mode select --zscore test/zscore_demo.tsv --panel test/demo_chr10 --ld_coord GRCh37 -o results_panel/zscore_Z_chunks_GRCH37 -v

# 参考面板 + GRCh38 blocks
Rscript run_pipeline.R --mode select --zscore test/zscore_demo.tsv --panel test/demo_chr10 --ld_coord GRCh38 -o results_panel/zscore_Z_chunks_GRCH38 -v

# 真实基因型 + GRCh37 blocks
Rscript run_pipeline.R --mode select --zscore test/zscore_demo.tsv --geno_plink test/demo_chr10 --panel test/demo_chr10 --ld_coord GRCh37 -o results_geno/zscore_Z_chunks_GRCH37 -v

# 真实基因型 + GRCh38 blocks
Rscript run_pipeline.R --mode select --zscore test/zscore_demo.tsv --geno_plink test/demo_chr10 --panel test/demo_chr10 --ld_coord GRCh38 -o results_geno/zscore_Z_chunks_GRCH38 -v
```

### 4) 局部遗传相关（示例）

```bash
Rscript run_pipeline.R --mode correl   --zscore zscore.txt   --panel g1000_eur/g1000_eur   --py_accel -v -o results
```

> 扩展到全基因组时，可对 1–22 号染色体循环运行（如需特定窗口，可继续使用 `scripts/build_info_from_plink.R` 手动生成 Info）。

## 输入格式与自动识别

- **变异信息**：CSV，列包含 `rsid, chr, pos_bp`。可使用 `scripts/build_info_from_plink.R` 从 PLINK 面板生成。
- **GWAS 汇总**：
  - 标准：`CHR, POS, Z`。
  - 仅 `rsid + 数值列`：通过 `normalize_gwas_table()` 自动转换；脚本 `scripts/convert_zscore_to_gwas.R` 可提前批量处理。
  - 多表型：`--multi_gwas <file> --zcols col1,col2`。
- **基因型 / 参考面板**：
  - `--geno_rds`：RDS matrix（列名为 `rsid` 或 `chr:pos`）。
  - `--geno_csv`：CSV，配合 `--geno_format samples_by_snps|snps_by_samples`。
  - `--geno_plink`：PLINK 前缀，使用 `snpStats::read.plink` 读取。
  - 未提供真实基因型时，自动启用 `--ref_plink`（默认 `g1000_eur/g1000_eur`）。

## 输出说明

- Selection：`results/selection/<info>_<pheno>_selection.csv`
  - 列包含 `id, pval.orginal, pval.knockoff*, W, Qvalue, selected, nearest_gene`（其余元信息在 chunk summary 中保留）。
  - 若按 LD blocks 切块（`--ld_coord`），会额外生成 `<prefix>_chunk_summary.csv` 与 `<prefix>_chunks/`，记录每个分块的范围、来源（ld_block/fallback）及选择数。示例输出保存在 `results_panel/zscore_Z_chunks_GRCH37|GRCH38` 与 `results_geno/zscore_Z_chunks_GRCH37|GRCH38` 目录中。
- Correlation：`results/correlation/<info>_<pheno1>__<pheno2>_bivariate.csv` 及 `_significant_windows.csv`。
