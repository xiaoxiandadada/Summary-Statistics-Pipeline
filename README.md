# Summary Statistics Pipeline

用于 GWAS 汇总统计（summary statistics）的端到端分析工具，入口脚本为 `run_pipeline.R`。

- **变量选择** (`--mode select`): GhostKnockoff + LAVAKnock FDR 控制，输出显著变异。
- **局部遗传相关** (`--mode correl`): LAVA-Knock 双表型窗口分析，输出窗口级局部相关与 knockoff 过滤结果。

> 提示：流水线默认通过 `accelerator.py` 调用 NumPy/Numba 执行相关矩阵与 LD 处理，无需额外开关。

## 快速开始

以下命令可直接在仓库根目录执行。若只提供 `zscore.txt` 与参考面板 `g1000_eur`，pipeline 会自动生成所需的 Info/GWAS 文件，再依据用户指定的 LD 坐标（`--ld_coord`，需显式给出 `GRCh37` 或 `GRCh38`）进行切块。

### 1) 安装依赖

```bash
Rscript install_packages.R
```

### 2) 变量选择（参考面板 + 自动切块）

```bash
Rscript run_pipeline.R --mode select   --zscore zscore_GRCh37.txt   --panel g1000_eur   --ld_coord GRCh37   -v -o results
```

- 自动解析 `zscore.txt` 中的 `chr:pos:ref:alt` 或 `CHR/POS` 列，与参考面板按 `CHR+POS` 匹配。
- `--panel` 可省略（默认 `g1000_eur/g1000_eur`）。
- `--ld_coord` 必须显式指定 LD block 坐标，可填 `GRCh37`/`GRCh38` 或自定义 block 文件路径。
- `--threads`（可选）指定 Python 分块线程数；缺省值为可用 CPU 数。

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
Rscript run_pipeline.R --mode select --zscore test/zscore_demo.tsv --geno_plink test/demo_chr10 --panel test/demo_chr10 --ld_coord GRCh37 --threads 4 -o results_geno/zscore_Z_chunks_GRCH37 -v

# 真实基因型 + 自定义 block 文件
Rscript run_pipeline.R --mode select --zscore test/zscore_demo.tsv --geno_plink test/demo_chr10 --panel test/demo_chr10 \
  --ld_coord /path/to/my_blocks.tsv --threads 4 -o results_geno/custom_blocks -v

# 真实基因型 + GRCh38 blocks
Rscript run_pipeline.R --mode select --zscore test/zscore_demo.tsv --geno_plink test/demo_chr10 --panel test/demo_chr10 --ld_coord GRCh38 -o results_geno/zscore_Z_chunks_GRCH38 -v
```

### 4) 局部遗传相关（示例）

```bash
Rscript run_pipeline.R --mode correl   --zscore zscore.txt   --panel g1000_eur/g1000_eur   -v -o results
```

> 扩展到全基因组时，可对 1–22 号染色体循环运行（如需特定窗口，可继续使用 `scripts/build_info_from_plink.R` 手动生成 Info）。

### 性能与内存提示

- 默认 `--threads` 会使用全部可用 CPU，可根据机器资源手动调低（如 `--threads 4`）以避免竞争。
- Pipeline 会按 LD block 切块并逐块写盘，处理超大 zscore 时仍建议按染色体拆分输入，或针对高峰区域自定义更小的 block 文件以降低单次内存峰值。
- 若需进一步控制内存，可将参考面板或真实基因型按染色体拆分，配合 `--panel`/`--geno_plink` 循环批处理，并开启 `-v` 观察 chunk 进度。

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

## 目录速览

```
g1000_eur/                    # 默认参考面板 (.bed/.bim/.fam/.synonyms)
test/demo_chr10.*             # demo 真实基因型（PLINK）
test/zscore_demo.tsv          # demo zscore
results_panel/                # 参考面板运行示例
  ├── zscore_Z_chunks_GRCH37/
  └── zscore_Z_chunks_GRCH38/
results_geno/                 # 真实基因型运行示例
  ├── zscore_Z_chunks_GRCH37/
  └── zscore_Z_chunks_GRCH38/
scripts/                      # 辅助脚本
README.md / WORKFLOW.md       # 使用指南、流程说明
```

## 常见问题

- **Error in `snpStats::read.plink(...): unrecognised snp selected`**  
  表示传入的 rsid / `chr:pos` 在 PLINK `.bim` 中不存在。脚本会忽略缺失 SNP 并继续，但建议检查：
  1. `.bim` 与 zscore 是否在同一参考版本；
  2. 自定义 `--ld_coord` 文件的 chr 格式（是否带 `chr` 前缀）与单位（bp）；
  3. 如需查看详细日志可添加 `-v`，并在生成前过滤缺失 SNP。
