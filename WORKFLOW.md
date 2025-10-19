# LAVA-Knock Pipeline 工作流梳理

本项目支持两种分析模式：
- select：变量选择（GhostKnockoff + LAVAKnock FDR 控制），针对单表型
- correl：局部遗传相关性（窗口双表型 + Knockoff 过滤），针对两表型

## 输入
- 变异信息 CSV：列包含 `rsid, chr, pos_bp`
- GWAS 汇总统计：列包含 `CHR, POS, Z`（一或两个表型文件）；亦支持 `rsid + 数值列`（如 `chr:pos:ref:alt + zscore`），入口自动解析生成 `CHR/POS`
- 真实基因型（可选，推荐）：RDS、CSV 或 PLINK（矩阵：样本×SNP；列名为 rsid 或 `chr:pos`）

## 核心步骤
1) 数据对齐与准备
   - 读取 Info 与 GWAS，按 `chr + pos` 对齐；若 GWAS 仅含 `rsid`，自动拆解 `chr:pos`
   - 加载真实基因型（RDS/CSV/PLINK）；若未提供，默认回退到 `--ref_plink`（内置 `g1000_eur`）参考面板
   - 将基因型列名规范化为 `chr:pos` 并与 GWAS 顺序对齐（同时记录 rsid 以便输出）

2) GhostKnockoff 生成（两模式均使用）
   - 质控：MAF/MAC/方差过滤；可选 LD 修剪（hclust 或 Python 加速）
   - 计算收缩 LD（`corpcor::cor.shrink`）并生成 M 个 knockoff Z 分数

3) 模式分流
   - select：将原始与 knockoff Z 转为双侧 p 值，调用 `LAVAKnock()` 做 FDR 控制；无论使用真实基因型还是参考面板，均依据 `--ld_coord` 指定的 LD block 坐标逐块执行 GhostKnockoff，输出 block summary（必要时自动 fallback 到染色体范围）
   - correl：按窗口切片，调用 `LAVAKnock_bivariate()` 计算窗口双表型局部相关性，再用 `LAVAKnock()` 做窗口级 FDR 过滤

4) 输出
   - select：`results/selection/<info>_selection.csv`
     - Demo 结果示例：`results_panel/zscore_Z_chunks_GRCH37|GRCH38`、`results_geno/zscore_Z_chunks_GRCH37|GRCH38`
   - correl：`results/correlation/<info>_bivariate.csv` 与 `<info>_significant_windows.csv`（如有显著窗口）

## Python 加速接入
- Python 加速：`accelerator.py`（可选）
  - R 侧通过 `reticulate::source_python('accelerator.py')` 载入
  - 主要函数：`fast_correlation_matrix`、`ld_pruning` 等（R 侧封装为 `py_fast_correlation`/`py_ld_pruning` 回退到 R 实现）

在 `summary_pipeline.R` 的 `generate_ghostknockoff_scores()` 中，通过 `use_python_accel` 参数控制是否使用 Python 加速的相关系数与 LD 修剪。

## 脚本入口（单入口 + 可选加速）
- `run_pipeline.R`
  - `--mode select|correl` 选择模式；`--info` 指定窗口信息
  - `--gwas1/--gwas2` 或 `--multi_gwas + --zcols` 指定 GWAS 输入，支持自动解析 `rsid` 型格式
  - `--geno_rds`、`--geno_csv`、`--geno_plink` 输入真实基因型；`--geno_format` 描述 CSV 方向
  - `--ref_plink` 指定参考面板前缀（默认 `g1000_eur/g1000_eur`），仅在缺少真实基因型时使用
  - `--ld_coord` 控制 selection 模式下使用的 LD block 坐标（须指定 `GRCh37` 或 `GRCh38`）
  - 其他参数：`--n`、`--win`、`--step`、`-k`、`--fdr`、`-v`、`--py_accel` 等

## 清理脚本
- `scripts/clean_workspace.R`：清理临时与过程文件（默认预览，使用 `-y` 执行）

## 目录速览

```
g1000_eur/                    # 默认参考面板 (.bed/.bim/.fam/.synonyms)
test/demo_chr10.*             # 示例真实基因型（PLINK）
test/zscore_demo.tsv          # 示例 GWAS 子集（与 demo_chr10 对齐）
zscore_GRCh37.txt             # 全量 GRCh37 坐标汇总统计（用户提供）
zscore_GRCh38.txt             # 全量 GRCh38 坐标汇总统计（用户提供）
results_panel/
  ├── zscore_Z_chunks_GRCH37/ # 参考面板 + GRCh37 blocks 的示例输出
  └── zscore_Z_chunks_GRCH38/ # 参考面板 + GRCh38 blocks 的示例输出
results_geno/
  ├── zscore_Z_chunks_GRCH37/ # 真实基因型 + GRCh37 blocks 示例输出
  └── zscore_Z_chunks_GRCH38/ # 真实基因型 + GRCh38 blocks 示例输出
scripts/                      # Info/GWAS 转换、demo 构建等辅助脚本
README.md                     # 使用说明（含 demo 命令）
WORKFLOW.md                   # 工作流原理与步骤梳理
```
