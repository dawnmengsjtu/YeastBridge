# GSE42528 酵母 deletion-expression 真值层

## 结论

GSE42528 可以补上“酵母端与高等生物端同为基因扰动后表达变化”的关键缺口。它不是药物响应，也不是 HIP/HOP：每条输入都是一个酿酒酵母基因删除株相对于同一次双通道芯片内 `refpool` 或 wild-type 对照的表达效应。

官方 SuperSeries 分成两个 processed SubSeries：

- GSE42526：784 个删除株，按“显著 mRNA 变化不超过 3 个”归为 non-responsive；
- GSE42527：700 个删除株，按“显著 mRNA 变化至少 4 个”归为 responsive。

二者合计 1,484 个删除株。这个 responsive/non-responsive 分类由表达结果本身产生，只能作为来源分层字段，绝不能作为模型的功能标签，否则会发生标签泄漏。

权威记录：[GSE42528](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42528)、[GSE42526](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42526)、[GSE42527](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42527)。对应研究为 Kemmeren 等的 2014 年 Cell 论文：[PMID 24766815](https://pubmed.ncbi.nlm.nih.gov/24766815/)。

## 为什么使用官方 final matrix

GSE42528 总 `series_matrix` 有 2,633 个 GSM，主要是同一删除株的染料交换阵列。两个 SubSeries 另外提供了 `final_GeneExpressionMatrix`：

- 每个删除株一列最终 M 值；
- 一列与之严格配对的 p-value；
- 一次整合一个或两个 source arrays；
- 原始最终表有 6,182 个 gene-probe rows；其中 12 个基因各有 2 个 probe，按 systematic ID 聚合后是 6,170 个 response genes。

因此训练矩阵使用官方 final matrix，而总 series matrix 只负责把最终列完整追溯到 GSM、突变/对照通道、平台和处理说明。

两侧样本标题并非逐字完全相同：742 条可直接匹配，另外 1,891 条在
Series header 中带有开头的提交者批次前缀（`[hs1990] ` 964 条、
`[hs1991] ` 927 条），而 final matrix 省略了该前缀。预处理只允许剥离
这两个 anchored 前缀；处理后 2,633 对形成严格一一映射，歧义、未匹配、
未使用样本和手工 override 均为 0。所有原文和理由冻结在 alignment TSV，
没有采用模糊字符串匹配。

## 对照与归一化口径

原实验为双通道 spotted oligonucleotide array。每个 GSM 内，一个通道是 deletion strain，另一个是 `refpool` 或 wild type；染料为 Cy5/Cy3。提交者协议记录了：

1. print-tip LOESS；
2. 不做 background subtraction；
3. 使用 wild-type/wild-type hybridization 估计的 gene-specific dye-bias correction；
4. limma 2.12.0 汇总；
5. p-value 采用 Benjamini-Hochberg FDR correction。

派生矩阵不再次减对照、不强行翻转、不把缺失填成零。`yeast_deletion_signature_matrix.tsv.gz` 保存提交者最终 signed M；配对的 `yeast_deletion_bh_fdr_matrix.tsv.gz` 保存官方 processed significance values。

## 输出契约

运行：

```bash
python scripts/fetch_gse42528.py \
  --manifest configs/gse42528.sources.json \
  --data-root /versioned/gse42528_processed_v2

python scripts/prepare_gse42528.py \
  --manifest configs/gse42528.sources.json \
  --raw-dir /versioned/gse42528_processed_v2/raw \
  --output-dir /versioned/gse42528_processed_v2/processed
```

主要输出：

| 文件 | 粒度与用途 |
|---|---|
| `yeast_deletion_signature_matrix.tsv.gz` | 1,484 signatures × 6,170 yeast response genes；signed M，`NA` 保留 |
| `yeast_deletion_bh_fdr_matrix.tsv.gz` | 同形状，配对的 BH-FDR 数值 |
| `yeast_deletion_metadata.tsv` | 删除靶点、GSM、染料交换、对照、target ID 解析和 benchmark blocker |
| `response_gene_metadata.tsv` | 6,170 个 gene-level 特征及其 probe 聚合规则 |
| `probe_to_response_gene.tsv` | 6,182 个 source probes 到 6,170 个 ORF 的完整 lineage；12 个重复基因的 M 取均值、FDR 取保守最大值 |
| `series_sample_metadata.tsv` | 2,633 个 GSM 的样本和通道 lineage |
| `sample_title_alignment.tsv` | 全部 2,633 个 final/Series 原始标题、GSM、匹配方法和理由 |
| `sample_title_prefix_alignment.tsv` | 仅 1,891 个非精确前缀对齐，便于逐条审计 |
| `QUALITY_REPORT.json` | 完整性、唯一性、范围、方向和 claim boundary 审计 |
| `SOURCE_MANIFEST.json` / `DERIVED_MANIFEST.json` | URL、字节数和 SHA-256 冻结 |

metadata 的 `function_label` 固定为 `UNASSIGNED_NEEDS_FROZEN_GO`，`split` 留空，`benchmark_ready=0`。只有接入冻结的 GO/保守模块标签并重建跨物种 split group 后，才能进入 human-only / yeast-only / joint 三臂盲评。

## 2026-09-01 冻结实跑状态

正式输入根目录为：

`/public/home/mengxl/dzy/yeastbridge_re_mvp_stage_20260901/gse42528_processed_v2/processed/`

- 状态：`READY_KO_EXPRESSION_ARM`；远端字段级复核 `PASS`；
- M 与 BH-FDR：均为 1,484 × 6,170，缺失值均为 0；
- M 范围：-6.86872 至 6.11245；BH-FDR 范围：0 至 1；
- 1,484 个 signature ID 和 stable deletion target 均唯一，未解析靶点为 0；
- 2,633 个 GSM 全部追溯，lineage fraction 为 1.0；
- 6,182 probes 到 6,170 genes；12 个重复基因共 24 probes；
- target self-effect 中位 M 为 -2.52738，99.26% 为负，支持当前 signed 方向；
- 本地完整测试：45 passed；远端重新验源、编译、矩阵/metadata/manifest 逐字段复核通过。

关键冻结哈希：

| 资产 | SHA-256 |
|---|---|
| `yeast_deletion_signature_matrix.tsv.gz` | `7dea42f99b5811243b0350e6324a95cd995a9407cab0586f2546ac2cc004493b` |
| `yeast_deletion_bh_fdr_matrix.tsv.gz` | `6869d0d01c073993f62633f3e70729d1c2c0a8d641b395e6a203d0e1736faca3` |
| `yeast_deletion_metadata.tsv` | `81baed8e69d25cf4d3c7aca93327253a27f07f681fc6b05b45a33c3a798b3e5a` |
| `sample_title_prefix_alignment.tsv` | `570e0637561ee60b0f751c2860df547999d9b777dca66d2592c6a0883a2c22be` |
| `DERIVED_MANIFEST.json` | `cfad65ecc28d962650fa442edc3ab22019254db2778024d53301f9c35d35ae54` |
| `RUN_LINEAGE.json` | `8c450bc216797e5597af48026d57718e6519a10a7f5e60449eb3fa513364123b` |

## 靶点名称漂移：为什么不能把两条 LUG1 记录强行合并

2014 数据及 GPL11232 的时代注释把 final title 中的 `LUG1` 对应到
`YCR087C-A`；同一 final 数据还单独包含 `YLR352W` 删除株。当前 SGD 则把
标准名 `LUG1` 赋给 `YLR352W`（[SGD LUG1/YLR352W](https://www.yeastgenome.org/locus/S000004344)），
而 NCBI 仍把 `YCR087C-A` 作为独立 locus tag、Gene ID 850449
（[NCBI Gene 850449](https://www.ncbi.nlm.nih.gov/gene/850449)）。这是名称随时间
迁移，不是两个实验样本重复。为保留真实 deletion genotype，本数据层使用
实验同期 GPL 的 systematic ID：`LUG1` 行保留为 `YCR087C-A`，`YLR352W`
行保留为 `YLR352W`；若仅按当前标准名回填，会错误合并两条物理不同的删除株。

## 数据质量边界

- 可用：真实 deletion-versus-control 同模态表达签名、连续有符号效应、完整 GSM lineage。
- 不可直接宣称：这些 1,484 个删除覆盖所有酵母基因；当前只测了约四分之一基因组，并且是旧式双通道 array。
- 不可混用：Lee 2014 HIP/HOP/compound response 只能在 bridge 通过后做药物优先级和速度评估，不能给 KO truth 打标签。
- 不可泄漏：GSE42526/7 的 response class 不能当 `function_label`，也不能用来调 test threshold。
- 尚需补齐：冻结的 GO-BP/GO-slim/ortholog projection、跨物种连通分量 split，以及 human-only / yeast-only / joint 的正式比较。
