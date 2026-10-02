# 数据来源与处理记录

| 数据 | 用途 | 获取入口 | 已有版本/时间记录 |
| --- | --- | --- | --- |
| Lee et al. 2014 HIP/HOP | 化合物响应读出 | [E-MTAB-2391](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-2391) | 已恢复SDRF/IDF、HTTP响应头、18包哈希与预处理记录，见 [数据准备](PREPROCESSING.md) |
| Replogle K562 GWPS | 人侧扰动与 EV1 查询 | [研究论文及数据入口](https://doi.org/10.1016/j.cell.2022.05.013) | 原记录为 Replogle 2022；需补实际使用文件的版本与获取时间 |
| Norman 2019 | 207查询基本结构基准 | [GSE133349](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE133349) | 冻结矩阵、元数据与哈希在库 |
| Kemmeren 酵母 KO | 酵母扰动基准 | [GSE42528](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42528) | 冻结矩阵、元数据与哈希在库 |
| GO / GAF | 功能标签和注释投影 | [GO 下载与发布说明](https://geneontology.org/docs/download-ontology/) | 本仓库保留投影表；原 ontology/GAF 发布日期待补 |
| OrthoDB / OMA / InParanoid | 训练同源对与 EV0 验证 | [OrthoDB](https://www.orthodb.org/)、[OMA](https://omabrowser.org/)、[InParanoid](https://inparanoid.sbc.su.se/) | 映射表位于 `raw/externalvalidation/mappings/`；准确发布版本待补 |
| GSE125162 | B2微调的酵母单细胞表达数据 | [GSE125162](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125162) | 已恢复38,225细胞计数、准备代码和重建划分；原下载时间及历史环境快照未留存 |
| Adamson | 备用外部验证与预处理记录 | [GSE90546](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE90546) | 原项目记录2026-08；未纳入最终候选筛选 |
| LINCS L1000 | 历史备用方案 | [GSE92742](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742) | 原项目记录2026-09；对应验证模块未交付 |

“公开可下载”不等同于指定的再分发许可。原仓库将部分来源笼统写为“论文开放数据”“各自学术许可”或推测 CC-BY，不能据此确认授权。参赛附件还需保留各实际下载文件对应的许可/使用条款快照；本说明不替这些缺失记录指定许可。

## 已保存的处理与划分

- 表达基准的矩阵和元数据在 `raw/externalvalidation/frozen_inputs/gse_go_three_arm_v3_open_set/inputs_v3_open_set/`，数据构建说明见 [K562 输入](../docs/cross_species_match/K562_GSE42528_BRIDGE_INPUT.md)和[酵母表达数据](../docs/cross_species_match/GSE42528_DELETION_EXPRESSION.md)。
- `results_a6000/v51_per_query_ranks_trackA.tsv` 保留753对同源配对的 train/test 和 fold 记录。后续表示可用性筛选会取子集，计数写入对应运行记录。
- 基本结构评估采用分组划分；分组和去重规则见 [V9 协议](../docs/cross_species_match/V9_ESM2_BASIC_STRUCTURE_PROTOCOL.md)。外部查询构建与排除规则见 [EV 协议](../docs/cross_species_match/EV_PROTOCOL.md)。这些记录说明已有措施，不代表与赛事未知隐藏集的重叠已经核验。
- 蛋白嵌入的索引和提取记录在 `raw/tier1_esm2/`。计算流程校验冻结输入的 SHA-256。
- HIP/HOP响应、剂量单位、原处理记录及脚本已恢复。全量处理表达矩阵重建的8个NPZ数组均与冻结输入一致，见 [数据准备](PREPROCESSING.md)。
- B2逐样本划分和完全重复行检查见 [训练说明](../models/TRAINING.md)。随机细胞划分未按实验批次分组，不能据此排除所有数据泄漏。

## 演示数据

[合成演示](demo/README.md)由本仓库脚本生成，无需下载生物数据。它只检验计算衔接和输出格式，不能用于支持科研结论。
