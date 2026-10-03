# 数据来源与版本

逐项记录见 [datasets.json](datasets.json)，许可与引用见 [THIRD_PARTY_DATA.md](THIRD_PARTY_DATA.md)。版本使用原发布编号或实际文件的SHA-256；对原下载时间未留存的材料，记录本次恢复并用于新运行的时间，不将文件修改时间当作获取时间。

| 数据 | 用途 | 本次核实的版本与记录 |
| --- | --- | --- |
| [GSE125162](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125162) | B2微调 | 38,225细胞×6,733基因的准备矩阵，原计数和准备输入均有哈希；2026-10-03恢复后用于新训练 |
| [UniProt](https://www.uniprot.org/) / [SGD](https://www.yeastgenome.org/) | 蛋白序列、基因标识与预计算表示 | 人810条、酵母6,742条及最终1,175靶点的索引和嵌入均有哈希；ESM-2 layer 33均值池化，调用记录在`raw/tier1_esm2/` |
| [E-MTAB-2391](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-2391) | HIP/HOP化合物响应 | SDRF/IDF、2026-08-16 HTTP记录、18个原始归档哈希；[预处理与重建复核](PREPROCESSING.md) |
| [OrthoDB](https://www.orthodb.org/) | 配对监督与划分 | 保留`current`接口取得的实际映射哈希；原接口未记录发布号，不按脚本注释猜测为v11。753对冻结划分可直接核验 |
| [Kemmeren GSE42528/GSE42526](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42528) | 酵母菌株范围与历史表达基准 | 2026-09-01来源清单，GPL11232；实际表达矩阵来自GSE42526补充文件，见[来源记录](provenance/recovered-20261003/gse42528-source.json) |
| [Replogle K562](https://doi.org/10.25452/figshare.plus.20029387.v1) | 历史跨物种基准与外部验证 | Figshare v1；normalized bulk文件35773217，MD5已核实；[输入处理清单](provenance/recovered-20261003/replogle-source.json) |
| [Norman GSE133349](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE133349) | 历史207查询基本结构基准 | 冻结表达矩阵、元数据及哈希；不参与本次B2训练 |
| [GO / GAF](https://geneontology.org/docs/download-ontology/) | 历史功能标签及注释投影 | 本体、GO-slim、人GAF：2026-07-26；酵母GAF：2026-08-02。获取时间、字节数和哈希见[原来源清单](provenance/recovered-20261003/go-sources.json) |
| [OMA](https://doi.org/10.5281/zenodo.20922901) / [InParanoiDB](https://inparanoidb.sbc.su.se/about/) | 历史同源初始化与外部验证 | OMA All.May2026；InParanoid9。原转换表保留版本字段，许可分别为CC BY 4.0和CC BY-SA 4.0 |
| [Adamson GSE90546](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE90546) / [LINCS GSE92742](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742) | 历史辅助验证/备用方案 | 列明原用途；不纳入本次B2训练或最终筛选，不打包原始大文件 |

scFoundation和ESM-2的上游预训练数据属于第三方模型来源，单独披露在[预训练数据说明](THIRD_PARTY_DATA.md#上游预训练数据)；本项目没有重新实施这部分预训练。

## 处理、划分与数据泄漏控制

- B2按固定细胞顺序与seed 42划分36,314训练细胞、1,911验证细胞，实际清单随每次训练输出。未设置独立测试集，未按实验批次划分，因此验证MSE仅作为同研究内的训练诊断。
- [去重核验](provenance/recovered-20261003/split-audit.json)未发现重复细胞编号、完全重复计数行或跨训练/验证集的相同计数行。可运行 `python scripts/audit_data_split.py --counts raw/training/corpus_counts.npy --output outputs/split-audit.json` 复核。
- 配对监督沿用753对冻结清单：train 511、validation 100、test 142；各划分没有重叠的人源基因或酵母ORF。表示可用性筛选后训练使用510对，最终清单写入任务导出记录。对齐参数只在训练对上拟合。
- GO本体和GAF的发布日期不完全一致，原构建器记录未解析术语；本次保留原投影，不以“更新数据库”悄悄改动历史标签。
- HIP/HOP用于筛选、提名和确认；确认置换重复剂量选择，BH范围是两方向合并的候选家族，未校正此前提名选择。结果是探索性关联，不是独立测试集上的发现率保证。
- 历史选型、开发数据与评估协议继续保留，尤其不把已经用于方法选择的数据称为新的外部测试。规则见[V9协议](../docs/cross_species_match/V9_ESM2_BASIC_STRUCTURE_PROTOCOL.md)和[外部验证协议](../docs/cross_species_match/EV_PROTOCOL.md)。
- 项目没有访问赛事隐藏评测集。上游模型及公开数据可能包含相同基因/蛋白，已披露来源以供赛事方核验；不宣称排除了未知隐藏集的所有重叠。

[合成演示](demo/README.md)由脚本生成，只用于检验计算和输出，不支持生物学结论。
