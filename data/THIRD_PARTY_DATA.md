# 第三方数据与模型使用说明

本项目的用途是非商业科研与赛事复核。下列条款按2026-10-03可核实的官方记录列出，数据处理后的矩阵和映射保留原始来源；不把第三方材料重新授权为本项目原创。下载核查的URL、时间与响应哈希见 [online-source-checks.json](provenance/recovered-20261003/online-source-checks.json)。

| 来源 | 适用条款与归属 | 本项目处理与分发方式 |
| --- | --- | --- |
| Gene Ontology Consortium | [CC BY 4.0与引用要求](https://geneontology.org/docs/go-citation-policy/) | 注释投影是派生数据；保留GO来源、发布日期和许可。本体为2026-07-26，人GAF对应同版本，酵母GAF对应2026-08-02，未把不同日期写成同一版 |
| UniProt Consortium | [数据库版权部分采用CC BY 4.0](https://www.uniprot.org/api-documentation/support-data) | 蛋白索引与嵌入标明UniProt来源、模型和文件哈希；均值池化为本项目计算 |
| Saccharomyces Genome Database (SGD) | [数据库作者的使用说明](https://pmc.ncbi.nlm.nih.gov/articles/PMC11912841/)：CC BY 4.0 | 基因标识与注释转换到项目gene master；保留来源，内容版本由恢复记录中的文件哈希固定 |
| OrthoDB，EM Zdobnov lab | 官方页面结构化元数据给出[CC BY 4.0](https://www.orthodb.org/)；[证据记录](provenance/recovered-20261003/orthodb-license-evidence.json) | 保留下载映射的内容哈希和处理规则；未把BUSCO数据集的CC BY-ND条款套用到OrthoDB映射 |
| OMA Browser，All.May2026 | [该次发布记录](https://doi.org/10.5281/zenodo.20922901)明确为CC BY 4.0；[保存的元数据](provenance/recovered-20261003/oma-zenodo-metadata.json) | 映射为OMA预测经本项目基因ID转换的派生表；注明修改。旧FAQ曾给出CC BY-SA 2.5，本表依据实际All.May2026发布记录 |
| InParanoiDB 9 | [官方关于页](https://inparanoidb.sbc.su.se/about/)：CC BY-SA 4.0 | `raw/externalvalidation/mappings/inparanoid_yeast_human.tsv`为转换后的映射，保留署名和相同许可；不将其并入无约束的项目许可 |
| Replogle et al. 2022 | [Figshare v1](https://doi.org/10.25452/figshare.plus.20029387.v1)：CC BY 4.0；[元数据](provenance/recovered-20261003/replogle-figshare-metadata.json) | 使用公开细胞系扰动表达数据，保留文件ID、MD5、预处理和划分清单；不分发原始人体身份信息 |
| GEO：GSE125162、GSE133349、GSE42528/GSE42526、GSE90546、GSE92742 | [NCBI分子数据库使用政策](https://www.ncbi.nlm.nih.gov/home/about/policies/)：NCBI不另加使用或分发限制，原提交者的第三方权利仍可能适用 | 按accession和原研究署名；不自行标为CC-BY或公有领域。提交训练所需的处理矩阵与处理方法；备用数据只提供来源和用途说明 |
| Lee et al. 2014，E-MTAB-2391 | [EMBL-EBI使用条款](https://www.ebi.ac.uk/about/terms-of-use/)：对原作者提交数据不增加使用或再分发限制，原作者权利继续适用 | 提供公开条目的来源、SDRF/IDF、原始归档哈希和处理方法；未发现条目单独声明的CC许可，因此不指定一个不存在的SPDX许可 |
| scFoundation，BioMap/Tsinghua团队 | 源码[Apache-2.0](https://github.com/biomap-research/scFoundation/blob/main/LICENSE)；权重单独适用[MODEL_LICENSE](https://github.com/biomap-research/scFoundation/blob/main/MODEL_LICENSE) | B2由此微调，保留作者、许可和修改说明，限定非商业研究用途；初始权重从官方入口获取 |
| ESM-2，Meta FAIR团队 | [MIT](https://github.com/facebookresearch/esm/blob/main/LICENSE) | 使用`esm2_t33_650M_UR50D`，未微调ESM-2；分发的是带模型版本和输入索引的预计算表示 |

GEO和EMBL-EBI的访问政策不是原作者向所有用途授予的独立许可。这里准确披露其公开科研数据政策及出处，不作额外商业使用授权。赛事代码使用本地模型，无商业API、提示词或收费平台调用。

## 上游预训练数据

- ESM-2采用上游模型表所列 **UR50/D 2021_04**。本项目未下载该完整预训练语料，实际取得并调用的是上述模型及本项目蛋白序列。来源：[官方模型表](https://github.com/facebookresearch/esm#available-models)。
- scFoundation上游模型预训练于超过5,000万人类单细胞表达谱；数据来自GEO、Single Cell Portal、HCA和EMBL-EBI。逐数据集清单由原论文Supplementary Data 1、2提供，见[原论文数据可用性](https://doi.org/10.1038/s41592-024-02305-7)。本项目未重新收集或执行这部分预训练；实际初始权重由SHA-256固定，取得本次提交副本的时间是2026-10-03。
- 上游预训练可能包含公开基因/蛋白及相关实验数据；不能把“未使用赛事隐藏集”扩大为“上游模型与任何公开评估完全无重叠”。本项目披露其范围，对本地训练、配对监督和外部评估分别保存划分；隐藏评测集重叠由赛事方核验。
