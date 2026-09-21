# 数据来源与许可

| 数据 | 用途 | 来源/版本 | 获取时间 | 许可 |
|---|---|---|---|---|
| Lee et al. 2014 HIP/HOP | 酵母化学-遗传学响应（执行层读出） | E-MTAB-2391（ArrayExpress）；官方 per-screen FD 恢复件 | 2026-09 | CC-BY（ArrayExpress 惯例，以页面为准） |
| Replogle K562 GWPS | 人侧 CRISPRi 签名 / EV1 外样本查询 | Replogle 2022 附属（h5ad） | 项目冻结时 | 论文开放数据 |
| Norman 2019 | 基本结构基准（207 查询） | GSE133349 | 项目冻结时 | GEO 开放 |
| Kemmeren GSE42528 | 酵母 KO 签名 | GSE42528 | 项目冻结时 | GEO 开放 |
| GO / GAF | 功能词表与注释 | GO 官方（tier0_gaf_obo） | 项目冻结时 | CC-BY 4.0 |
| OrthoDB / OMA / InParanoid | 同源对（训练对与 EV0 外样本） | 各官方发布（mappings/） | 2026-09 | 各自学术使用许可 |
| Adamson GSE90546 | EV1b 备用外样本（未进主链） | GEO | 2026-08 | GEO 开放 |
| LINCS L1000 GSE92742 | EV2 备用（模块未交付） | GEO | 2026-09 | GEO 开放 |

数据划分与泄漏防控：官方 train/test 切分（753 对）、split_group 分组防泄漏（复制/剂量/
跨物种同源共组）、负例冻结种子；全部种子见 configs*/（seed 字段）。大体积原始件不入库，
按 MANIFEST.sha256 与 raw/externalvalidation/ 内记录的哈希回填；受限许可数据未使用。
