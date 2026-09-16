# 高等生物 KO/CRISPRi 基座审计

## 结论

主选是 Replogle 2022 的 **K562 genome-wide CRISPRi Perturb-seq**，不是继续扩大 Norman CRISPRa。原因很直接：它给出 9,866 个被扰动的人类基因和约 200 万个过滤后细胞，并且作者已经公开了只有约 375 MB 的 pseudobulk；比赛版不需要先下载 65.8 GB 单细胞对象，更不需要重跑约 21 TB 的 SRA 原始数据。

官方来源为 [Figshare processed dataset（DOI 10.25452/figshare.plus.20029387.v1）](https://doi.org/10.25452/figshare.plus.20029387.v1)，许可为 CC BY 4.0。论文为 [Replogle et al., Cell 2022](https://doi.org/10.1016/j.cell.2022.05.013)，开放全文见 [PMC9380471](https://pmc.ncbi.nlm.nih.gov/articles/PMC9380471/)，原始测序归档见 [NCBI BioProject PRJNA831566](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA831566)。

这个选择把人类端从 Norman 的 100 个 single perturbations 提升到近一万个 CRISPRi 靶基因。它仍然不是“证明人类+酵母一定更好”；它只是把之前样本数太少、GO 测试组全零的结构性问题补上，使三臂盲评真正可做。

## 官方 processed 文件与成本

| 用途 | 官方文件 | 大小 | 官方/实测 MD5 |
|---|---|---:|---|
| 主模型输入，先过 finite gate | [`K562_gwps_normalized_bulk_01.h5ad`](https://ndownloader.figshare.com/files/35773217?download=1) | 374,587,922 B | `a3dfaa94ea8724217f5ecb1e14a5f0c8` |
| 审计及重处理后备 | [`K562_gwps_raw_bulk_01.h5ad`](https://ndownloader.figshare.com/files/35774443?download=1) | 374,587,922 B | `4570b53c9d62ff6df281e622f0350060` |
| 可选，不用于比赛 MVP | [`K562_gwps_normalized_singlecell_01.h5ad`](https://ndownloader.figshare.com/files/35774440?download=1) | 65,830,941,948 B | `6cd393e369506849ebf959175989d632` |
| 可选，不用于比赛 MVP | [`K562_gwps_raw_singlecell_01.h5ad`](https://ndownloader.figshare.com/files/35775507?download=1) | 65,830,941,948 B | `887e3e6a8c8df6eadf7a3030a53c9546` |

两个 bulk 文件已在隔离服务器目录完整下载并复算 MD5，和 Figshare 一致。最小方案只需约 357 MiB；预留 1 GiB 足够保存输入、特征白名单和映射 manifest。raw 与 normalized 两个 bulk 都保留时约 714 MiB。

raw 与 normalized **只是字节数巧合相同，不是同一个文件**：MD5 不同，`X` 的全部 92,855,984 个元素都不同；但 shape、`obs`、`var` 和两侧索引完全相同。

## h5ad 里实际有什么

字段级读取结果为 11,258 profiles × 8,248 Ensembl response features：

- 10,673 个 target profiles，覆盖 9,866 个唯一 target symbols；
- 9,608 个 P1/P1P2 primary-transcript profiles，覆盖 9,605 个唯一 symbols；
- 按 `P1/P1P2 + num_cells_filtered >= 10 + target Ensembl ID 非空` 冻结后，主分析有 9,580 profiles / 9,577 symbols；
- 585 个 non-targeting controls，其中作者标记的 `core_control=true` 为 514 个；
- 过滤后 target cells 合计 1,914,250，对照 cells 合计 75,328；
- 8,248 个 `var` 行对应 8,246 个唯一 gene symbols；`HSPA14`、`TBCE` 各有两个 Ensembl feature，不能只按 symbol 静默合并。

`obs` 索引为 `gene_transcript`，关键列为：

`UMI_count_unfiltered`、`num_cells_unfiltered`、`num_cells_filtered`、`control_expr`、`fold_expr`、`pct_expr`、`core_control`、`mean_leverage_score`、`std_leverage_score`、`energy_test_p_value`、`anderson_darling_counts`、`mann_whitney_counts`、`z_gemgroup_UMI`、`mitopercent`、`TE_ratio`、`cnv_score_z`。

`var` 索引为 Ensembl `gene_id`，列为：

`gene_name`、`mean`、`std`、`cv`、`in_matrix`、`gini`、`clean_mean`、`clean_std`、`clean_cv`。

## normalized bulk 的必要质量门

normalized bulk 不能不检查就直接喂模型。完整扫描发现：

- raw `X` 非有限值为 0；
- normalized `X` 有 6,881 个 `+Inf`；
- 涉及 2,341 个 profiles 和 73 个 response features；没有 `NaN` 或 `-Inf`。

因此输入协议冻结为：对 normalized bulk 按 Ensembl feature 做全列检查，只保留在全部 11,258 行均为有限值的列；73 个问题列整列剔除，不将无穷值填成 0。实跑后保留 8,175 个 Ensembl features / 8,173 个唯一 symbols。精确 feature manifest 必须在读 GO 标签、训练或看 score 之前落盘并哈希。raw bulk 作为无非有限值的审计/重处理后备。

## 与 GSE42528 的直接桥覆盖

以下只计算冻结 OrthoDB 表里的直接 human–yeast edges，不把 GO 功能相似冒充成同源。GSE42528 当前有 1,484 个 deletion targets × 6,170 个 response genes。

| 对齐面 | 人类全集 | 可直接映射的人类基因 | 映射到的酵母基因 | 直接边 |
|---|---:|---:|---:|---:|
| 9,866 个全部 K562 targets → GSE deletion targets | 9,866 | 722（7.32%） | 493 | 881 |
| 9,577 个主分析 eligible targets → GSE deletion targets | 9,577 | 708（7.39%） | 493 | 857 |
| finite gate 后 K562 response symbols → GSE response genes | 8,173 | 2,269（27.76%） | 2,106 | 3,115 |

最后一行已经按确定性的 finite gate 重算；8,175 个 Ensembl feature 中有两个 symbol 重复，所以唯一 symbol 为 8,173。

虽然 target 直接覆盖只有约 7.3%，绝对数量已经是 700 级，不再是 Norman 的个位数/十位数组。GO-BP/GO-slim 可以扩大“功能坐标”覆盖，但必须以 K562 target/response scope 重新构建；不能复制 Norman 的 GO 标签、开发/测试划分或任何调参决定。

## 一对多 ortholog 怎么处理

722 个人类 target 与 493 个 GSE deletion targets 形成 881 条边，不是一一对应。把每条边当独立样本会重复计权并产生泄漏。

预注册规则如下：

1. 人类和酵母每条 perturbation signature 仍各自保留一行，不挑“最佳 ortholog”，也不把 881 条边扩成 881 个 IID 配对样本。
2. 用冻结 OrthoDB 构建 human-target ↔ yeast-deletion 的二部图，以完整 connected component 为 split unit。
3. 一个 component 的全部人类基因、酵母基因和边只能进入 train、validation、test 中的一个分区。
4. 先冻结 component manifest，再计算 K562-scope GO 标签、训练或看任何分数。

全量图有 427 个 components：271 个为严格 1:1，156 个为一对多/多对多；最大 component 含 22 个人类基因、最多 5 个酵母基因、最多 24 条边。这 427 个独立分组足以做比旧 Norman 11 个 GO 测试组更稳定的 grouped split。

## assay 可比，但不能硬等价

两侧共同点是“基因功能降低后测全转录组”：K562 是 CRISPRi knockdown，论文报告 K562 中位 knockdown 85.5%；GSE42528 是酵母 deletion。二者都能产生有方向的表达响应，因此适合比较保守功能程序。

不能直接比较原始数值大小：K562 是 10x 单细胞后 pseudobulk，并在感染后第 8 天测量；酵母是稳态 deletion clone 的双通道 microarray M 值。正式桥必须在物种内标准化，只比较方向、排序或保守程序结构；不能说 `human z=1` 等于 `yeast M=1`，也不能用药物响应给 KO 真值打标签。

## 推荐、备选与 NO-GO

- **推荐**：K562 genome-wide normalized bulk，先执行整列 finite gate；raw bulk同时保留作审计。状态是 `RECOMMEND_WITH_REQUIRED_FINITE_FEATURE_GATE`。
- **备选 1**：K562 essential-scale raw bulk，79,766,954 B，论文规模 2,057 genes。适合快速 smoke test 和同细胞系复现，不能替代主选的广覆盖。
- **备选 2**：RPE1 essential-scale raw bulk，95,350,546 B，约 2,393 genes。它是非癌、近二倍体、p53 阳性的稳健性验证，不能作为唯一人类基座。
- **NO-GO**：继续用 Norman 100 个 single CRISPRa 承担主三臂结论；直接下载 21 TB SRA 重做比赛版；把 normalized bulk 的 `Inf` 填零；按表现挑 ortholog；复用 Norman test/GO 标签。

所有 URL、字节数、MD5、字段、覆盖统计和协议边界冻结在 `configs/human_ko_base.sources.json`。该配置只证明数据已经就位并可进入冻结 benchmark，不构成 joint superiority 结果。
