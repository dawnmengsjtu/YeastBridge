# K562 CRISPRi–酵母 deletion 比赛输入数据卡

## 一句话结论

这条人类 KO-like ↔ 酵母 KO 路线已经从“有数据”推进到了 **`READY_COMPETITION_INPUT`**：Replogle K562 genome-wide CRISPRi pseudobulk 已通过全列 finite 硬门，OrthoDB 一对多映射已被封装为 427 个不可拆分的二部图连通分量，train/validation/test 之间的 component 泄漏为 0。

这个状态只表示“可以正式训练和盲评”，不是模型成绩，也不是“人+酵母一定更好”的证明。

## 适配器实际做了什么

1. 先校验 `configs/human_ko_base.sources.json` 中的冻结规则，要求 Norman split/label 未复用、未读模型结果，且 split 必须以 connected component 为单元。
2. 对 11,258 × 8,248 的 normalized bulk `X` 做分块全矩阵扫描。任何一列只要在任意 profile 出现 `NaN` 或正负无穷，整个 Ensembl feature 就被剔除，不填 0、不插值。
3. 保留每个 Ensembl feature 和 symbol；`HSPA14` 和 `TBCE` 的重复 symbol 被显式列出，没有静默合并。
4. 从冻结 OrthoDB 表中只保留 K562 target symbol 与 GSE42528 deletion target 之间的精确边，不按表现选“最佳同源基因”。
5. 把完整 human–yeast 二部图的 connected component 作为唯一 split unit，按冻结 SHA256 seed 排序，再用 Hamilton largest-remainder 分配 component 数。

## 正式 dry run 结果

隔离服务器路径：

`/public/home/mengxl/dzy/yeastbridge_re_mvp_stage_20260901/k562_bridge_input_v1/official_run`

### 人类输入质控

| 项目 | 结果 |
|---|---:|
| profiles × response features | 11,258 × 8,248 |
| 全列 finite 后保留 | 8,175 |
| 含非有限值而整列剔除 | 73 |
| `+Inf` 元素 | 6,881 |
| `NaN` / `-Inf` | 0 / 0 |
| 涉及非有限值的 profiles | 2,341 |
| target profiles / unique symbols | 10,673 / 9,866 |
| 主分析合格 profiles / symbols | 9,580 / 9,577 |

主分析合格的定义是 `P1` 或 `P1P2` transcript profile、`num_cells_filtered >= 10`、并且 target Ensembl ID 存在。

### 同源图与 split

| 项目 | 结果 |
|---|---:|
| 精确 OrthoDB 边 | 881 |
| 可映射 K562 target symbols | 722 |
| 可映射 GSE42528 deletions | 493 |
| connected components | 427 |
| 严格 1:1 components | 271 |
| 一对多/多对多 components | 156 |
| 跨 partition components | **0** |

| split | components | bridge-eligible human profiles | yeast deletions | OrthoDB edges |
|---|---:|---:|---:|---:|
| train | 299 | 508 | 336 | 602 |
| validation | 64 | 91 | 74 | 110 |
| test | 64 | 109 | 83 | 169 |

表中 yeast deletion 数和 edge 数在 split 间不必严格按 70/15/15，因为优先级是不拆 component，而不是把一个同源家族为了配平数字强行切开。

## 关键输出

- `response_feature_finite_manifest.tsv`：8,248 个响应 feature 的全量 finite/剔除记录。
- `duplicate_response_symbols.tsv`：重复 symbol 及其对应 Ensembl IDs。
- `retained_response_ensembl_ids.txt`：8,175 个冻结响应 feature 白名单。
- `human_profile_manifest.tsv`：全部 human profiles、合格原因、component 和 split。
- `bridge_eligible_human_profiles.tsv`：708 条合格且可直接对齐酵母 deletion 的 human profiles。
- `orthodb_bipartite_edges.tsv`、`ortholog_components.tsv`、`split_summary.tsv`：图、分量和三分区。
- `SOURCE_INPUTS.json`、`QC_REPORT.json`、`DERIVED_MANIFEST.json`：来源、质控和派生文件哈希。

`DERIVED_MANIFEST.json` SHA256 为 `899e4b725894308fcfd20889f1079477058d9f7c2dd33fbc1b1f0e910c931347`；第二次独立 replay 的 manifest、QC 和 component 表与首次逐字节一致。完整哈希见 `configs/k562_gse42528_bridge_input.formal_status.json`。

## 边界与下一步

- 直接 target ortholog 覆盖仍只有 7% 左右；这是一条高置信度直达路由，不是全面覆盖人类靶点。
- K562 CRISPRi knockdown 与酵母稳态 deletion 不是同一种 assay；后续必须在物种内标准化，比较方向/程序结构，不直接比原始数值。
- 本次没有 GO label、没有训练、没有读成绩。下一步应在这个已封死的 split 上定义 human-only、yeast-only、joint 三臂的同一盲评终点，才能回答联合是否真的更好。
