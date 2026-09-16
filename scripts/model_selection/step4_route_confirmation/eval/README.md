# YeastBridge 评测 harness(T1–T5)

统一评测入口,任务定义以 `docs/plan_v1.0.md` 第四节为准。四路线(A/B/C/D)的模型产物以"特征"身份接入本 harness,保证同数据划分、同协议、同指标。

## 用法(项目根 `/public/home/mengxl/dzy/yeastbridge` 下)

```bash
PY=/public/home/mengxl/dzy/envs/yeastbridge/bin/python
$PY -m eval.run --list                                        # 各任务数据就绪状态
$PY -m eval.run --task t1 --feature routea_scgpt_ft           # T1 细胞状态(需先提取细胞 embedding)
$PY -m eval.run --task t2 --feature routed_scyeast --model logistic
$PY -m eval.run --task t3 --feature routed_scyeast --model ridge
$PY -m eval.run --task t3 --model mean_profile                # 平凡基线(必须一起报告)
$PY -m eval.run --task t4 --feature routed_scyeast            # T4 通路改造排序(ko_sensitivity)
$PY -m eval.run --task t5 --feature routed_scyeast            # 数据效率曲线(1%..100%)
```

- 可选特征:`esm2_mean`(ESM2 650M 均值池化,1280d)| `routea_scgpt_ft`(路线A 微调 scGPT,512d)| `routea_init`(同源初始化对照臂,512d)| `routed_scyeast`(路线D scYeast pos_embedding,200d,w_knowledge ckpt)
- 结果:`results/harness/<task>_<feature>_<model>_s<seed>.json` + `_detail.csv`
- 实验追踪:MLflow sqlite 后端 `results/mlflow.db`;`--no-mlflow` 可跳过

## 结构

- `eval/data.py` — 数据加载层(gene_master / 各路线特征 / T2 标签 / Kemmeren T3 / T4 通路文献 / 细胞 embedding),lru_cache
- `eval/tasks.py` — T1–T5 定义与指标;`REGISTRY` 注册表;新任务/新 baseline 在此加
- `eval/run.py` — CLI + MLflow 记录

## 全特征 × 全任务结果总表(2026-08-09,seed=42)

| 任务 | 指标 | esm2_mean | routea_init | routea_scgpt_ft | routeb_protein(蛋白桥) | routec_scgpt(从零) | routed_scyeast | 基线 |
|---|---|---|---|---|---|---|---|---|
| T1 细胞状态 | ARI/NMI/sil(语料内) | — | — | 0.347/0.552/0.156 | **0.488/0.617/0.292**(留出版) | 0.314/0.499/0.009† | 0.352/0.483/0.233 | **零样本留出子集 ARI/NMI:D 0.514/0.413 > A 0.318/0.292 > B 0.187/0.232** |
| T2 必需基因 | AUROC | **0.824** | 0.688 | 0.749 | 0.817 | 0.620 | 0.738 | prevalence 0.19 |
| T2 | AUPRC | **0.586** | 0.405 | 0.457 | 0.545 | 0.280 | 0.394 | — |
| T3 扰动响应 | Spearman(α=val自动) | **0.172** | 0.151 | 0.153 | **0.172** | 0.151 | 0.149 | mean-profile 0.151 |
| T3 | top100 召回 | 0.217 | 0.201 | 0.204 | 0.216 | 0.194 | **0.218** | mean-profile 0.193 |
| T4 通路排序 | hit@3/hit@5(α=val自动) | 0.321/0.464 | 0.214/0.464 | 0.179/0.464 | 0.071/0.143 | 0.107/0.214 | 0.214/0.464 | 随机 0.140/0.234 |
| T4 | 方向一致率 | 0.571 | 0.571 | 0.571 | 0.571 | 0.571 | 0.571 | 多数类先验 **0.857** |
| T5 数据效率 | 1%→100% Spearman(α=val自动) | 0.039→0.172 单调 | 0.029→0.151 单调 | 0.029→0.153 单调 | 0.030→0.172 单调 | 0.029→0.151 单调 | 0.043→0.149 近单调 | — |

†路线C 的 T1 为语料内指标(训练集即评测集,与路线A/B 全数据版同口径);零样本对照以留出版 vs D 为准。

## 重要基线结论

- **T3 与 T5 的"阴性/非单调"曾是正则化假象**(2026-08-09 修复):固定 α=1.0 严重欠正则;改为 val 划分自动选 α 后,ESM2 0.172 反超 mean-profile(0.151),routeA_ft 0.153 打平略超,且**四特征 T5 曲线全部恢复单调**。教训:任何 T3/T5 结论必须在 val 选参口径下报告,锚点 mean-profile 仍须同表。
- **T4 排序有信号、方向无信号,且对打分器正则强度敏感**:α=1.0 时 routeA_ft hit@3 0.357;α=val 自动(多选 100)后各特征 hit@5 齐平 0.464、hit@3 收敛 0.18-0.32——强收缩下逐基因预测趋同,排名差异塌缩。T4 v1 结论须待 Kemmeren 直接 KO 上限臂对照后才能定;方向判定(0.571)仍不如全猜 overexpress(0.857)。
- **迁移增益的定位(A vs C,核心科学问题的定量答案)**:T2 +0.13 AUROC(0.749 vs 0.620,增益最集中);T3 ≈0(0.153 vs 0.151);T1 语料内 +0.03 ARI 但零样本被 D 反超。**人类预训练主干对基因级扰动响应几乎无增益,增益集中在必需基因先验(T2)与细胞状态几何(T1 silhouette 0.156 vs 0.009)**。
- 各向异性:各 embedding 随机对余弦基线偏高(~0.9),相似度分析用去均值相关或秩次。

## 消融实验(plan §七 高风险项,2026-08-10)

**消融1:路线A 无同源基因退化(scripts/routeA/run_ablations.py;results/routeA/ablation_no_ortholog.json)**

| 特征 | T2 AUROC ortholog/random | T3 Spearman ortholog/random |
|---|---|---|
| routea_scgpt_ft | 0.776 / **0.607** | 0.186 / 0.132 |
| routeb_protein | 0.773 / 0.816 | 0.198 / 0.157 |
| esm2_mean | 0.760 / 0.843 | 0.200 / 0.155 |
| routed_scyeast | 0.739 / 0.732 | 0.181 / 0.130 |

结论:plan 风险坐实——路线A 在无同源基因上 T2 退化 0.17 AUROC(0.776→0.607);**路线B 天然免疫(0.773→0.816 反而略升)**,T3 各特征 ortholog 组均略高(保守基因扰动响应更易学)。B>A 的论断在无同源基因维度上再次独立成立。

**消融2:T3 混杂/批次校正(results/harness/ablation_t3_confound.json)**

- **mean-profile 残差化**(真值与预测同减训练均值谱):routeB/ESM2 残差信号 **0.102-0.103**,routeA 0.033,scYeast 0.059——B/ESM2 的 T3 优势有约六成是共享成分,但残差仍为正且 B 最高,不是纯均值伪影。
- **跨子系列**(GSE42526↔GSE42527 互训互测):全体跌至 0.03-0.08——**批次迁移是 T3 的硬约束**,报告时必须声明结论限于同批次分布内评估。
- **responsive 子集**:全体 ~0.27(vs non_responsive 0.05-0.09)——强响应突变体才是 T3 信号区间,建议报告同时给全量与子集两个数。

**T4 上限臂(kemmeren_direct)与 v1 结论**

实测 KO 谱直接打分:仅 9/28 记录可用(文献改造靶点富集必需基因,而 Kemmeren 是非必需基因敲除库,**必需酶靶点结构性缺失**);hit@3 0.333 < 随机 0.689,方向一致率 0.111。**KO 敏感性代理在数据上限处即失效**——v1 排序信号不可用于改造建议。v2 方向:覆盖必需基因的打分(SGA 遗传互作/敲低预测头/GEARS 式图模型),不再调此代理。按 plan §六条款,T4 留作湿实验首轮校准,不阻塞主线。
- 路线D T1 口径(2026-08-09 定版):门控融合输出**表达位池化**(value>0 位均值)+ TPM 化输入;silhouette 从 -0.025 修复到 +0.233——负值主因是零值位稀释(池化),非 CPM/TPM 口径(TPM 单改反而更差,已记录)。
- **T1 双口径结论(2026-08-09 定版,论文级 finding)**:语料内 B(0.488)>A(0.359)>D(0.352),零样本留出子集 **D(0.514)>>A(0.318)>B(0.187)**——排序完全翻转。迁移路线(A/B)对见过条件的拟合优势随 token 表达力增强而扩大(蛋白>同源),但对未见条件的泛化反向恶化:token 表达力越强,语料记忆越重。酵母原生预训练(D)的泛化优势在此数据规模(33k 细胞)下不可撼动。T1 结论一律以零样本口径为准。
- 路线A T1 公平性细节:留出 YPDDiauxic/Proline/MinimalEtOH(4,766 细胞,27% 条件数)重训(feature `routea_scgpt_ft_h` / `routeb_protein_h`);子集指标复用全数据聚类标签;保留意见:单一留出划分,子集 NMI 不可与 11 类全表直接比。
- 路线 B/C:模型产出后在 `data.py` 的 `_FEATURE_LOADERS` 注册即复用全部任务。
