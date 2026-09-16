# Demo 第 00 步：模型选用证据（两步竞争性选择）

## 第一步：基座选用 — Norman 扰动检索基准

问题：哪个单细胞基座的预训练基因表示能正确检索未见扰动的程序方向？
任务规模：6 程序 / 80 基因 / 随机基线 0.1667。

| 基座 | accuracy | permutation p |
|---|---|---|
| scfoundation（胜出） | 0.400 | 0.0000 |
| scgpt | 0.350 | 0.0004 |
| geneformer | 0.325 | 0.0018 |

Gate：PASS。声明边界：单数据集的程序结构，不是响应回归声明。
证据文件（bit 级复制自旧项目）：legacy/norman_foundation/retrieval/retrieval_result.json

## 第二步：迁移机制选用 — 五路线对照（同 scF 基座）

A2 正交初始化 / B2 ESM2 蛋白注入 / C2 随机 / D' scYeast / E' 图原生：
B2 胜出（T2 AUC 0.820 / T3 Spearman 0.168）；A2/C2 随机水平，证明只有
ESM2 蛋白嵌入注入能把人源蛋白信息组织进酵母基因表。
证据文件：legacy/feasibility/transfer_routes/RESULTS_v3.md + REGISTRATION_v3.md

## 第三步：本 demo 的恢复链验证（2026-09-02）

- 训练资产按冻结配置重建：6/6 文件与 2026-08-27 冻结副本 sha256 BIT_MATCH
  （A2_init/C2_init/B2_esm2_matrix/A2_init.tsv/corpus_cells/build_stats）
- B2 注入模型权重（final_model.pt）在旧项目中已被删除（目录仅存基因表与
  train_meta.json）；本 demo 按 train_meta.json 冻结超参（epochs=6, batch=32,
  lr=1e-4, mask_p=0.30, zero_mask_p=0.03, seed=42, n_train=36314）重训恢复。
- 重训模型基因表将与冻结的 gene_table_final.npy 做一致性对比后使用。
