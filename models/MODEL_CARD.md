# 模型说明（Model Card）

本项目最终链不自训练深度模型；使用的模型与权重如下（均为冻结资产，哈希在册）：

| 模型 | 用途 | 版本/来源 | 许可 | 获取 |
|---|---|---|---|---|
| ESM-2 650M (esm2_t33_650M_UR50D) | 人/酵母蛋白序列嵌入（主语言） | fair-esm, Meta AI | 模型权重随 fair-esm 发布（非商用研究用） | 官方下载后放 raw/tier1_models/esm2_650M/，sha256=ea9d0522b335a877…（README 哈希清单核对） |
| B2 注入投影 final_model.pt | ESM-2 -> route-B 空间线性注入（b2_verbatim 成员） | 本项目早前在 scFoundation 骨干上训练（训练协议见 docs/model_selection/REGISTRATION_v3.md 与 scripts/model_selection/step4_route_confirmation/） | 自有 | 放 raw/tier1_models/route_b/（.gitignore 排除，按哈希回填） |
| scFoundation / scGPT / Geneformer | 仅用于基座选型对比（未进最终链） | 各官方发布 | 各自许可 | 选型数字见 docs/model_selection/RESULTS_v3.md |

输入输出：ESM-2 输入蛋白序列（FASTA），输出 1280 维均值池化嵌入；B2 输入 1280 维嵌入，
输出 768 维 route-B 空间向量。已知局限：嵌入对家族内细分分辨率受序列相似度约束；
对无序列信息的靶点不可用。第三方调用时间与参数：全部为本地推理，无 API 调用。
