# YeastBridge 模型说明

## 方法与用途

YeastBridge 使用蛋白表示与酵母化学遗传学响应谱进行候选筛选，面向 GPCR、离子通道等靶点。输入包括人和酵母的蛋白嵌入、B2 投影参数、酵母任务表及化合物响应矩阵；输出为带相关系数、置换 p 值和证据分类的候选清单。

最终 EJ 方法以 `k=20` 的倒数排名融合整合去 PC1 的 ESM-2 相似度、配对监督对齐、B2 投影。EJ-dc 再去除任务谱共享分量。项目实现工作包括匹配与融合、共模去除、双向响应匹配、两阶段统计确认和家族分解。方法定义见 [冻结记录](../docs/cross_species_match/ESM2_JOINT_METHOD_FREEZE.md)。

## 模型组件与来源

| 组件 | 版本与使用方式 | 训练状态 | 来源与许可 |
| --- | --- | --- | --- |
| ESM-2 | `esm2_t33_650M_UR50D`，fair-esm 2.0.0，layer 33 均值池化，1280维；主流程读取预计算嵌入 | 使用第三方预训练模型 | [官方仓库](https://github.com/facebookresearch/esm)，[MIT 许可](https://github.com/facebookresearch/esm/blob/main/LICENSE) |
| B2 投影 | `final_model.pt` 中 `state_dict.pos_emb.proj.weight/bias`，将1280维嵌入映射到768维 route-B 空间 | 项目曾在 scFoundation 骨干上训练，推理阶段读取冻结参数 | [本项目训练记录](TRAINING.md)；上游 [scFoundation](https://github.com/biomap-research/scFoundation) 的代码许可与模型权重许可分别适用 |
| Ridge / identity-anchored / PLS | 冻结超参数为 λ=1、λ=1、k=8；在指定训练配对上拟合 | 运行时确定性拟合 | 实现位于 `scripts/esm2_joint_task_export.py` |

scFoundation 官方代码采用 Apache-2.0；模型权重采用单独的[非商业研究许可](https://github.com/biomap-research/scFoundation/blob/main/MODEL_LICENSE)。B2 是基于该上游模型训练的组件，不应将整份 checkpoint 笼统标为“自有、无第三方限制”。分发前需要保留上游许可并核对适用条件。

scGPT、Geneformer 和 scYeast 在历史选型或基准中使用；相关脚本、结果与勘误保留在 `scripts/model_selection/`、`results_model_selection/` 和 `docs/model_selection/`。这些模型不是所有运行模式的依赖，历史 scGPT 训练脚本不等同于 B2 训练入口。

## 资产、参数与调用记录

嵌入的索引、运行时间和模型信息见 `raw/tier1_esm2/*/run_info.json` 与 `human_810/build_manifest.tsv`。算法配置及种子见 `configs/`。主流程使用本地文件，不调用商业 API。

B2 权重和 HIP/HOP 响应矩阵已从原环境恢复，下载、安装与 SHA-256 见 [资产说明](../data/ASSETS.md)。运行 `python main.py --mode check` 可核对当前完整流程是否具备输入条件。

## 输出解释与限制

输出指标是酵母响应谱的统计关联，不是结合亲和力或已验证药效。筛选、提名和确认使用同一响应数据；确认家族内的多重检验校正没有涵盖前序选择。任务谱相近的家族成员可能难以区分；残差显著也不直接证明作用机制。

历史候选清单保留原数值和证据标签。新运行重新计算两个方向联合确认家族的 BH 值，并把未通过确认阈值的记录标为候选。家族阈值、残差检验与解释见 [家族特异性协议](../docs/cross_species_match/FAMILY_SPECIFICITY_PROTOCOL.md)。

合成演示使用人工生成的小维度表示与响应，不使用 ESM-2/B2 权重。演示输出始终标为“合成演示，非科研候选”。
