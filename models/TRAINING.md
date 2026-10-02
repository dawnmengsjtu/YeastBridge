# B2 训练溯源

最终 EJ 方法依赖项目早前训练得到的 B2 投影。推理不重新训练 B2，不等于项目未开展训练。

已保留的证据：

- [`train_meta.json`](../raw/tier1_models/route_b/train_meta.json)：route B2，6 epochs，batch 32，梯度累积1，学习率0.0001，mask_p 0.3，zero_mask_p 0.03，seed 42，训练样本36,314，验证样本1,911。
- [训练协议](../docs/model_selection/REGISTRATION_v3.md)：scFoundation 路线的输入处理、冻结层、训练目标、优化器和评估任务。
- [模型选型结果](../docs/model_selection/RESULTS_v3.md)及[勘误](../docs/model_selection/ERRATUM.md)。
- `raw/tier1_models/route_b/gene_table_final.npy`、`A2_init.tsv`：已发布的基因表示表和行顺序。

尚未在当前仓库找到可完整复核的材料：

1. 产生该 B2 checkpoint 的确切训练源代码、环境依赖与入口命令。
2. 对应的数据版本、逐样本划分清单和预处理运行记录。
3. 逐 epoch 训练/验证日志、硬件与耗时记录。
4. 冻结的 `final_model.pt` 文件或可核实的获取地址。

`scripts/model_selection/step4_route_confirmation/finetune_scgpt.py` 是历史 scGPT 路线，依赖原环境的 scGPT 源码和未提交的 `protein_inject` 模块；不能用它声称已复现 scFoundation/B2 训练。当前不提供会伪装成功的 `train.py` 占位入口。

提交附件5的训练材料前，需从原实验环境恢复以上文件，并通过权重哈希确认对应关系。合成演示中的小型 ridge 拟合仅检验流程，不替代这些材料。
