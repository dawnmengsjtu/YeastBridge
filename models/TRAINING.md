# B2 训练与溯源

默认入口采用`submission-v2`，从官方scFoundation初始权重重新训练；旧B2只作为历史版本保留。新版本修正逐行遮蔽和验证模式，完整训练记录完成后保存到`models/training/evidence/submission-v2/`。

B2 将 ESM-2 的1280维蛋白表示通过可训练投影接入 scFoundation，并用酵母单细胞表达数据微调。最终筛选读取冻结 B2 权重；完整筛选不重新训练这个组件。

## 复跑训练

建议在 Linux、NVIDIA GPU 上另建 Python 3.10 环境，避免与 CPU 分析环境的 PyTorch 版本混用：

```bash
python3.10 -m venv .venv-training
source .venv-training/bin/activate
python -m pip install -r requirements-training.lock.txt
python scripts/install_training_assets.py --checkpoint /path/to/models.ckpt
python train.py --smoke --device cuda --output outputs/b2-smoke
python train.py --device cuda --output outputs/b2-training
```

安装器下载并校验准备好的训练输入。`models.ckpt` 是 scFoundation 的上游初始权重，从[官方模型说明](https://github.com/biomap-research/scFoundation/blob/main/model/README.md)提供的下载入口取得；本项目使用版本的 SHA-256 为 `9f40bf324d3d0084c4b288d06f5af4fddd12206e2a3f022551d12e89e33a0ea9`。安装器拒绝其他版本。离线时可另外传入 `--inputs /path/to/b2-inputs.tar.gz`。

`requirements-training.lock.txt`锁定Linux x86_64训练环境的33项传递依赖；新训练另保存实际`pip freeze`、操作系统、GPU、驱动、输入/输出哈希、开始结束时间及源码哈希。历史环境和本次训练环境分别记录。

默认训练使用 [b2_training_v2.json](../configs/b2_training_v2.json) 的完整配置。`--smoke` 仅用2个细胞做一次前向、反向更新并保存基因表，用于检验环境；不会产出可替代正式模型的权重。输出目录保留日志、配置、逐样本划分和运行环境，必须使用新目录。

2026-10-03已在 RTX A6000 上通过 smoke：Python 3.10、PyTorch 2.6.0、CUDA运行时12.4、驱动550.54.15；运行证据见 [smoke-20261003](training/evidence/smoke-20261003/)。这次检查不是完整6轮重训。早期同配方日志记录 B2 训练约161.7分钟，可作资源估计；实际耗时和显存取决于设备与批大小。

## 数据与配置

- 输入：GSE125162，38,225个细胞、6,733个基因；6,075个基因在原始矩阵中有对应列，其余填零。按 gene master 行顺序组装 ESM-2 矩阵，6个缺失蛋白表示填零。
- 预处理：每细胞 log1p(CPM-10k)，随机抽取最多1,200个表达基因，附加两个分辨率标记。原始准备脚本保存在 [original/scf_build_yeast_assets.py](training/original/scf_build_yeast_assets.py)。准备后的矩阵以文件哈希固定，见 [recovery.json](training/evidence/recovery.json)。
- 划分：`np.random.default_rng(42).permutation(38225)`，前1,911行为验证集，其余36,314行为训练集，无单独训练测试集。逐样本清单为 [split-reconstructed.tsv](training/evidence/split-reconstructed.tsv)，由留存代码和细胞顺序重建，不是当时写出的原始清单。
- 本次检查没有发现重复细胞编号、完全重复的计数行或跨训练/验证集的完全相同计数行。该划分按细胞随机抽取，未按实验批次分组；不据此声称消除了批次泄漏或与赛事隐藏集的重叠。
- 训练：6 epochs、batch 32、AdamW学习率1e-4、梯度裁剪1.0、seed 42、bf16、梯度检查点。冻结编码器前10层，训练末2层、数值嵌入、蛋白投影与解码器。完整参数见 [b2_training.json](../configs/b2_training.json)。

## 恢复了哪些材料

| 材料 | 位置与性质 |
| --- | --- |
| 原始训练/准备代码及配置 | `training/original/`，按恢复文件保留 |
| 可移植入口 | `train.py`、`scripts/scf_finetune_yeast.py`；增加路径参数、哈希检查、日志和划分导出 |
| scFoundation模型实现 | `training/vendor/scfoundation/`，含上游许可及改动说明 |
| 冻结模型 | `raw/tier1_models/route_b/final_model.pt`，SHA-256 `112b39af…9f2be7`，完整值见资产清单 |
| 权重一致性 | [checkpoint-consistency.json](training/evidence/checkpoint-consistency.json)：蛋白矩阵逐值一致，投影生成的基因表在CPU/GPU浮点容差内一致 |
| 当前服务器环境 | `environment-recovered-20261003.txt` 与 `recovery.json`，是恢复时环境，不是历史训练时的锁定记录 |
| 早期训练日志 | `training/evidence/earlier-run/scf_train.log`，包含2026-08-27的A2/B2/C2训练 |

早期 B2 基因表的 SHA-256 是 `c1cb4353552dc30b1a8134210729e563b8f1673d37042ad04b38acc850dda175`；当前冻结表是 `d5cd2912f7f798f236a41a31e2d23509402480415f74129e66a2c865dbae9c24`。两者不是同一次训练，因此不能把早期日志当作当前 checkpoint 的逐轮日志。当前权重对应训练的日志及当时环境快照仍未找到。

## submission-v2 的计算规则

- 从哈希锁定的scFoundation初始权重开始，未接续旧B2 checkpoint；训练6轮，不按最终候选筛选结果选择轮数或超参数。
- 按真实gene ID逐行排除分辨率标记与填充位；解码器收到相应padding mask。每个细胞最多抽取1,200个表达基因，遮蔽其中30%。由于不抽入未表达基因，原配置`zero_mask_p`在此表示构造中没有实际的零基因遮蔽对象。
- 验证使用`eval()`和`no_grad()`；以独立固定种子重放验证基因抽样与遮蔽，不消耗训练随机数流。验证集不参与反向更新。
- 每轮同时输出按批均值和按遮蔽条目数加权的训练/验证MSE。最终checkpoint内保存训练配置及逐轮指标，并以哈希关联外部日志和基因表。
- 训练输入和划分沿用冻结配方。固定随机种子不保证不同GPU、CUDA或PyTorch版本的训练位级一致；正式筛选使用已发布的最终权重和基因表。

## 历史实现的边界

使用`--config configs/b2_training.json`可显式运行历史配方。原实现验证时仍使用训练模式，并重新抽取基因和遮蔽；变长批次按最后两列排除分辨率位，短行的分辨率标记和填充位置未逐行单独排除。日志中的验证MSE应按这一实际实现解释，不能写成固定验证掩码下的独立评测。因此本次修正单独命名为`submission-v2`并重新训练、筛选，保留原权重及历史评估，不把旧结果归到新模型名下。

`scripts/model_selection/step4_route_confirmation/finetune_scgpt.py` 属于历史 scGPT 路线，不是本 B2 入口。合成演示中的小型 ridge 拟合也不替代 B2 训练。

## 从GEO原计数重建B2输入

准备矩阵之外还提供单独的B2数据处理入口，不要求重建A2/C2初始化。取得[原始计数压缩文件](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE125nnn/GSE125162/suppl/GSE125162_ALL-fastqTomat0-Counts.tsv.gz)，解压为TSV后运行：

```bash
python scripts/prepare_training_data.py \
  --counts-tsv /path/to/GSE125162_ALL-fastqTomat0-Counts.tsv \
  --output outputs/rebuilt-b2-inputs
```

程序按冻结gene master顺序排列6,733个基因，原计数中缺列者填零；蛋白表示按酵母ESM-2索引重排，缺失者填零。它分别验证原始TSV、基因表、蛋白嵌入、输出计数矩阵、蛋白矩阵和细胞顺序的SHA-256，输出`preprocessing.json`。此步骤不对已经冻结的输入另加QC或标准化；CPM-10k和log1p在训练中逐批执行。
