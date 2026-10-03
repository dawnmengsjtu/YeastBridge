# B2 训练与溯源

当前提交版本 `B2-submission-v2-20261003` 从官方scFoundation初始权重完成6轮酵母表达微调。完整训练证据位于 [submission-v2/](training/evidence/submission-v2/)，权重、逐轮日志、样本划分、输入和实际环境由 [training_manifest.json](training_manifest.json)关联。

## 复跑训练

建议使用Linux、Python 3.10及支持bfloat16、至少24 GB显存的NVIDIA GPU。训练环境与CPU筛选环境分别安装：

```bash
python3.10 -m venv .venv-training
source .venv-training/bin/activate
python -m pip install -r requirements-training.lock.txt
python scripts/install_training_assets.py --checkpoint /path/to/models.ckpt
python train.py --smoke --device cuda --output outputs/b2-smoke
python train.py --device cuda --output outputs/b2-training
```

`models.ckpt`从[scFoundation官方模型说明](https://github.com/biomap-research/scFoundation/blob/main/model/README.md)给出的入口获取。本配方初始权重的SHA-256为 `9f40bf324d3d0084c4b288d06f5af4fddd12206e2a3f022551d12e89e33a0ea9`，安装器拒绝其他版本。离线训练输入可通过 `--inputs /path/to/b2-inputs.tar.gz` 安装。重训前阅读随仓库提供的上游 [MODEL_LICENSE](training/vendor/scfoundation/MODEL_LICENSE)。

训练锁定文件列出33项依赖。实测Python 3.10.20、PyTorch 2.6.0、CUDA运行时12.4、Linux 3.10.0-1160系列/glibc 2.17、驱动550.54.15、单张RTX A6000 48 GB；实际完整环境见 [requirements-observed.txt](training/evidence/submission-v2/requirements-observed.txt)。本次训练用时160.3分钟，观察到约15 GB显存占用。固定种子不保证不同硬件或CUDA版本下训练位级一致；筛选复核使用发布的最终权重与基因表。

`--smoke`只取2个细胞进行一次前向、反向更新及基因表保存，不产出可替代正式模型的checkpoint。每次输出目录必须新建；正式训练保存逐轮指标、日志、逐样本划分、环境、配置及所有产物哈希。

## 数据与训练规则

- GSE125162：38,225个细胞、6,733个基因；6,075列可对应原计数，其余填零。蛋白矩阵按gene master排序，6个缺失表示填零。
- 每细胞CPM-10k后log1p；均匀抽取最多1,200个表达基因，附两个分辨率标记。逐行依据真实gene ID排除padding和特殊位，对表达基因以0.3概率遮蔽。由于未抽入零表达基因，配置中的`zero_mask_p=0.03`在此构造中没有实际遮蔽对象。
- `np.random.default_rng(42).permutation(38225)`：前1,911行验证，其余36,314行训练。实际逐样本清单见 [split.tsv](training/evidence/submission-v2/split.tsv)。无单独训练测试集。
- 6轮、batch 32、AdamW lr=1e-4、梯度裁剪1.0、seed 42、bf16及梯度检查点。冻结12层编码器中的前10层，更新末2层、数值嵌入、蛋白投影和解码器，约34.6M可训练参数。完整参数见 [b2_training_v2.json](../configs/b2_training_v2.json)。
- 验证使用`eval()`/`no_grad()`；独立固定种子重放基因抽样与遮蔽，并保留训练随机数流。验证不反向更新。采用第6轮checkpoint，未根据最终候选结果挑选轮数或超参数。

## 本次逐轮结果

下表MSE按遮蔽条目数量加权；原始文件另保存按批均值及遮蔽条目数。每轮时间包含训练和验证。

| 轮次 | 训练MSE | 验证MSE | 分钟 |
| --- | --- | --- | --- |
| 1 | 0.247347 | 0.186616 | 26.63 |
| 2 | 0.182002 | 0.167740 | 26.61 |
| 3 | 0.163296 | 0.154122 | 26.67 |
| 4 | 0.152610 | 0.146252 | 26.75 |
| 5 | 0.145459 | 0.142013 | 26.70 |
| 6 | 0.141101 | 0.137655 | 26.63 |

训练日志使用服务器时钟，筛选日志使用本机时钟；本次观测到服务器快约352秒，见[时钟核对](../data/provenance/recovered-20261003/clock-offset.json)。原始时间戳未改写，跨机器依赖关系按成功完成状态和模型哈希核验。

完整日志见 [train.log](training/evidence/submission-v2/train.log)。[checkpoint-consistency.json](training/evidence/submission-v2/checkpoint-consistency.json)核对了checkpoint内蛋白矩阵与输入逐值一致，以及CPU重算投影表与训练时保存表在1e-5容差内一致。

## 从原始计数重建输入

从[GEO](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE125nnn/GSE125162/suppl/GSE125162_ALL-fastqTomat0-Counts.tsv.gz)取得并解压原始计数，然后运行：

```bash
python scripts/prepare_training_data.py \
  --counts-tsv /path/to/GSE125162_ALL-fastqTomat0-Counts.tsv \
  --output outputs/rebuilt-b2-inputs
```

处理器按冻结gene master排列基因，计数与蛋白表示缺失项填零，核验原始TSV、gene master、ESM-2输入及输出哈希。2026-10-03已完整执行，计数矩阵、蛋白矩阵和细胞清单三个输出的SHA-256全部与训练输入一致，见[重建证据](../data/provenance/recovered-20261003/b2-preprocessing-rebuild.json)。原gene master的CRLF版本与Git中的LF版本经逐字节换行比较一致，处理器只允许这两个已知哈希。

[划分审计](../data/provenance/recovered-20261003/split-audit.json)未发现重复细胞编号、完全重复计数行或跨集完全重复行；细胞随机划分没有按批次隔离。配对监督的753对按人/酵母实体分组，训练、验证、测试无实体交叉；新筛选只使用其中具备表示的510个训练配对。

## 将新权重接入筛选

在已安装完整分析依赖的环境中运行：

```bash
python scripts/configure_trained_model.py --training-run outputs/b2-training --output outputs/b2-config
python main.py --mode full --task-config outputs/b2-config/task.json --output outputs/b2-screen
```

配置器拒绝未完成训练及smoke输出，核验权重、基因表、划分和逐轮记录，生成可移植配置。最终结果保存在 `outputs/b2-screen/results.csv`，模型版本与本次训练对应。使用同一工作目录切换两个虚拟环境即可共享输入和产物。

## 历史实现

`configs/b2_training.json`保留`historical-v1`配方，9月权重和原始结果单独归档。原实现验证仍使用训练模式，并重新抽取基因和遮蔽；变长批次按最后两列排除分辨率位，短行的特殊位和padding未逐行排除。本次将修正命名为`submission-v2`，从上游初始权重重训并重新筛选。

早期8月日志与9月冻结基因表不同，不能当作9月checkpoint的日志。原代码、早期日志、恢复环境和划分重建记录继续保留在 `training/original/` 与 `training/evidence/`，但不充当本次训练证据。历史scGPT训练及合成ridge演示也不替代本B2训练。
