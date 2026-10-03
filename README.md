# YeastBridge

YeastBridge 将人和酵母的蛋白表示与酵母化学遗传学响应结合，用于 GPCR、离子通道等人源靶点的小分子候选筛选。仓库提供训练、任务设计、双向筛选和结果导出的完整流程。

[候选清单](results.csv) · [模型说明](models/MODEL_CARD.md) · [数据来源](data/DATA_SOURCES.md) · [完整提交包](https://github.com/dawnmengsjtu/YeastBridge/releases/download/submission-v2-20261003/YeastBridge-full.zip) · [附件5核对](docs/submission/STATUS.md)

## 方法

```text
ESM-2蛋白表示 + B2酵母表达微调
    → 三成员匹配与排名融合（EJ）
    → 去除任务谱共享分量（EJ-dc）
    → HIP/HOP双向筛选 → 置换确认 → 家族与残差分析 → 候选清单
```

EJ以 `k=20` 的倒数排名融合整合去PC1的ESM-2相似度、配对监督对齐和B2投影。筛选使用Spearman相关、置换检验和多重检验校正。输出是待实验检验的统计关联。

当前发布版本为 `B2-submission-v2-20261003`，使用从scFoundation初始权重微调6轮的B2模型。训练配置、模型组成和适用范围见[模型说明](models/MODEL_CARD.md)。

## 安装

使用Python 3.10。导出和小规模演示只需CPU：

```bash
git clone https://github.com/dawnmengsjtu/YeastBridge.git
cd YeastBridge
python3.10 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

上述命令适用于Linux和macOS。核心依赖为NumPy 1.26.4、pandas 2.3.3、SciPy 1.15.3。完整筛选另安装：

```bash
python -m pip install -r requirements-full.txt
```

完整环境包括PyTorch 2.9.0、anndata 0.11.4、h5py 3.12.1和fair-esm 2.0.0。本次CPU环境的全部传递依赖见 [核心锁定文件](requirements.lock.txt)和[完整分析锁定文件](requirements-full.lock.txt)。完整筛选读取预计算表示，可在CPU运行，无需CUDA；实测环境为macOS 26.2、Apple M1、8 GB内存，线程上限设为4；系统数值库见[运行环境](repro/submission-v2-20261003/runtime-environment.json)。建议预留10 GB磁盘。GPU训练使用单独环境，见下文。

## 快速开始

### 导出提交清单

```bash
python predict.py
# 等效入口
bash run.sh
```

生成UTF-8 `results.csv` 和溯源文件 `results.metadata.json`。默认读取当前提交版本的冻结确认记录，生成80条记录，涉及78个靶点、5个化合物。指定 `--output outputs/published/results.csv` 可改变输出位置。

### 运行小规模演示

```bash
python main.py --mode demo --output outputs/demo
```

固定种子的合成示例包含10个模拟靶点、160株模拟酵母和3个示例结构，完成ridge拟合、任务设计、双向筛选、确认、家族分析及导出，通常数秒至一分钟。无需下载模型；输出明确标为合成演示，见[演示说明](data/demo/README.md)。

### 从模型重新筛选

```bash
python scripts/install_assets.py --download
python main.py --mode check
python main.py --mode full --output outputs/full
```

安装器下载当前B2权重与HIP/HOP响应矩阵并校验SHA-256。完整流程对每方向3,818,750对假设筛选，每方向最多提名40对，执行100,000次置换确认、家族与残差检验，再输出 `outputs/full/results.csv`。本次14步CPU复跑用时约17.8分钟；实测配置、日志和结果见 [运行记录](repro/submission-v2-20261003/run.json)。

运行配置、日志、候选清单和置换不确定性区间都在输出目录中。默认将线程上限设为4，可通过 `--threads` 调整。每次运行须使用新目录。离线安装见[资产说明](data/ASSETS.md)。

## 训练与数据准备

B2将ESM-2的1280维蛋白表示投影到768维，并在scFoundation骨干上使用酵母表达数据微调。Linux/NVIDIA训练环境使用Python 3.10、PyTorch 2.6.0、CUDA运行时12.4；本次在RTX A6000上训练6轮，用时约160.3分钟，驱动550.54.15，运行中观察到约15 GB显存占用。训练GPU须支持bfloat16，建议至少24 GB显存；为训练环境、初始权重和产物预留20 GB磁盘。

```bash
# 在单独的Python 3.10训练环境中
python -m pip install -r requirements-training.lock.txt
python scripts/install_training_assets.py --checkpoint /path/to/models.ckpt
python train.py --device cuda --output outputs/b2-training
python scripts/configure_trained_model.py --training-run outputs/b2-training --output outputs/b2-config
python main.py --mode full --task-config outputs/b2-config/task.json --output outputs/b2-screen
```

最后一条筛选命令在上述完整分析环境中执行。scFoundation初始权重的官方获取入口、许可、完整训练参数、划分及逐轮指标见[训练说明](models/TRAINING.md)。安装器也支持离线训练输入。`--smoke`可先检查训练环境，不替代正式训练。

[数据准备](data/PREPROCESSING.md)说明GSE125162计数、ESM-2蛋白表示和HIP/HOP响应矩阵的重建方法，并提供输入版本、处理命令及重建核对记录。

## 输出与解释

候选表包含编号、赛道、靶点、InChIKey、PubChem CID、SMILES、方向、相关系数、经验p值、BH校正值、剂量及单位、证据等级、模型版本和备注。字段、精度、排序与缺失值规则见[输出说明](docs/submission/RESULTS_SCHEMA.md)。

`+z`为耐药富集方向，`-z`为超敏方向。两个方向的提名记录合并计算BH；确认复用筛选数据，q值不包含前序提名选择，不能解释为全库错误发现率（FDR）。置换区间量化有限置换的抽样误差，不代表亲和力、药效或作用机制的置信区间。

历史清单可用 `python predict.py --legacy --output outputs/legacy/results.csv` 导出。原始跨物种检索、基本结构基准、家族特异性及外部验证保存在 [docs/cross_species_match/](docs/cross_species_match/)；其中的9月评估值对应历史版本，不作为本次B2新模型的性能成绩。

## 验证与打包

```bash
python scripts/verify_submission.py
python -m unittest discover -s tests -v
python scripts/build_submission.py --profile full
```

`verify_submission.py`核对最终权重、逐轮日志、筛选记录和候选表。[完整提交包](https://github.com/dawnmengsjtu/YeastBridge/releases/download/submission-v2-20261003/YeastBridge-full.zip)包含源码、当前权重、筛选输入、训练输入归档及日志；重训所需scFoundation初始权重按其官方说明另行获取。[精简复核包](https://github.com/dawnmengsjtu/YeastBridge/releases/download/submission-v2-20261003/YeastBridge-review.zip)用于阅读、导出和运行合成演示。每个包都有逐文件SHA-256清单。

比赛材料与附件5的对应关系见[核对记录](docs/submission/STATUS.md)。候选表采用附件5允许的默认CSV，依据见[格式说明](docs/submission/TEMPLATE_BASIS.md)。

## 目录

| 路径 | 内容 |
| --- | --- |
| `main.py`、`run.sh`、`predict.py` | 统一运行入口和候选导出 |
| `train.py`、`models/` | B2训练、模型说明及对应训练证据 |
| `scripts/`、`configs/` | 核心分析、数据准备代码和配置 |
| `data/`、`raw/` | 数据来源、许可、冻结输入 |
| `repro/submission-v2-20261003/` | 当前提交版本的筛选日志和结果 |
| `tasks/`、`panels/`、`results_a6000/`、`results_model_selection/` | 历史任务、面板和评估记录 |
| `docs/`、`tests/` | 方法协议、结果解释和验证 |

引用见 [CITATION.cff](CITATION.cff)，第三方模型与数据许可见 [NOTICE.md](NOTICE.md)。复现问题可提交到 [Issues](https://github.com/dawnmengsjtu/YeastBridge/issues)。
