# YeastBridge

YeastBridge 将人和酵母的蛋白表示与酵母化学遗传学响应结合，用于 GPCR、离子通道等人源靶点的小分子候选筛选。仓库包含分析代码、冻结模型与数据的获取方式、候选清单和评估记录。

[候选清单](results.csv) · [模型说明](models/MODEL_CARD.md) · [数据来源](data/DATA_SOURCES.md) · [附件5核对记录](docs/submission/STATUS.md)

## 方法

```text
蛋白表示 → 三成员匹配与排名融合（EJ）→ 去除任务谱共享分量（EJ-dc）
        → HIP/HOP双向筛选 → 置换确认与家族分解 → 候选清单
```

EJ融合去PC1的ESM-2相似度、配对监督对齐和B2投影，使用 `k=20` 的倒数排名融合。筛选采用Spearman相关、置换检验和多重检验校正。输出是供后续实验检验的统计关联，不是已测得的亲和力或药效。

## 安装

使用 Python 3.10。导出和合成演示只需CPU：

```bash
git clone https://github.com/dawnmengsjtu/YeastBridge.git
cd YeastBridge
python3.10 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

Linux、macOS使用上述激活方式；Windows使用 `.venv\Scripts\activate`。核心依赖为NumPy 1.26.4、pandas 2.3.3、SciPy 1.15.3。完整分析另外需要anndata 0.11.4、h5py 3.12.1、PyTorch 2.9.0、fair-esm 2.0.0：

```bash
python -m pip install -r requirements-full.txt
```

完整筛选读取预计算嵌入，也可在CPU上运行，无需CUDA。传递依赖锁定在 [核心环境](requirements.lock.txt)和[完整分析环境](requirements-full.lock.txt)中。B2训练使用单独的GPU环境，见[训练说明](models/TRAINING.md)。

## 快速开始

### 导出已有候选清单

```bash
python predict.py
# 等效入口
bash run.sh
```

生成 `results.csv` 和 `results.metadata.json`。现有清单含80条记录、78个靶点、5个化合物；该命令读取已保存的确认结果，不重新筛选。也可指定 `--output outputs/published/results.csv`。

### 运行小规模演示

```bash
python main.py --mode demo --output outputs/demo
```

演示从固定种子的合成表示开始，完成ridge拟合、任务设计、双向筛选、提名、确认、家族分析和结果导出，通常数秒至一分钟。10个模拟靶点、160株模拟酵母、3个示例结构均由脚本提供，无需下载模型。输出明确标为合成演示，详见[演示说明](data/demo/README.md)。

### 运行完整筛选

```bash
python scripts/install_assets.py --download
python main.py --mode check
python main.py --mode full --output outputs/full
```

下载约513 MiB的B2权重和HIP/HOP响应矩阵并校验哈希，随后执行EJ导出、EJ-dc、双向全库筛选、每方向最多40对提名、100,000次置换确认、家族与残差检验。候选文件来自本次运行，保存在 `outputs/full/results.csv`；配置和日志保存在同一目录。

默认限制BLAS/OpenMP为4线程，可通过 `--threads` 调整。完整库每方向有3,818,750对假设，需按小时规划CPU运行时间，并预留约10 GB磁盘。缺失文件或哈希不符时在计算前退出。重复运行应换用新的输出目录。离线安装与附件哈希见[资产说明](data/ASSETS.md)，2026-10-03的完整CPU复核用时约33分钟（4线程，与其他检查并行），生成的80个候选对与历史清单一致；实测记录见[验证记录](docs/submission/validation.json)。

## 训练与数据准备

B2是在scFoundation骨干上训练的组件。已恢复原训练代码、计数矩阵、逐样本划分重建方法，并提供可移植入口：

```bash
# 在单独的Python 3.10训练环境中
python -m pip install -r requirements-training.lock.txt
python scripts/install_training_assets.py --checkpoint /path/to/models.ckpt
python train.py --smoke --device cuda --output outputs/b2-smoke
python train.py --device cuda --output outputs/b2-training
```

上游权重下载、模型许可、训练参数和历史记录的对应关系见[TRAINING.md](models/TRAINING.md)。单步训练检查已经通过；早期日志与当前冻结模型并非同一次训练，当前模型对应的逐轮日志仍需找回。

HIP/HOP可从已处理表达矩阵或原始CEL重建。全量表达矩阵重建后的8个数组字段已与冻结响应逐项核对一致，命令和R环境见[数据准备](data/PREPROCESSING.md)。

## 输出与结果

候选表包含编号、赛道、靶点、InChIKey、PubChem CID、SMILES、方向、相关系数、经验p值、BH校正值、剂量及单位、证据等级、模型信息和备注。现有80条记录的剂量单位均为micromolar（µM）。字段、精度、排序和缺失值规则见[输出说明](docs/submission/RESULTS_SCHEMA.md)。

`+z` 表示耐药富集方向，`-z` 表示超敏方向。默认导出保留历史统计值；新运行对两个方向合并后的确认家族计算BH。两者均未校正前序提名选择，不应把确认家族内的q值解释为全库发现FDR。

主要评估：

- [跨物种检索](docs/cross_species_match/V8_ESM2_FORMAL_RESULT.md)：去PC1的ESM-2 test AUC为0.9638，配对监督联合臂为0.8770。
- [基本结构基准](docs/cross_species_match/V9_ESM2_BASIC_RESULT.md)：联合方法MRR为0.403，human-only为0.365；两者比较 `p=0.0513`，未达0.05阈值。
- [两阶段确认](docs/cross_species_match/ESM2_JOINT_DC_ADDENDUM.md)、[家族特异性](docs/cross_species_match/FAMILY_SPECIFICITY_RESULT.md)和[外部验证](docs/cross_species_match/EV_PROTOCOL.md)。

## 验证与打包

```bash
python -m unittest discover -s tests -v
python scripts/build_submission.py --profile review
python scripts/build_submission.py --profile inference --output dist/YeastBridge-inference.zip
```

可直接下载[精简复核包](https://github.com/dawnmengsjtu/YeastBridge/releases/download/reproducibility-assets-20261003/YeastBridge-review.zip)或[离线推理包](https://github.com/dawnmengsjtu/YeastBridge/releases/download/reproducibility-assets-20261003/YeastBridge-inference.zip)。

`review`是代码、文档、导出与演示的精简包；`inference`额外包含正式筛选所需冻结输入，可离线运行分析。包内含SHA-256清单。完整参赛包的材料核验使用 `--profile full`，当前会提示尚缺与最终权重对应的训练日志。数据许可记录和赛道最终模板的核对状态见[附件5核对记录](docs/submission/STATUS.md)。

## 目录

| 路径 | 内容 |
| --- | --- |
| `main.py`、`run.sh`、`predict.py` | 统一入口与候选导出 |
| `train.py`、`models/` | B2训练入口、模型说明与训练溯源 |
| `scripts/`、`configs/` | 分析和数据准备代码、配置 |
| `data/`、`raw/` | 数据说明、来源记录与冻结输入 |
| `tasks/`、`panels/` | 历史任务表和候选面板 |
| `results_a6000/`、`results_model_selection/`、`repro/` | 历史结果和运行记录 |
| `docs/`、`tests/` | 协议、结果、勘误和验证 |

`v5/v7/v8`等模块仍被最终方法调用；历史协议和勘误用于溯源。打包时排除临时缓存和不需要的中间文件。

引用见 [CITATION.cff](CITATION.cff)，第三方许可见 [NOTICE.md](NOTICE.md)。复现问题可提交到 [Issues](https://github.com/dawnmengsjtu/YeastBridge/issues)。
