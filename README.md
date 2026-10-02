# YeastBridge

YeastBridge 是面向 GPCR、离子通道等人源靶点的小分子候选筛选项目。它将人和酵母的蛋白表示进行跨物种匹配，再结合酵母化学遗传学响应数据，为后续人源化酵母及其他实验验证筛选候选化合物。

仓库包含分析代码、配置、评估记录和候选清单。可以直接查看 [results.csv](results.csv)，也可以运行 `predict.py` 重新导出清单；完整分析需要另外准备模型权重和原始数据。

## 工作流程

```text
人源靶点与酵母蛋白表示
        ↓
跨物种匹配与排序融合（EJ）
        ↓
去除共享分量，生成靶点对应的酵母基因排序（EJ-dc）
        ↓
与 HIP/HOP 化合物响应谱匹配，进行双向筛选
        ↓
两阶段统计确认与家族特异性分析
        ↓
候选清单 results.csv
```

EJ 使用倒数排名融合（RRF，`k=20`）整合三种匹配结果：去除第一主成分后的 ESM-2 相似度、配对监督对齐，以及 B2 投影。EJ-dc 在此基础上去除不同靶点间共享的任务谱分量。筛选层使用 Spearman 相关、置换检验和多重检验校正评估靶点与化合物响应谱的关联。

## 快速开始

### 导出已有候选清单

使用 Python 3.10，在终端运行：

```bash
git clone https://github.com/dawnmengsjtu/YeastBridge.git
cd YeastBridge

python3.10 -m venv .venv
source .venv/bin/activate
python -m pip install numpy==1.26.4 pandas==2.3.3

python predict.py
```

脚本读取仓库内已保存的确认结果、家族特异性结果和化合物信息，生成或覆盖根目录的 `results.csv`。这一步不需要下载模型权重，也不会重新执行筛选。

当前输出为 80 条候选记录，涉及 78 个靶点和 5 个化合物，包含以下信息：

| 字段 | 内容 |
| --- | --- |
| 候选编号、靶点 | 候选 ID 与人源靶点名称 |
| 化合物_InChIKey、化合物_CID、SMILES | 化合物标识与结构 |
| 方向、剂量 | 酵母响应方向及原始记录中的剂量 |
| spearman_rho、emp_p、q_bh | 相关系数、经验 p 值和 BH 校正值 |
| 证据等级、模型与版本、备注 | 统计证据分类及方法说明 |

`+z` 表示耐药富集方向，`-z` 表示超敏方向。证据等级按残差检验和任务谱聚类结果划分，用于区分靶点级与家族级的统计关联；候选的结合活性和作用机制仍需实验验证。

### 运行完整分析

完整分析使用 Python 3.10，建议在 Linux 下运行。使用预计算嵌入时，分析链可在 CPU 上运行。除上述依赖外，还需安装：

```bash
python -m pip install scipy==1.15.3 anndata==0.11.4 h5py==3.12.1 \
    torch==2.9.0 fair-esm==2.0.0
```

运行前需要完成数据和路径配置：

1. 按 [数据说明](data/DATA_SOURCES.md) 和 [模型说明](models/MODEL_CARD.md) 准备嵌入、B2 投影权重、表达数据及 HIP/HOP 响应矩阵。大文件未随仓库发布，部分 `raw/` 条目是指向原归档的符号链接，需要替换为本机可访问的文件。
2. 对照配置中的 SHA-256 和 `MANIFEST.sha256` 核对所需输入。该清单记录的是原始运行环境中的文件，公开仓库不包含清单中的全部内容。
3. 将 `configs/` 中的 `stage_root` 和输入、输出路径改为本机路径，同时检查被引用的配置，例如 `esm2_joint_tasks.json` 中的 `method.v7_config`。`configs_a6000_frozen/` 保存原运行配置，供结果核对。
4. 将正、负方向执行配置中的 `task_dir` 指向下方生成的 `repro/tasks_ej_dc/`，并分别设置 `response_npz`、`compound_table` 和 `results_dir`。

完成配置后，依次运行：

```bash
# 生成 EJ 任务表
python scripts/esm2_joint_task_export.py \
    --config configs/esm2_joint_tasks.json \
    --output repro/tasks_ej

# 去除共享分量
python scripts/task_axis_double_center.py \
    --input repro/tasks_ej \
    --output repro/tasks_ej_dc

# 耐药富集与超敏两个方向的筛选
python scripts/product_execute_hiphop.py \
    --config configs/product_execute.esm2joint_dc.json \
    --results-suffix _screen_repro

python scripts/product_execute_hiphop.py \
    --config configs/product_execute.esm2joint_dc_neg.json \
    --results-suffix _neg_screen_repro
```

两阶段确认的提名规则、置换次数和参数见 [EJ-dc 分析说明](docs/cross_species_match/ESM2_JOINT_DC_ADDENDUM.md)；家族分解见 [家族特异性分析协议](docs/cross_species_match/FAMILY_SPECIFICITY_PROTOCOL.md)。导出新一轮候选前，需要完成这些分析，并相应更新 `predict.py` 的结果目录及配套面板、化合物信息。

仓库中的 `run.sh` 保留了原运行环境的串联命令，使用 `configs_a6000_frozen/` 下的配置。迁移到新机器时需要同步修改路径；其最后一步仍由 `predict.py` 读取已有确认结果。

## 评估与结果

| 评估任务 | 已保存结果 | 详细记录 |
| --- | --- | --- |
| 跨物种检索 | ESM-2 去 PC1 的 test AUC 为 0.9638；配对监督联合臂为 0.8770 | [V8 结果](docs/cross_species_match/V8_ESM2_FORMAL_RESULT.md) |
| 基本结构基准 | 联合方法 MRR 为 0.403，human-only 为 0.365，yeast-only 为 0.356 | [V9 结果](docs/cross_species_match/V9_ESM2_BASIC_RESULT.md) |
| 候选确认与特异性分析 | 80 条候选记录；残差检验、任务谱聚类和家族均值分析用于区分证据等级 | [确认结果](docs/cross_species_match/ESM2_JOINT_DC_ADDENDUM.md)、[特异性结果](docs/cross_species_match/FAMILY_SPECIFICITY_RESULT.md) |
| 外部验证 | EV0 跨物种检索与 EV1 基本结构评估 | [验证协议](docs/cross_species_match/EV_PROTOCOL.md)、[运行结果](results_a6000/externalvalidation/) |

V9 中联合方法相对 human-only 的比较为 `p=0.0513`，未达到 0.05 阈值。候选确认使用同一响应数据上的两阶段筛选，家族内 BH 校正未覆盖第一阶段的选择多重度，不能将 80 条候选解释为 80 个独立验证的药物靶点发现。详细统计口径和效应量见对应结果文档。

## 仓库结构

```text
predict.py               从已有确认结果导出候选清单
results.csv              已导出的候选清单
run.sh                   原环境下的分析链入口
scripts/                 匹配、筛选、验证与模型选型脚本
configs/                 分析配置，运行前需适配本机路径
configs_a6000_frozen/     原运行配置
data/                   数据来源说明与化合物信息
models/                  模型来源与权重说明
raw/                     输入数据、索引和归档链接
tasks/                  靶点对应的酵母基因任务表
panels/                  筛选面板、候选列表及构建记录
results_a6000/           主分析与外部验证的运行结果
results_model_selection/ 模型选型结果
docs/                   方法、分析协议与结果说明
repro/                   复现输出
MANIFEST.sha256          原始运行文件的哈希清单
```

## 文档与引用

- [跨物种匹配文档](docs/cross_species_match/README.md)
- [EJ 方法与参数](docs/cross_species_match/ESM2_JOINT_METHOD_FREEZE.md)
- [模型选型结果](docs/model_selection/RESULTS_v3.md)及[勘误](docs/model_selection/ERRATUM.md)
- [数据来源](data/DATA_SOURCES.md)与[模型说明](models/MODEL_CARD.md)

引用本项目请使用 [CITATION.cff](CITATION.cff)。问题或复现反馈可提交到 [Issues](https://github.com/dawnmengsjtu/YeastBridge/issues)。
