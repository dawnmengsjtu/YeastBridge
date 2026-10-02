# YeastBridge

YeastBridge 面向 GPCR、离子通道等人源靶点，将跨物种蛋白匹配与酵母化学遗传学响应结合，用于小分子候选筛选。仓库包含分析代码、冻结参数、候选清单和评估记录。

[候选清单](results.csv) · [模型说明](models/MODEL_CARD.md) · [数据来源](data/DATA_SOURCES.md) · [参赛提交状态](docs/submission/STATUS.md)

目前可直接运行已发表于本仓库的候选清单导出和合成演示。完整科研流程仍需补齐 B2 权重与 HIP/HOP 响应矩阵；B2 训练材料和部分数据许可记录也尚待恢复，详见提交状态。

## 方法

```text
人和酵母的蛋白表示
    → 三成员匹配与 RRF 融合（EJ）
    → 去除任务谱共享分量（EJ-dc）
    → HIP/HOP 双向响应筛选
    → 候选提名、置换确认与家族分解
    → results.csv
```

EJ 融合去 PC1 的 ESM-2 相似度、配对监督对齐和 B2 投影，RRF 参数为 `k=20`。筛选使用 Spearman 相关、置换检验和多重检验校正。结果用于安排后续实验，指标不代表已测得的结合亲和力或药效。

## 安装

使用 Python 3.10。导出和合成演示在 CPU 上运行，无需 CUDA 或模型下载。Linux 和 macOS 可使用以下命令；Windows 激活环境时使用 `.venv\Scripts\activate`。

```bash
git clone https://github.com/dawnmengsjtu/YeastBridge.git
cd YeastBridge
python3.10 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

核心依赖为 NumPy 1.26.4、pandas 2.3.3、SciPy 1.15.3。测试环境中的传递依赖见 [requirements.lock.txt](requirements.lock.txt)。正式 EJ 导出另需 anndata 0.11.4、h5py 3.12.1、torch 2.9.0 和 fair-esm 2.0.0：

```bash
python -m pip install -r requirements-full.txt
```

已验证的完整分析环境版本锁定见 [requirements-full.lock.txt](requirements-full.lock.txt)。

完整流程读取预计算嵌入，分析阶段不要求 GPU。历史模型训练的硬件、驱动、完整依赖和耗时记录尚待补齐；现有旧训练脚本不包含在已验证运行范围内。

## 运行

### 导出提交清单

```bash
python predict.py
# 等效主入口
bash run.sh
```

输出根目录的 `results.csv` 和 `results.metadata.json`。清单为80条记录，涉及78个靶点、5个化合物；元数据记录输入哈希、代码版本、排序和统计口径。此命令读取已有确认结果，不重新筛选。

可以指定其他输出文件：

```bash
python predict.py --output outputs/published/results.csv
```

### 运行小规模演示

```bash
python main.py --mode demo --output outputs/demo
```

演示从固定种子的合成表示开始，完成小型 ridge 拟合、任务设计、共模去除、双向筛选、提名、确认、家族分析和结果导出。输入、拟合参数、配置、日志和结果均保存到输出目录，通常为数秒至一分钟，实际耗时随机器与首次导入依赖而变化。

演示包含10个模拟靶点、160株模拟酵母和3个示例结构，不使用正式模型权重。所有候选均标为“合成演示，非科研候选”。详见 [演示说明](data/demo/README.md)。重复运行请使用新的输出目录，程序不会混用旧运行结果。

### 运行完整筛选

先按 [资产说明](data/ASSETS.md) 安装原始冻结文件，再检查输入：

```bash
python main.py --mode check
python main.py --mode full --output outputs/full
```

完整流程会执行 EJ 导出、EJ-dc、双向全库筛选、每方向最多40对提名、100,000次置换确认、家族和残差检验，最后从本次运行生成候选清单。缺失文件或哈希不符时在计算前退出，不会用已保存结果替代新计算。

完整库每个方向约382万对，运行量明显大于演示。原环境说明为单机 CPU 数小时；本次尚未在完整数据上验证该耗时和峰值内存，不应按演示资源估算正式任务。

## 输出与结果

清单包含候选编号、赛道、靶点、InChIKey、PubChem CID、SMILES、方向、Spearman 相关、经验 p 值、BH 校正值、剂量、证据等级、模型说明和备注。字段类型、精度、排序及缺失值规则见 [输出说明](docs/submission/RESULTS_SCHEMA.md)。

`+z` 为耐药富集方向，`-z` 为超敏方向。默认导出保留历史数值；新运行对两个方向合并后的确认家族计算 BH。两者均未校正前序提名选择，候选不等于独立验证的药物靶点发现。原始剂量单位尚待核实。

主要评估记录：

- [跨物种检索](docs/cross_species_match/V8_ESM2_FORMAL_RESULT.md)：去 PC1 的 ESM-2 test AUC 为0.9638，配对监督联合臂为0.8770。
- [基本结构基准](docs/cross_species_match/V9_ESM2_BASIC_RESULT.md)：联合方法 MRR 为0.403，human-only 为0.365，yeast-only 为0.356；与 human-only 比较的 `p=0.0513` 未达0.05阈值。
- [两阶段确认](docs/cross_species_match/ESM2_JOINT_DC_ADDENDUM.md)与[家族特异性](docs/cross_species_match/FAMILY_SPECIFICITY_RESULT.md)。
- [外部验证协议](docs/cross_species_match/EV_PROTOCOL.md)与[运行结果](results_a6000/externalvalidation/)。

## 验证与打包

```bash
python -m unittest discover -s tests -v
python scripts/build_submission.py --profile review
```

`dist/YeastBridge-review.zip` 是用于候选导出与合成演示的精简复核包，含实际文件的 SHA-256 清单。它省去大体积数据、重复嵌入、历史中间任务和绘图取数表，不作为完整科研提交包。

完整包使用 `--profile full`。在关键资产或 B2 训练证据不齐时，打包器会拒绝生成完整包。附件5各项现状与剩余材料见 [提交状态](docs/submission/STATUS.md)。

## 仓库结构

| 路径 | 内容 |
| --- | --- |
| `main.py`、`run.sh` | 统一运行入口 |
| `predict.py`、`results.csv` | 候选导出与已提交清单 |
| `scripts/` | 分析、校验、资产安装和打包代码 |
| `requirements*.txt`、`tests/` | 环境与验证 |
| `data/`、`models/` | 数据、演示和模型说明 |
| `configs/` | 主流程相对路径配置，以及历史基准配置 |
| `configs_a6000_frozen/` | 原环境配置，仅作溯源 |
| `raw/`、`tasks/`、`panels/` | 已保存输入、任务表与面板 |
| `results_a6000/`、`results_model_selection/`、`repro/` | 历史结果及复现记录 |
| `docs/` | 协议、结果、勘误和提交说明 |

`v5/v7/v8` 等模块仍被最终方法导入，不能按版本名称删除。历史协议、相反证据与勘误保留用于溯源；临时备份和绘图中间表可从 Git 历史恢复。

引用信息见 [CITATION.cff](CITATION.cff)，第三方许可说明见 [NOTICE.md](NOTICE.md)。问题与复现反馈请提交到 [Issues](https://github.com/dawnmengsjtu/YeastBridge/issues)。
