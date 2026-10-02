# HIP/HOP 数据准备与复核

主筛选读取冻结响应矩阵，不需要重新处理原始CEL。响应矩阵的来源链是：

```text
E-MTAB-2391的18个CEL压缩包 + 平台探针坐标表
  → log2、分位数归一化、每探针组median polish
  → 同批次vehicle中位数差值
  → 同菌株探针组均值
  → strain_response.npz
```

原实现没有背景卷积校正。NPZ键名虽然为 `z_score`，内容实际是条形码表达的vehicle对比值，不是除以标准差得到的标准z分数。批次匹配按数组文件名第一个下划线前的标记分组；该标记不能自动等同于完整实验批次。vehicle自身使用全局vehicle中位数。

## 从已处理表达矩阵重建

```bash
python scripts/prepare_hiphop.py --output outputs/hiphop-rebuild
```

命令下载带SHA-256的 `expr_full.tsv.gz`，结合仓库SDRF生成响应文件和溯源JSON。离线时使用 `--expression /path/to/expr_full.tsv.gz`。需要核心Python依赖。

2026-10-03已在原服务器执行重建，8个数组字段逐项完全一致，`z_score` 最大绝对误差为0。NPZ压缩容器的字节哈希可随压缩实现不同；数组内容比较见 [rebuild-comparison-20261003.json](provenance/hiphop/rebuild-comparison-20261003.json)。主链仍使用原冻结文件及其哈希。

## 从原始CEL重建

```bash
conda env create -f environment-preprocessing.yml
conda activate yeastbridge-hiphop
# 在已安装核心Python依赖的环境中调用此Rscript
python scripts/prepare_hiphop.py --from-cel --rscript /path/to/yeastbridge-hiphop/bin/Rscript \
    --output outputs/hiphop-from-cel
```

命令按 [原始归档清单](provenance/hiphop/hiphop_readiness_final.json)下载18个官方压缩包、逐包检查哈希，并按冻结条件表的顺序处理3,850个CEL。原始压缩包约2.8 GB，解压、R矩阵及输出需要额外磁盘和内存；建议至少16 GB可用内存、10 GB可用磁盘。完整CEL处理的本次耗时未测量；本次通过了12数组的R流程测试和全量表达矩阵到响应矩阵的复核。

R环境：4.2.3、affyio 1.68.0、preprocessCore 1.60.2。记录见 [R环境](provenance/hiphop/r-environment-recovered-20261003.txt)和 [CEL测试](provenance/hiphop/cel-smoke-summary-20261003.txt)。R脚本按原实现保留；Python移植版增加输出保护，并修复带显式行名列表头的兼容性，没有更改冻结输入的计算结果。

## 来源与条件信息

- 官方数据：[E-MTAB-2391](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-2391)，Lee et al., DOI [10.1126/science.1250217](https://doi.org/10.1126/science.1250217)。
- 原始记录：SDRF、IDF、HTTP响应头、归档哈希、处理记录在 `data/provenance/hiphop/`。HTTP响应头记录获取响应的时间；没有将文件修改时间当作下载时间。
- 3,850个条件包含494个vehicle、3,356个处理条件，涉及3,250个化合物、5,668株酵母。`data/response_conditions.tsv` 从冻结NPZ提取条件、化合物、剂量与单位，并保留源文件哈希。
- 原始库的剂量单位有micromolar、nanomolar、millimolar、picomolar及percent。现有80条候选均为micromolar；导出程序按化合物和剂量匹配单位，遇到缺失或歧义时停止。
- 冻结数据的具体再分发许可仍需结合来源条款核对；原始项目的许可记录为unknown，没有据公开可下载推定为CC-BY。
