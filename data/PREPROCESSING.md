# 数据准备与复核

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

# 保持前面安装了requirements.txt的Python环境，指定R环境中的Rscript路径
python scripts/prepare_hiphop.py --from-cel --rscript /path/to/yeastbridge-hiphop/bin/Rscript \
    --output outputs/hiphop-from-cel
```

命令按 [原始归档清单](provenance/hiphop/hiphop_readiness_final.json)下载18个官方压缩包、逐包检查哈希，并按冻结条件表的顺序处理3,850个CEL。原始压缩包约2.8 GB，解压、R矩阵及输出需要额外磁盘和内存；建议至少16 GB可用内存、10 GB可用磁盘。完整CEL处理的本次耗时未测量；本次通过了12数组的R流程测试和全量表达矩阵到响应矩阵的复核。

R环境：4.2.3、affyio 1.68.0、preprocessCore 1.60.2。记录见 [R环境](provenance/hiphop/r-environment-recovered-20261003.txt)和 [CEL测试](provenance/hiphop/cel-smoke-summary-20261003.txt)。R脚本按原实现保留；Python移植版增加输出保护，并修复带显式行名列表头的兼容性，没有更改冻结输入的计算结果。

## 来源与条件信息

- 官方数据：[E-MTAB-2391](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-2391)，Lee et al., DOI [10.1126/science.1250217](https://doi.org/10.1126/science.1250217)。
- 原始记录：SDRF、IDF、HTTP响应头、归档哈希、处理记录在 `data/provenance/hiphop/`。HTTP响应头记录获取响应的时间；没有将文件修改时间当作下载时间。
- 3,850个条件包含494个vehicle、3,356个处理条件，涉及3,250个化合物、5,668株酵母。`data/response_conditions.tsv` 从冻结NPZ提取条件、化合物、剂量与单位，并保留源文件哈希。
- 原始库的剂量单位有micromolar、nanomolar、millimolar、picomolar及percent。9月历史清单的80条候选均为micromolar；当前提交表按每条候选的原始条件记录剂量单位。导出程序按化合物和剂量匹配单位，遇到缺失或歧义时停止。
- EMBL-EBI公开数据使用政策及原作者权利说明见[第三方数据条款](THIRD_PARTY_DATA.md)。原始项目曾把许可记为unknown；当前披露保留平台条款和作者署名，不将未声明CC许可的数据另行标为CC-BY。


## ESM-2蛋白表示重建

当前流程使用的三组ESM-2输入已从原环境恢复：6,742条酵母蛋白、810条人源配对记录和1,175个筛选靶点。压缩FASTA、来源快照、调用参数及权重哈希见 [esm2/](provenance/esm2/)。2026-10-03使用fair-esm 2.0.0、PyTorch 2.6.0及RTX A6000重新提取，三组索引和嵌入矩阵均与冻结输入逐字节一致，见[重建比较](provenance/esm2/rebuild-comparison.json)和[运行记录](provenance/esm2/rebuild-run.json)。

使用单独的Linux/CUDA环境安装 `requirements-embeddings.txt`。从[ESM官方模型入口](https://github.com/facebookresearch/esm#available-models)取得 `esm2_t33_650M_UR50D.pt` 及其同目录的 `esm2_t33_650M_UR50D-contact-regression.pt`；模型与辅助文件的哈希记在重建运行记录中。模型采用MIT许可，蛋白序列保留UniProt署名与CC BY 4.0条款。

```bash
python -m pip install -r requirements-embeddings.txt
mkdir -p outputs/esm2-inputs
for name in human_810 universe_1175 yeast_6742; do
  gzip -dc data/provenance/esm2/${name}.fasta.gz > outputs/esm2-inputs/${name}.fasta
  python scripts/extract_esm2_candidates.py \
    --fasta outputs/esm2-inputs/${name}.fasta \
    --gene-master raw/externalvalidation/mappings/gene_master.tsv \
    --model /path/to/esm2_t33_650M_UR50D.pt \
    --outdir outputs/rebuilt-esm2/${name} --device cuda
done
```

使用layer 33，GPU以fp16运行，按残基取float32均值，排除BOS/EOS。输出同时包含 `.npy`、索引和调用参数。正式筛选默认读取已冻结的表示；跨硬件或软件版本重建可能有浮点差异，因此哈希核验以已发布输入为准。
