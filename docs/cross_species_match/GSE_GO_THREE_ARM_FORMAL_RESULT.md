# GSE42528 × GO × Norman 正式结果数据卡

结论先说：**数据和跨物种评估链已经打通，但这一版模型明确 `NO_GO`。**

这不等于“酵母不能筛药”。它说明的是：把人类 CRISPRa 和酵母 deletion-expression
投到 GO-slim，再用“主 GO 标签 + 最近质心”做跨物种功能分类，当前这个简单模型
一个 test 都没有精确分对，更不能证明“人 + 酵母”比任何单臂更好。

## 数据层已经补齐

- 人类端：Norman CRISPRa 223 条扰动签名，207 条在冻结 GO 空间有主标签；
- 酵母端：GSE42528 的 1,484 个删除株 × 6,170 个 systematic response gene，1,282 条
  在 v1 冻结 GO 范围可进入 benchmark；
- 数据语义：GSE42526/27 的 responsive/non-responsive 只是 provenance，没有当标签；
- 防泄漏：人类和酵母的 perturbation self-feature 分别 mask 了 346 和 1,279 个单元格；
  缺失值用 `NA`，从未填0；
- 分组：11 个人类 test 样本对应 11 个独立 component，跨 train/test 基因泄漏为0；
- 投影：16,583 条 GO 边、52 个双物种模块，所有签名的最低绝对效应覆盖率为
  0.6264，高于事先锁定的0.5。

因此，这次失败不是“没数据”、“酵母矩阵没读进来”或“泄漏被拦了”；失败点在
当前的标签/模型对这个任务没有足够辨别力。

## 三臂结果

主结果保留全部 11 条 test。训练空间是人类-train 与酵母-train 共同拥有的 19 类；
2 条不在训练空间的 test 没有被删除，三臂均按 accuracy/top-k/MRR=0 计入。

| 臂 | Accuracy | MRR | Top-3 accuracy |
|---|---:|---:|---:|
| human-only | 0.0000 | 0.0906 | 0.0909 |
| yeast-only | 0.0000 | 0.1688 | 0.1818 |
| joint | 0.0000 | 0.1451 | 0.1818 |

9 条 known-class 子集只是次要描述，不能影响 GO/NO_GO。它的三臂 accuracy 同样都是0；
MRR 分别为 human 0.1108、yeast 0.2063、joint 0.1774。

以 accuracy 为事先锁定的主指标：

- joint 相对 human-only：增益0，95% bootstrap CI `[0, 0]`，exact sign-flip `p=1`；
- joint 相对 yeast-only：增益0，95% bootstrap CI `[0, 0]`，exact sign-flip `p=1`；
- 两个比较都没越过事先冻结的2个百分点增益门槛。

最终状态是 `NO_GO_JOINT_ROUTE_SINGLE`。程序输出的 `yeast_only` 只是“三臂 accuracy
都为0时的预设 tie-break，且它的描述性 MRR 最高”；**不能把这句话包装成“酵母模型
已可用”**。

## 三次尝试的完整审计

1. v1：`BLOCKED_UNSEEN_TEST_CLASSES`。发现 2/11 test 类在某一训练臂不存在，没删 test，
   没计算模型分数。
2. v2：`ABORTED_PRE_MODEL_UNSUPPORTED_PROJECTION_EVIDENCE`。open-set 保留全 test 后，loader 在
   读投影文件时发现文件仍包含配置未允许的 ortholog 审计行，模型前立即停止，结果
   目录无文件。
3. v3：仅把 benchmark 投影文件实体过滤为配置已声明的 GO rows；数据、GO、seed、
   split、19类训练空间、11条 test、0.5 joint权重和2%门槛全部不变。预模型检查通过后
   只正式运行一次，得到上述 NO_GO。

可随代码包分发的机器可读摘要在 `configs/gse_go_three_arm.formal_status.json`；
完整轻量运行文件在本机 `result/gse_go_three_arm_v3_formal/`，由于 `result/` 默认排除于代码包，
也可按下述 SHA-256 从服务器独立取回。

## 精确输出与哈希

服务器根目录：

`/public/home/mengxl/dzy/yeastbridge_re_mvp_stage_20260901/gse_go_three_arm_v3_open_set`

主要 SHA-256：

- protocol freeze：`964bb4f24ee6edc1c264d6dbd7c0c595e96899af4439203fc11943fa6bd0601e`
- metrics：`14ecb967c8f6d605e7b6045eae5e94a55be6629f417d6bef72d45e12c7092583`
- predictions：`be0664536e049526e603a5a2cac3219e5f54b7d9516b5e61109c0fab1cf50801`
- protocol lock：`72bd9b4f540313f0701d8ad1e59af35c31b3659b1a11ded188d19884d87316e9`
- split lock：`7f92ad2adba779380eca167526d66a7e5e9abdf2398047c074a00d7a1e0c2fc7`
- integration report：`259d303f8fac03415c0e5790bb970c043a156ae600cb1ef5d0123c7044805f91`
- GSE v2 derived manifest：`cfad65ecc28d962650fa442edc3ab22019254db2778024d53301f9c35d35ae54`
- GSE matrix：`7dea42f99b5811243b0350e6324a95cd995a9407cab0586f2546ac2cc004493b`
- GO源投影：`f55eeeef5261bf10d7cf6dc40c1272ef2b0dc244aa77cab800750fb2f74925db`

轻量结果已复制到 `result/gse_go_three_arm_v3_formal/`，不包含 60 MB 的 GSE 大矩阵。

## 复现命令

按顺序运行：

```bash
python scripts/freeze_gse_go_three_arm_protocol.py \
  --protocol configs/gse42528_go_three_arm.v3_open_set.json \
  --adapter scripts/prepare_gse_go_three_arm.py \
  --benchmark-program scripts/cross_species_match.py \
  --go-projection /frozen/projection.shared.tsv \
  --norman /frozen/norman_task.npz \
  --modules /frozen/yeast_modules_fine.tsv \
  --orthologs /frozen/orthodb_yeast_human_s288c.tsv \
  --gse-derived-manifest /frozen/GSE42528/DERIVED_MANIFEST.json \
  --future-input-dir /new/inputs_v3_open_set \
  --output /new/gse_go_three_arm_v3_open_set.freeze.json

python scripts/prepare_gse_go_three_arm.py \
  --protocol configs/gse42528_go_three_arm.v3_open_set.json \
  --freeze-manifest /new/gse_go_three_arm_v3_open_set.freeze.json \
  --norman /frozen/norman_task.npz \
  --modules /frozen/yeast_modules_fine.tsv \
  --orthologs /frozen/orthodb_yeast_human_s288c.tsv \
  --go-projection /frozen/projection.shared.tsv \
  --gse-dir /frozen/GSE42528/processed \
  --output-dir /new/inputs_v3_open_set

python scripts/cross_species_match.py \
  --config /new/inputs_v3_open_set/benchmark_config.json
```

## 可以怎么对外说

可以说：“我们已把高等生物基因扰动签名、酵母删除表达签名、冻结功能投影、无泄漏
分组和三臂审计全部打通，并且诚实得到了一个 NO_GO；这个反证告诉我们，简单 GO
最近质心不足以支撑跨物种联合主张。”

不能说：“联合模型已优于单模型”、“酵母已预测人体药效”或“当前 yeast-only 模型已可以
直接筛新药”。HIP/HOP 的已知药物正对照筛选是另一个实验问题，不能用它替这次跨物种
三臂 NO_GO 翻案。
