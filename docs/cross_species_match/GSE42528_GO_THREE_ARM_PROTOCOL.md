# GSE42528 × 冻结 GO × Norman 正式三臂协议（v1）

状态：协议和适配器已冻结；数据完成后才允许运行正式模型。此文档不声称两物种
联合一定优于单物种。

## 这一层回答什么

人类端用 Norman CRISPRa 的基因扰动后表达变化，酵母端用 GSE42528 的基因删除后
表达变化。两边的响应基因都投影到同一个冻结 GO Biological Process / GO-slim
坐标，然后在同一批未见人类扰动上比较：

1. `human_only`：只用人类训练签名；
2. `yeast_only`：只用酵母训练签名；
3. `joint`：两组得分按事先冻结的 0.5 / 0.5 权重合并。

这是一个“联合是否真带来额外信息”的可证伪比较，不是靠故事默认联合有用。

## 标签和分组如何冻结

- 每个扰动基因的 `function_label` 只由
  `assets/go_bridge_v1/projection.shared.tsv` 决定；规则是
  `GO 边权重 × 两物种信息量` 最高的主标签。多标签 GO 边仍全部保留在投影表中；
  主标签只是为当前单类别评估而设。
- GSE42526 的 `non-responsive` 和 GSE42527 的 `responsive` 是根据表达变化数量事后
  分组的结果字段，只能作数据来源记录，绝不能当功能标签、分层条件或阈值
  调参依据。适配器要求原始 `function_label` 仍是
  `UNASSIGNED_NEEDS_FROZEN_GO`；如果学生预先填了 `responsive`，程序会直接拒绝。
- Norman 双基因扰动先按共用基因连成不可拆分的 component。如果一个酵母删除
  基因和 Norman 中的人类扰动是 OrthoDB 同源对，两者强制使用同一
  `split_group` 和同一 train/test 归属。
- 没有 Norman 对应同源基因的 GSE 删除株只进入酵母训练库，不冒充人类 test
  样本。类别空间只由锁定后两臂训练集的标签交集决定，不读 test 标签来删类。

## 防止“看答案”

两个物种都会把被扰动基因本身从该行响应特征中 mask 成 `NA`。否则，删除株
中的目标 mRNA 自身下降可能变成一个几乎直接的基因身份提示。源数据本来缺失的
单元格也保持 `NA`；数字0只能表示“确实观测到零效应”。

## 事先锁定的 GO / NO_GO 门槛

`configs/gse42528_go_three_arm.v1.json` 在正式结果之前已锁定：

- `joint` 必须同时严格超过 `human_only` 和 `yeast_only`；
- 对每个单臂的绝对 accuracy 增益都必须 **> 0.02**；
- 按独立 `split_group` 做配对 bootstrap，95% CI 下界也必须 **> 0.02**；
- 单侧配对置换检验 `p <= 0.05`，且 test 至少有 8 个独立 component。

任何一项不过都必须输出 `NO_GO_JOINT_ROUTE_SINGLE`。不允许在看到正式分数后改
0.02、0.5 权重、随机种子、类别空间或投影证据。

## 运行

```bash
python scripts/freeze_gse_go_three_arm_protocol.py \
  --protocol configs/gse42528_go_three_arm.v1.json \
  --adapter scripts/prepare_gse_go_three_arm.py \
  --benchmark-program scripts/cross_species_match.py \
  --norman /path/to/norman_task.npz \
  --modules /path/to/yeast_modules_fine.tsv \
  --orthologs /path/to/orthodb_yeast_human_s288c.tsv \
  --go-projection assets/go_bridge_v1/projection.shared.tsv \
  --gse-derived-manifest /versioned/GSE42528/processed/DERIVED_MANIFEST.json \
  --future-input-dir /new/versioned/gse_go_three_arm_inputs_v1 \
  --output /new/versioned/protocol_freeze/gse_go_three_arm_v1.freeze.json

python scripts/prepare_gse_go_three_arm.py \
  --protocol configs/gse42528_go_three_arm.v1.json \
  --freeze-manifest /new/versioned/protocol_freeze/gse_go_three_arm_v1.freeze.json \
  --norman /path/to/norman_task.npz \
  --modules /path/to/yeast_modules_fine.tsv \
  --orthologs /path/to/orthodb_yeast_human_s288c.tsv \
  --go-projection assets/go_bridge_v1/projection.shared.tsv \
  --gse-dir /versioned/GSE42528/processed \
  --output-dir /new/versioned/gse_go_three_arm_inputs_v1

python scripts/cross_species_match.py \
  --config /new/versioned/gse_go_three_arm_inputs_v1/benchmark_config.json
```

冻结器只读 config、代码、投影和 GSE 的派生 manifest，不打开两个表达矩阵；它会在
formal result 之前钉死时间、SHA-256、未来输入/结果目录和 2% 增益门槛。适配器随后
会先复核这个 freeze，再核对 Norman、OrthoDB、module 表、GO 投影和 GSE 派生文件的 SHA-256，
并拒绝覆盖旧输出目录。代码和模型输入可提交；GSE 大矩阵作为独立数据包用 manifest
对齐。

## 不能越界的结论

即使返回 `GO_JOINT`，也只说明“在这个锁定的跨物种基因扰动功能匹配任务上，
酵母表达签名提供了人类签名之外的信息”。它不等于药物有效、药物和靶点结合、
酵母可以代替人/鼠验证，也没有实测速度和成本。GSE 基因范围超出 v1 冻结 GO
映射的删除株会作为失败案例留在 `gse_inventory.tsv`，不会临时扩映射追求好结果。
