# 外样本验证（external validation）协议（EV，2026-09-16 冻结）

独立运行、单次登记。本协议先于任何 EV 运行写定；运行后不回改判据。

## 通用纪律（三模块共同）

1. **先冻结后运行**：每模块运行前，输入文件 SHA-256 与判据写入 run 记录；
   看过结果后改判据 = 整模块作废。
2. **算子 verbatim**：所有对齐器/超参/嵌入取自主链冻结产物，运行时断言一致
   （复用 esm2_joint_task_export.py 的 assert 模式）；EV 运行中禁止任何调参。
3. **运行登记**：判据先冻结后运行，产物登记于 results_a6000/externalvalidation/。
4. **阴性读取框架（预注册）**：EV 打平或失败时，默认解读为"主张未获外样本
   支持"，不得事后归因于平台差异；仅当运行前登记的已知技术差异（如 EV1 的
   测定平台差异）可作为注记，且不得改变结论措辞。
5. **产物等级**：全部 EV 结果为 candidate_only；进入对外材料需独立核对
   （数字与脚本 sha 复核）后才可引用。

## EV0 — 检索层外样本（第二采样确认）

- **输入**：`yeastbridge/data/mappings/{orthodb,oma,inparanoid}_yeast_human.tsv`
  （sha 冻结）；人/酵母 ESM-2 冻结资产（YeastBridge_0912/raw/tier1_esm2）；
  753 对冻结清单（用于排除）。
- **对集构造规则（冻结）**：三源并集 → 删除与 753 对完全相同的 (human,yeast)
  对（对级排除）；同一 human 基因多条新对全部保留；报告时另给"训练基因级
  排除"子集（human 基因出现在 753 训练对者剔除）作敏感性。
- **运行**：V8 `track_a` 冻结配方逐字重跑（臂 A 裸余弦 / 臂 B 去PC1 /
  臂 C 冻结超参三算子秩平均，`c_recipe` 用 V8 注册 picks，不重选）。
- **判据（冻结）**：主门 = 臂 C test AUC 的 MWU p<0.05 且 ΔAUC(C−A)
  bootstrap CI95 下界>0（与主链门 A 同构）；副指标 = OOF 秩、any-hit@10。
- **定位**：同源数据库生态内的第二采样，报告措辞为"second-sample
  confirmation"，不得称 fully independent。

## EV1 — 基本结构外样本（Adamson GSE90546）

- **输入**：GSE90546_RAW.tar（3/6 批次恢复件，sha 冻结）+ sidecar identities；
  Norman 207 元数据（标签冻结表）；V9 三臂机器与冻结 ESM-2/GO 资产。
- **预处理（ev_preprocess_adamson.py，规则冻结）**：good-coverage 细胞按
  guide 聚 pseudobulk；签名 = log1p(CP10K) 目标−对照 的跨批次均值，逐签名
  z 化；每目标每批次 ≥5 细胞；对照池 = 非基因符号 guide。
- **标签规则（冻结）**：Norman 207 单基因签名的 gene→function_label 查表
  直传；未命中基因保留于矩阵、标 label_unmapped，**不进门端点**。
- **运行**：V9 `basic_structure` 机器逐字复跑；查询 = 基因 ESM-2 嵌入均值
  （与 Norman 查询同式）；LOGO 组 = 单基因；嵌套 (family,k) 在新集组内照
  V9 程序选择（复现方法而非搬 pick）。
- **判据（冻结）**：主门 = joint MRR 同时 > human-only 与 yeast-only 且
  sign-flip p<0.05（与 G2 同构、同阈值）；G1 型对照（joint vs 组内最优单臂）
  一并登记。
- **预注册注记**：Adamson 为同细胞系（K562）不同实验/批次，测定平台差异
  属登记在案的技术差异，作注记不影响结论措辞。


## 目录

```
YeastBridge_0912（脚本在 scripts/，产物在 results_a6000/externalvalidation/，数据在 raw/externalvalidation/）
  EV_PROTOCOL.md          本协议（冻结）
  scripts/                EV 专用脚本（预处理/运行器，sha 随 run 记录）
  data/adamson/extracted/ 解包件（不删 tar 原件）
  runs/
    rehash/               A6000 全量重对报告（前置检查）
    ev0_retrieval/        EV0 运行输出
    ev1_adamson/          EV1 预处理与运行输出
  logs/                   后台任务日志
```

## 前置检查（已完成于 2026-09-16）

- - 输入完整性：MANIFEST.sha256 6275 项在 A6000 原位全量重对通过
  （见 results_a6000/externalvalidation/rehash/）。
- FD 恢复状态以 RECOVERY_LOG_v1.tsv 为准；静态站 2026-09-16 重试仍超时
  （HTTP 000/20s），已记录。

## EV1 修订记录（2026-09-16，预处理后、任何评分运行前）

预处理产出：Adamson 92 靶 × 32938 基因矩阵（runs/ev1_adamson/，含批次/细胞数报告）。
查询源设计：Adamson 靶集与 Norman 单基因标签表交集有限，据此将 EV1
分为两个子模块（判据不变）：

- **EV1a（主门）**：查询 = Replogle K562 GWPS 靶基因 ∩ Norman 标签表（79 基因）
  − 753 训练对的 human 基因（训练基因级排除，防对齐器循环）；标签 = 同一冻结
  查表直传；V9 三臂机器 verbatim。声明：与训练侧同测定平台（K562 GWPS），
  外样本性来自基因级 held-out + 标签从未用于任何轮次。
- **EV1b（次级，平台迁移）**：Adamson 92 靶；需先冻结标签规则 v2
  （GO-BP 类词集由 Norman 已标签基因的 ancestor 闭包派生，多数匹配指派），
  规则文档独立成册后方可运行；定位 = 探索性扩展。
- Adamson 预处理产物与 2 基因命中事实原样登记，不回改。
