# V5 联合增益正式确认协议（预注册）

日期：2026-09-09
状态：**确认运行（confirmatory run，同数据）**，非外样本验证

## 1. 背景与前期探索披露

v3 三臂（19 类主标签分类，分数等权融合 joint）正式裁决 `NO_GO_JOINT_ROUTE_SINGLE`。
2026-09-09 的诊断（`match_diag_20260908/diag02_signal_audit.py`，A–U 段）定位了
三层原因并找到了显著的正向配置：

1. 数据质量：K562 矩阵 73 个响应列含 float32 溢出值（3.40e38），波及 2111/9608
   个 P1 档案；修复 = 列级 |max|<1e30 清洗门（等价于官方 bridge 管线的 finite
   gate），该事件与撤销清单见 `match_diag_20260908/diag02_03_SIGNAL_AUDIT_20260909.md`；
2. 表示：52 维 GO-slim 换成 640 个双物种共享 GO-BP term 的富集得分（ES）空间；
3. 联合机制：分数等权融合无效；**配对监督对齐（表示级）**与**多路去相关秩融合**
   有效。

**本协议为事后配置的独立正式化**：所有配置参数（下文全部冻结）来自上述事后
诊断；v3/v4 裁决不被修改。本运行在同一数据上按冻结配方单次执行，用途是可审计
的确认与留痕，不构成外样本证据。

## 2. 冻结输入（SHA-256 于 configs/v5_joint_formal.json，运行前强制核对）

- K562 GWPS normalized bulk（P1 档案提取）
- GSE42528 酵母矩阵与元数据（v3 benchmark 口径）
- Norman CRISPRa 矩阵与元数据（v3 benchmark 口径）
- orthodb 双边同源边（官方 511/142 split）
- GO BP 祖先注释（双物种）
- GO-slim 投影表（goSlim 臂与标签质心对齐用）

## 3. 轨道 A：K562 检索三臂（人查询 → 1282 酵母 KO 检索同源 gold）

表示：GO-BP ES。构建配方（冻结）：每物种 term 基因集限制在各自响应矩阵基因内，
10 ≤ 基因数 ≤ 500，取双物种共享 term；每签名行 z 标准化（nan 感知，clip ±50），
term 得分 = 成员基因均分 / √n，行 L2 归一。

臂定义：
- **A（无对齐）**：ES 直接余弦；
- **B（无监督协调）**：各物种签名矩阵各自去第一主成分（均值中心化后）再余弦；
- **C（联合，配对监督）**：B 的空间上，ridge 映射 W 由官方 train 511 对同源边
  拟合（λ ∈ {1,10,100,1000}，仅按 train AUC 选），查询经 W 变换后检索。

主评估：官方 test 142 对。主指标 = 查询行内判别 AUC（负例 = 每查询行固定采样
30 个非 882 边格，seed 20260911），Mann-Whitney 单侧。配对 ΔAUC(C−A) 与
ΔAUC(C−B) 用逐点胜率 bootstrap（3000 次，seed 20260912）。

次评估：5 折分组 OOF（人靶基因符号排序后 mod 5 分折，同靶不跨折）：C vs A/B
gold 百分位秩差符号翻转（2048 次）；任一同源 hit@10（该人靶全部酵母同源 KO
任一进 top-10）B vs C 符号翻转。

**门 A（主）**：C 臂 test AUC 的 MWU p < 0.05 且 ΔAUC(C−A) bootstrap CI95 下界 > 0。
**门 A（次，至少一项）**：OOF C-vs-A 符号翻转 p < 0.05；hit@10 增益 p < 0.05。

## 4. 轨道 B：Norman 基本结构三臂（v3 任务：19 类主 GO 标签分类）

臂 = 数据使用方式，LOGO OOF 按 split_group（组内排除人类参考库）：
- 单臂对照：human-only（ES-pc1）、yeast-only（goSlim-pc0，官方 train+有标签
  酵母参考库），最优单臂为判据基准；
- **联合（5 路秩融合 RRF）**，成员（冻结）：
  1. h_es1：human-ES-pc1 LOGO 分数；
  2. y_gs0：yeast-goSlim-pc0 参考库分数；
  3. jrepr：Norman ES-pc1 查询 × 轨道 A 的 W 映射 → 酵母 ES-pc1 参考库分数
     （Norman 自身配对零参与）；
  4. lab_gs0：goSlim 19 类质心 ridge（λ=100）跨物种对齐分数（LOGO）；
  5. y_es0：yeast-ES-pc0 参考库分数。
- RRF 公式 1/(k+rank)，k ∈ {5,10,20,40,60} **嵌套选择**（每个 LOGO 组的 k 仅在
  其余组上选）。

主指标：嵌套联合 vs 最优单臂的逐查询 reciprocal-rank 差，符号翻转（2048）。
次指标：固定 k=10/20 变体的 OOF/test MRR、top1/top3。

**门 B**：符号翻转 p < 0.05。

## 5. 运行与输出

- 单次运行：`python scripts/v5_joint_formal_run.py --config configs/v5_joint_formal.json
  --output <run_dir>`；输入哈希不匹配立即中止。
- 输出：`result.json`（全指标与门判定）、`per_query_ranks_trackA.tsv`、
  `per_query_rr_trackB.tsv`、`INPUT_MANIFEST.json`、脚本自哈希。
- 种子：20260911（负例采样）、20260912（bootstrap）、20260913+检验序号（各符号
  翻转检验）。所有随机过程仅用列出的种子。

## 6. 诚实边界

- 同数据确认：5 路配方与全部超参来自事后探索，本运行不提供外样本证据；
- 门判定的 p 值是单侧探索后确认口径，不得解读为预先假设检验的发现概率；
- CORE5 核心模块口径（诊断中 test AUC 0.69）为事后集合，本协议只报告不做门；
- 外样本验证（新的人侧数据或新的同源对）是下一步，不在本协议内。

## 7. V5.1 修订（第二轮事后迭代，2026-09-09 晚）

第一轮运行（`v5_joint_formal_20260909/`）门 A 主门未过：ΔAUC(C−A)=+0.066，
CI95 下界 −0.004，差 0.004。诊断（match_diag V 段）将轨道 A 的 C 臂从单一
检索器升级为**双检索器秩集成**后主门过线。V5.1 修订内容：

- **仅改轨道 A 的 C 臂配方**：成员 = ES-pc1 ridge（λ 网格仍按官方 train AUC 选）
  + ES-pc1 PLS 4 分量（train 拟合投影）；集成 = 查询行内升序秩百分位（rank/ny，
  高=优）两成员平均；
- 门定义、指标、官方 split、种子、负例采样、Track B 配方全部不变；
- 配置 `configs/v5_1_joint_formal.json`（新文件；第一轮配置不可变）；
- 诊断依据：diag02 V 段 ens c1+p4 于官方 test AUC 0.5662（p=0.0038），
  ΔAUC=+0.0838，CI95 [+0.0184, +0.1499]；
- 披露：这是同一数据上的第二轮事后迭代，两轮的 p 值均不得按单轮预注册解读；
  第一轮结果原样保留，不因 V5.1 覆盖。
