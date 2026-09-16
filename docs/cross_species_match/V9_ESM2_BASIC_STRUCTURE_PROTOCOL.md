# V9 ESM2 轴基本结构三臂正式协议（human-only / yeast-only / joint，2026-09-12冻结）

配置：`configs/v9_esm2_basic.json`（inputs 与 trackA 与 V8 冻结配置逐字节同源）
脚本：`scripts/v9_esm2_basic_run.py`（自哈希入 result.json）
口径澄清（2026-09-12）：**三臂两门 = v3 基本结构语义**：
human-only / yeast-only / joint；交付标准 = 一个方法，**joint 显著最优**，
且**跨物种匹配显著**。此前 V7/V8 按 V5.1 门结构（无对齐/无监督/联合 +
融合门）执行，属语义误读，本轮以 V9 更正（V7/V8 产物照常冻结保留）。

## 1. 单一方法的两个记录面

**G3 跨物种匹配显著（检索面）**：V8 冻结配方逐字重跑（同数据应逐位复现，
兼作再现性自检）：臂 A 裸余弦 / 臂 B 去 PC1 / 臂 C 配对监督 pc 空间三成员
秩集成。判定：门 A 主门（C p<0.05 且 ΔAUC(C−A) CI95>0）+ 次门 + B/C 各自
test AUC p<0.05。

**G1/G2 joint 显著最优（基本结构面）**：Norman 207 → 19 类 GO-BP，
esm2 轴，V8 同查询/参考/LOGO 口径：
- human-only = h_es1（LOGO 人侧带标签参考，npc1）；
- yeast-only = y_es1（酵母带标签参考，npc1）；
- joint = 嵌套 (family,k) 零泄漏选择的 RRF 融合，家族集合**首次包含中间
  家族**：pair_hy{h,y}、triple_hyj9{h,y,j9}、quad_hyjb{h,y,j9,y_b2}、all5；
- j9 = V8 过门 A 的完整 C 配方算子（npc* 空间 ridge+identity-anchored+PLS
  秩平均，strain 级分数按类取 max），joint 中的配对监督跨物种成员。

## 2. 门（运行前冻结）

- **G1**：joint（嵌套）vs 比较基线（5 单臂家族组内选优，零泄漏）：
  Δ MRR 符号翻转 p<0.05（seed 20261016）；
- **G2（"joint 显著最优"的严格形式）**：joint vs 固定 human-only 与固定
  yeast-only 双双配对符号翻转 p<0.05（seeds 20261017/20261018）、两 Δ 均值
  >0、且 joint 点估计高于两单臂；
- **G3**：第 1 节检索面判定全部成立。

## 3. 迭代与披露纪律

- 第五轮同数据事后迭代（V5→V5.1→V7→V8→V9），p 值不得按单轮预注册解读；
- 设计动机全部来自已登记结果：V8 单臂 MRR（h 0.365 / y 0.356 / jrepr8
  0.328 / y_b2 0.276 / lab_gs0 0.231）显示 all5 等权被弱成员稀释，而强弱
  搭配的中间家族未测；家族集合扩展由此动机驱动，选择本身零泄漏；
- **若 G2 仍不过：停止配方搜索**，登记结论"esm2 轴强表示下 joint 无超越
  单臂的增量"，"joint 显著最优"的在册证据保持为 GO 签名轴 V5.1 轨道 B
  （0.348 vs 0.286，p=0.0024），并如实呈报两轴分工；
- V7/V8 冻结产物只读；输出仅写 `v9_esm2_basic_20260912/`。
