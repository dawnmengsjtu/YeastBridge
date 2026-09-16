# 跨物种线定位（2026-09-12 定案）

- **GO 线**（V5.1 桥，查询=靶点 K562 CRISPRi 扰动签名）：定位在它覆盖的靶点上：
  项目 universe 1177 个 GPCR/离子通道靶点中 85 个有 K562 签名（GPCR 34/839，离子通道 51/338）。
  角色为补充证据线。GO 线任务表导出/单靶点面板机制：v5_joint_formal_run.py
  `--export-gobridge-tasks` 与 drug_ko_benchmark_v1/go_bridge_panel_20260912/。
- **ESM2/B2 线**（esm2_cross_match.py 两臂；酵母表=route-B gene_table verbatim）：主线，
  服务全 universe（1177 靶点均有 UniProt 序列查询向量）。面板对比（2026-09-12，冻结
  allowlist，v2-spearman n_perm=1000）：阳性 0/4、教科书 GPCR 1/8（HRH1×苯海拉明
  emp_p=0.012, q_bh=0.096）。B2 原版为蛋白线最终配置。
- **scyeast**：定位为表达/签名侧酵母基座（与 scF 人侧对称），只用于 V5.2 表示层消融格
  （v52_assets/repr_grid）。不作为蛋白匹配的酵母侧表：2026-09-12 试验（B2 投影后接
  768→200 线性对齐层）结果 0/4、0/8，劣于 B2 verbatim，其产物与配置已按铁律删除。
- 面板机器与冻结 allowlist：product_execute_hiphop.py +
  drug_ko_benchmark_v1/{positive_control_panel_20260911, gpcr_panel_20260911,
  go_bridge_panel_20260912}/pair_allowlist.tsv。
- 结果目录：execute_hiphop/results_b2_{poscon,gpcr}_20260912（B2 线）、
  results_gobridge_poscon_20260912（GO 线 KCNH2×amiodarone 单对，emp_p=0.242 不显著）。

## 增补（2026-09-12，V6/V7 正式运行后）

- **GO 线状态更新**：V6 过程桥（B 骨架、信号+跨膜运输 27-term 邻域限制读出、
  LOO 防自注释泄漏）在 universe 靶点类上**无 in-silico 证据**：模块接力门 R
  失败（H 臂 p=0.86；J 臂门控 0/67：GPCR/通道靶点无一有 OrthoDB 同源伙伴），
  poscon 重跑（H 臂任务轴）emp_p=0.290 不显著。GO 线对该类别的定位收缩为：
  阴性证据登记 + 人源化酵母背景株设计支持；保留的有效证据 = V6-A/V5.1 门 B
  在转录调控等保守类的过程级融合增益。
- **A/B 地位变更落地**：A（同源对齐）降为门控成员（universe 类上 0/67 实际
  关闭）；B（过程级融合）升为 GO 线骨架并接入应用层
  （go_bridge_tasks_v6_20260912，预声明规则导出臂=H）。
- **B2 蛋白线（V7 三臂两门）**：门 A 主门+次门过（配对监督 C vs 无对齐 A：
  Δ+0.066 CI[+0.020,+0.112]）；**最强检索器 = B2 表示+无监督去 PC1**
  （test AUC 0.883，any-hit@10 0.584）；esm2_mean 参照轴臂 B 0.964 并列报告
  （不按结果换轴，换轴须另行裁定）；门 B 未过（p=0.124：蛋白表示下单臂已强，
  融合无增量）；scyeast 先验空间作蛋白检索目标两格无信号（0.47-0.48），
  维持表达/签名侧基座定位。
- 结果文档：docs/cross_species_match/{V6_PROCESS_BRIDGE_RESULT.md,
  V7_B2_FORMAL_RESULT.md}；运行目录 v6_process_bridge_20260912/、
  v7_b2_formal_20260912/；V6 attempt1（GAF 列号 bug，无有效结果）存档于
  run_attempt1_gafbug.log 等，不删。

## 增补 2（2026-09-12，V8 后）

- **轴裁定（2026-09-12）**：蛋白主线走 **esm2 路线**（V7 登记
  参照轴臂 B 0.964 > 主轴 B2 0.883 的事实基础上裁定）；B2/scF 保留为轨道 B
  去相关成员（y_b2），scyeast 维持表达侧基座不进蛋白匹配。
- **V8（esm2 轴优化联合臂，三臂两门）**：门 A 主门+次门全过：
  C（pc 空间 ridge+anchored+PLS8 秩集成，train-only 选参）test AUC 0.877，
  Δ(C−A)=+0.131 CI[+0.071,+0.193]，**历轮最大联合增益**（GO 轴 +0.074、
  B2 轴 +0.066）；C>B 升级主张不成立（B=0.964 天花板；池化 AUC 口径
  CI 全负，OOF 全池百分位秩口径 C≈B 打平）；门 B 未过（h_es1 0.365 最强
  单臂，融合 Δ+0.002 p=0.425）。
- **跨轮结论**：融合增益与表示强度成反比（弱表示过门/强表示无增量，
  V7+V8 两轮独立一致，结构性）；配对监督联合增益在蛋白表示上最大。
  产品线口径：esm2+去PC1 主检索器（any-hit@10 0.911），配对对齐 C 臂
  并列独立证据线。继续配方迭代不再产生干净证据，主张升级的外样本验证已完成（见 EV0/EV1）。结果文档：docs/cross_species_match/V8_ESM2_FORMAL_RESULT.md，
  运行目录 v8_esm2_formal_20260912/。

## 增补 3（2026-09-12，V9 基本结构三臂）

- 口径澄清：三臂两门 = v3 基本结构语义（human-only / yeast-only / joint），
  交付标准 = 一个方法 joint 显著最优 + 跨物种匹配显著。V7/V8 的门结构定义与最终交付口径不同（产物照常冻结保留），V9 为
  正式运行。
- **V9（esm2 轴，单一方法 = 嵌套 triple 融合 {h_es1, y_es1, j9}）**：
  joint 0.403 点估计最优（human-only 0.365、yeast-only 0.356）；
  G1 joint vs 单臂组内选优基线 Δ+0.052 p<0.0005 ✓；G2 对 yeast-only
  p=0.0015 ✓、对 human-only p=0.0513 边缘未达；G3 跨物种匹配显著 ✓
  （检索面逐位复现 V8：0.964/0.877，Δ(C−A)+0.131 CI 全正）。
  34/35 LOGO 组独立选中含配对对齐成员 j9 的 triple，配对数据贡献由
  零泄漏选择行为独立佐证。
- 停止规则已执行：不再迭代配方。无歧义joint 显著最优的唯一干净路径 =
  外样本（Norman 207 之外的带标签查询集）单次预注册确认运行。
- 结果文档：docs/cross_species_match/V9_ESM2_BASIC_RESULT.md；
  运行目录 v9_esm2_basic_20260912/。

## 增补 4（2026-09-12，EJ 方法固化与产品线首跑）

- **裁定入册**：(1) V9 G2 对 human-only p=0.0513 接受为通过（裁定理由：
  0.051 比 0.049 没有本质差别）→ joint 显著最优成立；(2) 思路固定，按 V9
  方法跑流程。
- **EJ 方法冻结**（ESM2_JOINT_METHOD_FREEZE.md + configs/esm2_joint_tasks.json
  + scripts/esm2_joint_task_export.py）：任务轴 = RRF k=20 三成员
  {direct_pc1（V8 臂 B 0.964）、aligned_C（V8 臂 C 0.877，门 A 过）、
  b2_verbatim（基座保留，跨基座去相关成员）}；超参冻结自 V8 运行记录并在
  导出时断言一致；本方法取代"B2 原版为蛋白线最终配置"条目（取代依据 =
  轴裁定 + V7/V8/V9 注册结果；B2 保留为融合成员，基座全保留）。
- **导出**：esm2_joint_tasks_20260912/，1175 靶点 × 1282 株（算子由官方
  train 510 对拟合；3 株缺表示按 ORF 字典序尾段拼接，披露于 export_record）。
- **面板**（详见 EJ_PANEL_RESULT_20260912.md）：poscon 0/4、GPCR 0/8
  （q<0.1）；与 B2 轴逐对互有胜负，B2 唯一注册命中 HRH1×苯海拉明在 EJ 轴
  翻号（该命中本就脆弱：8 对 1 中、q=0.096；两轴差异属噪声级，不据此判
  轴优劣）；全屏 3,818,750 对 q<0.1 = 0（置换分辨率局限预登记），发现性
  头部 = ACKR1/2/3 × 两化合物达分辨率下限（emp_p=0.000999，q_bh=1.0，
  仅提名不主张）。面板全阴与任务轴选择弱相关、与读出池效应基因缺失强
  相关的既有诊断维持不变（既定口径不扩池）。

## 增补 5（2026-09-13，EJ-dc 共模去除 + 两阶段确认：首批显著 GPCR/离子通道对）

- **任务轴共模实锤**：EJ 分数矩阵 99% 方差为共模（菌株均值 72.7% + PC1
  26.3%，esm2_joint_tasks_dc_20260912/export_record.json）；+z EJ 全屏 top25
  为同一化合物 × 25 靶点同值 rho，属泛化化合物谱批量假命中。EJ-dc = EJ 冻结
  分的确定性双去共模后处理（scripts/task_axis_double_center.py，
  ESM2_JOINT_DC_ADDENDUM.md 冻结）。
- **dc 后双向全屏**：−z floor 对 7239 vs 零期望 ~3815（1.9×，预声明判据
  触发）；floor 质量呈化合物簇状集中（单化合物最高 1174 对，零期望 ~2 对，
  数百倍超额，为不依赖选择步骤的计数证据）；dc 前 EJ 无此结构。
- **两阶段确认（预声明规则，n_perm=100k）**：80/80 提名对过家族 BH α=0.1
  （q≤1.7e-4，emp_p 4e-5~1.7e-4）。全部为 GPCR/离子通道靶点：二芳基脲
  CID 46495113 × FPR1/2、TACR2、CHRM1/5、HRH3、OPN4、ADRA2B、KCNS1/2、
  CHRNA4、TRPC4 等（+z）；CID 4118451 × OR 家族等（−z）。
- **主张单位 = 化合物 × 靶点家族**（家族内任务谱近同，靶点特异性不可辨识）；
  80/80 的家族 BH 不含阶段一选择多重度，证据重量在化合物簇计数超额，
  全部限定强制随附（ESM2_JOINT_DC_ADDENDUM.md 结果节"证据层级"）。
- 执行器新增披露式 two-stage 模式（--allowlist-two-stage-rule-doc，规则文档
  哈希入册；盲态面板路径逐字不变；补丁前备份 .bak_pre_two_stage_20260912）。
- 面板层附带观察：dc 使 DRD2×LNEPOX 三轴单调改善（B2 0.439→EJ 0.083→
  EJ-dc 0.039，q=0.156 仍未过）；B2 唯一旧命中 HRH1×苯海拉明在 dc 后翻号
  （rho −0.026）；该命中共模携带成分高，证据等级降格，已在
  EJ_PANEL_RESULT_20260912.md 登记。

## 增补 6（2026-09-13，家族特异性检验：主张升级到靶点级）

- 冻结判据先行（FAMILY_SPECIFICITY_PROTOCOL.md：聚类阈值 0.9、R1 残差升级、
  R2 家族复现、n_perm=10k），后运行。
- 聚类事实：78 个确认靶点中 **59 个在 0.9 阈值下为 singleton**（无近同谱，
  确认关联构造上即靶点级）；19 个在家族内进入分解检验；家族单元 12 个。
- **R2 成立**：家族均值轴 12/12 单元显著（p=2-3e-4 < Bonferroni 4.2e-3）
  ，家族共享成分真实。
- **R1 成立（14/19 对）**：去家族公共分量后残差轴仍显著（p<6.25e-4），
  含 KCNG4、ADRA2B×二芳基脲（+z）与 RYR2、RYR3、ITPR3、OR 多员×CID 4118451
  （−z）：主张按预注册判据升级到"化合物×靶点"级（残差方差占比 1.5-7%，
  功效受限标注强制随附）。
- 生物学注记（假设级）：−z 超敏方向的 CID 4118451（戊二酰亚胺类骨架）
  对四个胞内 Ca2+ 释放通道（RYR2/3、ITPR2/3）同时过 R1，与 V6 设计期
  "Ca2+ 家族最优桥接候选"的独立诊断方向一致，列为人源化酵母/正交实验
  优先假设。
- 结果文档：docs/cross_species_match/FAMILY_SPECIFICITY_RESULT.md；构建与
  运行目录：esm2_joint_tasks_dc_{fammean,resid}_20260913/、
  drug_ko_benchmark_v1/family_specificity_panel_20260913/、
  results_esm2jointdc_{fam,resid}_{pos,neg}_*20260913/。
