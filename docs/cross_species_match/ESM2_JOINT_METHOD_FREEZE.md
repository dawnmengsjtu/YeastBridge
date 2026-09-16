# ESM2-Joint 产品线方法冻结记录（2026-09-12）

## 裁定记录（本冻结的授权依据）

1. 轴裁定（2026-09-12）：蛋白主线走 **esm2 路线**；B2/scF 保留为成员
   （"尽可能保留基座模型"）；scyeast 维持表达侧基座不进蛋白匹配。
2. 结果接受（2026-09-12）：V9 基本结构三臂 G2 对 human-only 的
   p=0.0513 **接受为通过**（裁定理由：0.051 与 0.049 不构成实质差别）；
   据此 **joint 显著最优成立**（vs 最优单臂基线 p<0.0005；vs yeast-only
   p=0.0015；vs human-only p=0.0513 裁定接受）。
3. 方法固定："固定该思路，按照这里的方法跑流程"。

## 冻结方法（EJ 任务轴）

对每个 universe 靶点（1175，全部有 ESM2 查询向量），任务分 =
**RRF k=20 三成员融合**：

| 成员 | 定义 | 注册证据 |
|---|---|---|
| m1 direct_pc1 | cos(靶点 ESM2 去PC1, 池菌株 ESM2 去PC1) | V8 臂 B：检索 AUC 0.9638，any-hit@10 0.911 |
| m2 aligned_C | npc*=3 空间三算子秩平均（ridge λ=1 + identity-anchored λ=1 + PLS k=8，官方 train 510 对训练） | V8 臂 C：0.8770，门 A 主+次过，Δ(C−A)+0.131 CI[+0.071,+0.193] |
| m3 b2_verbatim | cos(ESM2→scF-B2 注入投影, route-B 酵母表)（esm2_cross_match route_b_top 逐字同式） | B2 原版蛋白线在册配置；面板唯一历史命中（HRH1×苯海拉明 q_bh=0.096）出自该成员 |

- 超参冻结自 V8 运行记录（`v8_esm2_formal_20260912/results/result.json`
  armC_recipe_info.picks：npc=3, ridge λ=1, anchored λ=1, PLS k=8），导出
  脚本运行时断言一致；RRF k=20 冻结自 V9 选择结果（34/35 LOGO 组选
  triple|k=20）。
- 结构同构于 V9 获胜 triple {h_es1, y_es1, j9}：强单物种视角 + 跨物种
  直连 + 配对监督对齐的去相关秩聚合；universe 任务无标签，故用冻结配方
  而非逐靶嵌套选择。
- PC 基拟合矩阵：fx = 810 蛋白冻结资产矩阵（v7 assets），fy = 池酵母
  ESM2 矩阵（1279 株）；transductive 约定与 V5.1/V7/V8 注册口径一致。
- 池 1282 株中 3 株缺 esm2/routeB 行：尾段按 ORF 字典序单调拼接披露。

## 取代关系与地位

- 本方法**取代** "B2 原版为蛋白线最终配置"（LINE_POSITIONING 20260912
  原条目）作为产品线任务轴，取代依据为轴裁定与 V7/V8/V9 注册结果，
  非结果挑选：m3 成员保留 B2 verbatim，基座全保留。
- GO 线（V6 H 臂导出）对该靶点类阴性已登记（V6_PROCESS_BRIDGE_RESULT），
  产品线不再使用 GO 线任务轴；GO 线证据保留为保守类过程级融合（V5.1
  轨道 B / V6-A）。
- 对外主张口径（V9_ESM2_BASIC_RESULT 判定总表 + 本文件裁定记录第 2 条）：
  "esm2 轴 joint（三成员融合）显著最优：vs 最优单臂 p<0.0005、vs
  yeast-only p=0.0015、vs human-only p=0.051（按项目裁定接受）；跨物种
  匹配显著（0.964/0.877，Δ(C−A) CI 全正）"。

## 流程（本冻结的执行范围）

1. `scripts/esm2_joint_task_export.py --config configs/esm2_joint_tasks.json`
   → `esm2_joint_tasks_20260912/`（1175 任务表 + export_record.json）；
2. 面板重跑（执行器配置 `configs/product_execute.esm2joint.json`，除
   task_dir/results_dir 外与冻结 gobridge 模板逐项一致：v2-spearman、
   sign=+1、n_perm=1000、seed=42、BH α=0.1）：
   - poscon：冻结 allowlist sha `59a19823656ef4d0550811a29a4a08c010c78a8c58cef6c97f39e7d58af08b93`（4 对）
   - gpcr：冻结 allowlist sha `d0ce81401be879c66d039a6c2deefb285340d3180e50913896c7423c885312ed`（8 对）
   - screen：全 universe × 3250 化合物（发现性，报告不设门；置换分辨率
     局限预登记）
   结果目录 `results_esm2joint_{poscon,gpcr,screen}_20260912`。
3. 对照基线（已在册，不重跑）：results_b2_{poscon,gpcr,screen}_20260912
   （B2 verbatim 单成员任务轴：poscon 0/4、gpcr 1/8）。

## 边界

- 面板仍受 1282 池效应基因缺失制约（既定口径不扩池）；面板阴性不能
  反证方法无效，面板阳性则为方法+池联合证据；此局限在所有面板结果
  文档中原样随附。
- 本冻结不改动任何 V5.1/V6/V7/V8/V9 产物；方法参数若需变更 = 新冻结
  文件 + 新目录。
