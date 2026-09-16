# V7 B2 轴联合增益正式协议（三臂两门，2026-09-12 冻结）

配置：`configs/v7_b2_joint_formal.json`（输入 SHA-256 冻结，运行前强制核对）
脚本：`scripts/v7_b2_formal_run.py`（自哈希入 result.json）
本轮任务（2026-09-12）：把 B2 路线（scf-esm2-scyeast）的三臂两门完整跑一遍。

## 1. 问题

V5.1 的三门证据建立在 GO-BP ES 签名表示上（test AUC 0.5687 天花板，
repr_grid 四格表示层消融均未突破）。本运行把**同一协议骨架**（三臂语义、
两道门、冻结 split、种子族、指标）搬到 B2 蛋白轴表示上：

- 人侧：ESM2-650M（layer33 mean-pooled，810 蛋白新推理，冻结哈希）→
  scF-B2 注入投影 `pos_emb.proj`（`esm2_cross_match.py route_b_top` 逐字同式）；
- 酵母侧：route-B gene_table verbatim（主轴）；yeast ESM2 表（参照轴）；
  scyeast 先验空间（附加格）。
- 轴纪律（2026-09-11 轴裁定）：B2 为主轴，esm2_mean 为注册参照臂并列
  报告；不按结果选轴。

## 2. 轨道 A（检索三臂，门 A）

任务语义与 V5.1 轨道 A 相同（人查询→池内酵母菌株检索，gold=直系同源），
表示层换蛋白嵌入。主分析对集 = V5.1 冻结 753 对清单 ∩ B2 表示可用
（同数据可比性）；扩展对集（全部可行 edges，含 validation）仅描述性。

- 臂 A 无对齐：cos(X_b2, Y_rb)。注意与 V5.1 不同，蛋白轴的无对齐臂是
  强基线（同源性在序列空间直接可测），门 A 因此是更严格的主张：
  **配对监督对齐是否还在 B2 强基线之上加得动**；
- 臂 B 无监督：各空间去 PC1；
- 臂 C 联合：ridge（λ 网格，train-pair AUC 选参）+ PLS4 秩平均集成
  （V5.1 armC_recipe 同构）。
- 门 A 主门：C test AUC MWU p<0.05 且 ΔAUC(C−A) bootstrap CI95 下界>0；
  次门（至少一项）：OOF C-vs-A 符号翻转 p<0.05，或 any-hit@10 增益 p<0.05。
- 参照轴（esm2_mean）同构三臂全部并列报告，不设项目门。
- 附加格（仅对齐器 test AUC，repr_grid 式）：X_b2→scyeast、X_b2→yESM2、
  X_es→routeB、X_es→scyeast。

## 3. 轨道 B（Norman 19 类，门 B）

查询表示 = 扰动基因（1-2 个）ESM2 嵌入均值→l2n（B2 轴再过注入投影）；
酵母参考 = V5.1 yeast_ref_mask 同配方 ∩ 表示表覆盖；五成员
{h_b2, y_b2, jrepr_b2, lab_gs0, y_es}，其中 lab_gs0 为 V5.1 逐字复用的
表示无关去相关成员。家族 = 5 单臂 + all5；嵌套 (family,k) 组内选组外评，
比较基线 = 单臂组内选优。门 B：Δ MRR 符号翻转 p<0.05。

## 4. 机械适配披露（与 V5.1 的全部差异）

1. 查询从签名向量换为蛋白嵌入 → 轨道 A 不再依赖 K562 签名存在性；
2. row_negs 在 V7 池（route_b∩yESM2 覆盖的池菌株）内以同 seed(20260911)
   同配方（每查询 30，排除全部 orthodb 伙伴）重抽；
3. 轨道 B 成员表示换代，任务/标签/LOGO/k 网格/嵌套机制不变；
4. 人嵌入 build（UniProt 抓取+ESM2 推理）在看任何结果前完成并哈希冻结
   （810/811 ok；NEDD8-MDP1 无 reviewed 条目剔除，涉及的 edges 计入披露）。

## 5. 纪律

- 输出：`v7_b2_formal_20260912/results/`（INPUT_MANIFEST / result.json /
  per_query_ranks_trackA.tsv / per_query_rr_trackB.tsv）。
- 单次运行，门与指标不重配；结果无论方向原样登记；V5.1 产物只读。
- GPU 纪律：推理仅用 GPU 1（CUDA_VISIBLE_DEVICES=1）。
