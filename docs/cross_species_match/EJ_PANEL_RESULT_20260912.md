# EJ 方法产品线首跑结果（任务导出 + 冻结面板，2026-09-12）

方法冻结：`ESM2_JOINT_METHOD_FREEZE.md`（轴裁定：esm2；V9 G2 对
human-only p=0.0513 接受；固定思路跑流程）
导出：`esm2_joint_task_export.py`（config sha 与超参断言记录于
`esm2_joint_tasks_20260912/export_record.json`；1175 靶点 × 1282 株，
三成员 RRF k=20：direct_pc1 + aligned_C + b2_verbatim；train 510 对）
面板配置：`configs/product_execute.esm2joint.json`（除 task_dir/results_dir
外与冻结 gobridge 模板逐项一致：v2-spearman、sign=+1、n_perm=1000、seed=42）

## 冻结 allowlist 面板（确证性，与 B2 轴在册结果对照）

**poscon 4 对（results_esm2joint_poscon_20260912）**：q<0.1 共 0。

| 对 | EJ rho / emp_p | B2 rho / emp_p（在册） |
|---|---|---|
| DRD2×LNEPOX | +0.025 / **0.083** | +0.011 / 0.439 |
| KCNH2×amiodarone | +0.003 / 0.416 | +0.008 / 0.279 |
| KCNQ1×amiodarone | +0.005 / 0.356 | +0.006 / 0.314 |
| TRPV1×capsaicin | +0.014 / **0.257** | −0.021 / 0.994 |

**GPCR 8 对（results_esm2joint_gpcr_20260912）**：q<0.1 共 0。

| 对 | EJ rho / emp_p | B2 rho / emp_p（在册） |
|---|---|---|
| **HRH1×苯海拉明(PCHPOR)** | **−0.023 / 0.962** | **+0.031 / 0.012（q_bh=0.096，B2 唯一命中）** |
| HRH1×ZPEIMT | +0.009 / 0.246 | −0.003 / 0.605 |
| HRH1×XXPDBL | −0.012 / 0.788 | +0.014 / 0.163 |
| HTR2A×JJCFRY | +0.003 / 0.416 | +0.016 / 0.107 |
| HTR2A×ZUXABO | −0.003 / 0.811 | −0.002 / 0.826 |
| DRD2×ZPEIMT | +0.010 / 0.232 | −0.001 / 0.550 |
| DRD2×ZUXABO | −0.001 / 0.767 | −0.005 / 0.890 |
| DRD4×ZUXABO | −0.003 / 0.813 | +0.013 / 0.311 |

## 读法（如实登记）

1. **两轴面板层面互有胜负、总体全阴**：EJ 改善 4 对（TRPV1、DRD2×LNEPOX、
   DRD2×ZPEIMT、HRH1×ZPEIMT），恶化 4 对，其中关键一条：B2 轴唯一注册
   命中 HRH1×苯海拉明在 EJ 轴上符号翻转（+0.031→−0.023）。
2. 该命中的在册边界本就脆弱（8 对中 1 对、q=0.096、方向依赖 +z 耐药口径）；
   EJ 对任务排序的重排足以使其翻号，该差异为噪声级，不据此比较两轴优劣。
3. 面板全阴与任务轴选择弱相关、与**读出矩阵瓶颈**强相关的既有诊断保持
   成立：1282 株池缺效应基因（PMC1/VCX1/ERG 系/TOK1/TRK1 等，V6 审计在案；
   既定口径不扩池），且置换分辨率不支持 3.8M 对全筛的 BH α=0.1（预登记）。
4. **方法的注册证据在检索/基本结构层**（V8 门 A、V9 G1/G3 + G2 按裁定接受），面板层不构成对方法的有效检验；此边界随所有面板数字强制
   随附。

## 全屏（发现性，报告不设门）

`results_esm2joint_screen_20260912/`：1175 靶点 × 3250 化合物 ≈ 3.82M 对，
n_perm=1000、seed=42、BH/BY/Gumbel 尾拟合 q 全列输出；top 名单在
per_target_top/。运行日志 `esm2_joint_tasks_20260912/screen_run.log`。

## 全屏完成数字（2026-09-12）

3,818,750 对：q_bh<0.1 = **0**，q_tail_bh<0.1 = **0**（Gumbel 尾口径亦全阴）；
置换分辨率不支持 BH α=0.1（预登记）。发现性 top 名单见 per_target_top/ 与
exec_matrix.tsv 按 emp_p 排序头部。
