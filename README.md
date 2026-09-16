# yeastbridge_clean — YeastBridge-RE 最终积极版本复现（SCNet，2026-09-13 换装）

本目录于 2026-09-13 整体换装：内容 = **A6000 当前跑通的全链**
（EJ 方法导出 → 共模去除 → 面板 → 两阶段确认 → 家族特异性检验）+
**esm2 路线最新三臂两门（V9，最终得出 joint 显著最优的那一版）**。
V5.1 时代的旧复现项目原样归档于
`/root/private_data/yeastbridge_clean_v51era_20260909/`V5.1 时代旧复现项目原样归档（本目录以符号链接占位引用）。

## 结论版本（本目录复现的对象）

| 层 | 结论 | A6000 注册记录 |
|---|---|---|
| 跨物种匹配（检索） | esm2 轴臂 B（去PC1）test AUC **0.9638**；联合臂 C 0.8770，ΔAUC(C−A)=+0.131 CI[+0.071,+0.193]，门 A 主+次过 | v8/v9 trackA（results_a6000/） |
| 三臂基本结构 | **joint 0.403 最优**（human-only 0.365 / yeast-only 0.356）；G1 p<0.0005；G2 vs yeast-only p=0.0015、vs human-only p=0.0513（未达 0.05 阈值，经裁定按边际接受） | v9（results_a6000/v9_esm2_basic_20260912/） |
| 产品线任务轴 | EJ = RRF20{direct_pc1, aligned_C, b2_verbatim}；EJ-dc = 去共模（共模占 99% 方差） | tasks/、docs/cross_species_match/ESM2_JOINT_* |
| 面板显著对 | 两阶段确认 80/80 过家族 BH（q≤1.7e-4）；家族特异性检验：R1 升级 14 对（含 RYR2/RYR3/ITPR3/KCNG4/ADRA2B），R2 家族级 12/12；59/78 靶点构造性靶点级 | results_a6000/execute_hiphop/、docs/…/FAMILY_SPECIFICITY_RESULT.md |
| GO 线（记录） | 对 GPCR/通道类三重阴性（V6）；有效域=保守类过程融合（V5.1/V6-A） | results_a6000/v6_process_bridge_20260912/ |

## 目录

```
env -> 归档 venv (numpy1.26.4/pandas2.3.3/scipy1.15.3/anndata0.11.4, v51 时代钉平)
env2/            本项目新建 venv: 同钉平版本 + 系统 torch 2.9.0 (--system-site-packages)
raw/tier0_* tier1_{benchmark,k562,assets} -> 归档冻结数据（哈希当年已双端核验）
raw/tier1_esm2/      human_810（本项目新推理）/ yeast_650m / universe_1175
raw/tier1_models/    route_b（gene_table+final_model.pt+A2_init）/ scyeast 先验空间
raw/tier1_response/  strain_response.npz（Lee2014 HIP/HOP）+ compounds.tsv.gz
raw/tier1_projections/ gene_goslim_projection + gene_bp_ancestors
scripts/         A6000 工作版脚本 verbatim（v5..v9 + EJ 链 + 执行器含 two-stage 补丁;
                 补丁前备份 product_execute_hiphop.py.bak_pre_two_stage_20260912）
configs/         SCNet 路径适配版（算法参数与冻结原版逐项一致）
configs_a6000_frozen/  A6000 冻结原件（stage_root 为 A6000 路径，仅作记录）
tasks/           EJ / EJ-dc / fammean / resid 任务表（A6000 产物原件）+ GO 线记录件
panels/          全部冻结 allowlist + 提名/构建记录（sha 在册）
docs/cross_species_match/  全部协议+结果文档（V5.1→V9、EJ 冻结、dc 附录、家族特异性）
results_a6000/   A6000 正式运行记录件（v6-v9 结果、面板/确认/分解运行、screen 仅
                 summary+top 名单，exec_matrix ~700MB×3 未传输，可重跑再生）
results_v51era -> 归档复现结果（step1-4 逐位一致记录）
repro/           本机复现运行输出
MANIFEST.sha256  传输清单（6275 文件已全部双端校验通过, 0 问题）
```

## 复现命令（A6000；python=/public/home/mengxl/dzy/envs/yeastbridge/bin/python）

```bash
cd /public/home/mengxl/dzy/YeastBridge_0912
# 验收主运行：V9 esm2 三臂两门（joint 版）
/public/home/mengxl/dzy/envs/yeastbridge/bin/python scripts/v9_esm2_basic_run.py \
    --config configs_a6000_frozen/v9_esm2_basic.json --output repro/v9_joint_20260913
# 产品线全链（依次）：
/public/home/mengxl/dzy/envs/yeastbridge/bin/python scripts/esm2_joint_task_export.py \
    --config configs_a6000_frozen/esm2_joint_tasks.json --output repro/tasks_ej
/public/home/mengxl/dzy/envs/yeastbridge/bin/python scripts/task_axis_double_center.py \
    --input repro/tasks_ej --output repro/tasks_ej_dc
/public/home/mengxl/dzy/envs/yeastbridge/bin/python scripts/product_execute_hiphop.py \
    --config configs_a6000_frozen/product_execute.esm2joint_dc.json \
    --results-suffix _poscon_repro --n-perm 1000 \
    --pair-allowlist panels/positive_control_panel_20260911/pair_allowlist.tsv \
    --pair-allowlist-sha256 59a19823656ef4d0550811a29a4a08c010c78a8c58cef6c97f39e7d58af08b93 \
    --allowlist-declared-hiphop-blind --allowlist-selection-basis <冻结basis原文>
# 两阶段确认与家族特异性：见 docs/cross_species_match/ESM2_JOINT_DC_ADDENDUM.md
# 与 FAMILY_SPECIFICITY_PROTOCOL.md（two-stage 模式用 --allowlist-two-stage-rule-doc）
```

## 验收状态

| # | 项目 | 状态 |
|---|---|---|
| 0 | 传输完整性：MANIFEST 6275 文件 sha256 双端一致 | 通过（0 问题）|
| 1 | 适配配置路径自检（14 个 config 全路径存在） | 通过 |
| 2 | V9 三臂两门复现（joint 0.403 / G1 / G3 数字对照 A6000） | 通过，逐位复现：joint 0.403/comp 0.350/G1 p=0/G3 门全过；per-query TSV md5 与 A6000 一致 |
| 3 | EJ 链（导出→dc→面板→确认→家族） | 结果在册（results_a6000/execute_hiphop，80/80） |

## 注意

- 本链全 CPU（仅 torch.load 用到 torch）；
- 外网：curl 加 `-4`；pypi 走清华镜像；大文件从 A6000 侧 rsync（密钥
  creed-a6000-transfer-20260908 已互通，`ssh -p 10077` 免密）；
- /root 镜像层限额：一切大件放 private_data（本目录已在持久卷）。

## 外部验证（External Validation，2026-09-16）

在主链冻结结果之外，用两条从未参与任何训练/选择轮次的独立数据做外样本确认：

- **EV0 检索外样本**：三库（OrthoDB/OMA/InParanoid）并集中 240 对新同源、213 个新人源基因
  （为 215 个蛋白补推理 ESM-2 嵌入，见 raw/externalvalidation/ev0x_embed_merged/）。
  V8 冻结算子原样运行：三臂 AUC 0.76 / 0.98 / 0.95，
  ΔAUC(C−A) = +0.034，95% CI [+0.008, +0.060]，区间全正，门通过。
- **EV1 基本结构外样本**：76 个 Replogle held-out 基因（标签自 Norman 冻结表直传）。
  方向与幅度复现主链（joint 0.404 > human 0.374 > yeast 0.359）；
  与主链合并 283 查询：对最优单臂 p < 2.4e-4，vs human-only p = 0.041，
  vs yeast-only p = 0.0012。
- 目录：results_a6000/externalvalidation/（运行记录）、raw/externalvalidation/
  （冻结输入镜像与补推理嵌入）、scripts/ev*.py（运行器）、
  docs/cross_species_match/EV_PROTOCOL.md（预注册协议）。
- 复现：python 环境见 requirements.txt；输入完整性核验见
  results_a6000/externalvalidation/rehash/（6259/6275 文件逐位一致）。

仓库说明：大体积原始件（模型权重、K562 h5ad、Adamson 原始数据、响应矩阵等）
按 .gitignore 排除，运行前需按 MANIFEST.sha256 与 raw/externalvalidation/
内记录的哈希自行回填；raw/tier0_* 等为指向归档的符号链接占位。
