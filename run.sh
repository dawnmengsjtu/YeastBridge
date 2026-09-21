#!/usr/bin/env bash
# 主运行入口（附件5 第二节）：完整链 = 任务轴导出 -> 共模去除 -> 酵母执行筛选 -> 标准候选清单
# 运行前提：按 README 配置环境；大体积输入（嵌入/模型/响应矩阵）按 data/DATA_SOURCES.md
# 的哈希清单回填到 raw/ 后执行。全链 CPU 即可；完整运行约数小时（3,818,750 对）。
set -e
PY=${PY:-python}
cd "$(dirname "$0")"
# 1) 任务轴导出（1,175 靶点 -> 酵母株排序）
$PY scripts/esm2_joint_task_export.py --config configs_a6000_frozen/esm2_joint_tasks.json --output repro/tasks_ej
# 2) 共模双去除（EJ-dc）
$PY scripts/task_axis_double_center.py --input repro/tasks_ej --output repro/tasks_ej_dc
# 3) 执行筛选（双向；化合物响应矩阵与 allowlist 见 configs_a6000_frozen/product_execute.esm2joint_dc*.json）
$PY scripts/product_execute_hiphop.py --config configs_a6000_frozen/product_execute.esm2joint_dc.json --results-suffix _screen_repro
$PY scripts/product_execute_hiphop.py --config configs_a6000_frozen/product_execute.esm2joint_dc_neg.json --results-suffix _neg_screen_repro
# 4) 标准候选清单（读取在册确认结果生成 results.csv；换为新跑结果时改 EH 路径）
$PY predict.py
echo "完成：results.csv"
