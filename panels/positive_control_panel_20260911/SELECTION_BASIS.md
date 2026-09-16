# 选对依据（冻结时盲态声明，2026-09-11）

本 allowlist 的选择依据**全部来自 HIP/HOP 数据集之外**：
- amiodarone/quinidine × 其已知通道靶（KCNH2/KCNQ1/KCNA5）：来自 target_screen.json
  注册基准药清单 + ChEMBL 已知靶注释；
- capsaicin×TRPV1、haloperidol×DRD2：demo 对子的已知靶（文献/ChEMBL）。

未查看 strain_response.npz 内容；exec_matrix 既往结果仅用于确认这些对子
**尚未被计算过**（它们不在 17655 行里），没有任何对子因其结果而被加入或移除。
