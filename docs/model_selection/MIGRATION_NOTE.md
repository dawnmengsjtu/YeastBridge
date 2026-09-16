# 迁移说明(2026-09-15,按项目惯例归位)

自 yeastbridge_combine / yeastbridge_re_mvp_stage_20260901 / yeastbridge / yeastbridge_re
迁移;仅复制脚本与结果,未重新运行。归位:

- 脚本:scripts/model_selection/step1_scF_selection/norman_retrieval.py;
  scripts/model_selection/step4_route_confirmation/{build_routeA_init,finetune_scgpt,
  extract_cell_embeddings,run_ablations}.py + eval/(原 harness 逐字)。
- 结果:results_model_selection/step1_scF_selection/retrieval_result.json;
  step2_scyeast_exclusion/repr_grid_result.json;
  step4_route_confirmation/five_route_results.json(v2 已作废,供 ERRATUM 对照)。
- 冻结输入:raw/model_selection/norman_assets/raw_oof_full.npz;
  raw/model_selection/norman_features/{scgpt,geneformer,scfoundation,*_target}.npz。
- 文档:docs/model_selection/{MODEL_SELECTION_EVIDENCE,RESULTS_v3,REGISTRATION_v3,ERRATUM}.md。

路径改写(已执行):
- norman_retrieval.py: ROOT/SRC/FEATS/OUT 全部指向本项目(raw/… 与 results_model_selection/step1/),
  资产已随迁,可直接复跑(需 sklearn 环境)。
- finetune_scgpt.py / build_routeA_init.py: 运行命令 cd 与 python 改为本项目 env2;
  MLFLOW_URI 改为本项目 results_model_selection/step4_route_confirmation/mlflow.db。
- 未随迁的大资产(酵母语料 h5ad、routeA init npy、finetune 权重)仍留原冻结绝对路径,
  脚本内以注释标明;重跑 step4 前需先镜像这些资产或改注释内路径。
- step2 生成脚本 = 本项目 scripts/v5_joint_formal_run.py 的 repr_grid()(2026-09-12 逐字抽出);
  其输入 scyeast_prior_space.npy(4.5M)原路径 v52_assets/,未随迁(结果已在册)。
- T1(cell state)为 REGISTERED DEFERRED,不在四步链上,不补跑。
