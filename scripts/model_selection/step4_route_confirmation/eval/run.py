"""评测 runner:统一入口 + MLflow 记录 + 结果落盘。

用法(项目根下):
  /public/home/mengxl/dzy/envs/yeastbridge/bin/python -m eval.run --task t2 --feature esm2_mean --model logistic
  /public/home/mengxl/dzy/envs/yeastbridge/bin/python -m eval.run --task t3 --feature esm2_mean --model ridge
  /public/home/mengxl/dzy/envs/yeastbridge/bin/python -m eval.run --task t3 --model mean_profile   # 平凡基线
  /public/home/mengxl/dzy/envs/yeastbridge/bin/python -m eval.run --task t5                        # 数据效率曲线
  /public/home/mengxl/dzy/envs/yeastbridge/bin/python -m eval.run --list                           # 各任务就绪状态
"""
from __future__ import annotations
import argparse, json, os, sys, time

import mlflow

PROJECT = "/public/home/mengxl/dzy/yeastbridge"
RESULTS = f"{PROJECT}/results/harness"
MLRUNS = f"sqlite:///{PROJECT}/results/mlflow.db"

sys.path.insert(0, PROJECT)
from eval.tasks import REGISTRY  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--task", choices=sorted(REGISTRY))
    ap.add_argument("--feature", default="esm2_mean")
    ap.add_argument("--model", default=None)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--alpha", default="auto", help="ridge 正则强度;auto=val 划分自动选择(默认)")
    ap.add_argument("--train-fraction", type=float, default=1.0)
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--no-mlflow", action="store_true")
    args = ap.parse_args()

    if args.list:
        for t, spec in sorted(REGISTRY.items()):
            ok, why = spec["ready"]()
            print(f"{t.upper()}: {'READY' if ok else 'BLOCKED'} - {why}")
        return
    if not args.task:
        ap.error("--task 或 --list 必填")

    spec = REGISTRY[args.task]
    ok, why = spec["ready"]()
    if not ok:
        raise SystemExit(f"任务 {args.task.upper()} 数据未就绪: {why}")
    model = args.model or {"t1": "kmeans", "t2": "logistic", "t3": "ridge", "t4": "ko_sensitivity", "t5": "ridge"}[args.task]

    t0 = time.time()
    metrics, detail = spec["run"](feature=args.feature, model_name=model, seed=args.seed,
                                  alpha=args.alpha, train_fraction=args.train_fraction)
    dur = time.time() - t0

    os.makedirs(RESULTS, exist_ok=True)
    tag = f"{args.task}_{args.feature}_{model}_s{args.seed}" + (
        f"_fr{args.train_fraction}" if args.train_fraction < 1.0 else "")
    detail_path = f"{RESULTS}/{tag}_detail.csv"
    detail.to_csv(detail_path, index=False)
    record = {"task": args.task, "feature": args.feature, "model": model, "seed": args.seed,
              "metrics": metrics, "seconds": round(dur, 1), "detail": detail_path,
              "time": time.strftime("%Y-%m-%d %H:%M:%S")}
    with open(f"{RESULTS}/{tag}.json", "w") as fh:
        json.dump(record, fh, indent=2, ensure_ascii=False)

    if not args.no_mlflow:
        mlflow.set_tracking_uri(MLRUNS)
        mlflow.set_experiment("yeastbridge_harness")
        with mlflow.start_run(run_name=tag):
            mlflow.log_params({"task": args.task, "feature": args.feature, "model": model,
                               "seed": args.seed, "alpha": args.alpha,
                               "train_fraction": args.train_fraction})
            mlflow.log_metrics({k: v for k, v in metrics.items() if isinstance(v, (int, float))})
            mlflow.log_artifact(detail_path)

    print(json.dumps(record, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
