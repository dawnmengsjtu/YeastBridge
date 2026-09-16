#!/usr/bin/env python3
"""run_esm2joint_chain — esm2-joint 线轴级总脚本（yeastbridge_agent.py --esm2-joint-chain 调用）。

链（全部注册脚本 + 冻结配置，无临时代码；用户指令 2026-09-13）：
  v9       三臂两门 esm2（V9，最终得出 joint 显著最优的版本；G3 含检索三臂）
  export   EJ 任务轴（1175 靶点 × 1282 株，RRF20 三成员）
  dc       EJ-dc 共模去除（列均值 + PC1）
  panels   EJ-dc 轴冻结 allowlist 面板（poscon 4 对 + gpcr 8 对）
  screens  EJ-dc 双向全屏（重；仅 steps 显式含 screens 时跑）
  nominate 两阶段机械提名（ej_dc_nominate.py，规则冻结于 DC 附录 §3）
  confirm  100k 置换双向确认（执行器 two-stage 披露模式）
  family   家族特异性（family_specificity_build + fam/resid 四运行）

行为：步骤产物已存在 → SKIP_EXISTS（审计模式）；--force 重跑到
{stage}/repro_chain/<step>/，绝不覆盖在册产物。
用法：
  python scripts/run_esm2joint_chain.py [--steps all|v9,export,...] [--stage-root DIR]
      [--py BIN] [--force]
"""
import argparse
import json
import subprocess
import time
from pathlib import Path

DEF_STAGE = Path("/public/home/mengxl/dzy/yeastbridge_re_mvp_stage_20260901")
DEF_PY = "/public/home/mengxl/dzy/envs/yeastbridge/bin/python"

# 冻结 allowlist 常量（sha 与在册记录一致）
POSCON_SHA = "59a19823656ef4d0550811a29a4a08c010c78a8c58cef6c97f39e7d58af08b93"
GPCR_SHA = "d0ce81401be879c66d039a6c2deefb285340d3180e50913896c7423c885312ed"
RULE_DOC = "docs/cross_species_match/ESM2_JOINT_DC_ADDENDUM.md"
FAM_RULE_DOC = "docs/cross_species_match/FAMILY_SPECIFICITY_PROTOCOL.md"
BLIND_BASIS = ("rerun of frozen panel allowlist (declared blind at freeze); task axis "
               "= EJ-dc export (frozen recipe RRF-20: direct_pc1 + aligned_C + "
               "b2_verbatim + common-mode removal; user ruling 2026-09-12); "
               "allowlist bytes identical to frozen sha")
STEP_NAMES = ["v9", "export", "dc", "panels", "screens", "nominate",
              "confirm", "family"]


class Ctx:
    def __init__(self, stage, py, force):
        self.stage = Path(stage)
        self.mvp = self.stage / "yeastbridge_re_mvp"
        self.scripts = self.mvp / "scripts"
        self.cfgs = self.mvp / "configs"
        self.eh = self.stage / "product_pipeline/product/execute_hiphop"
        self.panels = self.stage / "drug_ko_benchmark_v1"
        self.py, self.force = py, force
        self.repro = self.stage / "repro_chain"

    def out(self, canonical, step):
        """产物路径：正典目录（缺失时运行写入）或 --force 时的 repro 目录。"""
        c = self.stage / canonical
        if self.force:
            return self.repro / step
        return c


def run_cmd(cmd, timeout, log):
    print("  $ " + " ".join(str(c) for c in cmd[:8]) + (" ..." if len(cmd) > 8 else ""),
          flush=True)
    t0 = time.time()
    r = subprocess.run([str(c) for c in cmd], capture_output=True, text=True,
                       timeout=timeout)
    dt = time.time() - t0
    tail = (r.stdout or "").strip().splitlines()[-3:]
    for line in tail:
        print("    " + line[:300], flush=True)
    if r.returncode != 0:
        print("    STDERR: " + (r.stderr or "")[-800:], flush=True)
    log.append({"cmd": [str(c) for c in cmd], "rc": r.returncode,
                "seconds": round(dt, 1)})
    return r.returncode


def executor_call(ctx, config, suffix, allowlist=None, sha=None, blind=False,
                  rule_doc=None, basis=None, n_perm=None):
    cmd = [ctx.py, ctx.scripts / "product_execute_hiphop.py",
           "--config", ctx.cfgs / config, "--results-suffix", suffix]
    if n_perm:
        cmd += ["--n-perm", str(n_perm)]
    if allowlist:
        cmd += ["--pair-allowlist", allowlist, "--pair-allowlist-sha256", sha]
        if blind:
            cmd += ["--allowlist-declared-hiphop-blind"]
        if rule_doc:
            cmd += ["--allowlist-two-stage-rule-doc", ctx.mvp / rule_doc]
        cmd += ["--allowlist-selection-basis",
                basis or BLIND_BASIS]
    return cmd


def step_v9(ctx, log):
    out = ctx.out("v9_esm2_basic_20260912/results", "v9")
    if (out / "result.json").exists() and not ctx.force:
        return "SKIP_EXISTS"
    out.mkdir(parents=True, exist_ok=True)
    rc = run_cmd([ctx.py, ctx.scripts / "v9_esm2_basic_run.py",
                  "--config", ctx.cfgs / "v9_esm2_basic.json", "--output", out],
                 10800, log)
    return "OK" if rc == 0 else f"ERR:{rc}"


def step_export(ctx, log):
    out = ctx.out("esm2_joint_tasks_20260912", "export")
    if (out / "export_record.json").exists() and not ctx.force:
        return "SKIP_EXISTS"
    rc = run_cmd([ctx.py, ctx.scripts / "esm2_joint_task_export.py",
                  "--config", ctx.cfgs / "esm2_joint_tasks.json", "--output", out],
                 7200, log)
    return "OK" if rc == 0 else f"ERR:{rc}"


def step_dc(ctx, log):
    src = ctx.stage / "esm2_joint_tasks_20260912"
    out = ctx.out("esm2_joint_tasks_dc_20260912", "dc")
    if (out / "export_record.json").exists() and not ctx.force:
        return "SKIP_EXISTS"
    rc = run_cmd([ctx.py, ctx.scripts / "task_axis_double_center.py",
                  "--input", src, "--output", out], 3600, log)
    return "OK" if rc == 0 else f"ERR:{rc}"


def step_panels(ctx, log):
    statuses = []
    for name, cfg, sha, panel in (
            ("poscon", "product_execute.esm2joint_dc.json", POSCON_SHA,
             "positive_control_panel_20260911"),
            ("gpcr", "product_execute.esm2joint_dc.json", GPCR_SHA,
             "gpcr_panel_20260911")):
        resdir = ctx.eh / f"results_esm2jointdc_{name}_20260912"
        if (resdir / "execute_summary.json").exists() and not ctx.force:
            statuses.append(f"{name}:SKIP_EXISTS")
            continue
        suffix = f"_{name}_20260912" if not ctx.force else f"_{name}_repro"
        rc = run_cmd(executor_call(
            ctx, cfg, suffix,
            allowlist=ctx.panels / panel / "pair_allowlist.tsv",
            sha=sha, blind=True), 3600, log)
        statuses.append(f"{name}:{'OK' if rc == 0 else 'ERR'}")
    return ",".join(statuses)


def step_screens(ctx, log):
    statuses = []
    for name, cfg in (("pos", "product_execute.esm2joint_dc.json"),
                      ("neg", "product_execute.esm2joint_dc_neg.json")):
        canon = ctx.eh / (f"results_esm2jointdc_screen_20260912" if name == "pos"
                          else f"results_esm2jointdc_neg_screen_20260912")
        if (canon / "execute_summary.json").exists() and not ctx.force:
            statuses.append(f"{name}:SKIP_EXISTS")
            continue
        suffix = "_screen_20260912" if not ctx.force else "_screen_repro"
        rc = run_cmd(executor_call(ctx, cfg, suffix), 21600, log)
        statuses.append(f"{name}:{'OK' if rc == 0 else 'ERR'}")
    return ",".join(statuses)


def step_nominate(ctx, log):
    out = ctx.out("drug_ko_benchmark_v1/ej_dc_confirmation_panel_20260912", "nominate")
    if (out / "nomination_record.json").exists() and not ctx.force:
        return "SKIP_EXISTS"
    rc = run_cmd([ctx.py, ctx.scripts / "ej_dc_nominate.py",
                  "--pos-screen", ctx.eh / "results_esm2jointdc_screen_20260912",
                  "--neg-screen", ctx.eh / "results_esm2jointdc_neg_screen_20260912",
                  "--out", out], 1800, log)
    return "OK" if rc == 0 else f"ERR:{rc}"


def step_confirm(ctx, log):
    pan = ctx.stage / "drug_ko_benchmark_v1/ej_dc_confirmation_panel_20260912"
    rec = json.loads((pan / "nomination_record.json").read_text())
    statuses = []
    for tag, cfg, resdir in (
            ("pos", "product_execute.esm2joint_dc.json",
             "results_esm2jointdc_confirm100k_pos_20260912"),
            ("neg", "product_execute.esm2joint_dc_neg.json",
             "results_esm2jointdc_neg_confirm100k_neg_20260912")):
        if (ctx.eh / resdir / "execute_summary.json").exists() and not ctx.force:
            statuses.append(f"{tag}:SKIP_EXISTS")
            continue
        sha = rec.get(f"pair_allowlist_{tag}40.tsv.sha256")
        al = pan / f"pair_allowlist_{tag}40.tsv"
        if sha is None:
            import hashlib
            sha = hashlib.sha256(open(al, "rb").read()).hexdigest()
        suffix = (f"_confirm100k_{tag}_20260912" if not ctx.force
                  else f"_confirm100k_{tag}_repro")
        basis = (f"stage-2 confirmation family nominated by frozen mechanical rule "
                 f"(ESM2_JOINT_DC_ADDENDUM section 3) from stage-1 EJ-dc screens "
                 f"(direction {tag} top40 by z); outcome-derived by design; FDR "
                 f"scoped within 80-pair union family; structure in nomination_record")
        rc = run_cmd(executor_call(ctx, cfg, suffix, allowlist=al, sha=sha,
                                   blind=False, rule_doc=RULE_DOC, basis=basis,
                                   n_perm=100000), 21600, log)
        statuses.append(f"{tag}:{'OK' if rc == 0 else 'ERR'}")
    return ",".join(statuses)


def step_family(ctx, log):
    pan = ctx.stage / "drug_ko_benchmark_v1/family_specificity_panel_20260913"
    if (pan / "build_record.json").exists() and not ctx.force:
        build_st = "SKIP_EXISTS"
    else:
        rc = run_cmd([ctx.py, ctx.scripts / "family_specificity_build.py",
                      "--stage", ctx.stage], 7200, log)
        build_st = "OK" if rc == 0 else f"ERR:{rc}"
        if rc != 0:
            return build_st
    runs = []
    for tag, cfg, al, resdir in (
            ("fam_pos", "product_execute.esm2jointdc_fam.json",
             pan / "fam_allowlist_pos.tsv", "results_esm2jointdc_fam_pos_20260913"),
            ("fam_neg", "product_execute.esm2jointdc_fam_neg.json",
             pan / "fam_allowlist_neg.tsv", "results_esm2jointdc_fam_neg_neg_20260913"),
            ("resid_pos", "product_execute.esm2jointdc_resid.json",
             pan / "resid_allowlist_pos.tsv", "results_esm2jointdc_resid_pos_20260913"),
            ("resid_neg", "product_execute.esm2jointdc_resid_neg.json",
             pan / "resid_allowlist_neg.tsv", "results_esm2jointdc_resid_neg_neg_20260913")):
        if (ctx.eh / resdir / "execute_summary.json").exists() and not ctx.force:
            runs.append(f"{tag}:SKIP_EXISTS")
            continue
        import hashlib
        sha = hashlib.sha256(open(al, "rb").read()).hexdigest()
        suffix = ("_pos_20260913" if tag.endswith("pos") else "_neg_20260913") \
            if not ctx.force else f"_{tag}_repro"
        basis = (f"{tag} family derived from frozen confirmation family per "
                 f"FAMILY_SPECIFICITY_PROTOCOL (mechanical rule, outcome-derived "
                 f"by design, disclosed two-stage mode)")
        rc = run_cmd(executor_call(ctx, cfg, suffix, allowlist=al, sha=sha,
                                   blind=False, rule_doc=FAM_RULE_DOC,
                                   basis=basis, n_perm=10000), 7200, log)
        runs.append(f"{tag}:{'OK' if rc == 0 else 'ERR'}")
    return build_st + ";" + ",".join(runs)


STEPS = {"v9": step_v9, "export": step_export, "dc": step_dc,
         "panels": step_panels, "screens": step_screens, "nominate": step_nominate,
         "confirm": step_confirm, "family": step_family}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--steps", default="all")
    ap.add_argument("--stage-root", default=str(DEF_STAGE))
    ap.add_argument("--py", default=DEF_PY)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    ctx = Ctx(args.stage_root, args.py, args.force)
    names = STEP_NAMES if args.steps == "all" else \
        [s.strip() for s in args.steps.split(",") if s.strip()]
    # screens 重步骤：all 不含 screens（需显式点名），其余按序
    if args.steps == "all":
        names = [n for n in names if n != "screens"]
    log = []
    summary = {}
    for name in names:
        print(f"[chain] step {name} ...", flush=True)
        t0 = time.time()
        try:
            st = STEPS[name](ctx, log)
        except subprocess.TimeoutExpired:
            st = "TIMEOUT"
        except Exception as ex:  # noqa: BLE001
            st = f"EXC:{str(ex)[:80]}"
        summary[name] = {"status": st, "seconds": round(time.time() - t0, 1)}
        print(f"[chain] step {name}: {st}", flush=True)
    outdir = ctx.stage / "chain_runs"
    outdir.mkdir(parents=True, exist_ok=True)
    fp = outdir / f"chain_status_{time.strftime('%Y%m%d_%H%M%S')}.json"
    fp.write_text(json.dumps({"steps": summary, "commands": log,
                              "stage_root": str(ctx.stage), "force": ctx.force},
                             indent=1, ensure_ascii=False))
    print("[chain] status ->", fp)
    bad = [k for k, v in summary.items() if "ERR" in str(v) or "TIMEOUT" in str(v)
           or "EXC" in str(v)]
    print("[chain] " + ("ALL OK/SKIP" if not bad else f"FAILED STEPS: {bad}"))
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
