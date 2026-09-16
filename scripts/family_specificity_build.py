#!/usr/bin/env python3
"""family_specificity_build — 家族特异性检验的任务轴与 allowlist 构建。

协议：docs/cross_species_match/FAMILY_SPECIFICITY_PROTOCOL.md（判据先冻结）
输入（全部冻结产物，只读）：
  - esm2_joint_tasks_dc_20260912/（EJ-dc 任务表 + export_record.json）
  - drug_ko_benchmark_v1/ej_dc_confirmation_panel_20260912/pair_allowlist_{pos,neg}40.tsv
输出：
  - esm2_joint_tasks_dc_fammean_20260913/  （家族均值轴任务表 FAM_<id>）
  - esm2_joint_tasks_dc_resid_20260913/    （提名靶点残差轴任务表）
  - drug_ko_benchmark_v1/family_specificity_panel_20260913/
      fam_allowlist_{pos,neg}.tsv + build_record.json
"""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

ST = Path("/public/home/mengxl/dzy/yeastbridge_re_mvp_stage_20260901")
DC = ST / "esm2_joint_tasks_dc_20260912"
CONF = ST / "drug_ko_benchmark_v1/ej_dc_confirmation_panel_20260912"
FAMDIR = ST / "esm2_joint_tasks_dc_fammean_20260913"
RESIDDIR = ST / "esm2_joint_tasks_dc_resid_20260913"
PANEL = ST / "drug_ko_benchmark_v1/family_specificity_panel_20260913"
CORR_CUT = 0.9  # 协议冻结: average-linkage 距离阈值 0.1
for d in (FAMDIR, RESIDDIR, PANEL):
    d.mkdir(parents=True, exist_ok=True)


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


# ---- 1) 读 dc 分数矩阵 ----
rec = json.loads((DC / "export_record.json").read_text())
tail = set(rec["tail_strains"])
files = sorted(DC.glob("yeast_task_*.tsv"))
print(f"dc task files: {len(files)}", flush=True)
targets, rows, orfs = [], [], None
for fp in files:
    df = pd.read_csv(fp, sep="\t")
    sc = dict(zip(df.yeast_gene, df.ensemble_score))
    if orfs is None:
        orfs = list(df.yeast_gene)
    rows.append([sc[o] for o in orfs])
    targets.append(fp.name[len("yeast_task_"):-len(".tsv")])
S_full = np.asarray(rows, dtype=np.float64)
t_idx = {t: i for i, t in enumerate(targets)}
pool_cols = np.array([o not in tail for o in orfs])
S = S_full[:, pool_cols]
print(f"S: {S.shape} (pool {pool_cols.sum()})", flush=True)

# ---- 2) 聚类（协议冻结配方） ----
C = np.corrcoef(S)
D = 1.0 - C
np.fill_diagonal(D, 0.0)
D = (D + D.T) / 2
Z = linkage(squareform(D, checks=False), method="average")
lab = fcluster(Z, t=1.0 - CORR_CUT, criterion="distance")
lab7 = fcluster(Z, t=1.0 - 0.7, criterion="distance")  # 描述性敏感性
from collections import Counter
sizes = Counter(lab)
print(f"clusters@0.9: {len(sizes)} (size>=2: {sum(1 for v in sizes.values() if v>=2)}), "
      f"largest {max(sizes.values())}; clusters@0.7: {len(Counter(lab7))}", flush=True)

# ---- 3) 提名对 ----
noms = []
for d in ("pos", "neg"):
    a = pd.read_csv(CONF / f"pair_allowlist_{d}40.tsv", sep="\t", dtype=str)
    a["direction"] = "+z" if d == "pos" else "-z"
    noms.append(a)
nom = pd.concat(noms, ignore_index=True)
nom_targets = sorted(set(nom.target_id))
print(f"nominated pairs: {len(nom)}, targets: {len(nom_targets)}", flush=True)

fam_of = {t: int(lab[t_idx[t]]) for t in nom_targets}
fam_members = {}
for t, f in fam_of.items():
    fam_members.setdefault(f, []).append(t)
# 家族 = 含提名靶点的簇；成员 = 全体（含未提名）
fam_all = {}
for f in fam_members:
    fam_all[f] = [targets[i] for i in np.flatnonzero(lab == f)]
singletons = [t for t, f in fam_of.items() if len(fam_all[f]) < 2]
print(f"families with nominations: {len(fam_all)}; nominated singletons: {len(singletons)}", flush=True)

# ---- 4) 家族均值轴任务表 ----
orf_pos = {o: i for i, o in enumerate(orfs)}
tail_list = sorted(tail)


def write_table(dirpath, name, pool_scores):
    full = np.empty(len(orfs))
    full[pool_cols] = pool_scores
    if tail_list:
        lo = float(np.min(pool_scores))
        for a_, o in enumerate(tail_list):
            full[orf_pos[o]] = lo - 0.001 - a_ * 1e-9
    order = np.argsort(-full, kind="stable")
    fp = dirpath / f"yeast_task_{name}.tsv"
    with fp.open("w", encoding="utf-8", newline="") as fh:
        fh.write("rank\tyeast_gene\tensemble_score\n")
        for rk, j in enumerate(order, 1):
            fh.write(f"{rk}\t{orfs[j]}\t{full[j]:.9f}\n")
    return fp


fam_name = {}
for f in sorted(fam_all):
    fid = f"FAM{f:05d}"
    fam_name[f] = fid
    m = S[[t_idx[t] for t in fam_all[f]]].mean(axis=0)
    write_table(FAMDIR, fid, m)
print(f"fam tables: {len(fam_all)}", flush=True)

# ---- 5) 残差轴任务表 + 方差分解 ----
resid_var_frac = {}
for t in nom_targets:
    f = fam_of[t]
    if len(fam_all[f]) < 2:
        continue
    m = S[[t_idx[x] for x in fam_all[f]]].mean(axis=0)
    s = S[t_idx[t]]
    coef = float(s @ m) / float(m @ m)
    r = s - coef * m
    resid_var_frac[t] = float(r @ r) / float(s @ s) if float(s @ s) > 0 else np.nan
    write_table(RESIDDIR, t, r)
print(f"resid tables: {len(resid_var_frac)}", flush=True)

# ---- 6) 家族单元 allowlist（按方向） ----
for d, sym in (("pos", "+z"), ("neg", "-z")):
    sub = nom[nom.direction == sym]
    units = sorted({(fam_name[fam_of[t]], ik) for t, ik in
                    zip(sub.target_id, sub.inchikey)
                    if len(fam_all[fam_of[t]]) >= 2})
    pd.DataFrame(units, columns=["target_id", "inchikey"]).to_csv(
        PANEL / f"fam_allowlist_{d}.tsv", sep="\t", index=False)
    print(f"fam_allowlist_{d}: {len(units)} units", flush=True)

build_rec = {
    "protocol": "FAMILY_SPECIFICITY_PROTOCOL.md",
    "script_sha256": sha(Path(__file__).resolve()),
    "dc_export_record_sha256": sha(DC / "export_record.json"),
    "allowlist_sha256": {f"pair_allowlist_{d}40": sha(CONF / f"pair_allowlist_{d}40.tsv")
                         for d in ("pos", "neg")},
    "corr_cut": CORR_CUT,
    "n_clusters_0.9": len(sizes), "n_clusters_0.7_descriptive": len(Counter(lab7)),
    "families": {fam_name[f]: {"size": len(fam_all[f]),
                               "nominated": sorted(fam_members[f])}
                 for f in sorted(fam_all)},
    "nominated_singletons": singletons,
    "resid_var_frac": {t: round(v, 6) for t, v in resid_var_frac.items()},
}
(PANEL / "build_record.json").write_text(json.dumps(build_rec, indent=1))
for f in sorted(PANEL.iterdir()):
    print(f.name, sha(f))
print("[done]", flush=True)
