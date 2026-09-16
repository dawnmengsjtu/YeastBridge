"""任务层:T1-T5 任务定义(以 docs/plan_v1.0.md 第四节为准)。

每个任务暴露:
  - run(feature, model_name, seed, **kw) -> dict(metrics)  + 详情 DataFrame
  - data_ready() -> (bool, 原因)
当前实现:T2、T3(数据已齐);T1 需细胞编码器(路线A/C/D 产物),T4 缺通路注释,T5 复用 T3 的数据点扫描。
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import train_test_split

from . import data


# ---------------- T2 必需基因预测 ----------------
def t2_data_ready():
    return True, "essential_genes.tsv + esm2 embeddings 就绪"


def run_t2(feature: str = "esm2_mean", model_name: str = "logistic", seed: int = 42, **kw):
    lab = data.essentiality()
    lab = lab[lab["essentiality"].isin(["essential", "nonessential"])]  # 丢弃 unknown/conflicting
    X, kept = data.feature_matrix(feature, lab["systematic"].tolist())
    lab = lab[lab["systematic"].isin(kept)]
    y = (lab["essentiality"] == "essential").astype(int).to_numpy()
    # 保持 X 与 y 行序一致(lab 经 isin 过滤后顺序未变,因为 kept 来自 lab 顺序)
    Xtr, Xte, ytr, yte = train_test_split(X, y, test_size=0.2, stratify=y, random_state=seed)
    Xtr, Xva, ytr, yva = train_test_split(Xtr, ytr, test_size=0.125, stratify=ytr, random_state=seed)  # 0.8*0.125=0.1
    if model_name == "logistic":
        clf = LogisticRegression(max_iter=2000, class_weight="balanced")
    else:
        raise ValueError(model_name)
    clf.fit(Xtr, ytr)
    p = clf.predict_proba(Xte)[:, 1]
    metrics = {
        "auroc": float(roc_auc_score(yte, p)),
        "auprc": float(average_precision_score(yte, p)),
        "n_train": int(len(ytr)), "n_val": int(len(yva)), "n_test": int(len(yte)),
        "prevalence_test": float(yte.mean()),
    }
    detail = pd.DataFrame({"systematic": lab["systematic"].iloc[np.arange(len(lab))].to_numpy()}).iloc[0:0]  # 占位,详情略
    return metrics, detail


# ---------------- T3 扰动响应预测 ----------------
def t3_data_ready():
    return True, "kemmeren parquet + split_assignments + esm2 embeddings 就绪"


def run_t3(feature: str = "esm2_mean", model_name: str = "ridge", seed: int = 42,
           alpha="auto", train_fraction: float = 1.0, **kw):
    """特征=被敲基因的向量;标签=该株全基因组 log2FC;逐株 Spearman + top100 召回。
    附带 mean-profile 平凡基线(预测训练集均值)作锚点。
    alpha="auto": 在 val 划分上从 {0.1,1,10,100} 选 Spearman 最优,再评 test(避免测试集选参)。"""
    split, expr = data.kemmeren()
    mut_sys = split.set_index("mutant_name")["systematic"].to_dict()
    # 对齐:表达矩阵行(mutant 名) → 被敲基因 systematic → 特征
    emb_genes = data.feature_genes(feature)
    rows = []
    for m in expr["mutant"]:
        g = mut_sys.get(m, "")
        rows.append(g if g in emb_genes else None)
    expr = expr.assign(_ko=rows)
    expr = expr[expr["_ko"].notna()].reset_index(drop=True)
    X, kept = data.feature_matrix(feature, expr["_ko"].tolist())
    expr = expr[expr["_ko"].isin(kept)].reset_index(drop=True)
    gene_cols = [c for c in expr.columns if c not in ("mutant", "_ko")]
    Y = expr[gene_cols].to_numpy(np.float32)
    smap = split.set_index("mutant_name")["split"].to_dict()
    sl = expr["mutant"].map(smap)
    tr = sl.eq("train").to_numpy(); va = sl.eq("val").to_numpy(); te = sl.eq("test").to_numpy()
    rng = np.random.RandomState(seed)
    if train_fraction < 1.0:
        idx_tr = np.where(tr)[0]
        tr = np.zeros_like(tr)
        tr[rng.choice(idx_tr, max(1, int(len(idx_tr) * train_fraction)), replace=False)] = True
    if model_name == "ridge":
        if alpha == "auto":
            # val 划分上选 alpha,再评 test
            best = None
            for a in (0.1, 1.0, 10.0, 100.0):
                m_a = Ridge(alpha=a).fit(X[tr], Y[tr])
                Pv = m_a.predict(X[va])
                sv = np.nanmean([spearmanr(Pv[i], Y[va][i]).statistic for i in range(len(Pv))])
                if best is None or sv > best[1]:
                    best = (a, sv)
            mdl = Ridge(alpha=best[0]).fit(X[tr], Y[tr])
            chosen_alpha, val_spearman = best
        else:
            mdl = Ridge(alpha=float(alpha)).fit(X[tr], Y[tr])
            chosen_alpha, val_spearman = float(alpha), np.nan
    elif model_name == "mean_profile":  # 平凡基线
        mdl = None
        chosen_alpha, val_spearman = np.nan, np.nan
    else:
        raise ValueError(model_name)
    if mdl is not None:
        P = mdl.predict(X[te])
    else:
        P = np.repeat(Y[tr].mean(0, keepdims=True), te.sum(), axis=0)
    T = Y[te]
    spear = np.array([spearmanr(P[i], T[i]).statistic for i in range(len(P))])
    # top-100 |log2FC| 召回
    rec = []
    for i in range(len(P)):
        a = set(np.argsort(-np.abs(P[i]))[:100])
        b = set(np.argsort(-np.abs(T[i]))[:100])
        rec.append(len(a & b) / 100)
    metrics = {
        "spearman_mean": float(np.nanmean(spear)),
        "spearman_median": float(np.nanmedian(spear)),
        "top100_recall_mean": float(np.mean(rec)),
        "n_train": int(tr.sum()), "n_test": int(te.sum()),
        "train_fraction": train_fraction,
        "alpha": chosen_alpha, "val_spearman": val_spearman,
    }
    muts = expr["mutant"].to_numpy()[te]
    detail = pd.DataFrame({"mutant": muts, "ko_gene": expr["_ko"].to_numpy()[te],
                           "spearman": spear, "top100_recall": rec})
    return metrics, detail


# ---------------- T5 数据效率曲线 ----------------
def run_t5(feature: str = "esm2_mean", model_name: str = "ridge", seed: int = 42,
           fractions=(0.01, 0.02, 0.03, 0.05, 0.1, 0.25, 0.5, 1.0), **kw):
    """复用 T3,在不同训练数据比例下评估(plan T5:各点 Spearman)。"""
    rows = []
    for fr in fractions:
        m, _ = run_t3(feature=feature, model_name=model_name, seed=seed, train_fraction=fr)
        rows.append({"fraction": fr, **m})
    detail = pd.DataFrame(rows)
    return {"curve": "见 detail", "n_points": len(rows)}, detail


# ---------------- T1 细胞状态解析 ----------------
def t1_data_ready():
    import os
    have = [f for f, p in data._CELL_EMB_PATHS.items() if os.path.exists(p)]
    if have:
        return True, f"细胞 embedding 就绪: {','.join(have)}"
    return False, "缺细胞 embedding(路线A: scripts/routeA/extract_cell_embeddings.py;路线D: scripts/routeD/extract_scyeast_cell_embeddings.py)"


def run_t1(feature: str = "routea_scgpt_ft", model_name: str = "kmeans", seed: int = 42,
           label_key: str = "condition", **kw):
    """零样本细胞 embedding → kmeans 聚类 → 与条件标签比对。
    指标:ARI / NMI / silhouette(cosine, 5k 子采样; 标签作簇的条件可分性)。"""
    from sklearn.cluster import KMeans
    from sklearn.metrics import (adjusted_rand_score, normalized_mutual_info_score,
                                 silhouette_score)
    X, obs_names = data.cell_embeddings(feature)
    lab = data.sc_labels().loc[obs_names]
    y = lab[label_key].astype(str)
    mask = (y.notna() & (y != "")).to_numpy()
    X, y = X[mask], y[mask].to_numpy()
    if model_name != "kmeans":
        raise ValueError(model_name)
    k = int(pd.Series(y).nunique())
    km = KMeans(n_clusters=k, n_init=10, random_state=seed).fit(X)
    rng = np.random.RandomState(seed)
    sub = rng.choice(len(X), min(5000, len(X)), replace=False)
    metrics = {
        "ari": float(adjusted_rand_score(y, km.labels_)),
        "nmi": float(normalized_mutual_info_score(y, km.labels_)),
        "silhouette_cosine": float(silhouette_score(X[sub], y[sub], metric="cosine")),
        "n_cells": int(len(X)), "n_clusters": k, "label_key": label_key,
    }
    detail = pd.DataFrame({"cell": obs_names[mask], "label": y, "cluster": km.labels_})
    # 留出条件子集指标(路线A 零样本公平性臂): 复用全数据聚类标签, 只在留出条件细胞上算 ARI/NMI
    hj = data.PROJECT + "/data/routeA/holdout_conditions.json"
    if feature.endswith("_h"):
        import json as _json
        import os as _os
        if _os.path.exists(hj):
            hold = set(_json.load(open(hj))["conditions"])
            hm = np.array([c in hold for c in y])
            if hm.sum() > 0:
                metrics["ari_holdout"] = float(adjusted_rand_score(y[hm], km.labels_[hm]))
                metrics["nmi_holdout"] = float(normalized_mutual_info_score(y[hm], km.labels_[hm]))
                metrics["n_cells_holdout"] = int(hm.sum())
    return metrics, detail


# ---------------- T4 通路改造排序 ----------------
def t4_data_ready():
    return True, "data/pathway/ 通路清单 + 21 条文献记录就绪;ko_sensitivity 排序 v1"


def _t3_fit(feature, seed, alpha):
    """复用 T3 的对齐逻辑训练 ridge(被敲基因特征 -> 全基因组 log2FC),供 T4 打分。
    alpha="auto" 时在 val 划分上从 {0.1,1,10,100} 选 Spearman 最优(与 run_t3 同口径)。"""
    split, expr = data.kemmeren()
    mut_sys = split.set_index("mutant_name")["systematic"].to_dict()
    emb_genes = data.feature_genes(feature)
    rows = []
    for m in expr["mutant"]:
        g = mut_sys.get(m, "")
        rows.append(g if g in emb_genes else None)
    expr = expr.assign(_ko=rows)
    expr = expr[expr["_ko"].notna()].reset_index(drop=True)
    X, kept = data.feature_matrix(feature, expr["_ko"].tolist())
    expr = expr[expr["_ko"].isin(kept)].reset_index(drop=True)
    gene_cols = [c for c in expr.columns if c not in ("mutant", "_ko")]
    Y = expr[gene_cols].to_numpy(np.float32)
    smap = split.set_index("mutant_name")["split"].to_dict()
    sl = expr["mutant"].map(smap)
    tr = sl.eq("train").to_numpy(); va = sl.eq("val").to_numpy()
    if alpha == "auto":
        best = None
        for a in (0.1, 1.0, 10.0, 100.0):
            m_a = Ridge(alpha=a).fit(X[tr], Y[tr])
            Pv = m_a.predict(X[va])
            sv = np.nanmean([spearmanr(Pv[i], Y[va][i]).statistic for i in range(len(Pv))])
            if best is None or sv > best[1]:
                best = (a, sv)
        alpha = best[0]
    mdl = Ridge(alpha=float(alpha)).fit(X[tr], Y[tr])
    return mdl, gene_cols


def run_t4(feature: str = "esm2_mean", model_name: str = "ko_sensitivity", seed: int = 42,
           alpha: float = 1.0, topk=(1, 3, 5), **kw):
    """KO 敏感性代理排序 v1:
      对通路候选基因 g,用路线模型的 T3-ridge 预测"敲除 g 后通路成员基因的平均 log2FC" s(g)。
      s(g) 很负(敲除则通路塌)→ g 是通量依赖节点 → overexpress 候选(s 升序排前);
      s(g) 为正(敲除反升)→ g 是刹车/竞争支路 → knockdown 候选(s 降序排前)。
    指标:方向一致率 + hit@k(对照:解析随机基线 k/N 与多数类先验)。
    记录处理:多基因组合记录(systematic 含 '/')按成员展开(combo=True);
      混合方向('overexpress + knockdown')记录 v1 不拆,丢弃并在 detail 标注。
    model_name="kemmeren_direct"(上限臂):不经模型,直接用 g 的 Kemmeren 实测 KO 谱打分
      (g 必须是 1484 被测突变体之一,否则标 unusable),用于给模型打分器封顶。
    """
    if model_name not in ("ko_sensitivity", "kemmeren_direct"):
        raise ValueError(model_name)
    rec = data.engineering_records()
    pw = data.pathway_genes()
    rows = []
    for _, r in rec.iterrows():
        for g in r["systematic"].split("/"):
            g = g.strip()
            if not g:
                continue
            rows.append({"pathway_id": r["pathway_id"], "systematic": g,
                         "direction": r["direction"], "combo": "/" in r["systematic"],
                         "reference": r["reference"]})
    R = pd.DataFrame(rows)
    if model_name == "kemmeren_direct":
        split, expr = data.kemmeren()
        mut_sys = split.set_index("mutant_name")["systematic"].to_dict()
        kem_cols = [c for c in expr.columns if c != "mutant"]
        Yall = expr[kem_cols].to_numpy(np.float32)
        profiles = {}
        for i, m in enumerate(expr["mutant"]):
            g0 = mut_sys.get(m, "")
            if g0:
                profiles[g0] = Yall[i]
        mdl = None
    else:
        mdl, kem_cols = _t3_fit(feature, seed, alpha)
        profiles = None
    col_pos = {g: i for i, g in enumerate(kem_cols)}
    kem_set = set(kem_cols)

    detail_rows = []
    for pid, grp in R.groupby("pathway_id"):
        pw_genes = pw.loc[pw["pathway_id"] == pid, "systematic"].tolist()
        universe = sorted(set(pw_genes) | set(grp["systematic"]))
        out_genes = [g for g in pw_genes if g in kem_set]
        scores = {}
        if mdl is None:  # kemmeren_direct: 实测谱打分
            for g in universe:
                if g in profiles:
                    idx = [col_pos[og] for og in out_genes if og != g]
                    scores[g] = float(profiles[g][idx].mean()) if idx else np.nan
        else:
            Xu, kept_u = data.feature_matrix(feature, universe)
            if not kept_u:
                continue
            P = mdl.predict(Xu)
            for g, p in zip(kept_u, P):
                idx = [col_pos[og] for og in out_genes if og != g]
                scores[g] = float(p[idx].mean()) if idx else np.nan
        order_oe = sorted(scores, key=lambda g: scores[g])          # 越负越是过表达候选
        order_kd = sorted(scores, key=lambda g: -scores[g])         # 越正越是敲低候选
        for _, r in grp.iterrows():
            g = r["systematic"]
            d = r.to_dict()
            if r["direction"] not in ("overexpress", "knockdown"):
                d.update(usable=False, note="混合方向组合,v1 不拆")
            elif g not in scores or np.isnan(scores[g]):
                d.update(usable=False, note="该特征无此基因 embedding 或无通路输出列")
            else:
                order = order_oe if r["direction"] == "overexpress" else order_kd
                rank = order.index(g) + 1
                d.update(usable=True, score=scores[g], rank=rank, n_candidates=len(order),
                         pred_direction="overexpress" if scores[g] < 0 else "knockdown",
                         direction_match=(scores[g] < 0) == (r["direction"] == "overexpress"))
            detail_rows.append(d)
    D = pd.DataFrame(detail_rows)
    U = D[D["usable"] == True]  # noqa: E712
    metrics = {
        "n_records": int(len(D)), "n_usable": int(len(U)),
        "direction_consistency": float(U["direction_match"].mean()) if len(U) else np.nan,
        "chance_direction": float((U["direction"] == "overexpress").mean()) if len(U) else np.nan,
    }
    for k in topk:
        metrics[f"hit_at_{k}"] = float((U["rank"] <= k).mean()) if len(U) else np.nan
        metrics[f"chance_hit_at_{k}"] = float((np.minimum(k, U["n_candidates"]) / U["n_candidates"]).mean()) if len(U) else np.nan
    return metrics, D


# ---------------- T4 旧 stub 已移除(逻辑已实现于上) ----------------


REGISTRY = {
    "t1": {"ready": t1_data_ready, "run": run_t1},
    "t2": {"ready": t2_data_ready, "run": run_t2},
    "t3": {"ready": t3_data_ready, "run": run_t3},
    "t4": {"ready": t4_data_ready, "run": run_t4},
    "t5": {"ready": t3_data_ready, "run": run_t5},
}
