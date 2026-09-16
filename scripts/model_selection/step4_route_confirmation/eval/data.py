"""数据加载层:统一读取 YeastBridge 已就绪资产,基因 ID 以 SGD systematic name 对齐。

所有函数返回 pandas/numpy 对象,首次调用后缓存(进程内)。
路径硬编码到项目根,勿在项目外使用。
"""
from __future__ import annotations
import functools
import numpy as np
import pandas as pd

PROJECT = "/public/home/mengxl/dzy/yeastbridge"


@functools.lru_cache(maxsize=None)
def gene_master() -> pd.DataFrame:
    return pd.read_csv(f"{PROJECT}/data/mappings/gene_master.tsv", sep="\t", dtype=str).fillna("")


@functools.lru_cache(maxsize=None)
def esm2_embeddings() -> tuple[pd.DataFrame, np.ndarray]:
    """返回 (index_df, X):index_df 列 systematic/common/uniprot/seq_len,X 为 Nx1280 float32。"""
    idx = pd.read_csv(f"{PROJECT}/data/processed/esm2_650m/index.tsv", sep="\t", dtype=str).fillna("")
    X = np.load(f"{PROJECT}/data/processed/esm2_650m/esm2_mean_fp32.npy")
    assert len(idx) == len(X)
    return idx, X


@functools.lru_cache(maxsize=None)
def essentiality() -> pd.DataFrame:
    """T2 标签:systematic/common/essentiality/evidence/source。"""
    return pd.read_csv(f"{PROJECT}/data/essential/essential_genes.tsv", sep="\t", dtype=str).fillna("")


@functools.lru_cache(maxsize=None)
def kemmeren() -> tuple[pd.DataFrame, pd.DataFrame]:
    """T3 数据:返回 (meta, expr)。expr 行=突变体,列=systematic ID,值=log2FC;
    meta 含 mutant/split 列(按扰动基因 7:1:2, seed=42)。"""
    expr = pd.read_parquet(f"{PROJECT}/data/ko_kemmeren/kemmeren_2014/kemmeren_t3_expr_log2fc.parquet")
    split = pd.read_csv(f"{PROJECT}/data/ko_kemmeren/kemmeren_2014/split_assignments.tsv", sep="\t", dtype=str)
    return split, expr


def esm2_feature_matrix(genes: list[str]) -> tuple[np.ndarray, list[str]]:
    """按给定 systematic ID 列表取 ESM2 特征;缺失的基因被丢弃。
    返回 (X, kept_genes),行序与 kept_genes 一致。"""
    return feature_matrix("esm2_mean", genes)


# ---------------- 路线A 特征 ----------------
@functools.lru_cache(maxsize=None)
def routea_gene_embeddings(kind: str = "ft") -> tuple[pd.DataFrame, np.ndarray]:
    """路线A 基因 embedding。kind='ft': 微调后 encoder.embedding.weight;
    kind='init': 微调前的同源初始化表(对照臂, 分离"人类先验"与"酵母微调"的贡献)。
    返回 (index_df[systematic], X[N,512]),仅含真实基因(特殊符号除外)。"""
    import json
    vocab = json.load(open(f"{PROJECT}/data/routeA/routeA_vocab.json"))
    if kind == "ft":
        import torch
        sd = torch.load(f"{PROJECT}/models/routeA/scgpt_yeast_ft.pt", map_location="cpu", weights_only=True)
        W = sd["encoder.embedding.weight"].numpy().astype(np.float32)
    elif kind == "scratch":
        import torch
        sd = torch.load(f"{PROJECT}/models/routeC/scgpt_yeast_scratch.pt", map_location="cpu", weights_only=True)
        W = sd["encoder.embedding.weight"].numpy().astype(np.float32)
    elif kind == "init":
        W = np.load(f"{PROJECT}/data/routeA/routeA_init_embeddings.npy").astype(np.float32)
    else:
        raise ValueError(kind)
    order = sorted(((g, i) for g, i in vocab.items() if not g.startswith("<")), key=lambda x: x[1])
    names = [g for g, _ in order]
    X = W[[i for _, i in order]]
    return pd.DataFrame({"systematic": names}), X


_FEATURE_LOADERS = {
    "esm2_mean": lambda: esm2_embeddings(),
    "routea_scgpt_ft": lambda: routea_gene_embeddings("ft"),
    "routea_init": lambda: routea_gene_embeddings("init"),
    "routec_scgpt": lambda: routea_gene_embeddings("scratch"),
    "routeb_protein": lambda: routeb_gene_embeddings(),
    "routed_scyeast": lambda: scyeast_gene_embeddings(),
}


def feature_genes(feature: str) -> set:
    """该特征可用的 systematic 基因集合。"""
    if feature not in _FEATURE_LOADERS:
        raise ValueError(f"未知特征 {feature},可选 {sorted(_FEATURE_LOADERS)}")
    idx, _ = _FEATURE_LOADERS[feature]()
    return set(idx["systematic"])


def feature_matrix(feature: str, genes: list[str]) -> tuple[np.ndarray, list[str]]:
    """通用特征矩阵:按 systematic 列表取向量,缺失丢弃,返回 (X, kept)。"""
    if feature not in _FEATURE_LOADERS:
        raise ValueError(f"未知特征 {feature},可选 {sorted(_FEATURE_LOADERS)}")
    idx, Xall = _FEATURE_LOADERS[feature]()
    pos = {g: i for i, g in enumerate(idx["systematic"])}
    take = [(i, g) for i, g in enumerate(genes) if g in pos]
    kept = [g for _, g in take]
    dim = Xall.shape[1]
    X = np.stack([Xall[pos[g]] for _, g in take]) if take else np.zeros((0, dim), np.float32)
    return X, kept


# ---------------- T4 通路改造排序输入 ----------------
@functools.lru_cache(maxsize=None)
def engineering_records() -> pd.DataFrame:
    """文献改造记录(21 条 DOI 核验):pathway_id/gene/systematic/direction/method/effect/product/reference/doi。"""
    return pd.read_csv(f"{PROJECT}/data/pathway/engineering_records.tsv", sep="\t", dtype=str).fillna("")


@functools.lru_cache(maxsize=None)
def pathway_genes() -> pd.DataFrame:
    """KEGG 142 条通路的基因清单:pathway_id/pathway_name/systematic/common。"""
    return pd.read_csv(f"{PROJECT}/data/pathway/kegg_sce_pathway_genes.tsv", sep="\t", dtype=str).fillna("")


# ---------------- 路线D scYeast 特征 ----------------
@functools.lru_cache(maxsize=None)
def scyeast_gene_embeddings() -> tuple[pd.DataFrame, np.ndarray]:
    """scYeast 预训练的逐基因嵌入表 pos_embedding.emb.weight(5812×200, w_knowledge ckpt)。
    行序由 models/scyeast/gene_order_5812.txt 给出(由 ckpt 与 gene2vec 精确匹配恢复)。
    注意:这是静态基因表征;scYeast 的样本级表征需经表达门控融合,另行接入。"""
    import torch
    ck = f"{PROJECT}/src/external/scYeast/models/final_checkpoint_spearman0.5_w_knowledge_huber_v4.pth"
    sd = torch.load(ck, map_location="cpu", weights_only=True)["net"]
    W = sd["pos_embedding.emb.weight"].numpy().astype(np.float32)
    genes = [l.strip() for l in open(f"{PROJECT}/models/scyeast/gene_order_5812.txt") if l.strip()]
    assert len(genes) == W.shape[0], f"基因序 {len(genes)} 与嵌入表 {W.shape[0]} 不一致"
    return pd.DataFrame({"systematic": genes}), W


# ---------------- 路线B 蛋白桥接特征 ----------------
@functools.lru_cache(maxsize=None)
def routeb_gene_embeddings() -> tuple[pd.DataFrame, np.ndarray]:
    """路线B 微调后的基因 embedding(protein_inject 注入层展开, 6736×512, 含特殊符号行)。
    训练后由 finetune 脚本落盘 models/routeB/routeb_gene_embeddings.npy,行序=routeA_vocab 基因序。"""
    gm = gene_master()
    W = np.load(f"{PROJECT}/models/routeB/routeb_gene_embeddings.npy").astype(np.float32)
    return pd.DataFrame({"systematic": gm["systematic"]}), W[: len(gm)]


# ---------------- T1 细胞级输入 ----------------@functools.lru_cache(maxsize=None)
def sc_labels() -> pd.DataFrame:
    """GSE125162 细胞标签(genotype/condition),索引=obs_name。"""
    import anndata as ad
    a = ad.read_h5ad(f"{PROJECT}/data/processed/sc_gse125162/counts_raw.h5ad", backed="r")
    return a.obs[["genotype", "condition"]].copy()


_CELL_EMB_PATHS = {
    "routea_scgpt_ft": f"{PROJECT}/data/routeA/cell_emb_gse125162.npz",
    "routea_scgpt_ft_h": f"{PROJECT}/data/routeA/cell_emb_gse125162_holdout.npz",
    "routec_scgpt": f"{PROJECT}/data/routeC/cell_emb_gse125162.npz",
    "routeb_protein": f"{PROJECT}/data/routeB/cell_emb_gse125162.npz",
    "routeb_protein_h": f"{PROJECT}/data/routeB/cell_emb_gse125162_holdout.npz",
    "routed_scyeast": f"{PROJECT}/data/routeD/cell_emb_gse125162.npz",
    "routed_scyeast_tpm": f"{PROJECT}/data/routeD/cell_emb_gse125162_tpm.npz",
    "routed_scyeast_tpm_exprpool": f"{PROJECT}/data/routeD/cell_emb_gse125162_tpm_exprpool.npz",
}


@functools.lru_cache(maxsize=None)
def cell_embeddings(feature: str) -> tuple[np.ndarray, np.ndarray]:
    """细胞 embedding。返回 (X,N,D 的 X[N,D], obs_names[N])。支持 routea_scgpt_ft / routed_scyeast。"""
    if feature not in _CELL_EMB_PATHS:
        raise ValueError(f"特征 {feature} 无细胞级 embedding(可选 {sorted(_CELL_EMB_PATHS)})")
    d = np.load(_CELL_EMB_PATHS[feature], allow_pickle=True)
    return d["emb"].astype(np.float32), d["obs_names"].astype(str)
