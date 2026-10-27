#!/usr/bin/env python3
"""Lee 2014 per-screen FD 重算管线 v1（C 方案）。

从 E-MTAB-2391 原始 CEL（GenFlex Tag16K v2, CEL v3 文本格式）重算单个 screen 的
per-strain fitness defect，并与作者发布的 FD 文件（金标准）比对校准。

方法（v1 初版，以作者 FD 为目标实证调优）：
1. 解析 CEL v3 Intensity 段 -> cell -> 强度
2. probeset_map: probeset -> 2 cells, 取均值, log2
3. log2FC = log2(treated) - log2(batch DMSO median)
4. 实验内跨菌株 z 变换 -> recalc FD
5. 与作者 FD 比对: Pearson/Spearman/top-K 菌株重合

明确口径: 输出为"按公开方法重算"，非作者发布文件；作者文件的 expt.zyg.rank 等
全库排名字段不可复现。
"""
import argparse
import gzip
import sys
import json
import hashlib
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/public/home/mengxl/dzy/yeastbridge_vs/data/external/"
            "paired_compounds/hiphop_lee2014_e_mtab_2391")


def parse_cel_v3(path: Path) -> np.ndarray:
    """解析 CEL v3 文本格式的 Intensity 段，返回 (rows*cols,) 强度数组。"""
    opener = gzip.open if str(path).endswith(".gz") else open
    intens = None
    n = None
    with opener(path, "rt", errors="replace") as fh:
        section = None
        for line in fh:
            line = line.strip()
            if line.startswith("["):
                section = line
                continue
            if section == "[HEADER]":
                if line.startswith("Cols="):
                    cols = int(line.split("=")[1])
                elif line.startswith("Rows="):
                    rows = int(line.split("=")[1])
                    n = rows * cols
                    intens = np.full(n, np.nan)
            elif section == "[INTENSITY]" and n is not None:
                parts = line.split()
                if len(parts) in (3, 5):
                    try:
                        if len(parts) == 3:
                            x, y, v = int(parts[0]), int(parts[1]), float(parts[2])
                        else:
                            x, y, v = int(parts[0]), int(parts[1]), float(parts[2])
                    except ValueError:
                        continue
                    intens[y * cols + x] = v
    if intens is None or np.isnan(intens).all():
        raise ValueError(f"no intensity parsed from {path}")
    return intens


def load_probeset_map() -> pd.DataFrame:
    df = pd.read_csv(ROOT / "probeset_map.tsv", sep="\t", header=None,
                     names=["probeset", "cell"])
    return df


def barcode_array(args):
    """Read the existing CEL parser, then summarize five spatial tag replicates."""
    path, cells = args
    intensity = parse_cel_v3(Path(path))
    values = intensity[cells]
    if not np.isfinite(values).all() or (values <= 0).any():
        raise ValueError(f"non-positive or missing barcode intensity: {path}")
    # The trimmed mean reduces the influence of one damaged spatial replicate.
    return np.sort(values, axis=1)[:, 1:-1].mean(axis=1).astype(np.float32)


def build_barcode_response(args):
    """Rebuild the complete response using tag probes; never ORF-probe controls."""
    from scipy.stats import rankdata
    out = Path(args.out)
    if out.exists() and any(out.iterdir()):
        raise ValueError('Output must be a new or empty directory')
    out.mkdir(parents=True, exist_ok=True)
    data_root = Path(args.data_root)
    pmap = pd.read_csv(args.barcode_map, sep="\t")
    if pmap.cell.duplicated().any() or pmap.groupby("PROBE").orf.nunique().max() != 1:
        raise ValueError("ambiguous barcode-to-strain mapping")
    grouped = pmap.groupby("PROBE", sort=True)
    if not grouped.size().eq(5).all():
        raise ValueError("TAG4 map must provide exactly five cells per barcode")
    tags = grouped.first()
    cells = np.array(grouped.cell.apply(list).tolist(), dtype=int)
    meta = pd.read_csv(data_root / "raw/E-MTAB-2391.sdrf.txt", sep="\t", dtype=str).fillna("")
    meta = meta.apply(lambda s: s.str.strip()).sort_values("Array Data File").reset_index(drop=True)
    if meta["Array Data File"].duplicated().any():
        raise ValueError("duplicate CEL metadata")
    names = meta["Array Data File"].to_numpy()
    vehicle = meta["Factor Value[compound]"].str.upper().eq("DMSO").to_numpy()
    if not np.array_equal(vehicle, meta["Comment[inchikey]"].eq("")):
        raise ValueError("blank compound identifiers do not exactly identify DMSO")
    paths = [data_root / "cel_all" / name for name in names]
    if not all(p.is_file() for p in paths):
        raise ValueError("missing raw CEL")
    cache_dir = Path(args.cache_dir) if args.cache_dir else out
    cache = cache_dir / "barcode_intensities.npy"
    map_sha = hashlib.sha256(Path(args.barcode_map).read_bytes()).hexdigest()
    contract = {"map_sha256": map_sha, "arrays": names.tolist(),
                "summary": "five spatial replicates, trim one at each end, arithmetic mean"}
    if cache.exists():
        if json.loads((cache_dir / "barcode_cache.json").read_text()) != contract:
            raise ValueError("cached barcode intensities have a different input contract")
        intensity = np.load(cache)
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            vectors = []
            for i, vec in enumerate(pool.map(barcode_array, [(p, cells) for p in paths], chunksize=8), 1):
                vectors.append(vec)
                if i % 100 == 0:
                    print(f"barcode CEL {i}/{len(paths)}", flush=True)
        intensity = np.stack(vectors, axis=1)
        np.save(cache, intensity)
        (cache_dir / "barcode_cache.json").write_text(json.dumps(contract, indent=2))
    # Alternate barcodes for one deletion are resolved solely by vehicle signal,
    # before any target profile, drug annotation or retrieval result is read.
    tags["vehicle_median_intensity"] = np.median(intensity[:, vehicle], axis=1)
    chosen = tags.reset_index().sort_values(["vehicle_median_intensity", "PROBE"], ascending=[False, True])
    chosen = chosen.drop_duplicates(["orf", "tag"]).sort_values(["orf", "tag"])
    strain_universe = pd.read_csv(args.strain_universe, sep="\t")
    if strain_universe.orf.duplicated().any():
        raise ValueError("duplicate published HIP/HOP strain identifiers")
    chosen = chosen[chosen.orf.isin(strain_universe.orf)].copy()
    chosen["assay"] = chosen.orf.map(strain_universe.set_index("orf").assay)
    ix = tags.index.get_indexer(chosen.PROBE)
    expression = np.log2(intensity[ix].astype(float))
    # Normalize uptags and downtags separately, as specified for barcode assays.
    group = chosen.tag if not args.assay_aware else chosen.assay + "_" + chosen.tag
    for kind in sorted(group.unique()):
        take = np.flatnonzero(group.to_numpy() == kind)
        block = expression[take]
        target = np.sort(block, axis=0).mean(axis=1)
        for j in range(block.shape[1]):
            ranks = rankdata(block[:, j], method="average") - 1
            expression[take, j] = np.interp(ranks, np.arange(len(target)), target)
    dates = np.array(["_".join(n.split("_")[:3]) for n in names])
    global_reference = np.median(expression[:, vehicle], axis=1)
    response = np.empty_like(expression)
    sample_records = []
    for date in sorted(set(dates)):
        matched = vehicle & (dates == date)
        reference = np.median(expression[:, matched], axis=1) if matched.any() else global_reference
        cols = np.flatnonzero(dates == date)
        response[:, cols] = expression[:, cols] - reference[:, None]
        for j in cols:
            sample_records.append({"array_name": names[j], "date": date,
                                   "n_same_date_vehicles": int(matched.sum()),
                                   "reference": "same-date DMSO" if matched.any() else "global DMSO fallback"})
    strains = sorted(set(chosen.orf))
    if args.tag_audit_only:
        # Keep independent barcode measurements available without altering the
        # previously frozen mean-tag response or selecting on screening results.
        tag_arrays = {}
        for kind in ["up", "down"]:
            take = chosen.tag.to_numpy() == kind
            index = chosen.loc[take, "orf"].to_numpy()
            tag_arrays[kind] = pd.DataFrame(response[take], index=index).reindex(strains).to_numpy(dtype=np.float32)
            tag_arrays[kind + "_vehicle_intensity"] = chosen.loc[take].set_index("orf").vehicle_median_intensity.reindex(strains).to_numpy()
        np.savez_compressed(out / "tag_response.npz", strain_orfs=np.array(strains),
                            array_names=names.astype(str), is_vehicle=vehicle, **tag_arrays)
        print("Exported separate barcode measurements; frozen mean-tag response unchanged", flush=True)
        return
    if args.stable_tag:
        control = intensity[ix][:, vehicle].astype(float)
        center = np.median(control, axis=1)
        chosen["control_robust_cv"] = np.median(np.abs(control-center[:, None]), axis=1)/center
        best = chosen.reset_index(drop=True).sort_values(["control_robust_cv", "PROBE"], kind="stable")
        best_rows = best.drop_duplicates("orf").set_index("orf").reindex(strains)
        rows = chosen.reset_index(drop=True).reset_index().set_index("PROBE").loc[best_rows.PROBE, "index"].to_numpy()
        values = response[rows]
        chosen["selected_for_response"] = chosen.PROBE.isin(best_rows.PROBE)
    else:
        values = np.stack([response[chosen.orf.to_numpy() == g].mean(axis=0) for g in strains])
    if args.assay_aware:
        # HIP/HOP were grown as different pools. Standardize their contrasts
        # separately before combining strains (Barazandeh et al., 2022).
        assay = strain_universe.set_index("orf").assay.reindex(strains).to_numpy()
        for kind in sorted(set(assay)):
            take = np.flatnonzero(assay == kind)
            block = values[take]
            center = np.median(block, axis=0)
            scale = 1.4826*np.median(np.abs(block-center), axis=0)
            # A singleton same-day control is exactly zero by construction.
            if np.any((scale <= 0) & ~vehicle):raise ValueError("degenerate drug assay scale")
            values[take] = (block-center)/np.where(scale > 0, scale, 1)
    chosen.to_csv(out / "selected_barcode_probes.tsv", sep="\t", index=False)
    pd.DataFrame(sample_records).to_csv(out / "array_reference_audit.tsv", sep="\t", index=False)
    method = ("TAG4 barcode (not ORF-control) probes; five-replicate trimmed means; "
              "log2 quantile normalization separately for uptags/downtags; same-date "
              "DMSO median contrast with disclosed global fallback; mean of available tags; "
              "published HIP/HOP strain universe")
    if args.assay_aware or args.stable_tag:
        if args.stable_tag:method = method.replace("mean of available tags", "tag with minimum control robust CV")
        if args.assay_aware:method = method.replace("separately for uptags/downtags", "separately for HIP/HOP and uptags/downtags; per-assay robust scaling")
        method += ("; diagnostic published-protocol factors: assay_aware=" + str(args.assay_aware)
                   + ", stable_tag=" + str(args.stable_tag)
                   + "; quantile normalization retained, not a reproduction of median-polish batch correction")
    np.savez_compressed(out / "strain_response.npz", strain_orfs=np.array(strains),
                        array_names=names.astype(str), compound_inchikeys=meta["Comment[inchikey]"].to_numpy(dtype=str),
                        doses=meta["Factor Value[dose]"].to_numpy(dtype=str),
                        dose_units=meta["Unit[concentration unit]"].to_numpy(dtype=str),
                        is_vehicle=vehicle, z_score=values.astype(np.float32), method=np.array([method]))
    summary = {"method": method, "n_arrays": len(names), "n_vehicle": int(vehicle.sum()),
               "n_compounds": int(meta.loc[~vehicle, "Comment[inchikey]"].nunique()),
               "n_strains": len(strains), "n_tag_sequences": len(tags), "n_selected_tags": len(chosen),
               "n_arrays_without_same_date_vehicle": sum(r["n_same_date_vehicles"] == 0 for r in sample_records),
               "response_sign": "positive = relative barcode abundance; negative = depletion",
               "scale": ("per-assay median/MAD-standardized barcode contrast; not author-published FD" if args.assay_aware
                         else "log2 barcode contrast, not a z statistic or author-published FD"),
               "assay_aware": args.assay_aware, "stable_tag": args.stable_tag,
               "barcode_map_sha256": map_sha, "no_target_or_drug_labels_used": True}
    summary["strain_universe_path"] = args.strain_universe
    summary["strain_universe_sha256"] = hashlib.sha256(Path(args.strain_universe).read_bytes()).hexdigest()
    summary["n_published_strains"] = len(strain_universe)
    summary["published_strains_without_barcode_map"] = sorted(set(strain_universe.orf) - set(strains))
    (out / "response_rebuild_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2), flush=True)


def probeset_log2_matrix(cel_paths, pmap: pd.DataFrame) -> pd.DataFrame:
    """对多个 CEL 计算 probeset 级 log2 强度矩阵（列=样本）。"""
    rows = {}
    groups = pmap.groupby("probeset")["cell"].apply(list)
    cell_index = {c: i for i, c in enumerate(sorted(set(pmap["cell"])))}
    for p in cel_paths:
        intens = parse_cel_v3(Path(p))
        vec = {}
        for ps, cells in groups.items():
            vals = [intens[c] for c in cells if c < len(intens)
                    and np.isfinite(intens[c])]
            if vals:
                vec[ps] = np.log2(np.mean(vals))
        rows[Path(p).stem] = vec
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--build-response", action="store_true")
    ap.add_argument("--data-root", type=Path, default=ROOT,
                    help="Directory containing cel_all/ and raw/E-MTAB-2391.sdrf.txt")
    ap.add_argument("--tag-audit-only", action="store_true")
    ap.add_argument("--barcode-map")
    ap.add_argument("--strain-universe")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--cache-dir", help="Reuse the input-validated raw barcode cache")
    ap.add_argument("--assay-aware", action="store_true")
    ap.add_argument("--stable-tag", action="store_true")
    ap.add_argument("--treated-cel")
    ap.add_argument("--dmso-cels", nargs="+")
    ap.add_argument("--author-fd", help="作者 FD 文件（金标准）")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    if args.build_response:
        if not args.barcode_map or not args.strain_universe or args.workers < 1:
            ap.error("--build-response requires --barcode-map, --strain-universe and positive --workers")
        build_barcode_response(args)
        return
    if not all([args.treated_cel, args.dmso_cels, args.author_fd]):
        ap.error("single-screen comparison requires treated CEL, DMSO CELs and author FD")

    pmap = load_probeset_map()
    print(f"probeset map: {pmap['probeset'].nunique()} probesets", flush=True)

    mat = probeset_log2_matrix([args.treated_cel] + list(args.dmso_cels), pmap)
    print(f"matrix: {mat.shape}", flush=True)

    treated_col = Path(args.treated_cel).stem
    dmso_cols = [Path(c).stem for c in args.dmso_cels]
    dmso_med = mat[dmso_cols].median(axis=1)
    log2fc = mat[treated_col] - dmso_med

    # z 变换（跨菌株）
    fd = (log2fc - log2fc.mean()) / log2fc.std(ddof=0)
    recalc = pd.DataFrame({"probeset": fd.index, "recalc_fd": fd.values})
    recalc["orf"] = recalc["probeset"].str.split(".").str[0]

    author = pd.read_csv(args.author_fd, sep="\t")
    merged = author.merge(recalc[["orf", "recalc_fd"]], left_on="ORF",
                          right_on="orf", how="inner")
    n = len(merged)
    pearson = np.corrcoef(merged["FD"], merged["recalc_fd"])[0, 1]
    from scipy.stats import spearmanr
    spear = spearmanr(merged["FD"], merged["recalc_fd"]).statistic
    top_author = set(merged.nlargest(50, "FD")["orf"])
    top_recalc = set(merged.nlargest(50, "recalc_fd"))["orf"] if False else set(
        merged.nlargest(50, "recalc_fd")["orf"])
    overlap = len(top_author & top_recalc)

    summary = {
        "n_strains_matched": int(n),
        "pearson": round(float(pearson), 4),
        "spearman": round(float(spear), 4),
        "top50_overlap": overlap,
        "author_top5": merged.nlargest(5, "FD")[["ORF", "GENE", "FD"]].to_dict("records"),
        "recalc_top5": merged.nlargest(5, "recalc_fd")[["ORF", "GENE", "recalc_fd"]].to_dict("records"),
    }
    import json
    Path(args.out).write_text(json.dumps(summary, indent=1, default=str))
    print(json.dumps(summary, indent=1, default=str))


if __name__ == "__main__":
    sys.exit(main())
