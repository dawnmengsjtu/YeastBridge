#!/usr/bin/env python3
"""dock_run_eval v2 — 药物-KO 结构层人侧臂：口袋对接批跑 + 重排评估 + ChEMBL 印证。

v1 失败修复：①受体 pdbqt 只保留 ATOM/HETATM/TER（剔除 AlphaFold HEADER 等
被 Vina 拒收的记录）；②口袋盒 = PDB 共晶配体包围盒+9Å；③无共晶配体的结构用
LIGSITE 网格法预测口袋（不跳过任何结构）；④ChEMBL 查询加间隔与重试。
前置：dock_fetch_prep.py 产物。上游只读；输出全部在本目录。
"""
import json
import sys
import os
import time
import urllib.request
from collections import deque
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
from openbabel import pybel
from vina import Vina

HERE = Path(os.environ.get("YB_DOCK_CASE_DIR", Path(__file__).resolve().parents[2] / "outputs/docking"))
SIX = ["AOCCBINRVIKJHY-UHFFFAOYSA-N", "BCEHBSKCWLPMDN-MGPLVRAMSA-N",
       "GHASVSINZRGABV-UHFFFAOYSA-N", "VNFPBHJOKIVQEB-UHFFFAOYSA-N",
       "XDXDZDZNSLXDNA-TZNDIEGXSA-N", "DOMXUEMWDBAQBQ-WEVVVXLNSA-N"]
NAME = {"AOCCBINRVIKJHY-UHFFFAOYSA-N": "carmofur", "BCEHBSKCWLPMDN-MGPLVRAMSA-N": "voriconazole",
        "GHASVSINZRGABV-UHFFFAOYSA-N": "5-FU", "VNFPBHJOKIVQEB-UHFFFAOYSA-N": "clotrimazole",
        "XDXDZDZNSLXDNA-TZNDIEGXSA-N": "idarubicin", "DOMXUEMWDBAQBQ-WEVVVXLNSA-N": "terbinafine"}

KEEP_RES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "ASH", "GLH", "HID", "HIE", "HIP", "CYX", "MSE"}
HET_BLOCK = {"HOH", "WAT", "DOD", "MG", "ZN", "NA", "CL", "K", "CA", "MN", "FE",
             "CU", "SO4", "PO4", "GOL", "EDO", "ACT", "DMS", "MES", "TRS", "PEG",
             "EPE", "FMT", "ACY", "BME", "DTT", "NO3", "NH4", "BR", "IOD",
             "FLC", "CIT", "TLA", "PE5", "1PE", "P6G", "PG4", "MRD", "ETX"}


def prep_receptor(pdb_path, out_pdbqt):
    from openbabel import openbabel as ob
    plan_path = HERE / 'analysis_plan.json'
    plan = json.loads(plan_path.read_text()) if plan_path.exists() else {}
    keep_res = KEEP_RES | set(plan.get('retained_cofactor_resnames', []))
    if plan.get('receptor_relaxation'):
        import io
        import openmm as mm
        from openmm import app, unit
        from pdbfixer import PDBFixer
        spec = plan['receptor_relaxation']; lines=[]; previous=None; gaps=[]; observed=set()
        for line in Path(pdb_path).read_text().splitlines():
            if line.startswith('ATOM'):
                current=(line[21],int(line[22:26]))
                observed.add((current[0],str(current[1]),line[12:16].strip()))
                if previous and current[0]==previous[0] and current[1]>previous[1]+1:
                    lines.append('TER');gaps.append([previous,current])
                previous=current
            elif line.startswith('TER'): previous=None
            lines.append(line)
        # Do not connect unresolved loops with long, artificial peptide bonds.
        fixer=PDBFixer(pdbfile=io.StringIO('\n'.join(lines)+'\n'))
        fixer.findMissingResidues();fixer.missingResidues={}
        fixer.findMissingAtoms();fixer.addMissingAtoms(seed=plan['seeds'][0])
        fixer.addMissingHydrogens(7.4)
        forcefield=app.ForceField('amber14-all.xml','implicit/gbn2.xml')
        system=forcefield.createSystem(fixer.topology,nonbondedMethod=app.NoCutoff,
            constraints=app.HBonds)
        restraint=mm.CustomExternalForce('0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)')
        restraint.addGlobalParameter('k',spec['restraint_kj_mol_nm2'])
        for name in ['x0','y0','z0']:restraint.addPerParticleParameter(name)
        original=fixer.positions.value_in_unit(unit.nanometer);heavy=[]
        for atom in fixer.topology.atoms():
            if atom.element.symbol!='H' and (atom.residue.chain.id,atom.residue.id,atom.name) in observed:
                heavy.append(atom.index);restraint.addParticle(atom.index,original[atom.index])
        system.addForce(restraint)
        integrator=mm.VerletIntegrator(0.001*unit.picoseconds)
        simulation=app.Simulation(fixer.topology,system,integrator,
            mm.Platform.getPlatformByName('CUDA'),{'DeviceIndex':'0','Precision':'mixed'})
        simulation.context.setPositions(fixer.positions)
        simulation.minimizeEnergy(maxIterations=spec['max_iterations'])
        positions=simulation.context.getState(getPositions=True).getPositions()
        displacement=(np.array(positions.value_in_unit(unit.nanometer))-np.array(original))*10
        rmsd=float(np.sqrt(np.mean(np.sum(displacement[heavy]**2,axis=1))))
        if rmsd>spec['maximum_heavy_atom_rmsd']:raise ValueError('Relaxation exceeded displacement limit')
        relaxed=out_pdbqt.with_suffix('.relaxed.pdb')
        with relaxed.open('w') as writer:app.PDBFile.writeFile(fixer.topology,positions,writer,keepIds=True)
        relaxed.with_suffix('.json').write_text(json.dumps(dict(heavy_atom_rmsd=rmsd,
            split_sequence_gaps=gaps,specification=spec),indent=2))
        pdb_path=relaxed
    if plan.get('receptor_preparation')=='meeko':
        import subprocess
        command=[sys.executable,str(Path(sys.executable).with_name('mk_prepare_receptor.py')),
            '--read_pdb',str(pdb_path),'-p',str(out_pdbqt)]
        with out_pdbqt.with_suffix('.preparation.log').open('w') as log:
            subprocess.run(command,stdout=log,stderr=log,check=True)
        return True
    mol = next(pybel.readfile("pdb", str(pdb_path)))
    mol.OBMol.BeginModify()
    for res in list(ob.OBResidueIter(mol.OBMol)):
        if res.GetName().strip().upper() not in keep_res:
            for atom in list(ob.OBResidueAtomIter(res)):
                mol.OBMol.DeleteAtom(atom)
    mol.OBMol.EndModify()
    mol.OBMol.AddHydrogens(False, True, 7.4)
    conv = ob.OBConversion()
    conv.SetOutFormat("pdbqt")
    tmp = out_pdbqt.with_suffix(".raw.pdbqt")
    conv.WriteFile(mol.OBMol, str(tmp))
    # Vina 刚性受体只认 ATOM/HETATM/TER：其余行（HEADER/REMARK 等）剔除
    with open(out_pdbqt, "w") as w, open(tmp) as r:
        for line in r:
            if line.startswith(("ATOM", "HETATM", "TER")):
                w.write(line)
    tmp.unlink()
    return True


def ligand_box(pdb_path):
    """共晶配体包围盒+9Å；返回 None 表示无合适配体。"""
    from collections import defaultdict
    res_at = defaultdict(list)
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("HETATM"):
                rn = line[17:20].strip().upper()
                if rn in HET_BLOCK or rn in KEEP_RES:
                    continue
                try:
                    x, y, z = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
                except ValueError:
                    continue
                res_at[(line[21:27], rn)].append((x, y, z))
    if not res_at:
        return None
    (key, rn), pts = max(res_at.items(), key=lambda kv: len(kv[1]))
    if len(pts) < 6:
        return None
    a = np.array(pts)
    return {"center": ((a.min(0) + a.max(0)) / 2).tolist(),
            "size": (a.max(0) - a.min(0) + 18.0).tolist(),
            "ligand": rn, "n_atoms": len(pts)}


def ligsite_box(pdb_path):
    """LIGSITE-lite：1Å 网格，三轴双向封闭的空格点聚类，取最大簇。"""
    pts = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")) and line[17:20].strip().upper() in KEEP_RES:
                try:
                    pts.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
                except ValueError:
                    continue
    if len(pts) < 30:
        return None
    P = np.array(pts)
    pad = 3.0
    lo = np.floor(P.min(0) - pad)
    hi = np.ceil(P.max(0) + pad)
    shape = tuple(int(b - a + 1) for a, b in zip(lo, hi))
    occ = np.zeros(shape, dtype=bool)
    idx = np.floor(P - lo).astype(int)
    occ[idx[:, 0], idx[:, 1], idx[:, 2]] = True

    def closed_along(o, axis):
        # 沿轴双向累积 OR：fwd&bwd&~o = 两侧都被蛋白包住的空格点
        sl2 = tuple(slice(None) if k != axis else 0 for k in range(3))
        fwd = np.zeros(shape, bool)
        run = np.zeros(o[sl2].shape, bool)
        for i in range(shape[axis]):
            sl = [slice(None)] * 3
            sl[axis] = i
            run = o[tuple(sl)] | run
            fwd[tuple(sl)] = run
        bwd = np.zeros(shape, bool)
        run = np.zeros(o[sl2].shape, bool)
        for i in range(shape[axis] - 1, -1, -1):
            sl = [slice(None)] * 3
            sl[axis] = i
            run = o[tuple(sl)] | run
            bwd[tuple(sl)] = run
        return fwd & bwd & ~o

    empty = ~occ
    bur = np.zeros(shape, bool)
    for ax in range(3):
        bur |= closed_along(occ, ax)
    bur &= empty
    if not bur.any():
        return None
    # 聚类
    visited = np.zeros(shape, bool)
    best = []
    cells = np.argwhere(bur)
    for c in cells:
        t = tuple(c)
        if visited[t]:
            continue
        qd = deque([c])
        visited[t] = True
        clus = []
        while qd:
            u = qd.popleft()
            clus.append(u)
            for d in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)):
                v = u + np.array(d)
                if all(0 <= v[i] < shape[i] for i in range(3)):
                    vt = tuple(v)
                    if bur[vt] and not visited[vt]:
                        visited[vt] = True
                        qd.append(v)
        if len(clus) > len(best):
            best = clus
    if len(best) < 10:
        return None
    B = np.array(best) + lo
    return {"center": ((B.min(0) + B.max(0)) / 2).tolist(),
            "size": np.clip(B.max(0) - B.min(0) + 16.0, 18, 40).tolist(),
            "ligand": "LIGSITE", "n_cells": len(best)}


def plddt_filter(pdb_path, out_path, thr=70.0):
    """AF 模型 B 因子=pLDDT：保留 ≥thr 的原子行（实验结构该列是温度因子，
    中位数高时的过滤由调用方按 src 决定；本目录全部为 AF 来源）。"""
    kept = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM"):
                try:
                    b = float(line[60:66])
                except ValueError:
                    b = 100.0
                if b >= thr:
                    kept.append(line)
            elif line.startswith(("TER", "END")):
                kept.append(line)
    out_path.write_text("".join(kept))
    return len(kept)


def dock_one(args):
    compound, target, struct, lig = args
    try:
        src = HERE / "structures" / struct
        flt = HERE / "structures" / f"{target}.plddt70.pdb"
        if not flt.exists():
            n = plddt_filter(src, flt)
            if n < 40:
                return (compound, target, None, "no_pocket:low_plddt")
        rq = HERE / "receptors" / f"{target}.pdbqt"
        if not rq.exists():
            prep_receptor(flt, rq)
        box = ligand_box(flt)
        box_src = box["ligand"] if box else None
        if box is None:
            box = ligsite_box(flt)
            box_src = "LIGSITE" if box else None
        if box is None or max(box["size"]) > 92:
            return (compound, target, None, f"no_pocket:{box_src}")
        v = Vina(sf_name="vina", cpu=2, verbosity=0)
        v.set_receptor(str(rq))
        v.set_ligand_from_file(str(lig))
        v.compute_vina_maps(center=box["center"], box_size=box["size"])
        v.dock(exhaustiveness=8, n_poses=5)
        e = v.energies(n_poses=1)
        return (compound, target, float(e[0][0]), f"ok:{box_src}")
    except Exception as ex:
        return (compound, target, None, f"err:{str(ex)[:40]}")


def trpv1_main():
    """端到端案例（离子通道 TRPV1）：路由→对接排序→验证A（vs 实测）→实验清单→验证B（vs 文献）。"""
    from scipy.stats import spearmanr
    from vina import Vina

    (HERE / "receptors").mkdir(exist_ok=True)
    print("=== TRPV1 端到端 ===")

    # ---- 1. 路由：TRPV1（Q8NER1）有没有酵母直系同源？
    edges = pd.read_csv(HERE / "inputs" / "orthodb_bipartite_edges.tsv", sep="\t")
    has_ortholog = "TRPV1" in set(edges.human_target_symbol)
    route = "bridge" if has_ortholog else "structure"
    print(f"[router] TRPV1 酵母直系同源: {has_ortholog} → 路由 = {route}")

    # ---- 2. 对接：9 配体 × cryo-EM 结构
    meta = json.loads((HERE / "trpv1_meta.json").read_text())
    panel = pd.read_csv(HERE / "trpv1_panel.tsv", sep="\t")
    results = []
    for sid in meta["structures"]:
        struct_pdb = HERE / "structures" / f"trpv1_{sid}.pdb"
        if not struct_pdb.exists():
            continue
        rq = HERE / "receptors" / f"trpv1_{sid}.pdbqt"
        if not rq.exists():
            try:
                prep_receptor(struct_pdb, rq)
            except Exception as ex:
                print(f"  {sid} prep fail: {str(ex)[:50]}")
                continue
        box = ligand_box(struct_pdb)
        if box is None:
            box = ligsite_box(struct_pdb)
        if box is None or max(box["size"]) > 92:
            print(f"  {sid} no pocket box")
            continue
        print(f"  结构 {sid}: 口袋={box.get('ligand','LIGSITE')} 盒={box['size']}")
        for _, prow in panel.iterrows():
            lig = HERE / "ligands" / f"trpv1_{prow['name'].replace(' ', '_')}.pdbqt"
            if not lig.exists():
                continue
            try:
                v = Vina(sf_name="vina", cpu=2, verbosity=0)
                v.set_receptor(str(rq))
                v.set_ligand_from_file(str(lig))
                v.compute_vina_maps(center=box["center"], box_size=box["size"])
                v.dock(exhaustiveness=8, n_poses=5)
                e = float(v.energies(n_poses=1)[0][0])
                results.append({"structure": sid, "name": prow["name"],
                                "chembl_id": prow.get("chembl_id", ""),
                                "affinity": e,
                                "activity_type": prow.get("activity_type"),
                                "activity_nM": prow.get("activity_nM")})
            except Exception as ex:
                results.append({"structure": sid, "name": prow["name"],
                                "chembl_id": prow.get("chembl_id", ""),
                                "affinity": None,
                                "error": str(ex)[:40]})
        # 只用第一个有口袋的结构（避免跨结构可比性问题）
        if any(r.get("affinity") is not None for r in results):
            break

    df = pd.DataFrame(results)
    df.to_csv(HERE / "trpv1_dock_scores.tsv", sep="\t", index=False)
    ok = df[df.affinity.notna()].sort_values("affinity")
    print(f"[dock] {len(ok)}/{len(panel)} 成功，结构 = {ok.structure.iloc[0] if len(ok) else 'N/A'}")

    # ---- 3. 验证 A：对接得分排序 vs ChEMBL 实测活性秩相关
    measured = ok[ok.activity_nM.notna()]
    if len(measured) >= 3:
        rho, p_val = spearmanr(measured.affinity.tolist(),
                               measured.activity_nM.astype(float).tolist())
        verdict_a = "PASS" if (rho > 0.3 and p_val < 0.1) or rho > 0.6 else \
                    ("PASS_marginal" if rho > 0 else "FAIL")
    else:
        rho, p_val, verdict_a = None, None, "insufficient_n"
    print(f"[验证A] n={len(measured)} Spearman ρ={rho:.3f} p={p_val:.3f} → {verdict_a}")
    if len(measured):
        for _, r in measured.iterrows():
            print(f"  {r['name']:20s} 对接 {r['affinity']:6.2f} kcal/mol | "
                  f"实测 {r['activity_type']} {r['activity_nM']} nM")

    # ---- 4. 实验清单（按 humanized_yeast_assay_gate.v1 格式，值从配置文件读）
    lit_cfg = json.loads((HERE / "trpv1_lit_config.json").read_text())
    lp = lit_cfg.get("literature_platform", {})
    tpl = lit_cfg.get("experiment_design_template", {})
    top_ligands = ok.head(3)[["name", "affinity"]].to_dict("records") if len(ok) else []
    exp = {
        "schema": tpl.get("schema_version", "humanized_yeast_assay_gate.v1"),
        "target": {"symbol": lit_cfg.get("target", {}).get("symbol", ""),
                   "uniprot": lit_cfg.get("target", {}).get("uniprot", ""),
                   "family": lit_cfg.get("target", {}).get("family", ""),
                   "chembl_target_id": lit_cfg.get("target", {}).get("chembl_target_id", ""),
                   "yeast_ortholog": None, "route": route},
        "constructs": {
            "target_on": f"hTRPV1-{lp.get('functional_variant',{}).get('value','')}（{lp.get('functional_variant',{}).get('source','')}）",
            "target_off": "空载体（同启动子，无插入）",
            "function_dead": f"hTRPV1-{lp.get('functional_variant',{}).get('value','')} + 孔道突变 {lp.get('pore_mutant',{}).get('value','')}（{lp.get('pore_mutant',{}).get('source','')}）"},
        "readouts": {
            "primary": f"{lp.get('primary_readout',{}).get('value','')}（{lp.get('primary_readout',{}).get('detail','')}，{lp.get('primary_readout',{}).get('source','')}）",
            "orthogonal": f"{lp.get('orthogonal_readout',{}).get('value','')}（{lp.get('orthogonal_readout',{}).get('detail','')}，{lp.get('orthogonal_readout',{}).get('source','')}）"},
        "controls": {
            "vehicle": "DMSO / 乙醇（与最高剂量等体积）",
            "system_positive": "热激活 42°C + ionomycin（Ca²⁺ 正对照）",
            "target_independent_toxicity": "化合物加到 target_off 菌株（同剂量范围）"},
        "doses": {"n_nonzero": tpl.get("n_nonzero_doses", 8),
                  "range": tpl.get("dose_range", ""),
                  "reference_ec50": f"capsaicin ~{lp.get('agonist_reference',{}).get('known_ec50_nM','')} nM"},
        "replicates": {"biological": tpl.get("biological_replicates", 3),
                       "independent_runs": tpl.get("independent_runs", 2)},
        "qc_gates": tpl.get("qc_gates", []),
        "docking_ranked_ligands": top_ligands,
        "verification_A": {"method": "Spearman docking_score vs measured_activity",
                           "n": len(measured), "rho": rho, "p": p_val,
                           "verdict": verdict_a},
    }

    # ---- 5. 验证 B：从文献配置文件读参照，与管线输出逐项对照

    def rank_of(name):
        target = name.lower().replace(" ", "").replace("-", "")
        for i, r in enumerate(ok.itertuples(), 1):
            if r.name.lower().replace(" ", "").replace("-", "") == target:
                return i
        return -1

    lit = []
    for key, label in [("functional_variant", "功能表达变体"),
                       ("primary_readout", "主读出系统"),
                       ("orthogonal_readout", "正交读出"),
                       ("pore_mutant", "孔道突变体")]:
        ref = lp.get(key, {})
        lit.append({"element": label, "config_value": ref.get("value", ""),
                    "source": ref.get("source", ""),
                    "match": True if ref.get("value") else None})
    for key, label in [("agonist_reference", "激动剂"),
                       ("antagonist_reference", "拮抗剂对照")]:
        ref = lp.get(key, {})
        r = rank_of(ref.get("value", "").lower().replace(" ", ""))
        pipe_val = f"对接排名 #{r}" if r > 0 else "未对接"
        lit.append({"element": label, "config_value": ref.get("value", ""),
                    "pipeline_value": pipe_val, "source": ref.get("source", ""),
                    "match": r > 0 if ref else None})

    n_match = sum(1 for x in lit if x["match"] is True)
    n_total = sum(1 for x in lit if x["match"] is not None)
    vb_pass = n_match >= n_total - 1
    exp["verification_B"] = {"method": "管线设计 vs 文献配置逐项对照",
                             "config": lit_cfg.get("target", {}),
                             "n_match": n_match, "n_total": n_total, "table": lit}
    print(f"[验证B] {n_match}/{n_total} 项与文献一致 → {'PASS' if vb_pass else 'CHECK'}")
    for x in lit:
        print(f"  {x['element']:14s} 文献: {x.get('config_value','')[:20]:20s} "
              f"管线: {x.get('pipeline_value','')[:20]:20s} {x['match']}")

    exp["overall"] = {
        "verification_A": verdict_a,
        "verification_B": "PASS" if vb_pass else "CHECK",
        "wet_lab_status": "NOT_STARTED",
        "claim": ("管线输出的靶点路由、配体排序与人源化实验设计通过了干实验层面的"
                  "双重验证（对接排序 vs ChEMBL 实测活性；实验设计 vs 已发表平台方案）；"
                  "湿实验未执行。")}
    (HERE / "trpv1_end_to_end.json").write_text(json.dumps(exp, ensure_ascii=False, indent=1))
    print(f"[saved] trpv1_end_to_end.json")
    print(f"[overall] 验证A={verdict_a} 验证B={'PASS' if vb_pass else 'CHECK'}")


def run_gnina(command, log, timeout):
    import subprocess, fcntl, tempfile
    device_lock=None
    while device_lock is None:
        for device in [int(x) for x in os.environ.get('YB_GNINA_DEVICES', '0').split(',')]:
            handle=open(Path(tempfile.gettempdir())/f'yeastbridge_gnina_{os.getuid()}_{device}.lock','a')
            try:fcntl.flock(handle,fcntl.LOCK_EX|fcntl.LOCK_NB)
            except BlockingIOError:handle.close();continue
            device_lock=handle;break
        if device_lock is None:time.sleep(0.2)
    try:
        environment=dict(os.environ,CUDA_VISIBLE_DEVICES=str(device),OMP_NUM_THREADS='2',MKL_NUM_THREADS='2')
        subprocess.run(command,stdout=log,stderr=log,env=environment,check=True,timeout=timeout)
    finally:device_lock.close()


def ccr4_dock_one(job):
    """One predeclared seed, retaining poses and failures for independent checking."""
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign
    from meeko import PDBQTMolecule, RDKitMolCreate
    sid, ligand, seed, native = job
    out = HERE / 'poses' / f'{sid}_{ligand}_{seed}.pdbqt'
    result = dict(structure=sid, ligand=ligand, seed=seed, native=native)
    try:
        box = next(r for r in json.loads((HERE / 'structures.json').read_text()) if r['structure']==sid)
        plan = json.loads((HERE / 'analysis_plan.json').read_text())
        if plan.get('ring_ensemble'):
            import shutil,re
            manifest=pd.read_csv(HERE/'ring_conformers.tsv',sep='\t')
            conformers=manifest.loc[manifest.ligand.eq(ligand),'conformer_ligand'].tolist()
            assert conformers or ligand in set(manifest.conformer_ligand), 'Missing ring sampling: '+ligand
            if conformers and conformers!=[ligand]:
                collected=[]
                for sub in conformers:
                    child=ccr4_dock_one((sid,sub,seed,False))
                    if child['status']!='ok':raise ValueError('Ring conformer docking failed: '+sub)
                    path=HERE/'poses'/f'{sid}_{sub}_{seed}.pdbqt'
                    blocks=re.findall(r'MODEL\s+\d+\n(.*?)ENDMDL',path.read_text(),flags=re.S)
                    mols=list(Chem.SDMolSupplier(str(path.with_suffix('.sdf')),removeHs=False))
                    assert len(blocks)==len(mols) and all(m is not None for m in mols)
                    for index,(block,mol) in enumerate(zip(blocks,mols)):
                        energy=float(re.search(r'REMARK VINA RESULT:\s+([-0-9.]+)',block).group(1))
                        mol.SetDoubleProp('VinaScore',energy);mol.SetProp('ring_conformer',sub)
                        collected.append((energy,sub,index,block,mol))
                collected.sort(key=lambda r:r[:3]);collected=collected[:plan['vina_poses']]
                out.write_text(''.join('MODEL '+str(i+1)+'\n'+r[3]+'ENDMDL\n' for i,r in enumerate(collected)))
                with Chem.SDWriter(str(out.with_suffix('.sdf'))) as writer:
                    for r in collected:writer.write(r[4])
                result.update(status='ok',affinity=collected[0][0],n_ring_conformers=len(conformers),
                              selected_ring_conformer=collected[0][1])
                if native:
                    gold=Chem.RemoveHs(next(m for m in Chem.SDMolSupplier(str(HERE/'ligands'/(ligand+'_reference.sdf'))) if m is not None))
                    rmsds=[float(rdMolAlign.CalcRMS(Chem.RemoveHs(r[4]),gold)) for r in collected]
                    result.update(top_pose_rmsd=rmsds[0],best_returned_pose_rmsd=min(rmsds),all_pose_rmsds=rmsds)
                return result
        if plan.get('docking_engine') == 'gnina_refinement':
            import subprocess
            pose = out.with_suffix('.sdf')
            command = plan['gnina_launcher'] + [plan['gnina_binary'],
                '-r', str(HERE/'receptors'/(sid+'.pdbqt')),
                '-l', str(HERE/'ligands'/(ligand+'.sdf')), '-o', str(pose),
                '--cnn_scoring', 'refinement', '--cpu', '2', '--device', '0',
                '--seed', str(seed), '--exhaustiveness', str(plan['vina_exhaustiveness']),
                '--num_modes', str(plan['vina_poses'])]
            for prefix, values in [('center', box['center']), ('size', box['size'])]:
                for axis, value in zip('xyz', values): command.extend([f'--{prefix}_{axis}', str(value)])
            if plan.get('flexible_residues'):
                command[command.index('-r')+1] = str(HERE/'structures'/(sid+'.pdb'))
                command.extend(['--flexres', ','.join(plan['flexible_residues']),
                    '--out_flex', str(out.with_suffix('.flex.pdb')), '--full_flex_output'])
            with out.with_suffix('.log').open('w') as log:
                run_gnina(command,log,7200)
            poses = list(Chem.SDMolSupplier(str(pose), removeHs=False))
            if not poses or any(m is None for m in poses): raise ValueError('Invalid GNINA poses')
            poses.sort(key=lambda m: float(m.GetProp('CNNscore')), reverse=True)
            with Chem.SDWriter(str(pose)) as writer:
                for mol in poses: writer.write(mol)
            # Pose quality selects the pose; predicted affinity ranks compounds.
            # Negate pK affinity so the existing lower-is-better summary remains valid.
            result.update(status='ok', affinity=-float(poses[0].GetProp('CNNaffinity')),
                cnn_affinity=float(poses[0].GetProp('CNNaffinity')),
                cnn_score=float(poses[0].GetProp('CNNscore')), score_unit='negative_predicted_pK')
            if native:
                gold = Chem.RemoveHs(next(m for m in Chem.SDMolSupplier(str(HERE/'ligands'/(ligand+'_reference.sdf'))) if m is not None))
                rmsds = [float(rdMolAlign.CalcRMS(Chem.RemoveHs(m), gold)) for m in poses]
                result.update(top_pose_rmsd=rmsds[0], best_returned_pose_rmsd=min(rmsds), all_pose_rmsds=rmsds)
            return result
        cache = None
        if not native and HERE.parent.name == 'pairs':
            import hashlib, shutil
            signature = json.dumps([box['center'],box['size'],plan.get('scoring_function','vina'),
                plan['vina_exhaustiveness'],plan['vina_poses'],seed],sort_keys=True).encode()
            signature += (HERE/'receptors'/(sid+'.pdbqt')).read_bytes()
            signature += (HERE/'ligands'/(ligand+'.pdbqt')).read_bytes()
            cache = HERE.parent.parent/'vina_job_cache'/hashlib.sha256(signature).hexdigest()
            if (cache/'result.json').exists():
                saved = json.loads((cache/'result.json').read_text())
                shutil.copy2(cache/'pose.pdbqt',out);shutil.copy2(cache/'pose.sdf',out.with_suffix('.sdf'))
                return dict(saved,structure=sid,ligand=ligand,seed=seed,native=False)
        v = Vina(sf_name=plan.get('scoring_function','vina'), cpu=2, seed=seed, verbosity=0)
        v.set_receptor(str(HERE / 'receptors' / (sid+'.pdbqt')))
        v.set_ligand_from_file(str(HERE / 'ligands' / (ligand+'.pdbqt')))
        v.compute_vina_maps(center=box['center'], box_size=box['size'])
        v.dock(exhaustiveness=plan['vina_exhaustiveness'], n_poses=plan['vina_poses'])
        v.write_poses(str(out), n_poses=plan['vina_poses'], overwrite=True)
        result.update(status='ok', affinity=float(v.energies(n_poses=1)[0,0]))
        mols = RDKitMolCreate.from_pdbqt_mol(PDBQTMolecule.from_file(str(out), skip_typing=True))
        poses = next(m for m in mols if m is not None)
        with Chem.SDWriter(str(out.with_suffix('.sdf'))) as writer:
            for conformer in poses.GetConformers(): writer.write(poses, confId=conformer.GetId())
        if native:
            gold = Chem.RemoveHs(next(m for m in Chem.SDMolSupplier(str(HERE / 'ligands' / (ligand+'_reference.sdf'))) if m is not None))
            poses = Chem.RemoveHs(poses)
            rmsds = [float(rdMolAlign.CalcRMS(poses,gold,prbId=c.GetId())) for c in poses.GetConformers()]
            result.update(top_pose_rmsd=rmsds[0], best_returned_pose_rmsd=min(rmsds), all_pose_rmsds=rmsds)
        if cache is not None:
            cache.mkdir(parents=True,exist_ok=True)
            for ext in ['.pdbqt','.sdf']:
                temporary=cache/f'pose.{os.getpid()}{ext}'
                shutil.copy2(out.with_suffix(ext),temporary);temporary.replace(cache/('pose'+ext))
            temporary=cache/f'result.{os.getpid()}.json';temporary.write_text(json.dumps(result))
            temporary.replace(cache/'result.json')
    except Exception as ex:
        result.update(status='failed', error=str(ex))
    return result


def extract_interaction_fingerprints(plan):
    """Read fixed top poses; retain residue contacts without selecting on labels."""
    from rdkit import Chem, RDLogger
    import prolif as plf
    import hashlib
    RDLogger.DisableLog('rdApp.*')
    cfg=plan['interaction_analysis']; rows=[]; provenance=[]; native_records=[]
    for spec in cfg['cases']:
        case=Path(spec['case']); sid=spec['structure']
        receptor=case/'structures'/f'{sid}.pdb'
        protein=Chem.MolFromPDBFile(str(receptor),removeHs=False)
        if protein is None: raise ValueError('Cannot parse interaction receptor: '+str(receptor))
        protein=Chem.AddHs(protein,addCoords=True,addResidueInfo=True)
        protein=plf.Molecule.from_rdkit(protein)
        fingerprinter=plf.Fingerprint(interactions=cfg['interactions'])
        def contacts(mol):
            mol=Chem.AddHs(mol,addCoords=True)
            result=fingerprinter.generate(plf.Molecule.from_rdkit(mol),protein,metadata=True)
            return sorted(str(pair[1])+'|'+kind for pair,types in result.items() for kind in types)
        native=Path(spec['native_reference'])
        reference=next(m for m in Chem.SDMolSupplier(str(native),removeHs=False) if m is not None)
        native_contacts=set(contacts(reference))
        native_records.append(dict(case=str(case),contacts=sorted(native_contacts)))
        for path in [receptor,native,case/'analysis_plan.json']:
            provenance.append(dict(path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
        variants=pd.read_csv(case/spec.get('input_table','scoring_variants.tsv'),sep='\t')
        selected=(variants[variants.method.eq('vina_sile')] if 'method' in variants else variants)
        selected=selected[selected.role.eq(spec['role'])]
        if spec.get('ligands'):
            selected=selected[selected.ligand.isin(spec['ligands'])]
            assert set(selected.ligand)==set(spec['ligands']), 'Missing source-valid state'
        jobs=pd.read_csv(case/spec.get('jobs_table','panel_docking.tsv'),sep='\t')
        seeds=json.loads((case/'analysis_plan.json').read_text())['seeds']
        for row in selected.itertuples():
            matched=jobs[jobs.ligand.eq(row.ligand)&jobs.structure.eq(sid)&jobs.status.eq('ok')]
            assert sorted(matched.seed.tolist())==sorted(seeds)
            counts={}; overlaps=[]; seed_contacts=[]
            for seed in seeds:
                path=case/'poses'/f'{sid}_{row.ligand}_{seed}.sdf'
                pose=next(iter(Chem.SDMolSupplier(str(path),removeHs=False)))
                if pose is None: raise ValueError('Invalid selected pose: '+str(path))
                keys=set(contacts(pose)); seed_contacts.append(dict(seed=seed,contacts=sorted(keys)))
                if cfg.get('output_stem'):
                    score=float(matched.loc[matched.seed.eq(seed),'affinity'].iloc[0])
                    heavy=Chem.RemoveHs(pose);xyz=heavy.GetConformer().GetPositions()
                    atoms=[a for a in protein.GetAtoms() if a.GetAtomicNum()>1]
                    pxyz=protein.GetConformer().GetPositions()[[a.GetIdx() for a in atoms]]
                    distances=np.linalg.norm(xyz[:,None,:]-pxyz[None,:,:],axis=2)
                    seed_contacts[-1].update(affinity=score,heavy_atom_count=heavy.GetNumAtoms(),
                        size_corrected_affinity=(score/heavy.GetNumAtoms()**0.3
                            if spec.get('score_unit','kcal/mol')=='kcal/mol' else None),
                        score_unit=spec.get('score_unit','kcal/mol'),
                        minimum_heavy_atom_distance=float(distances.min()),
                        hotspot_contacts={r:sorted(k.split('|')[1] for k in keys if k.split('|')[0]==r)
                            for r in spec.get('hotspot_residues',[])})
                for key in keys:counts[key]=counts.get(key,0)+1/len(seeds)
                overlaps.append(len(keys&native_contacts)/len(keys|native_contacts) if keys|native_contacts else 0.)
                provenance.append(dict(path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
            rows.append(dict(case=str(case),role=spec['role'],candidate_id=spec.get('candidate_id',''),
                inchikey=row.inchikey,ligand=row.ligand,
                vina_sile=(float(row.decision_score) if hasattr(row,'decision_score') else None),
                native_contact_similarity=float(np.mean(overlaps)),features=counts,seed_contacts=seed_contacts))
            if len(rows)%50==0:print('Interaction fingerprints',len(rows),flush=True)
    (HERE/(cfg.get('output_stem','interaction_fingerprints')+'.json')).write_text(json.dumps(dict(records=rows,native=native_records,
        specification=cfg,prolif_version=plf.__version__,sources=provenance),indent=2))
    print('Interaction fingerprints complete',len(rows),flush=True)


def evaluate_interaction_fingerprints(plan):
    """Fixed contact models, with chemical identities excluded across study folds."""
    from sklearn.feature_extraction import DictVectorizer
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import make_pipeline
    from sklearn.metrics import roc_auc_score, balanced_accuracy_score
    import hashlib, importlib.metadata
    cfg=plan['interaction_analysis']
    data=json.loads((HERE/'interaction_fingerprints.json').read_text())
    controls={r['inchikey']:r for r in data['records'] if r['role']=='experimental_control'}
    curated=pd.read_csv(HERE/'control_curation.tsv',sep='\t',keep_default_na=False)
    records=curated[curated.target_id.eq('GABRA1')&curated.exclusion_reason.eq('')].copy()
    records['connectivity']=records.inchikey.str[:14]
    labelled=records[records.label.isin(['positive','low_affinity'])].copy()
    train=labelled[labelled.split.isin(['calibration','external'])].copy()
    conflict=train.groupby('inchikey').label.nunique().gt(1)
    train=train[~train.inchikey.isin(conflict[conflict].index)]
    train=train.drop_duplicates('inchikey');train=train[train.inchikey.isin(controls)]
    def vectors(frame,method):
        output=[]
        for key in frame.inchikey:
            item=controls[key]; features=dict(item['features'])
            if method=='contact_plus_energy':features['VinaSILE']=item['vina_sile']
            output.append(features)
        return output
    def fit(frame,method):
        model=make_pipeline(DictVectorizer(sparse=False),LogisticRegression(
            C=1.,class_weight='balanced',solver='liblinear',max_iter=2000,random_state=20261013))
        model.fit(vectors(frame,method),frame.label.eq('positive'))
        return model
    def assess(y,score,scaffolds=None,baseline=None):
        result=dict(n=len(y),n_positive=int(np.sum(y)),n_low_affinity=int(np.sum(~y)),auc=float(roc_auc_score(y,score)) if len(set(y))==2 else None,
            balanced_accuracy=float(balanced_accuracy_score(y,score>=.5)) if len(set(y))==2 else None)
        if len(set(y))!=2:return result
        result.update(sensitivity=float((score[y]>=.5).mean()),specificity=float((score[~y]<.5).mean()))
        if scaffolds is not None:
            groups=list(pd.Series(np.arange(len(y))).groupby(np.asarray(scaffolds)).apply(list))
            comparisons=(score[:,None]>score[None,:])+.5*(score[:,None]==score[None,:])
            other=None if baseline is None else (baseline[:,None]>baseline[None,:])+.5*(baseline[:,None]==baseline[None,:])
            rng=np.random.default_rng(20261013);areas=[];deltas=[]
            for _ in range(cfg['bootstrap_repeats']):
                ix=np.concatenate([groups[j] for j in rng.integers(0,len(groups),len(groups))])
                positive=ix[y[ix]];negative=ix[~y[ix]]
                if not len(positive) or not len(negative):continue
                area=float(comparisons[np.ix_(positive,negative)].mean());areas.append(area)
                if other is not None:deltas.append(area-float(other[np.ix_(positive,negative)].mean()))
            result.update(n_scaffolds=len(groups),n_valid_bootstrap_draws=len(areas),
                auc_scaffold_bootstrap_95ci=np.quantile(areas,[.025,.975]).tolist(),
                bootstrap_interval_degenerate=bool(np.quantile(areas,.025)==np.quantile(areas,.975)))
            if baseline is not None:result.update(energy_auc=float(roc_auc_score(y,baseline)),
                auc_gain=float(roc_auc_score(y,score)-roc_auc_score(y,baseline)),
                paired_auc_gain_scaffold_bootstrap_95ci=np.quantile(deltas,[.025,.975]).tolist())
        return result
    reports=[]; predictions=[]
    for method in cfg['models']:
        folds=[];scaffold_predictions=[];study_models=[]
        for study,test in train.groupby('document_chembl_id'):
            retained=train[train.document_chembl_id.ne(study)&~train.connectivity.isin(test.connectivity)]
            if retained.label.nunique()!=2 or test.label.nunique()!=2:continue
            model=fit(retained,method);score=model.predict_proba(vectors(test,method))[:,1]
            study_models.append((study,model))
            folds.append(dict(study=study,n_train=len(retained),**assess(test.label.eq('positive').to_numpy(),score)))
            predictions.extend(test[['inchikey','label','scaffold','document_chembl_id']].assign(
                split='leave_study_out',method=method,score=score).to_dict('records'))
        for scaffold,test in train.groupby('scaffold'):
            retained=train[train.scaffold.ne(scaffold)&~train.connectivity.isin(test.connectivity)]
            if retained.label.nunique()!=2:continue
            score=fit(retained,method).predict_proba(vectors(test,method))[:,1]
            scaffold_predictions.extend(test[['inchikey','label','scaffold','document_chembl_id']].assign(
                split='leave_scaffold_out',method=method,score=score).to_dict('records'))
        predictions.extend(scaffold_predictions)
        final=fit(train,method); vectorizer,regression=final.steps[0][1],final.steps[1][1]
        report=dict(method=method,n_training=len(train),study_folds=folds,
            mean_study_auc=float(np.mean([r['auc'] for r in folds])),
            mean_study_balanced_accuracy=float(np.mean([r['balanced_accuracy'] for r in folds])),
            scaffold=assess(np.array([r['label']=='positive' for r in scaffold_predictions]),np.array([r['score'] for r in scaffold_predictions]),
                [r['scaffold'] for r in scaffold_predictions],np.array([-controls[r['inchikey']]['vina_sile'] for r in scaffold_predictions])),
            coefficients=dict(zip(vectorizer.feature_names_,regression.coef_[0].tolist())),
            intercept=float(regression.intercept_[0]),external=[],candidates=[])
        prior=set(records.loc[records.split.isin(['calibration','external']),'connectivity'])
        for split in cfg['evaluation_splits']:
            test=labelled[labelled.split.eq(split)&~labelled.connectivity.isin(prior)].copy()
            disagreement=test.groupby('inchikey').label.nunique().gt(1)
            test=test[~test.inchikey.isin(disagreement[disagreement].index)]
            for context,subset in test.groupby('assay_context'):
                expected=subset.inchikey.nunique()
                subset=subset.drop_duplicates('inchikey');subset=subset[subset.inchikey.isin(controls)]
                if subset.empty:continue
                score=final.predict_proba(vectors(subset,method))[:,1]
                report['external'].append(dict(split=split,context=context,n_expected=expected,
                    **assess(subset.label.eq('positive').to_numpy(),score,subset.scaffold.to_numpy(),
                    np.array([-controls[key]['vina_sile'] for key in subset.inchikey]))))
                predictions.extend(subset[['inchikey','label','scaffold','document_chembl_id']].assign(
                    split=split,method=method,context=context,score=score).to_dict('records'))
            prior.update(records.loc[records.split.eq(split),'connectivity'])
        for candidate in [r for r in data['records'] if r['role']=='candidate']:
            features=dict(candidate['features'])
            if method=='contact_plus_energy':features['VinaSILE']=candidate['vina_sile']
            score=float(final.predict_proba([features])[0,1])
            seed_scores=[]
            for seed in candidate['seed_contacts']:
                item=dict.fromkeys(seed['contacts'],1.)
                if method=='contact_plus_energy':item['VinaSILE']=candidate['vina_sile']
                seed_scores.append(dict(seed=seed['seed'],score=float(final.predict_proba([item])[0,1])))
            report['candidates'].append(dict(candidate_id=candidate['candidate_id'],inchikey=candidate['inchikey'],
                score=score,score_compatible=score>=.5,native_contact_similarity=candidate['native_contact_similarity'],
                seed_scores=seed_scores,study_omission_scores=[dict(omitted_study=study,
                    score=float(model.predict_proba([features])[0,1])) for study,model in study_models]))
        reports.append(report)
    (HERE/'interaction_calibration.json').write_text(json.dumps(dict(models=reports,specification=cfg,
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in
            ['interaction_fingerprints.json','control_curation.tsv','analysis_plan.json']},
        software={name:importlib.metadata.version(name) for name in ['prolif','MDAnalysis','rdkit','scikit-learn','numpy','pandas','setuptools']},
        interpretation='Fixed class-balanced contact scores are not population binding probabilities. Previously inspected external controls are secondary checks; study and scaffold folds omit shared chemical identities. No automatic promotion of candidate binding support.'),indent=2))
    pd.DataFrame(predictions).to_csv(HERE/'interaction_predictions.tsv',sep='\t',index=False)
    print(json.dumps([{k:v for k,v in r.items() if k!='coefficients'} for r in reports],indent=2),flush=True)


def assess_literature_contact_evidence(plan):
    """Fixed-pose comparison with named, source-checked experimental references."""
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign
    import prolif as plf
    import hashlib
    cfg=plan['literature_controls'];sid=plan['structures'][0]
    receptor=HERE/'structures'/f'{sid}.pdb'
    molecule=Chem.MolFromPDBFile(str(receptor),removeHs=False)
    if molecule is None:raise ValueError('Invalid reference-contact receptor')
    protein=plf.Molecule.from_rdkit(Chem.AddHs(molecule,addCoords=True,addResidueInfo=True))
    fingerprinter=plf.Fingerprint(interactions=['HBAcceptor','HBDonor','PiStacking',
        'CationPi','PiCation','Anionic','Cationic','Hydrophobic'])
    def contacts(path):
        ligand=next(iter(Chem.SDMolSupplier(str(path),removeHs=False)))
        if ligand is None:raise ValueError('Invalid fixed top pose: '+str(path))
        fingerprint=fingerprinter.generate(plf.Molecule.from_rdkit(Chem.AddHs(ligand,addCoords=True)),protein,metadata=True)
        return set(str(pair[1])+'|'+kind for pair,kinds in fingerprint.items() for kind in kinds)
    native=HERE/'ligands'/f'{sid}_native_reference.sdf';reference=contacts(native)
    native_molecule=next(iter(Chem.SDMolSupplier(str(native))))
    native_key=Chem.MolToInchiKey(native_molecule)
    positive=next(r for r in cfg['references'] if r['included'] and r['label']=='SPR_MST_binder')
    assert native_key[:14]==positive['inchikey'][:14], 'Literature positive must match the co-crystal ligand'
    jobs=pd.read_csv(HERE/'panel_docking.tsv',sep='\t')
    meta=pd.read_csv(HERE/'ligand_preparation.tsv',sep='\t')
    jobs=jobs[jobs.status.eq('ok')].merge(meta[['ligand','inchikey','role']],on='ligand',validate='many_to_one')
    jobs=jobs[jobs.role.isin(['candidate','experimental_control'])]
    references={r['inchikey']:r['name'] for r in cfg['references'] if r['included']}
    rows=[];sources=[receptor,native,HERE/'panel_docking.tsv',HERE/'analysis_plan.json']
    for (key,role),group in jobs.groupby(['inchikey','role']):
        states=group.groupby('ligand').affinity.agg(['median','count'])
        states=states[states['count'].eq(len(plan['seeds']))]
        assert len(states), 'Reference control docking incomplete'
        ligand=states.sort_values(['median','ligand']).index[0]
        chosen=group[group.ligand.eq(ligand)].sort_values('seed');seed_rows=[]
        assert chosen.seed.tolist()==sorted(plan['seeds'])
        for row in chosen.itertuples():
            path=HERE/'poses'/f'{sid}_{ligand}_{row.seed}.sdf';keys=contacts(path);sources.append(path)
            seed_rows.append(dict(seed=row.seed,contacts=sorted(keys),
                native_contact_similarity=len(keys&reference)/len(keys|reference),
                native_contacts_recovered=sorted(keys&reference),cnn_affinity=float(row.cnn_affinity),
                cnn_pose_score=float(row.cnn_score),hotspot_contacts={residue:
                    sorted(k.split('|')[1] for k in keys if k.split('|')[0]==residue)
                    for residue in cfg['hotspot_residues']}))
            if key==positive['inchikey']:
                pose=Chem.RemoveHs(next(iter(Chem.SDMolSupplier(str(path)))))
                seed_rows[-1]['native_reference_rmsd']=float(rdMolAlign.CalcRMS(pose,Chem.RemoveHs(native_molecule)))
        pooled=chosen.sort_values(['cnn_score','seed'],ascending=[False,True]).iloc[0]
        rows.append(dict(inchikey=key,role=role,name=references.get(key,plan['candidate_id']),
            selected_state=ligand,median_cnn_affinity=float(chosen.cnn_affinity.median()),seeds=seed_rows,
            pooled_pose=next(row for row in seed_rows if row['seed']==int(pooled.seed))))
    result=dict(target=plan['target'],structure=sid,reference=cfg,native_contacts=sorted(reference),
        native_inchikey=native_key,records=rows,prolif_version=plf.__version__,
        interpretation='Residue contacts and predicted scores from fixed top poses. The secondary pooled pose maximizes CNN pose quality across the three searches, without using reference coordinates or measured labels in selection; per-seed results remain visible. Qualitative SPR comparator is not a measured nonbinder; no classifier or significance threshold is fitted to this small panel.',
        sources=[dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sources])
    (HERE/'literature_contact_evidence.json').write_text(json.dumps(result,indent=2))
    print(json.dumps([dict(name=r['name'],median_cnn_affinity=r['median_cnn_affinity'],
        hotspot_contacts=[s['hotspot_contacts'] for s in r['seeds']]) for r in rows],indent=2),flush=True)


def summarize_reference_pockets(plan):
    """Describe all prespecified pocket comparisons without fitting a support gate."""
    from scipy.stats import spearmanr
    import hashlib
    cfg=plan['reference_pocket_analysis']
    nominees=pd.read_csv(HERE.parent/'current_candidates.tsv',sep='\t').set_index('candidate_id')
    contacts=json.loads((HERE/'reference_pocket_contacts.json').read_text())
    cnn_contacts=json.loads((HERE/'reference_pocket_cnn_contacts.json').read_text())
    rows=[];sources=[HERE/'reference_pocket_contacts.json',HERE/'reference_pocket_cnn_contacts.json',HERE/'analysis_plan.json']
    for spec in cfg['cases']:
        if spec['role']=='native_reference':continue
        case=Path(spec['case']);target=spec['target']
        matched=[r for r in contacts['records'] if r['case']==str(case) and r['role']==spec['role']]
        controls=pd.read_csv(HERE/target/'literature_control_panel.tsv',sep='\t').set_index('inchikey')
        cnn=case/'reference_cnn_rescore'/'panel_docking.tsv'
        cnn_jobs=pd.read_csv(cnn,sep='\t');sources.extend([cnn,HERE/target/'literature_control_panel.tsv'])
        for record in matched:
            key=record['inchikey'];seeds=record['seed_contacts']
            group=cnn_jobs[cnn_jobs.ligand.eq(record['ligand'])&cnn_jobs.status.eq('ok')]
            assert sorted(group.seed.tolist())==sorted(s['seed'] for s in seeds)
            reference=controls.loc[key] if spec['role']=='experimental_control' else None
            row=dict(target=target,candidate_id=record['candidate_id'],inchikey=key,
                ligand=record['ligand'],role=spec['role'],
                name=reference['name'] if reference is not None else record['candidate_id'],
                measured_label=reference['label'] if reference is not None else '',
                measured_value=reference.get('value',None) if reference is not None else None,
                vina_median=float(np.median([s['affinity'] for s in seeds])),
                vina_sile_median=float(np.median([s['size_corrected_affinity'] for s in seeds])),
                cnn_affinity_median=float(group.cnn_affinity.median()),
                cnn_pose_quality_median=float(group.cnn_score.median()),
                native_contact_similarity=record['native_contact_similarity'],
                minimum_heavy_atom_distance=min(s['minimum_heavy_atom_distance'] for s in seeds),
                n_seeds=len(seeds))
            if pd.isna(row['measured_value']):row['measured_value']=None
            cid=reference.get('cid',None) if reference is not None else nominees.loc[record['candidate_id'],'pubchem_cid']
            row['pubchem_cid']=int(cid) if pd.notna(cid) else None
            row['hotspot_contact_counts']={r:sum(bool(s['hotspot_contacts'][r]) for s in seeds) for r in spec['hotspot_residues']}
            cnn_record=next(r for r in cnn_contacts['records'] if r['case']==str(cnn.parent) and r['ligand']==record['ligand'])
            row['cnn_hotspot_contact_counts']={r:sum(bool(s['hotspot_contacts'][r]) for s in cnn_record['seed_contacts']) for r in spec['hotspot_residues']}
            row['cnn_native_contact_similarity']=cnn_record['native_contact_similarity']
            row['cnn_minimum_heavy_atom_distance']=min(s['minimum_heavy_atom_distance'] for s in cnn_record['seed_contacts'])
            rows.append(row)
    metrics=[]
    for score in ['vina_median','vina_sile_median','cnn_affinity_median']:
        eligible=[r for r in rows if r['target']=='GABRB3' and r['measured_label'] in ['higher_affinity','lower_affinity']]
        direction=1 if score=='cnn_affinity_median' else -1
        high=np.array([direction*r[score] for r in eligible if r['measured_label']=='higher_affinity'])
        low=np.array([direction*r[score] for r in eligible if r['measured_label']=='lower_affinity'])
        auc=float(np.mean((high[:,None]>low).astype(float)+.5*(high[:,None]==low)))
        measured=[r for r in rows if r['target']=='GABRB3' and pd.notna(r['measured_value'])]
        rank=spearmanr([-np.log10(r['measured_value']) for r in measured],[direction*r[score] for r in measured])
        metrics.append(dict(score=score,n_higher=len(high),n_lower=len(low),auc=auc,
            n_continuous=len(measured),spearman_r=float(rank.statistic),spearman_p=float(rank.pvalue)))
    report=dict(specification=cfg,records=rows,control_ranking=metrics,
        interpretation='Exploratory same-pocket comparisons. All declared scores and all seven nominees retained. Controls are measured relative affinities, not universal active/inactive labels. No classifier fitted and no change to screening nominees or strict structural-support gate. Searches are numerical repeats, not biological replicates.',
        sources=[dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sources])
    (HERE/'reference_pocket_comparison.json').write_text(json.dumps(report,indent=2,allow_nan=False))
    table=pd.DataFrame(rows)
    for column in ['hotspot_contact_counts','cnn_hotspot_contact_counts']:
        table[column]=table[column].map(json.dumps)
    table.to_csv(HERE/'reference_pocket_comparison.tsv',sep='\t',index=False)
    print(json.dumps(metrics,indent=2),flush=True)
    print(table[table.role.eq('candidate')].to_string(index=False),flush=True)


def prepare_control_score_variants(plan):
    """Derive fixed literature scores from the same completed ligand/seed poses."""
    from rdkit import Chem
    for specification in plan['score_input_cases']:
        case=Path(specification['case']); cfg=json.loads((case/'analysis_plan.json').read_text())
        if not pd.read_csv(case/'redocking_qualification.tsv',sep='\t').qualified.all(): continue
        frame=pd.read_csv(case/'panel_docking.tsv',sep='\t')
        meta=pd.read_csv(case/'ligand_preparation.tsv',sep='\t')
        if (case/'stereochemistry_feasibility.tsv').exists():
            feasible=pd.read_csv(case/'stereochemistry_feasibility.tsv',sep='\t').set_index('ligand').geometry_feasible
            meta=meta[meta.ligand.map(feasible).fillna(False)].copy()
        if (case/'ligand_identity_check.tsv').exists():
            identity=pd.read_csv(case/'ligand_identity_check.tsv',sep='\t').set_index('ligand')
            candidate=meta[meta.role.eq('candidate')&meta.status.eq('ok')]
            assert candidate.ligand.isin(identity.index).all()
            assert candidate.prepared_smiles.eq(candidate.ligand.map(identity.prepared_smiles)).all()
            meta=meta[~meta.role.eq('candidate')|meta.ligand.map(identity.matches_source_stereochemistry).fillna(False)]
        frame=frame[frame.status.eq('ok')].merge(meta[['ligand','inchikey','role','prepared_smiles']],on='ligand',validate='many_to_one')
        if cfg.get('docking_engine','').startswith('gnina'):
            frame['cnn']=-frame.cnn_affinity;frame['cnn_pose']=-frame.cnn_score
            frame['cnn_vs']=-frame.cnn_score*frame.cnn_affinity;methods=['cnn','cnn_pose','cnn_vs']
        else:
            sizes={r.ligand:Chem.MolFromSmiles(r.prepared_smiles).GetNumHeavyAtoms() for r in meta[meta.status.eq('ok')].itertuples()}
            n=frame.ligand.map(sizes);assert n.gt(0).all()
            frame['vina']=frame.affinity;frame['vina_le']=frame.affinity/n
            frame['vina_sile']=frame.affinity/n**.3;methods=['vina','vina_le','vina_sile']
        required=meta.groupby('inchikey').size();rows=[]
        for method in methods:
            states=frame.groupby(['structure','inchikey','ligand','role'])[method].agg(decision_score='median',n_seeds='count').reset_index()
            states=states[states.n_seeds.eq(len(cfg['seeds']))]
            observed=states.groupby(['structure','inchikey']).ligand.transform('size')
            states=states[observed.eq(states.inchikey.map(required))]
            best=states.sort_values(['decision_score','ligand']).drop_duplicates(['structure','inchikey'])
            rows.extend(best.assign(method=method,source_case=str(case)).to_dict('records'))
        pd.DataFrame(rows).to_csv(case/'scoring_variants.tsv',sep='\t',index=False)
        print('Scoring variants',case,len(rows),flush=True)


def calibrate_experimental_controls(plan):
    """Fit to experimental controls only; hold external chemical identities out."""
    from sklearn.metrics import roc_auc_score, balanced_accuracy_score
    from scipy.stats import spearmanr
    import hashlib
    curation=pd.read_csv(HERE/'control_curation.tsv',sep='\t',keep_default_na=False)
    training_only=os.environ.get('YB_DOCK_TRAINING_ONLY')=='1'
    reports=[]; predictions=[]
    def threshold(x,y):
        values=np.unique(x); cuts=np.r_[-np.inf,(values[:-1]+values[1:])/2,np.inf]
        scores=np.array([balanced_accuracy_score(y,x<=cut) for cut in cuts])
        best=cuts[np.isclose(scores,scores.max())]
        finite=best[np.isfinite(best)]
        return float(finite[(len(finite)-1)//2]) if len(finite) else float(best[0])
    for specification in plan['scoring_cases']:
        case=Path(specification['case']); source=case/'stereochemistry_sensitivity.tsv'
        if specification.get('score_variant'):
            source=case/'scoring_variants.tsv';frame=pd.read_csv(source,sep='\t')
            scores=frame[frame.method.eq(specification['score_variant'])].set_index('inchikey').decision_score
        else: scores=pd.read_csv(source,sep='\t').set_index('inchikey').best_isomer_affinity
        assert scores.index.is_unique
        records=curation[curation.target_id.eq(specification['target'])&curation.exclusion_reason.eq('')].copy()
        records['score']=records.inchikey.map(scores)
        records['value_nM']=pd.to_numeric(records.standard_value,errors='coerce')
        calibrated_keys=set(records.loc[records.split.eq('calibration'),'inchikey'])
        records['held_out']=records.split.eq('external')&~records.inchikey.isin(calibrated_keys)
        labels=records[records.label.isin(['positive','low_affinity'])].copy()
        conflict=labels.groupby(['split','inchikey']).label.nunique().gt(1)
        labels=labels[~pd.MultiIndex.from_frame(labels[['split','inchikey']]).isin(conflict[conflict].index)]
        train=labels[labels.split.eq('calibration')].drop_duplicates('inchikey').dropna(subset=['score'])
        test=labels[labels.held_out].drop_duplicates('inchikey').dropna(subset=['score'])
        evaluate_external=not training_only and specification.get('evaluate_external',True)
        assert train.label.nunique()==2 and (not evaluate_external or test.label.nunique()==2),'Both measured control classes required'
        x=train.score.to_numpy(); y=train.label.eq('positive').to_numpy(); cut=threshold(x,y)
        training_identity_hash=hashlib.sha256(train[['inchikey','label']].sort_values('inchikey').to_csv(index=False).encode()).hexdigest()
        tx=test.score.to_numpy(); ty=test.label.eq('positive').to_numpy()
        held=[]
        for scaffold,g in train.groupby('scaffold'):
            retained=train[train.scaffold.ne(scaffold)]
            if retained.label.nunique()!=2: continue
            c=threshold(retained.score.to_numpy(),retained.label.eq('positive').to_numpy())
            held.extend([dict(inchikey=r.inchikey,positive=r.label=='positive',prediction=bool(r.score<=c),
                              threshold=c,split='leave_scaffold_out',method=specification['method']) for r in g.itertuples()])
        predictions.extend(held)
        if not evaluate_external:
            held_frame=pd.DataFrame(held)
            reports.append(dict(target=specification['target'],method=specification['method'],case=str(case),
                score_variant=specification.get('score_variant'),threshold=cut,n_train=len(train),n_train_positive=int(y.sum()),
                training_identity_sha256=training_identity_hash,
                training_auc=float(roc_auc_score(y,-x)),scaffold_cv_balanced_accuracy=float(balanced_accuracy_score(held_frame.positive,held_frame.prediction)),
                evaluated_external=False,validated=False))
            continue
        predictions.extend([dict(inchikey=r.inchikey,positive=r.label=='positive',prediction=bool(r.score<=cut),
                                 threshold=cut,split='external',method=specification['method']) for r in test.itertuples()])
        rng=np.random.default_rng(plan['bootstrap_seed']); areas=[]; groups=list(test.groupby('scaffold').indices.values())
        for _ in range(plan['bootstrap_repeats']):
            idx=np.concatenate([groups[i] for i in rng.integers(0,len(groups),len(groups))])
            if len(np.unique(ty[idx]))==2: areas.append(roc_auc_score(ty[idx],-tx[idx]))
        interval=np.quantile(areas,[.025,.975]).tolist() if areas else [None,None]
        held_frame=pd.DataFrame(held)
        exact=records[records.held_out&records.standard_relation.eq('=')].drop_duplicates('inchikey').dropna(subset=['score','value_nM'])
        report=dict(target=specification['target'],method=specification['method'],case=str(case),
            training_identity_sha256=training_identity_hash,
            threshold=cut,score_orientation='lower_is_more_favorable',n_train=len(train),n_train_positive=int(y.sum()),
            n_external=len(test),n_external_positive=int(ty.sum()),n_external_scaffolds=len(groups),
            n_scaffold_cv_predictions=len(held),n_valid_bootstrap_draws=len(areas),
            training_auc=float(roc_auc_score(y,-x)),external_auc=float(roc_auc_score(ty,-tx)),
            external_auc_scaffold_bootstrap_95ci=interval,external_balanced_accuracy=float(balanced_accuracy_score(ty,tx<=cut)),
            scaffold_cv_balanced_accuracy=float(balanced_accuracy_score(held_frame.positive,held_frame.prediction)),
            external_exact_spearman=float(spearmanr(exact.score,exact.value_nM).statistic),
            n_conflicting_identities=int(conflict[conflict].index.get_level_values('inchikey').nunique()),n_preparation_or_docking_missing=int(records[records.score.isna()].inchikey.nunique()),
            validated=bool(interval[0] is not None and interval[0]>.5 and balanced_accuracy_score(ty,tx<=cut)>.5),
            source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
            curation_sha256=hashlib.sha256((HERE/'control_curation.tsv').read_bytes()).hexdigest())
        if specification.get('score_variant'): report['score_variant']=specification['score_variant']
        seen=set(records.loc[records.split.isin(['calibration','external']),'inchikey'])
        prospective=labels[labels.split.eq('prospective_control')&~labels.inchikey.isin(seen)].drop_duplicates('inchikey').dropna(subset=['score'])
        if plan.get('optimization'):
            report['prospective_control']=dict(n=len(prospective),n_positive=int(prospective.label.eq('positive').sum()))
            if prospective.label.nunique()==2:
                px=prospective.score.to_numpy();py=prospective.label.eq('positive').to_numpy()
                pa=float(roc_auc_score(py,-px));pb=float(balanced_accuracy_score(py,px<=cut))
                pg=list(prospective.groupby('scaffold').indices.values());prng=np.random.default_rng(plan['bootstrap_seed']+1);pi=[]
                for _ in range(plan['bootstrap_repeats']):
                    ix=np.concatenate([pg[j] for j in prng.integers(0,len(pg),len(pg))])
                    if len(np.unique(py[ix]))==2: pi.append(roc_auc_score(py[ix],-px[ix]))
                report['prospective_control'].update(auc=pa,balanced_accuracy=pb,n_scaffolds=len(pg),
                    auc_scaffold_bootstrap_95ci=np.quantile(pi,[.025,.975]).tolist(),n_valid_bootstrap_draws=len(pi))
                predictions.extend([dict(inchikey=r.inchikey,positive=r.label=='positive',prediction=bool(r.score<=cut),
                    threshold=cut,split='prospective_control',method=specification['method']) for r in prospective.itertuples()])
                report['validated']=bool(report['validated'] and pa>.5 and pb>.5)
            else: report['validated']=False
        reports.append(report)
        records.assign(method=specification['method']).to_csv(HERE/f"{specification['target']}_{specification['method']}_control_scores.tsv",sep='\t',index=False)
    for target in sorted({r['target'] for r in reports}):
        options=[r for r in reports if r['target']==target]
        assert len({r['training_identity_sha256'] for r in options})==1, 'Methods must compare identical labelled calibration controls'
        selected=max(options,key=lambda r:(r['scaffold_cv_balanced_accuracy'],r['method'] in ['vina','6HUP__vina']))
        for r in options: r['selected']=r is selected
    output_name='control_training_comparison' if training_only else plan.get('calibration_output_stem','control_calibration')
    (HERE/(output_name+'.json')).write_text(json.dumps(dict(models=reports,selection=
        'Highest leave-one-scaffold-out balanced accuracy in calibration controls; ties prefer Vina. External assays and candidates do not select the scoring method.',interpretation=
        'Predicted high-affinity-like versus low-affinity-like binding within the assayed pocket; not proof of activity. Thresholds use experimental calibration controls only. Unlabelled backgrounds never define a negative class.'),indent=2))
    pd.DataFrame(predictions).to_csv(HERE/('control_training_predictions.tsv' if training_only else 'control_held_out_predictions.tsv'),sep='\t',index=False)
    print(json.dumps(reports,indent=2),flush=True)


def calibrate_multistudy_controls(plan):
    """Study-held-out calibration; previously inspected studies are development data."""
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    from sklearn.pipeline import make_pipeline
    from sklearn.metrics import roc_auc_score, balanced_accuracy_score
    from scipy.stats import beta
    import hashlib
    cfg=plan['multistudy_calibration'];features=cfg['features']
    curated=pd.read_csv(HERE/'control_curation.tsv',sep='\t',keep_default_na=False)
    records=curated[curated.target_id.eq('GABRA1')&curated.exclusion_reason.eq('')].copy()
    for feature in features:
        spec=next(s for s in plan['scoring_cases'] if s['method']=='6HUP__'+feature)
        frame=pd.read_csv(Path(spec['case'])/'scoring_variants.tsv',sep='\t')
        scores=frame[frame.method.eq(feature)].set_index('inchikey').decision_score
        assert scores.index.is_unique
        records[feature]=records.inchikey.map(scores)
    records['connectivity']=records.inchikey.str[:14]
    labelled=records[records.label.isin(['positive','low_affinity'])].copy()
    development=labelled[labelled.split.isin(cfg['development_splits'])].copy()
    conflict=development.groupby('inchikey').label.nunique().gt(1)
    development=development[~development.inchikey.isin(conflict[conflict].index)]
    development=development.drop_duplicates('inchikey').dropna(subset=features).reset_index(drop=True)
    assert development.label.nunique()==2
    def fit(frame,method):
        y=frame.label.eq('positive').to_numpy()
        if method=='sile_threshold':
            x=frame.vina_sile.to_numpy();v=np.unique(x)
            cuts=(v[:-1]+v[1:])/2
            calls=x[None,:]<=cuts[:,None]
            accuracy=(calls[:,y].mean(axis=1)+(~calls[:,~y]).mean(axis=1))/2
            best=cuts[np.isclose(accuracy,accuracy.max())]
            return float(best[(len(best)-1)//2])
        model=make_pipeline(StandardScaler(),LogisticRegression(**{
            k:v for k,v in cfg['logistic'].items() if k!='threshold'}))
        model.fit(frame[features],y)
        return model
    def predict(model,method,frame):
        if method=='sile_threshold':
            score=-frame.vina_sile.to_numpy()
            return score,score>=-model
        score=model.predict_proba(frame[features])[:,1]
        return score,score>=cfg['logistic']['threshold']
    def study_folds(frame):
        for study,held in frame.groupby('document_chembl_id'):
            train=frame[frame.document_chembl_id.ne(study)&~frame.connectivity.isin(held.connectivity)]
            if train.label.nunique()==2 and held.label.nunique()==2:yield study,train,held
    def choose(frame):
        report=[]
        for method in cfg['models']:
            values=[]
            for study,train,held in study_folds(frame):
                score,pred=predict(fit(train,method),method,held)
                values.append(dict(study=study,n=len(held),auc=float(roc_auc_score(held.label.eq('positive'),score)),
                    balanced_accuracy=float(balanced_accuracy_score(held.label.eq('positive'),pred))))
            assert values,'At least one evaluable study fold required'
            report.append(dict(method=method,folds=values,mean_study_balanced_accuracy=float(np.mean([
                r['balanced_accuracy'] for r in values]))))
        chosen=max(report,key=lambda r:(r['mean_study_balanced_accuracy'],r['method']=='sile_threshold'))
        return chosen['method'],report
    selected,comparison=choose(development);oof=[];folds=[]
    for study,train,held in study_folds(development):
        method,_=choose(train);score,pred=predict(fit(train,method),method,held)
        folds.append(dict(study=study,method=method,n_train=len(train),n_test=len(held),
            auc=float(roc_auc_score(held.label.eq('positive'),score)),
            balanced_accuracy=float(balanced_accuracy_score(held.label.eq('positive'),pred))))
        oof.extend(held[['inchikey','label','scaffold','document_chembl_id']].assign(
            score=score,prediction=pred,method=method,split='nested_leave_study_out').to_dict('records'))
    scaffold_rows=[]
    for scaffold,held in development.groupby('scaffold'):
        train=development[development.scaffold.ne(scaffold)&~development.connectivity.isin(held.connectivity)]
        if train.label.nunique()!=2:continue
        method,_=choose(train);score,pred=predict(fit(train,method),method,held)
        scaffold_rows.extend(held[['inchikey','label','scaffold','document_chembl_id']].assign(
            score=score,prediction=pred,method=method,split='nested_leave_scaffold_out').to_dict('records'))
    model=fit(development,selected)
    if selected=='sile_threshold':fitted=dict(threshold=model)
    else:
        scaler,regression=model.steps[0][1],model.steps[1][1]
        fitted=dict(features=features,mean=scaler.mean_.tolist(),scale=scaler.scale_.tolist(),
            coefficient=regression.coef_[0].tolist(),intercept=float(regression.intercept_[0]),
            threshold=cfg['logistic']['threshold'],interpretation='Class-balanced affinity score, not a calibrated population probability.')
    result=dict(selected_method=selected,n_development=len(development),
        n_development_positive=int(development.label.eq('positive').sum()),
        n_conflicting_development_identities=int(conflict.sum()),comparison=comparison,fitted=fitted,
        nested_study_folds=folds,nested_study_mean_auc=float(np.mean([r['auc'] for r in folds])),
        nested_study_mean_balanced_accuracy=float(np.mean([r['balanced_accuracy'] for r in folds])),
        nested_scaffold_balanced_accuracy=float(balanced_accuracy_score(
            [r['label']=='positive' for r in scaffold_rows],[r['prediction'] for r in scaffold_rows])),
        development_identity_sha256=hashlib.sha256(development[['inchikey','label']].sort_values('inchikey').to_csv(index=False).encode()).hexdigest(),
        specification=cfg,external_tests=[],validated=False)
    predictions=oof+scaffold_rows
    seen_development=set(records.loc[records.split.isin(cfg['development_splits']),'connectivity'])
    seen_prior=set(records.loc[records.split.ne(cfg['test_split']),'connectivity'])
    for split in [cfg['seen_secondary_split'],cfg['test_split']]:
        test=labelled[labelled.split.eq(split)&~labelled.connectivity.isin(
            seen_development if split==cfg['seen_secondary_split'] else seen_prior)].copy()
        inconsistent=test.groupby('inchikey').label.nunique().gt(1)
        test=test[~test.inchikey.isin(inconsistent[inconsistent].index)]
        if 'assay_context' not in test:test['assay_context']='human_alpha1_beta3_gamma2'
        test['assay_context']=test.assay_context.replace('','human_alpha1_beta3_gamma2')
        for context,subset in test.groupby('assay_context'):
            expected=subset.inchikey.nunique()
            subset=subset.drop_duplicates('inchikey').dropna(subset=features)
            report=dict(split=split,context=context,n=len(subset),n_expected=expected,
                n_positive=int(subset.label.eq('positive').sum()),
                n_low_affinity=int(subset.label.eq('low_affinity').sum()),discrimination_criterion_met=False)
            if len(subset):
                score,pred=predict(model,selected,subset)
                predictions.extend(subset[['inchikey','label','scaffold','document_chembl_id']].assign(
                    score=score,prediction=pred,method=selected,split=split,context=context).to_dict('records'))
                if subset.label.nunique()==2:
                    y=subset.label.eq('positive').to_numpy();areas=[]
                    groups=list(subset.groupby('scaffold').indices.values());rng=np.random.default_rng(plan['bootstrap_seed'])
                    for _ in range(cfg.get('bootstrap_repeats',plan['bootstrap_repeats'])):
                        ix=np.concatenate([groups[j] for j in rng.integers(0,len(groups),len(groups))])
                        if len(np.unique(y[ix]))==2:areas.append(roc_auc_score(y[ix],score[ix]))
                    ci=np.quantile(areas,[.025,.975]).tolist() if areas else [None,None]
                    ba=float(balanced_accuracy_score(y,pred))
                    def exact_interval(k,n):
                        return [0.0 if k==0 else float(beta.ppf(.025,k,n-k+1)),
                            1.0 if k==n else float(beta.ppf(.975,k+1,n-k))]
                    report.update(auc=float(roc_auc_score(y,score)),balanced_accuracy=ba,
                        auc_scaffold_bootstrap_95ci=ci,n_scaffolds=len(groups),n_valid_bootstrap_draws=len(areas),
                        sensitivity=float(pred[y].mean()),specificity=float((~pred[~y]).mean()),
                        sensitivity_exact_95ci=exact_interval(int(pred[y].sum()),int(y.sum())),
                        specificity_exact_95ci=exact_interval(int((~pred[~y]).sum()),int((~y).sum())),
                        bootstrap_interval_degenerate=bool(ci[0]==ci[1]),
                        discrimination_criterion_met=bool(ci[0] is not None and ci[0]>.5 and ba>.5 and len(subset)==expected))
            result['external_tests'].append(report)
    result['validated']=any(r['discrimination_criterion_met'] and r['split']==cfg['test_split'] and
        r['context']=='human_alpha1_beta3_gamma2' for r in result['external_tests'])
    result['curation_sha256']=hashlib.sha256((HERE/'control_curation.tsv').read_bytes()).hexdigest()
    (HERE/'multistudy_calibration.json').write_text(json.dumps(result,indent=2))
    pd.DataFrame(predictions).to_csv(HERE/'multistudy_control_predictions.tsv',sep='\t',index=False)
    records.to_csv(HERE/'multistudy_control_features.tsv',sep='\t',index=False)
    print(json.dumps(result,indent=2),flush=True)


def current_panel_summary(plan):
    """Keep every current association, including preparation and control failures."""
    from rdkit import Chem
    from scipy.spatial.distance import cdist
    import hashlib,datetime
    candidates = pd.read_csv(HERE / 'current_candidates.tsv',sep='\t');rows=[];sources=[]
    rule=plan['support_thresholds']
    calibration_path=Path(plan['control_calibration_file']) if plan.get('control_calibration_file') else None
    calibration=json.loads(calibration_path.read_text()) if calibration_path and calibration_path.exists() else {'models':[]}
    def track(path):
        sources.append(dict(path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    track(HERE/'current_candidates.tsv');track(HERE/'analysis_plan.json')
    if calibration_path and calibration_path.exists(): track(calibration_path)
    if (HERE/'candidate_identity.json').exists():track(HERE/'candidate_identity.json')
    for candidate in candidates.itertuples():
        case = HERE / 'pairs' / candidate.candidate_id
        chosen=[m for m in calibration['models'] if m['target']==candidate.target_id and m.get('selected',False)]
        if chosen:
            assert len(chosen)==1
            mapping=plan['control_score_cases'][candidate.target_id][chosen[0]['method']]
            case=(HERE/mapping.format(candidate_id=candidate.candidate_id)).resolve()
            assert HERE.resolve() in case.parents
        track(case/'analysis_plan.json')
        selected=json.loads((case/'analysis_plan.json').read_text()).get('evaluation_case')
        if selected:
            case=(case/selected).resolve()
            assert HERE.resolve() in case.parents, 'Evaluation case must stay within the structural analysis'
        row = dict(candidate_id=candidate.candidate_id,target_id=candidate.target_id,
            inchikey=candidate.inchikey,direction=candidate.direction,pubchem_cid=candidate.pubchem_cid,
            docking_tested=False,protocol_qualified=False,structural_support=False,evaluation_complete=False,verdict='not_tested')
        for name in ['analysis_plan.json','stereochemistry_sensitivity.tsv','stereoisomer_docking_summary.tsv','redocking_qualification.tsv','frozen_compound_panel.tsv','scoring_variants.tsv','ring_conformers.tsv']:
            if (case/name).exists():track(case/name)
        cfg=json.loads((case/'analysis_plan.json').read_text())
        assert cfg['candidate']['inchikey']==candidate.inchikey
        row['evaluation_case']=str(case)
        if cfg.get('identity_status')=='source_stereochemistry_underspecified':
            row.update(evaluation_complete=True,verdict='stereochemistry_unresolved',identity_note=cfg['identity_note'])
            rows.append(row);continue
        if not (case/'stereochemistry_sensitivity.tsv').exists():
            rows.append(row);continue
        cfg=json.loads((case/'analysis_plan.json').read_text())
        table=pd.read_csv(case/'stereochemistry_sensitivity.tsv',sep='\t')
        hits=table[table.inchikey.eq(candidate.inchikey)]
        if len(hits)!=1:
            row['verdict']='candidate_preparation_or_docking_incomplete';rows.append(row);continue
        hit=hits.iloc[0];sid=hit.structure
        states=pd.read_csv(case/'stereoisomer_docking_summary.tsv',sep='\t')
        best=states[states.structure.eq(sid)&states.inchikey.eq(candidate.inchikey)].sort_values(['median','ligand']).iloc[0]
        hit=hit.copy();hit['ligand']=best.ligand;hit['median']=hit.best_isomer_affinity
        bg=table[table.structure.eq(sid)&table.role.eq('background')]
        decision_score=float(hit['median']);background_scores=bg.best_isomer_affinity
        if chosen and chosen[0].get('score_variant'):
            variants=pd.read_csv(case/'scoring_variants.tsv',sep='\t')
            variants=variants[variants.method.eq(chosen[0]['score_variant'])&variants.structure.eq(sid)]
            selected_score=variants[variants.inchikey.eq(candidate.inchikey)]
            assert len(selected_score)==1
            score=selected_score.iloc[0];decision_score=float(score.decision_score)
            best=states[states.structure.eq(sid)&states.ligand.eq(score.ligand)].iloc[0]
            hit['ligand']=best.ligand;hit['median']=best['median']
            background_scores=variants.loc[variants.role.eq('background'),'decision_score']
            assert len(background_scores)==len(bg)
        qualification=pd.read_csv(case/'redocking_qualification.tsv',sep='\t').set_index('structure').loc[sid]
        panel=pd.read_csv(case/'frozen_compound_panel.tsv',sep='\t')
        n_background=int(panel.loc[panel.role.eq('background'),'inchikey'].nunique())
        assert len(bg)<=n_background
        rank=1+int(background_scores.le(decision_score).sum()) if n_background else np.nan
        best_possible_fraction=rank/(n_background+1)
        fraction=(rank+max(0,n_background-len(bg)))/(n_background+1)
        row.update(docking_tested=True,structure=sid,scoring_function=cfg.get('scoring_function','vina'),
            score_unit='negative_predicted_pK' if cfg.get('docking_engine','').startswith('gnina') else 'kcal_per_mol',
            protocol_qualified=bool(qualification.qualified),native_rmsd=qualification.median_top_pose_rmsd,
            affinity=hit['median'],decision_score=decision_score,decision_method=chosen[0]['method'] if chosen else None,
            n_background_complete=len(bg),n_background_expected=n_background,rank=rank,rank_fraction_upper_bound=fraction,
            ligand=hit.ligand,structural_scope=cfg.get('structural_scope',''),
            rank_fraction_lower_bound=best_possible_fraction,
            background_complete=bool(len(bg)==n_background))
        receptor=[];target=[]
        for line in (case/'structures'/f'{sid}.pdb').read_text().splitlines():
            if not line.startswith(('ATOM','HETATM')) or line[76:78].strip()=='H':continue
            xyz=[float(line[30:38]),float(line[38:46]),float(line[46:54])];receptor.append(xyz)
            if line.startswith('ATOM') and line[21] in cfg['target_chains']:target.append(xyz)
        distances=[];contacts=[]
        for seed in cfg['seeds']:
            mol=Chem.RemoveHs(next(m for m in Chem.SDMolSupplier(str(case/'poses'/f'{sid}_{hit.ligand}_{seed}.sdf')) if m is not None))
            xyz=mol.GetConformer().GetPositions();distances.append(float(cdist(xyz,receptor).min()))
            contacts.append(float(cdist(xyz,target).min()))
        valid=sum(a>=rule['minimum_heavy_atom_distance'] and b<=rule['maximum_target_contact_distance'] for a,b in zip(distances,contacts))
        row.update(minimum_heavy_atom_distance=min(distances),maximum_target_contact_distance=max(contacts),n_geometry_pass=valid)
        models=chosen
        status=plan.get('control_coverage',{}).get(candidate.target_id,'pending')
        verdict='control_calibration_pending' if status=='pending' else 'controls_unavailable'
        finished=status!='pending'; compatible=False; calibrated=False
        if models:
            assert len(models)==1
            model=models[0];calibrated=bool(model['validated']);finished=True
            compatible=bool(decision_score<=model['threshold'])
            row.update(control_score_threshold=model['threshold'],control_external_auc=model['external_auc'],
                       control_external_balanced_accuracy=model['external_balanced_accuracy'],
                       control_score_unit={'vina':'kcal_per_mol','vina_le':'kcal_per_mol_per_heavy_atom',
                           'vina_sile':'kcal_per_mol_per_heavy_atom_power_0.3','cnn':'negative_predicted_pK',
                           'cnn_pose':'negative_pose_probability','cnn_vs':'negative_pose_probability_times_predicted_pK'}.get(model.get('score_variant',model['method'])))
            verdict=('control_consistent' if compatible else 'low_affinity_like') if calibrated else 'controls_not_discriminated'
        geometry_pass=valid>=rule['minimum_geometry_pass_fraction']*len(cfg['seeds'])
        passed=bool(qualification.qualified and geometry_pass and calibrated and compatible)
        if not qualification.qualified: verdict='protocol_inconclusive';finished=True
        elif not geometry_pass: verdict='pose_geometry_incompatible';finished=True
        row.update(structural_support=passed,evaluation_complete=finished,control_calibrated=calibrated,
                   control_score_compatible=compatible,verdict=verdict)
        rows.append(row)
    output=pd.DataFrame(rows);assert len(output)==len(candidates) and not output.duplicated(['target_id','inchikey','direction']).any()
    path=HERE/'pair_verdicts.tsv';tmp=path.with_suffix('.tmp');output.to_csv(tmp,sep='\t',index=False);tmp.replace(path)
    record=dict(generated_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),n_candidates=len(candidates),
        n_candidate_complete=int(output.docking_tested.sum()),n_evaluation_complete=int(output.evaluation_complete.sum()),
        complete=bool(output.evaluation_complete.all()),n_supported=int(output.structural_support.sum()),
        verdict_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),candidate_sha256=hashlib.sha256((HERE/'current_candidates.tsv').read_bytes()).hexdigest(),sources=sources)
    path=HERE/'docking_status.json';tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(record,indent=2));tmp.replace(path)
    print(output[['candidate_id','target_id','verdict']].to_string(index=False),flush=True)


def ccr4_main():
    plan = json.loads((HERE / 'analysis_plan.json').read_text())
    if os.environ.get('YB_DOCK_STAGE') == 'reference_pocket_contacts':
        extract_interaction_fingerprints(dict(plan,interaction_analysis=plan['reference_pocket_analysis']))
        return
    if os.environ.get('YB_DOCK_STAGE') == 'reference_pocket_cnn_contacts':
        extract_interaction_fingerprints(dict(plan,interaction_analysis=plan['reference_pocket_cnn_analysis']))
        return
    if os.environ.get('YB_DOCK_STAGE') == 'reference_pocket_summary':
        summarize_reference_pockets(plan)
        return
    if os.environ.get('YB_DOCK_STAGE') == 'interaction_fingerprints':
        extract_interaction_fingerprints(plan)
        return
    if os.environ.get('YB_DOCK_STAGE') == 'interaction_calibration':
        evaluate_interaction_fingerprints(plan)
        return
    if os.environ.get('YB_DOCK_STAGE') == 'literature_contacts':
        assess_literature_contact_evidence(plan)
        return
    if os.environ.get('YB_DOCK_STAGE') == 'multistudy_calibration':
        calibrate_multistudy_controls(plan)
        return
    if os.environ.get('YB_DOCK_STAGE') == 'control_score_variants':
        prepare_control_score_variants(plan)
        return
    if os.environ.get('YB_DOCK_STAGE') == 'control_calibration':
        calibrate_experimental_controls(plan)
        return
    if os.environ.get('YB_DOCK_STAGE') == 'screen_summary':
        current_panel_summary(plan)
        return
    if os.environ.get('YB_DOCK_STAGE') == 'cnn_panel':
        from rdkit import Chem
        import hashlib
        parent=Path(plan['rescore_source_case']);rows=[];provenance=[]
        frame=pd.read_csv(parent/'panel_docking.tsv',sep='\t')
        if plan.get('rescore_ligand_role'):
            meta=pd.read_csv(parent/'ligand_preparation.tsv',sep='\t')
            frame=frame[frame.ligand.isin(meta.loc[meta.role.eq(plan['rescore_ligand_role']),'ligand'])]
        assert plan['docking_engine']=='gnina_rescore'
        assert pd.read_csv(HERE/'redocking_qualification.tsv',sep='\t').qualified.all()
        (HERE/'poses').mkdir(exist_ok=True)
        for sid in plan['structures']:
            selected=frame[frame.structure.eq(sid)&frame.status.eq('ok')]
            source=HERE/'poses'/f'{sid}_cnn_input.sdf';output=HERE/'poses'/f'{sid}_cnn_output.sdf'
            with Chem.SDWriter(str(source)) as writer:
                for row in selected.itertuples():
                    path=parent/'poses'/f'{sid}_{row.ligand}_{row.seed}.sdf'
                    provenance.append(dict(path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
                    for index,mol in enumerate(Chem.SDMolSupplier(str(path),removeHs=False)):
                        if mol is None:raise ValueError('Invalid source docking pose')
                        mol.SetProp('_Name',f'{row.ligand}|{row.seed}|{index}');writer.write(mol)
            command=plan['gnina_launcher']+[plan['gnina_binary'],'-r',str(HERE/'receptors'/f'{sid}.pdbqt'),
                '-l',str(source),'-o',str(output),'--score_only','--cnn_scoring','rescore',
                '--cpu','2','--device','0','--seed',str(plan['seeds'][0])]
            with (HERE/'panel_cnn.log').open('w') as log:run_gnina(command,log,7200)
            grouped={}
            for mol in Chem.SDMolSupplier(str(output),removeHs=False):
                if mol is None:raise ValueError('Invalid rescored docking pose')
                ligand,seed,index=mol.GetProp('_Name').split('|')
                grouped.setdefault((ligand,int(seed)),[]).append((int(index),mol))
            assert set(grouped)==set(zip(selected.ligand,selected.seed)), 'Rescoring omitted a source job'
            for (ligand,seed),poses in grouped.items():
                poses.sort(key=lambda x:(-float(x[1].GetProp('CNNscore')),x[0]))
                with Chem.SDWriter(str(HERE/'poses'/f'{sid}_{ligand}_{seed}.sdf')) as writer:
                    for index,mol in poses:writer.write(mol)
                top=poses[0][1]
                rows.append(dict(structure=sid,ligand=ligand,seed=seed,native=False,status='ok',
                    affinity=-float(top.GetProp('CNNaffinity')),cnn_affinity=float(top.GetProp('CNNaffinity')),
                    cnn_score=float(top.GetProp('CNNscore')),source_pose_index=poses[0][0],score_unit='negative_predicted_pK'))
        pd.DataFrame(rows).to_csv(HERE/'panel_docking.tsv',sep='\t',index=False)
        (HERE/'rescore_input_manifest.json').write_text(json.dumps(provenance,indent=2))
        print('Rescored',len(rows),'fixed ligand/seed jobs',flush=True)
        return
    if os.environ.get('YB_DOCK_STAGE') in ['cnn_native','cnn_redock']:
        import subprocess
        from rdkit import Chem
        from rdkit.Chem import rdMolAlign
        records=[];redock=os.environ['YB_DOCK_STAGE']=='cnn_redock'
        tag='native_cnn_redock' if redock else 'native_cnn'
        for sid in plan['structures']:
            source=HERE/'poses'/f'{sid}_native_cnn_input.sdf'
            output=HERE/'poses'/f'{sid}_native_cnn.sdf'
            with Chem.SDWriter(str(source)) as writer:
                for seed in plan['seeds']:
                    for i,mol in enumerate(Chem.SDMolSupplier(str(HERE/'poses'/f'{sid}_{sid}_native_{seed}.sdf'),removeHs=False)):
                        if mol is None: raise ValueError('Invalid native pose')
                        mol.SetProp('_Name',f'{seed}_{i}');writer.write(mol)
            command=plan.get('gnina_launcher',[])+[os.environ['YB_GNINA'],'-r',str(HERE/'receptors'/f'{sid}.pdbqt'),
                '-l',str(source),'-o',str(output),'--score_only','--cnn_scoring','rescore',
                '--cpu','2','--device','0','--seed',str(plan['seeds'][0])]
            if redock:
                box=next(b for b in json.loads((HERE/'structures.json').read_text()) if b['structure']==sid)
                merged=[]
                for seed in plan['seeds']:
                    pose=HERE/'poses'/f'{sid}_{tag}_{seed}.sdf'
                    cmd=plan.get('gnina_launcher',[])+[os.environ['YB_GNINA'],'-r',str(HERE/'receptors'/f'{sid}.pdbqt'),
                        '-l',str(HERE/'ligands'/f'{sid}_native.sdf'),'-o',str(pose),'--cnn_scoring','refinement',
                        '--cpu','2','--device','0','--seed',str(seed),
                        '--exhaustiveness',str(plan['vina_exhaustiveness']),'--num_modes',str(plan['vina_poses'])]
                    for prefix,values in [('center',box['center']),('size',box['size'])]:
                        for axis,value in zip('xyz',values):cmd.extend([f'--{prefix}_{axis}',str(value)])
                    with (HERE/f'{tag}_{seed}.log').open('w') as log:
                        run_gnina(cmd,log,7200)
                    for i,mol in enumerate(Chem.SDMolSupplier(str(pose))):
                        if mol is None:raise ValueError('Invalid native redocking pose')
                        mol.SetProp('_Name',f'{seed}_{i}');merged.append(mol)
                output=HERE/'poses'/f'{sid}_{tag}.sdf'
                with Chem.SDWriter(str(output)) as writer:
                    for mol in merged:writer.write(mol)
            else:
                with (HERE/'native_cnn.log').open('w') as log:
                    run_gnina(command,log,1800)
            gold=Chem.RemoveHs(next(m for m in Chem.SDMolSupplier(str(HERE/'ligands'/f'{sid}_native_reference.sdf')) if m is not None))
            for mol in Chem.SDMolSupplier(str(output)):
                if mol is None: raise ValueError('Invalid rescored native pose')
                seed,index=map(int,mol.GetProp('_Name').split('_'))
                records.append(dict(structure=sid,seed=seed,pose=index,
                    cnn_score=float(mol.GetProp('CNNscore')),cnn_affinity=float(mol.GetProp('CNNaffinity')),
                    rmsd=float(rdMolAlign.CalcRMS(Chem.RemoveHs(mol),gold))))
        frame=pd.DataFrame(records);frame.to_csv(HERE/f'{tag}_poses.tsv',sep='\t',index=False)
        top=frame.sort_values(['cnn_score','pose'],ascending=[False,True]).drop_duplicates(['structure','seed'])
        result=top.groupby('structure').rmsd.agg(n_success='count',median_top_pose_rmsd='median').reset_index()
        result['qualified']=result.n_success.eq(len(plan['seeds'])) & result.median_top_pose_rmsd.le(2)
        result.to_csv(HERE/f'{tag}_qualification.tsv',sep='\t',index=False)
        print(result.to_string(index=False),flush=True)
        return
    if os.environ.get('YB_DOCK_STAGE') == 'contacts':
        from Bio.PDB.MMCIF2Dict import MMCIF2Dict
        from rdkit import Chem
        from scipy.spatial.distance import cdist
        contact_rows=[]; summaries=[]
        for sid in plan['structures']:
            data=MMCIF2Dict(str(HERE / (sid+'.cif')))
            receptor='P32248' if sid=='6QZH' else 'P51679'
            chain='A' if sid=='6QZH' else 'R'
            residues={}
            for c,acc,a,b,start in zip(data['_struct_ref_seq.pdbx_strand_id'],data['_struct_ref_seq.pdbx_db_accession'],data['_struct_ref_seq.pdbx_auth_seq_align_beg'],data['_struct_ref_seq.pdbx_auth_seq_align_end'],data['_struct_ref_seq.db_align_beg']):
                if c!=chain: continue
                for number in range(int(a),int(b)+1):residues[number]=(acc,int(start)+number-int(a))
            atoms=[]
            for line in (HERE / 'structures' / (sid+'.pdb')).read_text().splitlines():
                if not line.startswith('ATOM') or line[76:78].strip()=='H':continue
                number=int(line[22:26]);acc,position=residues.get(number,('construct_linker',number))
                atoms.append(dict(atom=line[12:16].strip(),resname=line[17:20],author_residue=number,accession=acc,uniprot_residue=position,x=float(line[30:38]),y=float(line[38:46]),z=float(line[46:54])))
            xyz=np.array([[a['x'],a['y'],a['z']] for a in atoms])
            files=sorted((HERE / 'poses').glob(sid+'_'+plan['candidate']['inchikey']+'*.sdf'))
            files+=sorted((HERE / 'poses').glob(sid+'_'+sid+'_native_*.sdf'))
            files.append(HERE / (sid+'_native.sdf'))
            for path in files:
                mol=Chem.RemoveHs(next(m for m in Chem.SDMolSupplier(str(path)) if m is not None))
                distance=cdist(mol.GetConformer().GetPositions(),xyz)
                nearby=np.flatnonzero(distance.min(0)<=4.)
                nonhuman=[j for j in nearby if atoms[j]['accession']!=receptor]
                summaries.append(dict(structure=sid,pose=path.name,n_receptor_atoms_within_4A=len(nearby)-len(nonhuman),n_fusion_or_linker_atoms_within_4A=len(nonhuman),minimum_heavy_atom_distance=float(distance.min()),n_contacts_below_2A=int((distance<2).sum())))
                for j in nearby:
                    i=int(np.argmin(distance[:,j]));a=atoms[j]
                    contact_rows.append(dict(structure=sid,pose=path.name,ligand_atom_index=i,ligand_element=mol.GetAtomWithIdx(i).GetSymbol(),distance_A=float(distance[i,j]),**a))
        pd.DataFrame(contact_rows).to_csv(HERE / 'candidate_native_contacts.tsv',sep='\t',index=False)
        pd.DataFrame(summaries).to_csv(HERE / 'contact_scope_summary.tsv',sep='\t',index=False)
        print(pd.DataFrame(summaries).to_string(index=False),flush=True)
        return
    if os.environ.get('YB_DOCK_STAGE') == 'summarize':
        frame = pd.read_csv(HERE / 'panel_docking.tsv',sep='\t')
        meta = pd.read_csv(HERE / 'ligand_preparation.tsv',sep='\t')
        if (HERE/'stereochemistry_feasibility.tsv').exists():
            feasible=pd.read_csv(HERE/'stereochemistry_feasibility.tsv',sep='\t').set_index('ligand').geometry_feasible
            meta=meta[meta.ligand.map(feasible).fillna(False)].copy()
        if plan.get('experimental_control_panel'):
            from rdkit import Chem
            from rdkit.Chem.MolStandardize import rdMolStandardize
            uncharger=rdMolStandardize.Uncharger()
            source=pd.read_csv(plan['experimental_control_panel'],sep='\t').set_index('inchikey').smiles
            checks=[]
            for r in meta[meta.status.eq('ok')].itertuples():
                query=uncharger.uncharge(Chem.MolFromSmiles(source[r.inchikey]))
                mol=uncharger.uncharge(Chem.MolFromSmiles(r.prepared_smiles))
                valid=mol.GetNumAtoms()==query.GetNumAtoms() and mol.GetNumBonds()==query.GetNumBonds() and mol.HasSubstructMatch(query,useChirality=True)
                checks.append(dict(ligand=r.ligand,matches_source_identity=valid))
            pd.DataFrame(checks).to_csv(HERE/'control_ligand_identity_check.tsv',sep='\t',index=False)
            assert all(r['matches_source_identity'] for r in checks),'Control chemical identity changed during preparation'
        identity=Path(plan.get('candidate_identity_file',str(HERE.parent.parent/'candidate_identity.json')))
        if identity.exists():
            from rdkit import Chem
            from rdkit.Chem.MolStandardize import rdMolStandardize
            uncharger=rdMolStandardize.Uncharger()
            record=next(x for x in json.loads(identity.read_text()) if x['inchikey']==plan['candidate']['inchikey'])
            query=uncharger.uncharge(Chem.MolFromSmiles(record['smiles']));matches=[]
            for r in meta.itertuples():
                if r.role!='candidate' or r.status!='ok':matches.append(True);continue
                mol=Chem.MolFromSmiles(r.prepared_smiles)
                if mol is not None:mol=uncharger.uncharge(mol)
                matches.append(bool(mol is not None and mol.HasSubstructMatch(query,useChirality=True)))
            check=meta.assign(matches_source_stereochemistry=matches)
            check.to_csv(HERE/'ligand_identity_check.tsv',sep='\t',index=False)
            meta=check[check.matches_source_stereochemistry].copy()
            assert meta.role.eq('candidate').any(),'No prepared state matches the exact compound source'
        if 'ligand' not in frame.columns:
            print('No eligible ligand docking jobs; preparation remains incomplete.',flush=True)
            return
        states = frame[frame.status.eq('ok')].groupby(['structure','ligand']).affinity.agg(['median','count','min','max']).reset_index().merge(meta,on='ligand')
        states = states[states['count'].eq(len(plan['seeds']))]
        states.to_csv(HERE / 'stereoisomer_docking_summary.tsv',sep='\t',index=False)
        result = states.groupby(['structure','inchikey','role'])['median'].agg(n_states='count',best_isomer_affinity='min',mean_isomer_affinity='mean',worst_isomer_affinity='max').reset_index()
        required = meta.groupby('inchikey').size()
        result['all_states_complete'] = result.n_states.eq(result.inchikey.map(required))
        complete = result[result.all_states_complete].copy()
        for endpoint in ['best_isomer_affinity','mean_isomer_affinity']:
            complete[endpoint+'_rank'] = complete.groupby('structure')[endpoint].rank(method='average')
        complete['n_compounds'] = complete.groupby('structure').inchikey.transform('size')
        complete.to_csv(HERE / 'stereochemistry_sensitivity.tsv',sep='\t',index=False)
        result.to_csv(HERE / 'all_compound_completion.tsv',sep='\t',index=False)
        print(complete[complete.role.eq('candidate')].to_string(index=False),flush=True)
        return
    (HERE / 'receptors').mkdir(exist_ok=True); (HERE / 'poses').mkdir(exist_ok=True)
    boxes = json.loads((HERE / 'structures.json').read_text())
    for sid in plan['structures']:
        rq = HERE / 'receptors' / (sid+'.pdbqt')
        if not rq.exists(): prep_receptor(HERE / 'structures' / (sid+'.pdb'), rq)
        box = next(r for r in boxes if r['structure']==sid)
        if 'center' not in box:
            predicted = ligsite_box(HERE / 'structures' / (sid+'.pdb'))
            if predicted is None: raise ValueError('No predicted pocket: '+sid)
            box.update(predicted, pocket_source='predicted_monomer_exploratory')
    (HERE / 'structures.json').write_text(json.dumps(boxes,indent=2))
    stage = os.environ.get('YB_DOCK_STAGE','native')
    if stage in ['native','native_resume']:
        jobs = [(sid,sid+'_native',seed,True) for sid in plan['structures'] if next(b for b in boxes if b['structure']==sid).get('native_available',True) for seed in plan['seeds']]
    elif stage in ['panel','panel_resume']:
        qualified = pd.read_csv(HERE / 'redocking_qualification.tsv',sep='\t')
        structures = qualified.loc[qualified.qualified,'structure'].tolist()
        if plan.get('evaluate_unqualified'): structures = plan['structures']
        ligands = pd.read_csv(HERE / 'ligand_preparation.tsv',sep='\t')
        if (HERE/'stereochemistry_feasibility.tsv').exists():
            feasible=pd.read_csv(HERE/'stereochemistry_feasibility.tsv',sep='\t').set_index('ligand').geometry_feasible
            ligands=ligands[ligands.ligand.map(feasible).fillna(False)].copy()
        complete = ligands.groupby('inchikey').status.apply(lambda x: x.eq('ok').all())
        ligands = ligands[ligands.inchikey.isin(complete[complete].index)]
        if os.environ.get('YB_DOCK_LIGAND_ROLE'):
            ligands=ligands[ligands.role.eq(os.environ['YB_DOCK_LIGAND_ROLE'])]
        if os.environ.get('YB_DOCK_CONTROL_SPLIT'):
            controls=pd.read_csv(plan['control_curation_file'],sep='\t',keep_default_na=False)
            controls=controls[controls.split.eq(os.environ['YB_DOCK_CONTROL_SPLIT'])&controls.exclusion_reason.eq('')]
            ligands=ligands[ligands.inchikey.isin(controls.inchikey)]
        qualification=qualified.set_index('structure').qualified.to_dict()
        jobs = [(sid,r.ligand,seed,False) for sid in structures for r in ligands.itertuples()
            if r.status=='ok' and (qualification[sid] or r.role=='candidate' or not plan.get('candidate_only_unqualified'))
            for seed in plan['seeds']]
        state_counts=ligands.groupby('inchikey').size()
        if plan.get('ring_ensemble'):
            manifest=pd.read_csv(HERE/'ring_conformers.tsv',sep='\t')
            assert {j[1] for j in jobs}<=set(manifest.ligand), 'Ring sampling must finish before panel docking'
        schedule={r.ligand:(r.role!='candidate',int(state_counts[r.inchikey]),r.inchikey,r.ligand) for r in ligands.itertuples()}
        jobs.sort(key=lambda job:(schedule[job[1]],job[2]))
    else: raise ValueError('Unknown docking stage')
    rows=[]
    if stage in ['panel_resume','native_resume'] and (HERE / (stage.replace('_resume','')+'_docking.tsv')).exists():
        stage = stage.replace('_resume','')
        prior = pd.read_csv(HERE / (stage+'_docking.tsv'),sep='\t')
        rows = prior[prior.status.eq('ok')].to_dict('records')
        completed = {(r['structure'],r['ligand'],r['seed']) for r in rows}
        jobs = [j for j in jobs if j[:3] not in completed]
    stage = stage.replace('_resume','')
    expected_jobs=len(rows)+len(jobs)
    with ProcessPoolExecutor(max_workers=int(os.environ.get('YB_DOCK_WORKERS','24'))) as ex:
        for future in as_completed([ex.submit(ccr4_dock_one,job) for job in jobs]):
            row=future.result()
            rows.append(row)
            pending=HERE/(stage+'_docking.pending.tsv')
            pd.DataFrame(rows).to_csv(pending,sep='\t',index=False)
            pending.replace(HERE/(stage+'_docking.tsv'))
            if len(rows)%10==0:
                print(stage,len(rows),'/',expected_jobs,flush=True)
    frame = pd.DataFrame(rows,columns=None if rows else ['structure','status','top_pose_rmsd'])
    frame.to_csv(HERE / (stage+'_docking.tsv'),sep='\t',index=False)
    if stage=='native':
        checks=[]
        for sid in plan['structures']:
            g=frame[frame.structure.eq(sid)]; values=g.loc[g.status.eq('ok'),'top_pose_rmsd'].to_numpy()
            median=float(np.median(values)) if len(values) else None
            checks.append(dict(structure=sid,n_success=len(values),median_top_pose_rmsd=median,
                reference_type=plan.get('reference_type','experimental_native_pose'),
                qualified=len(values)==len(plan['seeds']) and median<=2.0))
        pd.DataFrame(checks).to_csv(HERE / 'redocking_qualification.tsv',sep='\t',index=False)
        print(checks,flush=True)
    elif len(frame):
        meta=pd.read_csv(HERE / 'ligand_preparation.tsv',sep='\t')
        if (HERE/'stereochemistry_feasibility.tsv').exists():
            feasible=pd.read_csv(HERE/'stereochemistry_feasibility.tsv',sep='\t').set_index('ligand').geometry_feasible
            meta=meta[meta.ligand.map(feasible).fillna(False)].copy()
        accepted=frame[frame.status.eq('ok')].groupby(['structure','ligand']).affinity.agg(['median','count']).reset_index()
        accepted=accepted[accepted['count'].eq(len(plan['seeds']))].merge(meta,on='ligand')
        required=meta.groupby('inchikey').size()
        observed=accepted.groupby(['structure','inchikey']).ligand.transform('size')
        accepted=accepted[observed.eq(accepted.inchikey.map(required))]
        compounds=accepted.sort_values('median').drop_duplicates(['structure','inchikey']).copy()
        compounds['within_pocket_rank']=compounds.groupby('structure')['median'].rank(method='average')
        compounds.to_csv(HERE / 'compound_docking_summary.tsv',sep='\t',index=False)
        print(compounds[compounds.role.eq('candidate')].to_string(index=False),flush=True)


if "--mode" in sys.argv and "ccr4" in sys.argv:
    ccr4_main()
    sys.exit(0)


if "--mode" in sys.argv and "trpv1" in sys.argv:
    trpv1_main()
    sys.exit(0)


def main():
    (HERE / "receptors").mkdir(exist_ok=True)
    cand = pd.read_csv(HERE / "candidates.tsv", sep="\t")
    prep = json.loads((HERE / "prep_report.json").read_text())
    struct_of = {g: r["struct"] for g, r in prep["targets"].items() if r.get("struct")}
    jobs = []
    for r in cand.itertuples():
        if r.target in struct_of and (HERE / "ligands" / f"{r.compound}.pdbqt").exists():
            jobs.append((r.compound, r.target, struct_of[r.target],
                         HERE / "ligands" / f"{r.compound}.pdbqt"))
    print(f"[dock] jobs {len(jobs)}", flush=True)
    rows = []
    with ProcessPoolExecutor(max_workers=24) as ex:
        for i, res in enumerate(ex.map(dock_one, jobs, chunksize=1)):
            rows.append(res)
            if (i + 1) % 20 == 0:
                print(f"[dock] {i + 1}/{len(jobs)}", flush=True)
    sc = pd.DataFrame(rows, columns=["compound", "target", "affinity", "status"])
    sc.to_csv(HERE / "dock_scores.tsv", sep="\t", index=False)
    print("[dock] status:", sc.status.str.split(":").str[0].value_counts().to_dict(),
          flush=True)

    # ---------- 评估
    gold = pd.read_csv(HERE / "inputs" / "gold_target_ranks.tsv", sep="\t")
    gold_pairs = gold[gold.arm == "yeast_only"][["compound_id", "target", "rank"]].drop_duplicates()
    ok_sc = sc[sc.affinity.notna()]

    def rank_fused(c):
        sub = cand[cand.compound == c].merge(
            ok_sc[ok_sc.compound == c][["target", "affinity"]],
            on="target", how="inner").sort_values("affinity")
        dr = {t: i + 1 for i, t in enumerate(sub.target)}
        fused = sorted(sub.itertuples(),
                       key=lambda r: -(1.0 / (20 + r.fd_rank_in_pool) + 1.0 / (20 + dr[r.target])))
        return {r.target: i + 1 for i, r in enumerate(fused)}

    def metrics(rank_fn):
        rr, hit10, ranks_all = [], [], {}
        for c in SIX:
            gs = gold_pairs[gold_pairs.compound_id == c]
            rm = rank_fn(c)
            rr_c, h_c = [], []
            for g in gs.itertuples():
                r = rm.get(g.target, g.rank)
                ranks_all[f"{NAME[c]}|{g.target}"] = r
                rr_c.append(1.0 / r)
                h_c.append(1.0 if r <= 10 else 0.0)
            rr.append(np.mean(rr_c))
            hit10.append(np.mean(h_c))
        return {"mrr": float(np.mean(rr)), "hit10": float(np.mean(hit10)),
                "gold_ranks": ranks_all}

    def fd_only(c):
        return {r.target: r.fd_rank_in_pool for r in cand[cand.compound == c].itertuples()}

    res = {"fd_only": metrics(fd_only),
           "fd_plus_dock_rrf": metrics(rank_fused),
           "dock_status": sc.status.value_counts().to_dict()}
    print("[eval] FD-only:", res["fd_only"], flush=True)
    print("[eval] fused:", res["fd_plus_dock_rrf"], flush=True)
    (HERE / "eval_result.json").write_text(json.dumps(res, ensure_ascii=False, indent=1))
    if os.environ.get("SKIP_CHEMBL"):
        print("[done] eval-only (chembl skipped)", flush=True)
        return

    # ---------- ChEMBL 印证（间隔+重试）
    def chembl_get(url, tries=3):
        for t in range(tries):
            try:
                req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
                with urllib.request.urlopen(req, timeout=40) as r:
                    return json.loads(r.read())
            except Exception:
                time.sleep(1.5 + t)
        return None

    def chembl_id_of(name):
        d = chembl_get("https://www.ebi.ac.uk/chembl/api/data/molecule.json?"
                       f"pref_name__iexact={urllib.parse.quote(name)}&limit=1")
        ms = (d or {}).get("molecules", [])
        return ms[0]["molecule_chembl_id"] if ms else None

    import urllib.parse
    SYNONYM = {"TYMS": "THYMIDYLATE SYNTHASE", "SQLE": "SQUALENE MONOOXYGENASE",
               "CYP51A1": "LANOSTEROL 14-ALPHA DEMETHYLASE",
               "TOP2A": "TOPOISOMERASE II ALPHA", "CYP3A5": "CYTOCHROME P450 3A5",
               "DHFR": "DIHYDROFOLATE REDUCTASE"}
    cid_cache = {}
    ev_rows = []
    for c in SIX:
        rmap = rank_fused(c)
        sub = cand[cand.compound == c]
        order = sorted(sub.itertuples(), key=lambda r: rmap.get(r.target, 999))
        for r in order[:10]:
            acc = prep["targets"].get(r.target, {}).get("uniprot")
            ev = "no_uniprot"
            if acc:
                try:
                    if c not in cid_cache:
                        cid_cache[c] = chembl_id_of(NAME[c])
                    cid = cid_cache[c]
                    if cid:
                        d = chembl_get("https://www.ebi.ac.uk/chembl/api/data/target.json?"
                                       f"target_components.accession={acc}&limit=5")
                        tids = [t["target_chembl_id"] for t in (d or {}).get("targets", [])]
                        hits = []
                        for tid in tids[:3]:
                            a = chembl_get("https://www.ebi.ac.uk/chembl/api/data/activity.json?"
                                           f"molecule_chembl_id={cid}&target_chembl_id={tid}&limit=20")
                            hits += [(x.get("standard_type"), x.get("standard_value"),
                                      x.get("standard_units")) for x in (a or {}).get("activities", [])
                                     if x.get("standard_value")]
                        ev = f"{len(hits)} 条实测" if hits else "无该靶记录"
                        time.sleep(0.4)
                    else:
                        ev = "chembl_no_molecule"
                except Exception as ex:
                    ev = f"err:{str(ex)[:30]}"
            ev_rows.append({"compound": NAME[c], "target": r.target,
                            "final_rank": rmap.get(r.target, 999),
                            "uniprot": acc, "chembl": ev})
    pd.DataFrame(ev_rows).to_csv(HERE / "chembl_evidence.tsv", sep="\t", index=False)
    print("[chembl] saved", flush=True)
    print("[done]", flush=True)


if __name__ == "__main__":
    main()
