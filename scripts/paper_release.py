"""Portable reconstruction of the frozen ion-channel manuscript analysis."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import tempfile

ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT / "data/paper"


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_json(path, value):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n")


def manifest():
    result = json.loads((REFERENCE / "assets.json").read_text())
    extension = REFERENCE / 'extended_assets.json'
    if extension.exists():
        extra = json.loads(extension.read_text())
        if set(result) & set(extra):
            raise ValueError('Extended asset groups shadow the core manifest')
        result.update(extra)
    return result


def inventory(level="reference"):
    groups = ["reference", "inference", "numerical"] if level == "all" else [level]
    result = []
    for group in groups:
        for entry in manifest()[group]:
            p = ROOT / entry["path"]
            actual = digest(p) if p.is_file() else None
            result.append({"group": group, "path": entry["path"],
                           "status": "missing" if actual is None else "ok" if actual == entry["sha256"] else "hash-mismatch",
                           "expected_sha256": entry["sha256"], "actual_sha256": actual})
    return result


def require(level):
    bad = [x for x in inventory(level) if x["status"] != "ok"]
    if bad:
        raise ValueError("Paper assets are missing or changed:\n" + "\n".join(
            f"  {x['status']}: {x['path']}" for x in bad) +
            "\nInstall the matching paper archive; see docs/REPRODUCIBILITY.md.")


def install(archive):
    entries = {x["path"]: x for group in manifest().values() for x in group}
    installed = []
    with tempfile.TemporaryDirectory(prefix="yeastbridge-assets-") as td:
        staging = Path(td)
        with tarfile.open(archive, "r:gz") as t:
            members = t.getmembers()
            names = [m.name for m in members]
            if len(names) != len(set(names)) or any(not m.isfile() or m.name not in entries for m in members):
                raise ValueError("Archive contains unregistered, duplicate, or non-file members")
            # Validate the entire archive before installing any member.
            for m in members:
                p = staging / m.name
                p.parent.mkdir(parents=True, exist_ok=True)
                with t.extractfile(m) as source, p.open("wb") as dest:
                    shutil.copyfileobj(source, dest)
                if p.stat().st_size != entries[m.name]["size"] or digest(p) != entries[m.name]["sha256"]:
                    raise ValueError(f"Archive checksum mismatch: {m.name}")
                target = ROOT / m.name
                if not target.resolve().is_relative_to(ROOT.resolve()):
                    raise ValueError(f"Asset destination escapes checkout: {m.name}")
                if target.exists() and digest(target) != entries[m.name]["sha256"]:
                    raise ValueError(f"Refusing to replace a different local asset: {m.name}")
            for m in members:
                p = ROOT / m.name
                p.parent.mkdir(parents=True, exist_ok=True)
                if not p.exists():
                    temporary = p.with_name(p.name + ".installing")
                    try:
                        shutil.copyfile(staging / m.name, temporary)
                        temporary.replace(p)
                    finally:
                        temporary.unlink(missing_ok=True)
                installed.append(m.name)
    return {"verified_files": len(installed), "archive": Path(archive).name}


def export(output=None):
    require("reference")
    output = Path(output or ROOT / "results.csv").resolve()
    source = REFERENCE / "results.csv"
    if output == source.resolve():
        raise ValueError("Choose an output outside the frozen reference directory")
    output.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, output)
    write_json(output.with_suffix(".metadata.json"), {
        "mode": "published-export", "release": "paper-ion-channel-20261027",
        "rows": 25, "targets": 7, "compounds": 18,
        "sources": {"data/paper/results.csv": digest(source)},
        "result_sha256": digest(output), "computation": "Export of the frozen manuscript candidate table"})
    return output


def new_output(output):
    output = Path(output).resolve()
    if output.exists() and any(output.iterdir()):
        raise ValueError(f"Output must be a new or empty directory: {output}")
    output.mkdir(parents=True, exist_ok=True)
    return output


def screenwide_p(scores, means, standard_deviations, maxima):
    import numpy as np
    scores, means, standard_deviations = map(np.asarray, (scores, means, standard_deviations))
    maxima = np.asarray(maxima)
    if not len(maxima) or not np.isfinite(maxima).all():
        raise ValueError("Invalid screen-wide permutation maxima")
    if not np.isfinite(scores).all() or not np.isfinite(means).all() or not np.isfinite(standard_deviations).all() or (standard_deviations <= 0).any():
        raise ValueError("Invalid observed scores or permutation calibration")
    z = (scores - means) / standard_deviations
    return (1 + len(maxima) - np.searchsorted(np.sort(maxima), z, side="left")) / (len(maxima) + 1)


def reproduce(output, null_root=None, input_root=None):
    """Recompute observed correlations and decisions from explicit input roots."""
    import numpy as np
    import pandas as pd
    from scripts.product_execute_hiphop import compound_groups, max_over_compounds
    require("reference")
    require("numerical")
    output = new_output(output)
    records, tables = {}, {}
    for layer, expected in [("match", 338), ("residual", 95)]:
        base = (Path(input_root) if input_root is not None else ROOT / "raw/paper") / layer
        with np.load(base / "rank_inputs.npz", allow_pickle=False) as z:
            d = dict(z)
        with np.load(base / "screen_arrays.npz", allow_pickle=False) as z:
            saved = dict(z)
        if d["task"].shape != (expected, 5865) or len(d["compounds"]) != 3250:
            raise ValueError(f"Unexpected {layer} scope")
        if not np.isfinite(d["task"]).all() or not np.isfinite(d["response"]).all():
            raise ValueError("Non-finite rank inputs")
        groups = compound_groups(d["condition_compounds"])
        if not np.array_equal(groups["ids"], d["compounds"]):
            raise ValueError("Compound ordering mismatch")
        # Float64 independently checks the archived float32 CPU/GPU computation.
        observed = d["task"].astype(np.float64) @ d["response"].astype(np.float64).T / 5865
        null_base = Path(null_root) / layer if null_root is not None else base
        plan = json.loads((null_base / "global_tail/analysis_plan.json").read_text())
        null = pd.read_csv(null_base / "global_tail/screen_max_null.tsv", sep="\t").both.to_numpy()
        if len(null) != 100000 or plan["n"] != len(null) or plan["seed"] != 20261024:
            raise ValueError("Permutation count or seed differs from the manuscript")
        errors, frames = {}, []
        for direction, sign in [("pos", 1), ("neg", -1)]:
            recalculated = max_over_compounds(sign * observed, groups)
            original = saved[direction + "_spearman_rho"]
            errors[direction] = float(np.max(np.abs(recalculated - original)))
            if errors[direction] > 5e-7:
                raise ValueError(f"Observed {layer}/{direction} correlations differ: {errors[direction]}")
            # Use the explicitly selected score/calibration arrays after the independent
            # observed-score check, avoiding tie changes from the float64 audit.
            p = screenwide_p(original, saved[direction + "_null_mean"], saved[direction + "_null_sd"], null)
            frames.append(pd.DataFrame({"target_id": np.repeat(d["targets"], 3250),
                "inchikey": np.tile(d["compounds"], expected), "direction": direction,
                "spearman_rho_" + layer: original.ravel(), "p_" + layer + "_screenwide": p.ravel(),
                "source_condition_index_" + layer: saved[direction + "_source_condition_index"].ravel()}))
        tables[layer] = pd.concat(frames, ignore_index=True)
        records[layer] = {"targets": expected, "strains": 5865, "compounds": 3250,
            "directional_hypotheses": expected * 3250 * 2, "max_absolute_correlation_error": errors,
            "permutations": len(null), "seed": plan["seed"]}
    keys = ["target_id", "inchikey", "direction"]
    matches = tables["match"].query("p_match_screenwide <= 0.05")
    combined = matches.merge(tables["residual"], on=keys, how="left", validate="one_to_one")
    selected = combined.query("p_residual_screenwide <= 0.05").sort_values(keys).reset_index(drop=True)
    selected.insert(0, "candidate_id", [f"YB-IC-{i:03d}" for i in range(1, len(selected) + 1)])
    reference = pd.read_csv(REFERENCE / "screenwide_residual_supported.tsv", sep="\t").sort_values(keys).reset_index(drop=True)
    if selected[keys].to_dict("records") != reference[keys].to_dict("records"):
        raise ValueError("Reconstructed candidate identities differ from manuscript")
    for a, b in [("p_match_screenwide", "p_match_screenwide"), ("p_residual_screenwide", "p_residual_screenwide"),
                 ("spearman_rho_match", "spearman_rho"), ("spearman_rho_residual", "spearman_rho_residual")]:
        tolerance = 5e-7 if input_root is not None and a.startswith('spearman_rho') else 1e-12
        if not np.allclose(selected[a], reference[b], rtol=0, atol=tolerance):
            raise ValueError(f"Reconstructed field differs: {a}")
    if (len(matches), len(selected), selected.target_id.nunique(), selected.inchikey.nunique()) != (5216, 25, 7, 18):
        raise ValueError("Manuscript result counts differ")
    combined.to_csv(output / "screenwide_associations.tsv", sep="\t", index=False)
    selected.to_csv(output / "reconstructed_candidates.tsv", sep="\t", index=False)
    report = {"status": "passed", "mode": "numerical-reconstruction", "layers": records,
        "match_associations": len(matches), "candidates": len(selected), "targets": 7, "compounds": 18,
        "recomputed": ["all observed correlations", "condition maximization", "screen-wide probabilities", "candidate intersection"],
        "frozen_inputs": (["rank-standardized profiles", "1000-permutation calibration"] if input_root is None else []) +
                         (["100000-permutation maxima"] if null_root is None else []),
        "fresh_permutations": null_root is not None, "fresh_calibration": input_root is not None,
        "fresh_training": False}
    write_json(output / "verification.json", report)
    return report


def profiles(output, threads=4):
    """Rebuild sequence-to-yeast task profiles with the manuscript checkpoint."""
    import pandas as pd
    require("reference")
    require("inference")
    output = new_output(output)
    env = dict(os.environ, OPENBLAS_NUM_THREADS=str(threads), OMP_NUM_THREADS=str(threads), MKL_NUM_THREADS=str(threads))
    commands = [[sys.executable, str(ROOT / "scripts/esm2_joint_task_export.py"), "--config",
                 str(ROOT / "configs/paper_tasks.json"), "--response-universe", str(ROOT / "raw/paper/strain_response.npz"),
                 "--output", str(output / "all_profiles")]]
    logs = output / "logs"; logs.mkdir()
    with (logs / "profile_export.log").open("w") as f:
        subprocess.run(commands[0], cwd=ROOT, env=env, stdout=f, stderr=subprocess.STDOUT, check=True)
    source = output / "all_profiles"
    scoped = output / "ion_channel_profiles"; scoped.mkdir()
    targets = pd.read_csv(REFERENCE / "target_scope.tsv", sep="\t").target_id.tolist()
    record = json.loads((source / "export_record.json").read_text())
    for target in targets:
        shutil.copyfile(source / f"yeast_task_{target}.tsv", scoped / f"yeast_task_{target}.tsv")
    record["targets"] = {t: record["targets"][t] for t in targets}
    record["n_targets"] = len(targets)
    write_json(scoped / "export_record.json", record)
    commands.append([sys.executable, str(ROOT / "scripts/task_axis_double_center.py"), "--input", str(scoped),
                     "--output", str(output / "corrected_profiles")])
    with (logs / "common_mode.log").open("w") as f:
        subprocess.run(commands[1], cwd=ROOT, env=env, stdout=f, stderr=subprocess.STDOUT, check=True)
    corrected = json.loads((output / "corrected_profiles/export_record.json").read_text())
    if corrected["n_targets"] != 338 or corrected["pool"] != 5865 or record["train_pairs"] != 510:
        raise ValueError("Rebuilt profile dimensions or training partition differ")
    reference = json.loads((REFERENCE / "dc_record.json").read_text())
    for name, value in reference["common_mode_variance"].items():
        if abs(corrected["common_mode_variance"][name] - value) > 1e-8:
            raise ValueError(f"Common-mode reconstruction differs: {name}")
    rank_check = None
    numeric = ROOT / "raw/paper/match/rank_inputs.npz"
    if numeric.exists():
        import numpy as np
        from scripts.product_execute_hiphop import load_task_relevance, rank_standardize_rows
        require("numerical")
        with np.load(numeric) as z:
            positions = {g: i for i, g in enumerate(z["orfs"])}
            rows = [load_task_relevance(output / "corrected_profiles" / f"yeast_task_{t}.tsv", positions, 5865)[0]
                    for t in z["targets"]]
            rebuilt = rank_standardize_rows(np.array(rows)).astype(np.float32)
            rank_check = float(np.max(np.abs(rebuilt - z["task"])))
        if rank_check > 1e-7:
            raise ValueError(f"Rebuilt strain ranks differ from the frozen screen: {rank_check}")
    report = {"status": "passed", "mode": "frozen-model-profile-reconstruction", "targets": 338, "strains": 5865,
              "training_pairs": 510, "members": record["members"], "rrf_k": record["rrf_k"],
              "common_mode_variance": corrected["common_mode_variance"], "commands": commands, "fresh_training": False,
              "rank_matrix_max_absolute_error": rank_check}
    write_json(output / "verification.json", report)
    return report
