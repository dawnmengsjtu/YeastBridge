#!/usr/bin/env python3
"""Rebuild rank inputs and 1,000-permutation calibration from computed profiles."""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def residual_profiles(profiles, destination):
    """Apply the manuscript's fixed profile families and original projection."""
    import numpy as np
    import pandas as pd
    from scripts.paper_release import REFERENCE, new_output, write_json
    destination = new_output(destination)
    membership = pd.read_csv(REFERENCE / 'family_membership.tsv', sep='\t').set_index('target_id').profile_family
    targets = sorted(membership.index)
    record = json.loads((profiles / 'export_record.json').read_text())
    if record['tail_strains']:
        raise ValueError('The manuscript profiles must cover exactly the represented strains')
    orfs, rows = None, []
    for target in targets:
        table = pd.read_csv(profiles / f'yeast_task_{target}.tsv', sep='\t').set_index('yeast_gene').ensemble_score
        if orfs is None:
            orfs = table.index.to_numpy()
        if set(table.index) != set(orfs):
            raise ValueError('Inconsistent profile strain universe')
        rows.append(table.reindex(orfs).to_numpy())
    scores = np.asarray(rows, dtype=np.float64)
    non_estimable = []
    for i, target in enumerate(targets):
        members = np.flatnonzero(membership.reindex(targets).to_numpy() == membership[target])
        if len(members) < 2:
            continue
        mean = scores[members].mean(axis=0)
        score = scores[i]
        residual = score - float(score @ mean) / float(mean @ mean) * mean
        fraction = float(residual @ residual) / float(score @ score) if float(score @ score) > 0 else np.nan
        if not np.isfinite(fraction) or fraction <= 1e-12:
            residual[:] = 0
            non_estimable.append(target)
        order = np.argsort(-residual, kind='stable')
        with (destination / f'yeast_task_{target}.tsv').open('w') as handle:
            handle.write('rank\tyeast_gene\tensemble_score\n')
            for rank, j in enumerate(order, 1):
                handle.write(f'{rank}\t{orfs[j]}\t{residual[j]:.9f}\n')
    expected = json.loads((REFERENCE / 'residual_export_record.json').read_text())['non_estimable_targets']
    if non_estimable != expected:
        raise ValueError('Residual estimability differs from the manuscript')
    write_json(destination / 'export_record.json', dict(tail_strains=[], non_estimable_targets=non_estimable,
                                                      source_task_directory=str(profiles)))
    return non_estimable


def run(profiles, response_path, output, device, threads):
    import numpy as np
    import pandas as pd
    import torch
    from scripts import paper_release as paper
    from scripts.product_execute_hiphop import load_task_relevance, rank_standardize_rows, compound_groups, run_v2_spearman
    paper.require('all')
    if paper.digest(response_path) != paper.digest(ROOT / 'raw/paper/strain_response.npz'):
        raise ValueError('Response differs from the manuscript measurement snapshot')
    output = paper.new_output(output)
    non_estimable = residual_profiles(profiles, output / 'residual_profiles')
    torch.set_num_threads(threads)
    os.environ['YB_MATRIX_DEVICE'] = device
    os.environ['YB_PERMUTATION_BATCH'] = '1'
    os.environ['YB_STRAIN_STRATA'] = str(ROOT / 'raw/paper/published_strain_universe.tsv')
    with np.load(response_path, allow_pickle=True) as z:
        response = dict(z)
    orfs = response['strain_orfs'].astype(str)
    positions = {g:i for i,g in enumerate(orfs)}
    keys = response['compound_inchikeys'].astype(str)
    conditions = np.flatnonzero((~response['is_vehicle'].astype(bool)) & (keys != ''))
    keys = keys[conditions]
    doses = response['doses'].astype(str)[conditions]
    sensitivity = np.asarray(response['z_score'], dtype=np.float32)[:, conditions]
    groups = compound_groups(keys)
    reports = {}
    for layer, folder in [('match', profiles), ('residual', output / 'residual_profiles')]:
        dest = paper.new_output(output / layer)
        files = sorted(folder.glob('yeast_task_*.tsv'))
        if layer == 'residual':
            files = [p for p in files if p.stem.removeprefix('yeast_task_') not in non_estimable]
        targets = np.array([p.stem.removeprefix('yeast_task_') for p in files])
        values, supports = [], []
        for path in files:
            row, count = load_task_relevance(path, positions, len(orfs))
            values.append(row)
            supports.append(np.ones(len(orfs), dtype=bool) if count == len(orfs) else row > row.min())
        if not all(np.array_equal(supports[0], s) for s in supports):
            raise ValueError('Target coverage differs')
        support = np.flatnonzero(supports[0])
        rank_data = dict(task=rank_standardize_rows(np.stack(values)[:,support]).astype(np.float32),
                         response=rank_standardize_rows(sensitivity[support].T).astype(np.float32),
                         targets=targets, compounds=groups['ids'], orfs=orfs[support],
                         conditions=conditions, condition_compounds=keys,
                         support=np.arange(len(support)), group_order=groups['order'],
                         group_starts=groups['starts'], group_counts=groups['counts'])
        with np.load(ROOT / 'raw/paper' / layer / 'rank_inputs.npz') as frozen:
            checks = {key:bool(np.array_equal(value, frozen[key])) for key,value in rank_data.items()}
        if not all(checks.values()):
            paper.write_json(dest / 'rank_verification.json', checks)
            raise ValueError(f'{layer} rank inputs differ: {checks}')
        np.savez_compressed(dest / 'rank_inputs.npz', **rank_data)
        arrays, errors = {}, {}
        index = pd.MultiIndex.from_product([targets,groups['ids']], names=['target_id','inchikey'])
        with np.load(ROOT / 'raw/paper' / layer / 'screen_arrays.npz') as z:
            reference = dict(z)
        for direction, sign in [('pos',1),('neg',-1)]:
            table, audit = run_v2_spearman(files, positions, sign*sensitivity, keys, doses, conditions,
                                          1000, np.random.default_rng(42), target_batch_size=64)
            table = table.set_index(['target_id','inchikey']).reindex(index)
            for column in ['spearman_rho','null_mean','null_sd','source_condition_index','emp_p']:
                name = direction + '_' + column
                # Original arrays were serialized via a pandas TSV reader.
                # Preserve values computed here and quantify serialization roundoff.
                value = table[column].to_numpy().reshape(len(targets),len(groups['ids']))
                value = value.astype(np.int32 if column=='source_condition_index' else np.float64)
                arrays[name] = value
                errors[name] = float(np.max(np.abs(value-reference[name])))
                tolerance = (0 if column=='source_condition_index' else
                             5e-7 if column=='spearman_rho' else 1e-12)
                if errors[name] > tolerance:
                    raise ValueError(f'{layer}/{name} calibration differs: {errors[name]}')
            paper.write_json(dest / f'{direction}_execution.json', audit)
        np.savez_compressed(dest / 'screen_arrays.npz', **arrays)
        reports[layer] = dict(status='passed', rank_arrays_equal=checks, max_absolute_errors=errors,
                              permutations=1000, seed=42, fresh_calibration=True)
        paper.write_json(dest / 'verification.json', reports[layer])
    reports['inputs'] = dict(response_sha256=paper.digest(response_path),
                            script_sha256=paper.digest(Path(__file__)), profile_directory=str(profiles))
    paper.write_json(output / 'verification.json', reports)
    print(json.dumps(reports, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--profiles', type=Path, required=True, help='Regenerated corrected_profiles directory')
    parser.add_argument('--response', type=Path, default=ROOT / 'raw/paper/strain_response.npz')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--device', default='cuda:0')
    parser.add_argument('--threads', type=int, default=4)
    args = parser.parse_args()
    if args.threads < 1:
        parser.error('Thread count must be positive')
    for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:
        os.environ[key] = str(args.threads)
    run(args.profiles, args.response, args.output, args.device, args.threads)

if __name__ == '__main__':
    main()
