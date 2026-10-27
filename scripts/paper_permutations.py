#!/usr/bin/env python3
"""Regenerate the manuscript's assay-stratified, screen-wide maximum nulls.

The permutation and reduction operations preserve those in
strengthening_mapping_controls_20261011.py. Calibration arrays are checked
separately; no stored maximum-null values enter the computation.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def permute_within_assays(rng, labels):
    import numpy as np
    permutation = np.arange(len(labels))
    for label in sorted(set(labels)):
        indices = np.flatnonzero(labels == label)
        permutation[indices] = rng.permutation(indices)
    return permutation


def generate(layer, output, permutations, seed, device, threads, calibration=None):
    import numpy as np
    import pandas as pd
    import torch
    from scripts import paper_release as paper
    from scripts.product_execute_hiphop import compound_groups

    torch.set_num_threads(threads)
    torch.backends.cuda.matmul.allow_tf32 = False
    torch.backends.cudnn.allow_tf32 = False
    base = (Path(calibration) if calibration is not None else ROOT / 'raw/paper') / layer
    with np.load(base / 'rank_inputs.npz', allow_pickle=False) as z:
        d = dict(z)
    with np.load(base / 'screen_arrays.npz', allow_pickle=False) as z:
        s = dict(z)
    strata = ROOT / 'raw/paper/published_strain_universe.tsv'
    frame = pd.read_csv(strata, sep='\t').set_index('orf')
    if not frame.index.is_unique:
        raise ValueError('Duplicate strain annotation')
    labels = frame.assay.reindex(d['orfs'])
    if labels.isna().any():
        raise ValueError('Missing assay membership')
    labels = labels.to_numpy()
    x = torch.tensor(d['task'], device=device)
    y = torch.tensor(d['response'].T, device=device)
    groups = compound_groups(d['condition_compounds'])
    if not np.array_equal(groups['ids'], d['compounds']):
        raise ValueError('Compound ordering mismatch')
    order = torch.tensor(groups['order'], device=device)
    starts = torch.tensor(groups['starts'], device=device)
    multiple = [(i, int(a), int(n)) for i, (a, n) in enumerate(
        zip(groups['starts'], groups['counts'])) if n > 1]

    def reduce_scores(a):
        a = a[:, order]
        pos, neg = a[:, starts].clone(), -a[:, starts].clone()
        for i, start, count in multiple:
            pos[:, i] = a[:, start:start + count].amax(1)
            neg[:, i] = (-a[:, start:start + count]).amax(1)
        return pos, neg

    means = {di: torch.tensor(s[di + '_null_mean'], dtype=torch.float32, device=device)
             for di in ('pos', 'neg')}
    sds = {di: torch.tensor(s[di + '_null_sd'], dtype=torch.float32, device=device)
           for di in ('pos', 'neg')}
    if any(not torch.all(torch.isfinite(sds[di]) & (sds[di] > 0)) for di in sds):
        raise ValueError('Invalid calibration scale')
    strains = x.shape[1]
    out = paper.new_output(output / layer / 'global_tail')
    plan = dict(n=permutations, seed=seed, n_targets=len(d['targets']),
                n_compounds=len(d['compounds']), strata_sha256=paper.digest(strata),
                algorithm='Shared strain permutation within HIP/HOP; dose maximization; maximum across both directions',
                fresh_permutations=True, device=str(device), torch=torch.__version__,
                numpy=np.__version__, tf32=False, threads=threads,
                fresh_calibration=calibration is not None,
                rank_inputs_sha256=paper.digest(base / 'rank_inputs.npz'),
                calibration_sha256=paper.digest(base / 'screen_arrays.npz'))
    paper.write_json(out / 'analysis_plan.json', plan)
    rng = np.random.default_rng(seed)
    values = []
    started = time.monotonic()
    with (out / 'screen_max_null.tsv').open('w') as handle, torch.no_grad():
        handle.write('permutation\tpos\tneg\tboth\n')
        for i in range(permutations):
            perm = torch.tensor(permute_within_assays(rng, labels), device=device)
            maxima = [float(((a - means[di]) / sds[di]).max()) for di, a in
                      zip(('pos', 'neg'), reduce_scores(x[:, perm] @ y / strains))]
            values.append([*maxima, max(maxima)])
            handle.write(f'{i + 1}\t{maxima[0]}\t{maxima[1]}\t{max(maxima)}\n')
            if (i + 1) % 1000 == 0:
                handle.flush()
                print(json.dumps(dict(layer=layer, completed=i + 1, total=permutations,
                                      seconds=round(time.monotonic() - started, 1))), flush=True)
    reference = pd.read_csv(ROOT / 'raw/paper' / layer / 'global_tail/screen_max_null.tsv', sep='\t')
    error = None
    if seed == 20261024 and permutations <= len(reference):
        error = float(np.max(np.abs(np.asarray(values) - reference[['pos', 'neg', 'both']].iloc[:permutations].to_numpy())))
    report = dict(**plan, seconds=time.monotonic() - started,
                  archived_null_max_absolute_error=error,
                  float64_value_match=error == 0.0,
                  float32_values_equal=bool(np.array_equal(np.asarray(values,dtype=np.float32),
                      reference[['pos','neg','both']].iloc[:permutations].to_numpy(dtype=np.float32)))
                      if seed==20261024 and permutations<=len(reference) else None,
                  null_sha256=paper.digest(out / 'screen_max_null.tsv'),
                  script_sha256=paper.digest(Path(__file__)))
    paper.write_json(out / 'verification.json', report)
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--layer', choices=['match', 'residual', 'both'], default='both')
    p.add_argument('--permutations', type=int, default=100000)
    p.add_argument('--seed', type=int, default=20261024)
    p.add_argument('--device', default='cuda:0')
    p.add_argument('--threads', type=int, default=4)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--calibration', type=Path, help='Fresh output from paper_calibration.py')
    a = p.parse_args()
    if a.permutations < 1 or a.threads < 1:
        p.error('Permutation and thread counts must be positive')
    for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
        os.environ[key] = str(a.threads)
    from scripts import paper_release as paper
    paper.require('all')
    if a.calibration:
        verification = json.loads((a.calibration / 'verification.json').read_text())
        for layer in ['match','residual']:
            if verification[layer]['status'] != 'passed' or not verification[layer]['fresh_calibration']:
                raise ValueError('Calibration has not passed verification')
    output = paper.new_output(a.output)
    layers = ['match', 'residual'] if a.layer == 'both' else [a.layer]
    results = {layer: generate(layer, output, a.permutations, a.seed, a.device, a.threads,a.calibration)
               for layer in layers}
    if a.layer == 'both' and a.permutations == 100000 and a.seed == 20261024:
        results['screening'] = paper.reproduce(output / 'screening', null_root=output,input_root=a.calibration)
    paper.write_json(output / 'verification.json', results)
    print(json.dumps(results, indent=2))


if __name__ == '__main__':
    main()
