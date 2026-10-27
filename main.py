#!/usr/bin/env python3
"""YeastBridge manuscript export, asset checks and numerical reconstruction."""
import argparse
import json
import os
from pathlib import Path
import sys

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--mode', choices=['export', 'check', 'install', 'reproduce', 'profiles'], default='export')
    ap.add_argument('--output', type=Path)
    ap.add_argument('--level', choices=['reference', 'inference', 'numerical', 'all'], default='reference')
    ap.add_argument('--archive', type=Path, action='append')
    ap.add_argument('--threads', type=int, default=4)
    args = ap.parse_args()
    if args.threads < 1: ap.error('--threads must be positive')
    for name in ['OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS']:
        os.environ[name] = str(args.threads)
    from scripts import paper_release as paper
    if args.mode == 'export':
        print(paper.export(args.output))
    elif args.mode == 'check':
        rows = paper.inventory(args.level)
        report = {'ready': all(x['status'] == 'ok' for x in rows), 'level': args.level, 'assets': rows}
        if args.output: paper.write_json(args.output, report)
        print(json.dumps(report, indent=2))
        return 0 if report['ready'] else 2
    elif args.mode == 'install':
        if not args.archive: ap.error('--archive is required')
        for archive in args.archive: print(json.dumps(paper.install(archive)))
    elif args.mode == 'reproduce':
        print(json.dumps(paper.reproduce(args.output or paper.ROOT / 'outputs/numerical'), indent=2))
    elif args.mode == 'profiles':
        print(json.dumps(paper.profiles(args.output or paper.ROOT / 'outputs/profiles', args.threads), indent=2))
    return 0

if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except (ValueError, FileNotFoundError) as exc:
        print(str(exc), file=sys.stderr)
        raise SystemExit(2)
