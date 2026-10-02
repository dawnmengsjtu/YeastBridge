#!/usr/bin/env python3
"""Check final-pipeline assets against the frozen checksums without running analysis."""
import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.project_io import ROOT, sha256, write_json


def inventory(mode="full"):
    entries = {}
    configs = ["configs/esm2_joint_tasks.json", "configs/v7_b2_joint_formal.json"]
    for relative in configs:
        cfg = json.loads((ROOT / relative).read_text())
        root = (ROOT / cfg.get("stage_root", ".")).resolve()
        for value in cfg.get("inputs", {}).values():
            if not isinstance(value, dict):
                continue
            for pk, hk in [
                ("path", "sha256"),
                ("path_abs", "sha256"),
                ("index_path", "index_sha256"),
                ("index_path_abs", "index_sha256"),
                ("gene_order_path", "gene_order_sha256"),
            ]:
                if pk in value:
                    entries[str((root / value[pk]).resolve())] = value.get(hk)
    # This matrix is a processed asset: its checksum cannot be recovered from an accession alone.
    for line in (ROOT / "MANIFEST.sha256").read_text().splitlines():
        digest, name = line.split(maxsplit=1)
        if name == "raw/tier1_response/strain_response.npz":
            entries[str(ROOT / name)] = digest
    entries[str(ROOT / "raw/tier1_response/compounds.tsv.gz")] = None
    results = []
    for name, expected in sorted(entries.items()):
        p = Path(name)
        actual = sha256(p) if p.is_file() else None
        status = (
            "missing"
            if actual is None
            else "hash-mismatch" if expected and actual != expected else "ok"
        )
        results.append(
            {
                "path": str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p),
                "status": status,
                "expected_sha256": expected,
                "actual_sha256": actual,
            }
        )
    return results


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()
    items = inventory()
    report = {"ready": all(x["status"] == "ok" for x in items), "assets": items}
    if args.output:
        write_json(args.output, report)
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0 if report["ready"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
