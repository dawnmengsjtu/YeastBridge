#!/usr/bin/env python3
"""Run the recovered B2 training recipe and preserve stdout, inputs and runtime metadata."""
import argparse
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

from scripts.project_io import ROOT, revision, sha256, write_json


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--config", type=Path, default=ROOT / "configs/b2_training.json")
    p.add_argument("--output", required=True, type=Path)
    p.add_argument("--device", default="cuda", choices=["cuda", "cpu"])
    p.add_argument("--smoke", action="store_true")
    p.add_argument("--smoke-batch", type=int, default=2)
    args = p.parse_args()
    output = args.output.resolve()
    if output.exists() and any(output.iterdir()):
        p.error("Output must be a new or empty directory")
    cfg = json.loads(args.config.read_text())
    for key, asset in cfg["required_assets"].items():
        path = ROOT / asset["path"]
        if not path.is_file() or sha256(path) != asset["sha256"]:
            p.error(f"Missing or changed asset: {key}. See models/TRAINING.md")
    output.mkdir(parents=True, exist_ok=True)
    script = ROOT / "scripts/scf_finetune_yeast.py"
    command = [
        sys.executable,
        "-u",
        str(script),
        "--config",
        str(args.config.resolve()),
        "--output",
        str(output / "model"),
        "--device",
        args.device,
    ]
    if args.smoke:
        command += ["--smoke", "--smoke-batch", str(args.smoke_batch)]
    record = {
        "mode": "training-smoke" if args.smoke else "training",
        "status": "running",
        "code_revision": revision(),
        "source_sha256": sha256(script),
        "seed": cfg["seed"],
        "python": sys.version,
        "platform": platform.platform(),
        "packages": {
            x: importlib.metadata.version(x)
            for x in ["torch", "numpy", "einops", "local-attention"]
        },
        "required_assets": cfg["required_assets"],
        "command": command,
        "cuda_visible_devices": os.environ.get("CUDA_VISIBLE_DEVICES"),
    }
    try:
        record["gpu"] = subprocess.check_output(
            [
                "nvidia-smi",
                "--query-gpu=name,driver_version,memory.total",
                "--format=csv,noheader",
            ],
            text=True,
        ).splitlines()
    except (FileNotFoundError, subprocess.CalledProcessError):
        record["gpu"] = []
    write_json(output / "run.json", record)
    start = time.time()
    with (output / "train.log").open("w") as log:
        proc = subprocess.Popen(
            command,
            cwd=ROOT,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        for line in proc.stdout:
            print(line, end="", flush=True)
            log.write(line)
            log.flush()
        code = proc.wait()
    record.update(
        status="complete" if code == 0 else "failed",
        returncode=code,
        seconds=round(time.time() - start, 3),
        log_sha256=sha256(output / "train.log"),
    )
    write_json(output / "run.json", record)
    raise SystemExit(code)


if __name__ == "__main__":
    main()
