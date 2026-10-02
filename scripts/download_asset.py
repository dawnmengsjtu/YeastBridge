"""Checksum-verified downloads from the project's versioned asset manifest."""

import json
import os
from pathlib import Path
import tempfile
import urllib.request

from scripts.project_io import ROOT, sha256


def download(name):
    manifest = json.loads((ROOT / "data/asset_downloads.json").read_text())
    item = manifest["assets"][name]
    destination = (ROOT / item["path"]).resolve()
    if not destination.is_relative_to(ROOT):
        raise ValueError("Asset destination escapes repository")
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and sha256(destination) == item["sha256"]:
        print(f'Already verified: {item["path"]}', flush=True)
        return destination
    descriptor, temporary = tempfile.mkstemp(
        prefix=".download-", dir=destination.parent
    )
    try:
        request = urllib.request.Request(
            item["url"], headers={"User-Agent": "YeastBridge-asset-installer"}
        )
        with os.fdopen(descriptor, "wb") as f, urllib.request.urlopen(
            request, timeout=120
        ) as response:
            while chunk := response.read(8 * 1024 * 1024):
                f.write(chunk)
        if sha256(temporary) != item["sha256"]:
            raise ValueError(
                f"{name}: checksum mismatch; downloaded asset was not installed"
            )
        os.replace(temporary, destination)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)
    print(f'Downloaded and verified: {item["path"]}', flush=True)
    return destination
