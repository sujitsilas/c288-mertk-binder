#!/usr/bin/env python3
"""Download pinned ESM3 training weights; no accelerator or ESM install needed."""

import argparse
import hashlib
import json
from pathlib import Path

from huggingface_hub import HfApi, snapshot_download


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path("configs/weights.json"))
    parser.add_argument("--output", type=Path, default=Path("weights/esm3"))
    parser.add_argument(
        "--dry-run", action="store_true", help="List sizes without downloading weights"
    )
    args = parser.parse_args()
    config = json.loads(args.config.read_text())
    info = HfApi().model_info(config["repo_id"], revision=config["revision"], files_metadata=True)
    selected = [f for f in info.siblings if f.rfilename in config["allow_patterns"]]
    if len(selected) != len(config["allow_patterns"]):
        raise RuntimeError("The pinned repository is missing a requested file")
    total = sum(f.size or 0 for f in selected)
    print(f"{info.id}@{info.sha}: {total / 1e9:.2f} GB")
    for file in selected:
        print(f"  {file.rfilename}: {(file.size or 0) / 1e6:.1f} MB")
    if args.dry_run:
        return
    args.output.mkdir(parents=True, exist_ok=True)
    snapshot_download(
        repo_id=config["repo_id"],
        revision=info.sha,
        allow_patterns=config["allow_patterns"],
        local_dir=args.output,
        max_workers=2,
    )
    manifest = {"repo_id": info.id, "revision": info.sha, "files": []}
    for file in selected:
        path = args.output / file.rfilename
        with path.open("rb") as handle:
            digest = hashlib.file_digest(handle, "sha256").hexdigest()
        expected = file.lfs.sha256 if file.lfs is not None else None
        if path.stat().st_size != file.size or (expected and digest != expected):
            raise RuntimeError(f"Checksum or size mismatch: {path}")
        manifest["files"].append({"path": file.rfilename, "size": file.size, "sha256": digest})
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Verified weights saved to {args.output}; keep this directory out of git.")


if __name__ == "__main__":
    main()
