"""Prepare reproducible BindCraft jobs from an audited target construct."""

import json
import re
import shlex
import shutil
import subprocess
from pathlib import Path

from mertk.provenance import file_hash
from mertk.structures import read_design_chain


def validate_hotspots(hotspots, chain_id, mapping):
    if hotspots is None:
        return None
    if any(r["insertion_code"] for r in mapping):
        raise ValueError("Hotspot syntax does not support insertion codes; map the construct first")
    positions = {r["author_number"] for r in mapping}
    for token in hotspots.split(","):
        match = re.fullmatch(r"([A-Za-z]?)(\d+)(?:-(\d+))?", token.strip())
        if not match or (match[1] and match[1] != chain_id):
            raise ValueError(
                "Hotspots must be residue numbers/ranges on the selected chain, e.g. A25,A30-33"
            )
        start, end = int(match[2]), int(match[3] or match[2])
        if start > end or not set(range(start, end + 1)) <= positions:
            raise ValueError(f"Hotspot residue is absent: {token}")
    return hotspots


def prepare_job(
    target_pdb,
    target_chain,
    output,
    lengths=(65, 120),
    designs=10,
    hotspots=None,
    name="MerTK",
):
    output = Path(output).resolve()
    target_pdb = Path(target_pdb).resolve()
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(f"Choose an empty job directory: {output}")
    if not re.fullmatch(r"[A-Za-z0-9_-]+", name):
        raise ValueError("Use an alphanumeric binder name")
    if not 1 <= lengths[0] <= lengths[1] or designs < 1:
        raise ValueError("Invalid binder length range or design count")
    target, mapping = read_design_chain(target_pdb, target_chain)
    hotspots = validate_hotspots(hotspots, target_chain, mapping)
    spec = json.loads(Path("configs/external_tools.json").read_text())["bindcraft"]
    root = Path(spec["directory"]).resolve()
    if not root.exists():
        raise FileNotFoundError("Run scripts/setup_external.py --tool bindcraft")
    revision = subprocess.check_output(
        ["git", "-C", str(root), "rev-parse", "HEAD"], text=True
    ).strip()
    dirty = subprocess.check_output(
        ["git", "-C", str(root), "status", "--porcelain", "--untracked-files=no"], text=True
    )
    if revision != spec["revision"] or dirty:
        raise ValueError("BindCraft checkout differs from the pinned revision")
    output.mkdir(parents=True, exist_ok=True)
    shutil.copy2(target_pdb, output / "target.pdb")
    filters = root / "settings_filters/default_filters.json"
    advanced = root / "settings_advanced/default_4stage_multimer.json"
    shutil.copy2(filters, output / "filters.json")
    shutil.copy2(advanced, output / "advanced.json")
    settings = {
        "design_path": str(output / "designs"),
        "binder_name": name,
        "starting_pdb": str(output / "target.pdb"),
        "chains": target_chain,
        "target_hotspot_residues": hotspots,
        "lengths": list(lengths),
        "number_of_final_designs": designs,
    }
    (output / "target.json").write_text(json.dumps(settings, indent=2) + "\n")
    command = [
        "conda",
        "run",
        "--no-capture-output",
        "-n",
        "BindCraft",
        "python",
        "-u",
        str(root / "bindcraft.py"),
        "--settings",
        str(output / "target.json"),
        "--filters",
        str(output / "filters.json"),
        "--advanced",
        str(output / "advanced.json"),
    ]
    (output / "run.sh").write_text(
        "#!/usr/bin/env bash\nset -euo pipefail\ncd "
        + shlex.quote(str(root))
        + "\nexec "
        + shlex.join(command)
        + "\n"
    )
    provenance = {
        "tool": spec,
        "target_sha256": file_hash(output / "target.pdb"),
        "target_chain": target_chain,
        "target_sequence": target.sequence,
        "residue_mapping": mapping,
        "settings": settings,
        "filters_sha256": file_hash(output / "filters.json"),
        "advanced_sha256": file_hash(output / "advanced.json"),
        "launched": False,
    }
    (output / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    return output
