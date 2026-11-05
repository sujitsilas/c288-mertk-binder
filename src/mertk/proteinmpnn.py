"""Pinned standard ProteinMPNN baseline with an explicit fixed-position design mask."""

import importlib.util
import json
import subprocess
from pathlib import Path

import numpy as np
import torch

ALPHABET = "ACDEFGHIKLMNPQRSTVWYX"


def load_proteinmpnn(device):
    spec = json.loads(Path("configs/external_tools.json").read_text())["proteinmpnn"]
    root = Path(spec["directory"]).resolve()
    if not root.exists():
        raise FileNotFoundError("Run scripts/setup_external.py --tool proteinmpnn first")
    revision = subprocess.check_output(
        ["git", "-C", str(root), "rev-parse", "HEAD"], text=True
    ).strip()
    dirty = subprocess.check_output(
        ["git", "-C", str(root), "status", "--porcelain", "--untracked-files=no"], text=True
    )
    if revision != spec["revision"] or dirty:
        raise ValueError("ProteinMPNN checkout differs from the pinned source")
    module_spec = importlib.util.spec_from_file_location(
        "pinned_proteinmpnn", root / "protein_mpnn_utils.py"
    )
    module = importlib.util.module_from_spec(module_spec)
    module_spec.loader.exec_module(module)
    checkpoint = torch.load(root / spec["checkpoint"], map_location="cpu", weights_only=True)
    model = module.ProteinMPNN(
        num_letters=21,
        node_features=128,
        edge_features=128,
        hidden_dim=128,
        num_encoder_layers=3,
        num_decoder_layers=3,
        augment_eps=0,
        k_neighbors=checkpoint["num_edges"],
    )
    model.load_state_dict(checkpoint["model_state_dict"], strict=True)
    return model.to(device).eval().requires_grad_(False), spec


def mpnn_inputs(example, device):
    binder, target = str(example["binder_sequence"]), str(example["target_sequence"])
    coordinates = []
    for chain in ("binder", "target"):
        if f"{chain}_oxygen" not in example:
            raise ValueError(
                "ProteinMPNN needs backbone O atoms. Regenerate examples with scripts/prepare_data.py, then cache structure tokens."
            )
        coordinates.append(
            np.concatenate(
                (example[f"{chain}_backbone"], example[f"{chain}_oxygen"][:, None, :]), axis=1
            )
        )
    xyz = torch.tensor(np.concatenate(coordinates), device=device, dtype=torch.float32)[None]
    valid = torch.isfinite(xyz).all(dim=(-1, -2))
    design = torch.tensor(
        np.concatenate((example["design_mask"], np.zeros(len(target), dtype=bool))), device=device
    )[None]
    if not design.any() or not valid[design].all():
        raise ValueError("Every design position must have N/CA/C/O coordinates")
    sequence = torch.tensor([ALPHABET.index(aa) for aa in binder + target], device=device)[None]
    # Deliberately remove every design label even though upstream sample only reads fixed labels.
    sequence[design] = ALPHABET.index("X")
    chain_ids = torch.tensor([1] * len(binder) + [2] * len(target), device=device)[None]
    residue_index = torch.tensor(
        list(range(len(binder))) + list(range(len(binder) + 100, len(binder) + 100 + len(target))),
        device=device,
    )[None]
    return {
        "X": torch.nan_to_num(xyz),
        "S_true": sequence,
        "mask": valid.float(),
        "chain_mask": design.float(),
        "chain_encoding_all": chain_ids,
        "residue_idx": residue_index,
        "design": design,
    }


@torch.no_grad()
def sample_proteinmpnn(model, example, device, seed, temperature):
    if temperature <= 0:
        raise ValueError("temperature must be positive")
    inputs = mpnn_inputs(example, device)
    design = inputs.pop("design")
    torch.manual_seed(seed)
    randn = torch.randn(design.shape, generator=torch.Generator().manual_seed(seed)).to(device)
    output = model.sample(
        **inputs,
        randn=randn,
        temperature=temperature,
        omit_AAs_np=np.array([0] * 20 + [1], dtype=np.float32),
        bias_AAs_np=np.zeros(21, dtype=np.float32),
        chain_M_pos=torch.ones_like(inputs["chain_mask"]),
        bias_by_res=torch.zeros((*design.shape, 21), device=device),
    )
    sequence = output["S"]
    if not torch.equal(sequence[~design], inputs["S_true"][~design]):
        raise RuntimeError("ProteinMPNN changed a fixed residue")
    length = len(str(example["binder_sequence"]))
    return "".join(ALPHABET[index] for index in sequence[0, :length].cpu().tolist())
