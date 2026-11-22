"""Experimental, differentiable glycan excluded-volume loss for pinned BindCraft.

Sugars are transported by each attached Asn's N/CA/C frame. This is a steric
restraint, not glycan cofolding, an occupancy prediction, or an affinity score.
The ordinary filters and native-context checks remain mandatory.
"""

import json
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser, PDBParser

from mertk.provenance import file_hash

SUGARS = {"NAG", "BMA", "MAN", "FUC", "GAL", "NDG", "A2G", "SIA"}
MODE = "asn_frame_excluded_volume_v1"


def frames(backbone, xp=np):
    """Return orthonormal row bases from [..., N/CA/C, xyz]."""
    n, ca, c = backbone[..., 0, :], backbone[..., 1, :], backbone[..., 2, :]

    def norm(v):
        return v / xp.sqrt(xp.sum(v * v, axis=-1, keepdims=True) + 1e-8)

    x = norm(c - ca)
    z = norm(xp.cross(x, n - ca))
    y = xp.cross(z, x)
    return xp.stack([x, y, z], axis=-2)


def prepare_constraints(target_pdb, native_cif):
    target = PDBParser(QUIET=True).get_structure("crop", target_pdb)[0]["A"]
    residues = [r for r in target if r.id[0] == " "]
    numbers = [r.id[1] for r in residues]
    if numbers != list(range(93, 286)):
        raise ValueError("This guard has been reviewed only for the native A93-A285 crop")
    native = MMCIFParser(QUIET=True).get_structure("native", native_cif)[0]
    protein = native["A"]
    sites, local_atoms, atom_sites = [], [], []
    attachments = []
    for chain in native:
        sugars = [r for r in chain if r.resname in SUGARS]
        if not sugars:
            continue
        candidates = sorted(
            (float(np.linalg.norm(s["C1"].coord - r["ND2"].coord)), r.id[1])
            for s in sugars
            if s.resname == "NAG" and "C1" in s
            for r in protein
            if r.resname == "ASN" and "ND2" in r
        )
        if not candidates or candidates[0][0] > 2.2:
            raise ValueError(f"Cannot establish glycan attachment for chain {chain.id}")
        distance, number = candidates[0]
        if number not in numbers:
            continue
        if number in sites:
            raise ValueError("Ambiguous duplicate glycan attachment")
        atoms = np.array([a.coord for r in sugars for a in r if a.element != "H"])
        triad = np.array([protein[number][name].coord for name in ("N", "CA", "C")])
        target_triad = np.array([target[number][name].coord for name in ("N", "CA", "C")])
        if not np.allclose(triad, target_triad, rtol=0, atol=0.0006):
            raise ValueError("Native and crop frames differ")
        # PDB rounds to 0.001 A; retain exact native sugar positions in the crop frame.
        triad = target_triad
        local_atoms.extend(((atoms - triad[1]) @ frames(triad).T).tolist())
        atom_sites.extend([len(sites)] * len(atoms))
        sites.append(number)
        attachments.append(
            {
                "native_asn": number,
                "glycan_chain": chain.id,
                "covalent_distance_A": distance,
                "heavy_atoms": len(atoms),
            }
        )
    sequons = [
        r.id[1]
        for i, r in enumerate(residues[:-2])
        if r.resname == "ASN"
        and residues[i + 1].resname != "PRO"
        and residues[i + 2].resname in ("SER", "THR")
    ]
    if sorted(sites) != [114, 170, 207, 215, 234] or sorted(sequons) != sorted(sites):
        raise ValueError("Unexpected native crop glycosylation mapping")
    return {
        "mode": MODE,
        "target_sha256": file_hash(target_pdb),
        "native_cif_sha256": file_hash(native_cif),
        "target_length": len(residues),
        "native_numbers": numbers,
        "attachments": attachments,
        "anchor_indices": [numbers.index(n) for n in sites],
        "sequon_indices": [numbers.index(n) for n in sequons],
        "glycan_local_xyz": local_atoms,
        "glycan_atom_sites": atom_sites,
        "heavy_atom_clearance_A": 3.5,
        "sequon_ca_exclusion_A": 15.0,
        "weights": {"glycan_volume": 10.0, "sequon_volume": 2.0},
        "limitations": "Soft experimental restraints; one observed insect-cell glycoform and conservative sequon spheres. No full glycan ensemble, no glycans cofolded by AF2. Post-run native checks required.",
    }


def losses(positions, mask, spec, xp=np):
    """AF atom37 order: N=0, CA=1, C=2. Mask absent atoms before reduction."""
    nt = spec["target_length"]
    if positions.shape[0] <= nt or positions.shape[1:] != (37, 3):
        raise ValueError("Expected target followed by a binder in atom37 order")
    anchor = positions[xp.asarray(spec["anchor_indices"]), :3]
    basis = frames(anchor, xp)
    indices = xp.asarray(spec["glycan_atom_sites"])
    local = xp.asarray(spec["glycan_local_xyz"])
    glycan = xp.einsum("ai,aij->aj", local, basis[indices]) + anchor[indices, 1]
    binder = positions[nt:]
    distances = xp.sqrt(xp.sum((binder[:, :, None] - glycan[None, None]) ** 2, axis=-1) + 1e-8)
    radius = spec["heavy_atom_clearance_A"]
    penalty = xp.maximum(0, 1 - distances.min(axis=-1) / radius) ** 2
    per_residue = xp.max(xp.where(mask[nt:] > 0.5, penalty, 0), axis=-1)
    # Averaging the worst fifth reduces dilution as binder length increases.
    count = max(1, (binder.shape[0] + 4) // 5)
    observed = xp.sort(per_residue)[-count:].mean()
    sequons = positions[xp.asarray(spec["sequon_indices"]), 1]
    ca_distance = xp.sqrt(xp.sum((binder[:, 1, None] - sequons[None]) ** 2, axis=-1) + 1e-8)
    sequon_penalty = (
        xp.maximum(0, 1 - ca_distance.min(axis=-1) / spec["sequon_ca_exclusion_A"]) ** 2
    )
    sequon_penalty = xp.where(mask[nt:, 1] > 0.5, sequon_penalty, 0)
    result = {"glycan_volume": observed, "sequon_volume": xp.sort(sequon_penalty)[-count:].mean()}
    if "membrane" in spec:
        membrane = spec["membrane"]
        outer = positions[xp.asarray(membrane["outer_indices"]), 1].mean(axis=0)
        inner = positions[xp.asarray(membrane["inner_indices"]), 1].mean(axis=0)
        axis = inner - outer
        axis = axis / xp.sqrt(xp.sum(axis**2) + 1e-8)
        origin = positions[membrane["origin_index"], 1]
        depth = xp.sum((binder - origin) * axis, axis=-1)
        margin = membrane["extracellular_margin_A"]
        overlap = xp.maximum(0, (depth + margin) / margin) ** 2
        overlap = xp.max(xp.where(mask[nt:] > 0.5, overlap, 0), axis=-1)
        result["membrane_volume"] = xp.sort(overlap)[-count:].mean()
    return result


def attach(model, spec, target_pdb):
    import jax.numpy as jnp

    if spec["mode"] != MODE or file_hash(target_pdb) != spec["target_sha256"]:
        raise ValueError("Glycan guard target or method mismatch")
    if model._target_len != spec["target_length"]:
        raise ValueError("ColabDesign target residue mapping mismatch")

    def callback(inputs, outputs):
        sm = outputs["structure_module"]
        return losses(sm["final_atom_positions"], sm["final_atom_mask"], spec, jnp)

    model._callbacks["model"]["loss"].append(callback)
    model.opt["weights"].update(spec["weights"])
    print(
        "PROJECT_GLYCAN_GUARD_ATTACHED",
        json.dumps(
            {
                "mode": MODE,
                "sites": [a["native_asn"] for a in spec["attachments"]],
                "weights": spec["weights"],
            }
        ),
        flush=True,
    )
    return callback


def install(functions, constraints_path):
    """Wrap only hallucination; upstream source and final validation are unchanged."""
    spec = json.loads(Path(constraints_path).read_text())
    original = functions.binder_hallucination
    namespace = original.__globals__
    factory = namespace["mk_afdesign_model"]

    def guarded_factory(*args, **kwargs):
        model = factory(*args, **kwargs)
        prep = model.prep_inputs

        def guarded_prep(*args, **kwargs):
            result = prep(*args, **kwargs)
            attach(model, spec, kwargs["pdb_filename"])
            return result

        model.prep_inputs = guarded_prep
        return model

    def guarded_hallucination(*args, **kwargs):
        namespace["mk_afdesign_model"] = guarded_factory
        try:
            return original(*args, **kwargs)
        finally:
            namespace["mk_afdesign_model"] = factory

    functions.binder_hallucination = guarded_hallucination
