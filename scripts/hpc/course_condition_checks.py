"""Check native identities and conditions independently of a parser's residue namespace."""

import numpy as np


def audit_arrays(*, chain_id, res_id, res_name, atom_name, coord, known, sequence, hot, reference):
    target = chain_id == "A"
    ca = np.flatnonzero(target & (atom_name == "CA"))
    residues = reference["residues"]
    expected_ids = [r["author_number"] for r in residues]
    if len(ca) != len(residues) or len(set(res_id[ca])) != len(ca):
        raise ValueError("Wrong target residue extent or duplicate CA residue identifiers")
    if list(res_name[ca]) != [r["residue_name"] for r in residues]:
        raise ValueError("Target residue identities/order changed")
    # A materialized contig can be renumbered. Map ordered, verified native residues
    # to author numbers; never assume author IDs equal the installed parser's IDs.
    to_author = dict(zip(map(int, res_id[ca]), expected_ids, strict=True))
    if set(map(int, res_id[target])) != set(to_author):
        raise ValueError("Unexpected target residue identifiers")
    backbone = target & np.isin(atom_name, ["N", "CA", "C", "O"])
    expected_atoms = {
        (r["author_number"], name): np.asarray(point)
        for r in residues
        for name, point in r["backbone"].items()
    }
    indices = np.flatnonzero(backbone)
    keys = [(to_author[int(res_id[i])], str(atom_name[i])) for i in indices]
    if len(keys) != len(expected_atoms) or set(keys) != set(expected_atoms):
        raise ValueError("Wrong target backbone atom extent")
    actual = coord[indices]
    native = np.asarray([expected_atoms[key] for key in keys])
    if not np.isfinite(actual).all():
        raise ValueError("Non-finite target backbone coordinates")
    x, y = actual - actual.mean(0), native - native.mean(0)
    u, _, vt = np.linalg.svd(x.T @ y)
    sign = np.linalg.det(u @ vt)
    rotation = u @ np.diag([1, 1, sign]) @ vt
    deviation = np.linalg.norm(x @ rotation - y, axis=1)
    if deviation.max() > 0.01:
        raise ValueError("Native target backbone geometry changed")
    if not np.array_equal(known, backbone):
        raise ValueError("Internal coordinate mask changed")
    if not np.array_equal(sequence, target):
        raise ValueError("Target-only sequence conditioning changed")
    if (hot & ~target).any() or not np.isfinite(coord[hot]).all():
        raise ValueError("Hotspot outside target or with non-finite coordinates")
    observed_hot = [(to_author[int(res_id[i])], str(atom_name[i])) for i in np.flatnonzero(hot)]
    expected_hot = {tuple(item) for item in reference["hotspot_atoms"]}
    if len(observed_hot) != len(expected_hot) or set(observed_hot) != expected_hot:
        raise ValueError("Hotspot atom identities changed")
    return {
        "target_residues": len(ca),
        "conditioned_backbone_atoms": len(indices),
        "hotspot_residues": sorted({r for r, _ in expected_hot}),
        "hotspot_atoms": len(expected_hot),
        "parsed_to_author_residue_map": to_author,
        "max_backbone_deviation_A": float(deviation.max()),
        "glycans_present_in_generation": False,
        "reference_context_review_required": True,
    }
