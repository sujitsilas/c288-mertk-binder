#!/usr/bin/env python3
"""Author fixed-native-context RFD4 repair inputs and screen actual generated outputs."""

import argparse
import hashlib
import json
from pathlib import Path


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def author(suite):
    import numpy as np
    import rfproteina.conditions  # noqa: F401
    from atomworks.io import parse
    from atomworks.io.transforms.atom_array import (
        add_id_and_entity_annotations,
        add_iid_annotations,
    )
    from atomworks.io.utils.annotator import ensure_annotations
    from atomworks.io.utils.atom_array_plus import as_atom_array_plus, concatenate_atom_array_plus
    from atomworks.io.utils.io_utils import to_cif_file
    from atomworks.ml.conditions import C_CRD, C_SEQ
    from rfproteina.conditions import C_HOT
    from rfproteina.conditions_reporting import (
        check_index_group_molecule_consistency,
        check_index_group_res_id_gaps,
        run_pipeline_validation,
        sync_index_groups,
    )
    from rfproteina.utils.parse_presets import ANNOTATED_CIF, DISTILLATION

    record = json.loads((suite / "manifest.json").read_text())
    native = as_atom_array_plus(
        parse(str(suite / "inputs/37NT.cif"), config=DISTILLATION)["assemblies"]["1"][0]
    )
    ensure_annotations(native, "is_protein", "is_ligand")
    sugar = np.isin(native.res_name, ["NAG", "BMA", "MAN", "FUC", "GAL", "NDG", "A2G", "SIA"])
    native = native[((native.chain_id == "A") & native.is_protein) | sugar]
    assert np.isfinite(native.coord).all(), "Missing native coordinates require explicit review"
    assert ((native.atom_name == "CA") & native.is_protein).sum() == 392
    assert sugar.sum() > 0
    audits = {}
    for key, spec in record["experiments"].items():
        binder = as_atom_array_plus(
            parse(str(suite / spec["binder_input"]), config=DISTILLATION)["assemblies"]["1"][0]
        )
        ensure_annotations(binder, "is_protein", "is_ligand")
        binder = binder[binder.is_protein]
        binder.chain_id[:] = "Z"
        for aa in (native, binder):
            # Both source arrays are unannotated; initialize matching conditions before concatenation.
            C_CRD.set_annotation_from_ground_truth(aa, np.zeros(len(aa), dtype=bool))
            C_SEQ.set_annotation_from_ground_truth(aa, np.zeros(len(aa), dtype=bool))
        aa = concatenate_atom_array_plus([native.copy(), binder])
        aa = add_id_and_entity_annotations(aa, overwrite=True)
        aa = add_iid_annotations(aa, overwrite=True)
        # Source label_entity_id survives concatenation independently of chain_entity.
        # Give every chain its own entity so serialization cannot reuse receptor sequence for Z.
        entities = {chain: str(i + 1) for i, chain in enumerate(np.unique(aa.chain_id))}
        aa.set_annotation("label_entity_id", np.array([entities[c] for c in aa.chain_id]))
        fixed = aa.chain_id != "Z"
        C_CRD.set_annotation_from_ground_truth(aa, fixed)
        C_SEQ.set_annotation_from_ground_truth(aa, fixed)
        hot = (
            (aa.chain_id == "A")
            & np.isin(aa.res_id, spec["hotspots"])
            & ~np.isin(aa.atom_name, ["N", "CA", "C", "O"])
        )
        assert hot.any(), "Hotspots missing after parse"
        C_HOT.set_annotation(aa, hot)
        sync_index_groups(aa, preserve_manual_edits=True)
        checks = [
            run_pipeline_validation,
            check_index_group_res_id_gaps,
            check_index_group_molecule_consistency,
        ]
        violations = [msg for fn in checks for msg in fn(aa)]
        if violations:
            raise ValueError(violations)
        output = suite / spec["annotated_input"]
        output.parent.mkdir(parents=True, exist_ok=True)
        to_cif_file(
            aa, str(output), save_standard_annotations=True, chain_disambiguation="chain_iid"
        )
        back = as_atom_array_plus(parse(str(output), config=ANNOTATED_CIF)["assemblies"]["1"][0])
        assert len(back) == len(aa)

        # Serialization may reorder atoms; compare by complete atom identity.
        def values(arr):
            coords = C_CRD.annotation(arr, default="generate")
            seq = C_SEQ.annotation(arr, default="generate")
            hot = C_HOT.mask(arr)
            return {
                (str(c), int(r), str(n), str(a)): (xyz, v, bool(h))
                for c, r, n, a, xyz, v, h in zip(
                    arr.chain_id, arr.res_id, arr.res_name, arr.atom_name, coords, seq, hot
                )
            }

        before, after = values(aa), values(back)
        assert before.keys() == after.keys()
        for ident, (xyz, seq, hot) in before.items():
            x2, s2, h2 = after[ident]
            assert np.allclose(xyz, x2, atol=0.002, equal_nan=True) and seq == s2 and hot == h2

        def chemistry(arr):
            ids = [
                (str(c), int(r), str(n), str(a))
                for c, r, n, a in zip(arr.chain_id, arr.res_id, arr.res_name, arr.atom_name)
            ]
            bonds = sorted(
                (min(ids[i], ids[j]), max(ids[i], ids[j]), int(order))
                for i, j, order in arr.bonds.as_array()
            )
            charges = {ident: int(charge) for ident, charge in zip(ids, arr.charge)}
            return bonds, charges

        assert chemistry(aa) == chemistry(back), "Bond graph/orders or parsed charges changed"
        violations = [msg for fn in checks for msg in fn(back)]
        if violations:
            raise ValueError(violations)
        audits[key] = dict(
            input_sha256=digest(output),
            fixed_atoms=int(fixed.sum()),
            binder_ca=int(((aa.chain_id == "Z") & (aa.atom_name == "CA")).sum()),
            glycan_atoms=int(
                np.isin(aa.res_name, ["NAG", "BMA", "MAN", "FUC", "GAL", "NDG", "A2G", "SIA"]).sum()
            ),
            hotspot_atoms=int(hot.sum()) if hasattr(hot, "sum") else int(C_HOT.mask(aa).sum()),
            roundtrip_conditions_exact=True,
            pipeline_passed=True,
            bond_orders_and_parsed_charges_roundtrip_exact=True,
        )
    (suite / "input_validation.json").write_text(json.dumps(audits, indent=2) + "\n")
    print("RFD4_REPAIR_INPUTS_VALIDATED", json.dumps(audits), flush=True)


def quality(suite):
    from check_course_outputs import collect

    record = json.loads((suite / "manifest.json").read_text())
    reports = {}
    for key, spec in record["experiments"].items():
        reports[key] = collect(
            sorted((suite / key / "inference_samples").rglob("*_generated.cif")),
            suite / "inputs/target.pdb",
            suite / "inputs/37NT.cif",
            tuple(spec["hotspots"]),
        )
    (suite / "quality.json").write_text(
        json.dumps(
            {"records": reports, "qualified": False, "independent_refolding_performed": False},
            indent=2,
        )
        + "\n"
    )
    print(
        json.dumps(
            {
                key: {
                    "outputs": len(rows),
                    "geometry_passes": sum(r["preliminary_geometry_pass"] for r in rows),
                }
                for key, rows in reports.items()
            }
        )
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("stage", choices=["author", "quality"])
    p.add_argument("suite", type=Path)
    a = p.parse_args()
    (author if a.stage == "author" else quality)(a.suite)
