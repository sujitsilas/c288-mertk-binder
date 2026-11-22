"""Build an auditable molecular review from real structures and pending controls."""

import copy
import csv
import html as html_module
import io
import json
import re
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from Bio.PDB import PDBIO, MMCIFParser, PDBParser, Select
from Bio.SeqUtils import seq1

from mertk.protector_screen import bond_distances, rank_key, saved_poses
from mertk.provenance import file_hash
from mertk.structural import rigid_fit, rmsd

SUGARS = {"NAG", "BMA", "MAN", "FUC", "GAL", "NDG", "A2G", "SIA"}


class HeavyAtoms(Select):
    def accept_atom(self, atom):
        return atom.element not in {"H", "D"}


def sequence(chain):
    return "".join(seq1(r.resname, undef_code="X") for r in chain if r.id[0] == " ")


def native_numbering(chain, canonical):
    """Recover native numbering from exact sequence, not generator residue labels."""
    seq = sequence(chain)
    start = canonical.find(seq)
    if len(seq) < 30 or start < 0 or canonical.find(seq, start + 1) >= 0:
        raise ValueError("Target must match a unique contiguous native MERTK sequence")
    residues = [r for r in chain if r.id[0] == " "]
    return {r.id: start + index + 1 for index, r in enumerate(residues)}


def pdb_chain(chain, mapping=None, transform=None):
    result = copy.deepcopy(chain)
    # Detach children before renumbering, so old/new IDs cannot collide.
    residues = list(result)
    for r in residues:
        result.detach_child(r.id)
    for r in residues:
        if mapping is not None and r.id in mapping:
            r.id = (" ", mapping[r.id], " ")
        if transform is not None:
            rotation, translation = transform
            for atom in r:
                atom.coord = atom.coord @ rotation + translation
        result.add(r)
    writer = PDBIO()
    writer.set_structure(result)
    stream = io.StringIO()
    writer.save(stream, select=HeavyAtoms())
    return stream.getvalue()


def nearest_points(x, y):
    distances, indices = [], []
    for block in np.array_split(x, max(1, int(np.ceil(len(x) / 128)))):
        d = np.linalg.norm(block[:, None] - y[None], axis=-1)
        distances.extend(d.min(axis=1).tolist())
        indices.extend(d.argmin(axis=1).tolist())
    return np.array(distances), np.array(indices)


def inspect_pose(path, native, canonical, expected_target=None, expected_partner=None):
    parser = MMCIFParser(QUIET=True) if path.suffix == ".cif" else PDBParser(QUIET=True)
    model = parser.get_structure("prediction", path)[0]
    proteins = [c for c in model if any(r.id[0] == " " for r in c)]
    matches = []
    for chain in proteins:
        try:
            mapping = native_numbering(chain, canonical)
        except (ValueError, KeyError):
            continue
        matches.append((chain, mapping))
    if len(proteins) != 2 or len(matches) != 1:
        raise ValueError("Need exactly two protein chains and one exact native MERTK match")
    target, mapping = matches[0]
    partner = next(c for c in proteins if c.id != target.id)
    if expected_target is not None and sequence(target) != expected_target:
        raise ValueError("Prediction does not match the declared target input sequence")
    if expected_partner is not None and sequence(partner) != expected_partner:
        raise ValueError("Prediction does not match the declared partner input sequence")
    native_res = {r.id[1]: r for r in native["A"] if r.id[0] == " "}
    mapped = {mapping[r.id]: r for r in target if r.id in mapping}
    atoms = [a for r in partner if r.id[0] == " " for a in r if a.element != "H"]
    xyz = np.array([a.coord for a in atoms])
    gly = [a for r in native.get_residues() if r.resname in SUGARS for a in r if a.element != "H"]
    glyxyz = np.array([a.coord for a in gly])
    proteinxyz = np.array([a.coord for r in native_res.values() for a in r if a.element != "H"])
    contacts = []
    for n, r in mapped.items():
        if (
            n in native_res
            and nearest_points(np.array([a.coord for a in r if a.element != "H"]), xyz)[0].min()
            <= 4.5
        ):
            contacts.append(n)
    alignments = {}
    for key, limits in [
        ("Ig1", (93, 185)),
        ("Ig2", (186, 285)),
        ("Ig_pair", (93, 285)),
        ("Stalk", (485, 498)),
        ("Proximal_context", (386, 484)),
    ]:
        positions = sorted(n for n in mapped if n in native_res and limits[0] <= n <= limits[1])
        if len(positions) < (10 if key == "Stalk" else 20):
            continue
        x = np.array([mapped[n]["CA"].coord for n in positions])
        y = np.array([native_res[n]["CA"].coord for n in positions])
        rotation, translation = rigid_fit(x, y)
        moved = xyz @ rotation + translation
        gd, gi = nearest_points(moved, glyxyz)
        pd, _ = nearest_points(moved, proteinxyz)
        closest = int(gd.argmin())
        alignments[key] = {
            "target_ca_rmsd_A": rmsd(x @ rotation + translation, y),
            "glycan_min_A": float(gd.min()),
            "protein_min_A": float(pd.min()),
            "glycan_clashing_atoms": int((gd < 2.5).sum()),
            "closest_glycan_chain": gly[gi[closest]].parent.parent.id,
            "closest_binder_residue": atoms[closest].parent.id[1],
            "partner_pdb": pdb_chain(partner, transform=(rotation, translation)),
            "predicted_target_pdb": pdb_chain(target, mapping, (rotation, translation)),
        }
    return {
        "cleavage_atom_min_distances_A": bond_distances(xyz, mapped)
        if 490 in mapped and 491 in mapped
        else None,
        "source_sha256": file_hash(path),
        "target_sequence": sequence(target),
        "partner_sequence": sequence(partner),
        "raw_target_pdb": pdb_chain(target, mapping),
        "raw_partner_pdb": pdb_chain(partner),
        "contacts": sorted(contacts),
        "alignments": alignments,
        "target_chain": target.id,
        "partner_chain": partner.id,
    }


def safe_json(value):
    return json.dumps(value, separators=(",", ":"), allow_nan=False).replace("<", "\\u003c")


def build_review(suite, config_path):
    config = json.loads(config_path.read_text())
    record = json.loads((suite / "manifest.json").read_text())
    native_path = Path(config["reference"])
    native = MMCIFParser(QUIET=True).get_structure("native", native_path)[0]
    canonical_path = Path(config["canonical_sequence"])
    canonical = json.loads(canonical_path.read_text())["sequence"]["value"]
    cases, errors, sources = (
        [],
        [],
        {str(native_path): file_hash(native_path), str(canonical_path): file_hash(canonical_path)},
    )
    for source in [Path(__file__), config_path, Path(config["template"]), Path(config["library"])]:
        sources[str(source)] = file_hash(source)
    campaign_specs = [(suite, key, spec) for key, spec in record["experiments"].items()]
    for extra in config.get("additional_suites", []):
        extra = Path(extra)
        if extra.resolve() == suite.resolve():
            continue
        extra_record = json.loads((extra / "manifest.json").read_text())
        campaign_specs.extend(
            (extra, key, spec) for key, spec in extra_record["experiments"].items()
        )
    sequence_proposals = 0
    for source_suite, key, spec in campaign_specs:
        protector = config.get("include_protector") and spec.get("epitope") == "Stalk_P490S491"
        if not protector and not str(spec.get("epitope", "")).startswith(("Ig1", "Ig2")):
            continue
        campaign = source_suite / "campaigns" / key
        if protector:
            designs = campaign / "bindcraft/designs"
            sequence_proposals += sum(1 for _ in (designs / "MPNN/Sequences").glob("*.fasta"))
            for table in [designs / "trajectory_stats.csv", designs / "final_design_stats.csv"]:
                if table.is_file():
                    sources[str(table)] = file_hash(table)
            for path, meta in saved_poses(designs):
                try:
                    if "error" in meta:
                        raise ValueError(meta["error"])
                    pose = inspect_pose(
                        path,
                        native,
                        canonical,
                        expected_target=canonical[385:526],
                        expected_partner=meta["expected_sequence"],
                    )
                    fit = pose["alignments"]["Stalk"]
                except (ValueError, KeyError, IndexError) as error:
                    errors.append(f"{path.name}: {error}")
                    continue
                is_refold = meta["stage"] == "Sequence refold"
                cases.append(
                    {
                        **meta,
                        **pose,
                        "id": path.stem,
                        "domain": "Stalk",
                        "family": meta["name"].split("_mpnn")[0],
                        "kind": "design" if is_refold else "trajectory",
                        "engine": "BindCraft / AF2" if is_refold else "BindCraft design trajectory",
                        "status": meta["stage"],
                        "qualified": False,
                        "native_clash": fit["protein_min_A"] < 2 or fit["glycan_min_A"] < 2.5,
                        "source": str(path),
                    }
                )
                sources[str(path)] = pose["source_sha256"]
            continue
        stats_path = campaign / "bindcraft/designs/final_design_stats.csv"
        if not stats_path.exists():
            continue
        sources[str(stats_path)] = file_hash(stats_path)
        stats = {r["Design"]: r for r in csv.DictReader(stats_path.open())}
        domain = "Ig2" if str(spec["epitope"]).startswith("Ig2") else "Ig1"
        for path in sorted((campaign / "bindcraft/designs/Accepted").glob("*.pdb")):
            design = re.sub(r"_model\d+$", "", path.stem)
            if design not in stats:
                errors.append(f"No final-statistics row for {path.name}; not displayed")
                continue
            try:
                pose = inspect_pose(path, native, canonical)
            except (ValueError, KeyError) as error:
                errors.append(f"{path.name}: {error}")
                continue
            stat = stats[design]
            alignment = pose["alignments"][domain]
            clash = alignment["glycan_min_A"] < 2.5 or alignment["protein_min_A"] < 2.0
            cases.append(
                {
                    "id": design,
                    "name": design,
                    "domain": domain,
                    "family": design.split("_mpnn")[0],
                    "kind": "design",
                    "engine": "BindCraft / AF2",
                    "status": "Native-context clash" if clash else "Native review required",
                    "native_clash": clash,
                    "qualified": False,
                    "iptm": float(stat["Average_i_pTM"]),
                    "plddt": 100 * float(stat["Average_pLDDT"]),
                    "metric_note": "Mean AF2 values from BindCraft; not directly calibrated to control AF3 scores.",
                    "source": str(path),
                    **pose,
                }
            )
            sources[str(path)] = file_hash(path)
    refold_paths = dict(config.get("parent_refolds", {}))
    for repair_suite in config.get("rfd4_suites", []):
        repair_suite = Path(repair_suite)
        repair_manifest = repair_suite / "manifest.json"
        repair = json.loads(repair_manifest.read_text())
        sources[str(repair_manifest)] = file_hash(repair_manifest)
        for key, spec in repair["experiments"].items():
            for path in sorted((repair_suite / key / "inference_samples").rglob("*_generated.cif")):
                try:
                    pose = inspect_pose(path, native, canonical)
                    alignment = pose["alignments"]["Ig1"]
                except (ValueError, KeyError) as error:
                    errors.append(f"{path.name}: {error}")
                    continue
                clash = alignment["glycan_min_A"] < 2.5 or alignment["protein_min_A"] < 2
                design_id = f"{repair_suite.name}_{key}_{path.stem}"
                if config.get("repair_refold_directory"):
                    metrics_path = (
                        Path(config["repair_refold_directory"])
                        / key
                        / pose["source_sha256"]
                        / "metrics.json"
                    )
                    refold_paths[design_id] = metrics_path
                cases.append(
                    {
                        "id": design_id,
                        "name": f"RFD4 {path.stem.removesuffix('_generated')}",
                        "domain": "Ig1",
                        "family": Path(spec["parent"]).stem.split("_mpnn")[0],
                        "kind": "design",
                        "engine": "RFD4 partial diffusion",
                        "status": "Native-context clash" if clash else "Repair; refolding required",
                        "native_clash": clash,
                        "qualified": False,
                        "iptm": None,
                        "plddt": None,
                        "metric_note": "Generated repair geometry; independent AF2/AF3 confidence unavailable.",
                        "source": str(path),
                        **pose,
                    }
                )
                sources[str(path)] = file_hash(path)
    refolds = []
    parent_designs = {case["id"]: case for case in cases}
    for parent_id, metrics_file in refold_paths.items():
        metrics_path = Path(metrics_file)
        if not metrics_path.is_file():
            continue
        path = metrics_path.with_name("refold.pdb")
        parent = parent_designs.get(parent_id)
        if parent is None or not path.is_file():
            errors.append(f"Missing parent design or refold for {parent_id}")
            continue
        metrics = json.loads(metrics_path.read_text())
        audit = metrics["native_context"]
        if file_hash(path) != audit["sha256"]:
            errors.append(f"Refold coordinates differ from the scored file for {parent_id}")
            continue
        pose = inspect_pose(
            path, native, canonical, sequence(native["A"]), parent["partner_sequence"]
        )
        contact_pae = metrics["mean_contact_pae_A"]
        contact_text = f"{contact_pae:.1f} Å" if contact_pae is not None else "no contacts"
        refolds.append(
            {
                "id": f"{parent_id}_refold",
                "name": f"AF2 refold · {parent['name']}",
                "domain": parent["domain"],
                "family": parent["family"],
                "kind": "design",
                "engine": "AlphaFold2-Multimer",
                "status": "Native geometry passed"
                if audit["preliminary_geometry_pass"]
                else "Native geometry failed",
                "native_clash": audit["intact_protein_min_heavy_atom_distance_A"] < 2.0
                or audit["observed_glycan_clearance_A"] < 2.5,
                "qualified": False,
                "iptm": metrics["ipTM"],
                "plddt": None,
                "metric_note": (
                    f"Single-seed refold with a receptor-only template. Contact PAE {contact_text}. "
                    f"Binder pose RMSD {metrics['target_aligned_binder_ca_rmsd_A']:.2f} Å. "
                    f"Closest protein contact after whole-receptor alignment "
                    f"{audit['intact_protein_min_heavy_atom_distance_A']:.2f} Å."
                ),
                "source": str(path),
                **pose,
            }
        )
        if parent["engine"] == "RFD4 partial diffusion":
            parent["status"] = "Generated pose with AF2 refold"
            parent["metric_note"] = (
                "Select the corresponding AF2 refold to inspect the independently predicted interface."
            )
        sources[str(path)] = pose["source_sha256"]
        sources[str(metrics_path)] = file_hash(metrics_path)
    cases = refolds + cases
    controls_path = Path(config["controls_manifest"])
    controls = json.loads(controls_path.read_text())
    sources[str(controls_path)] = file_hash(controls_path)
    comparisons = []
    for case in controls["cases"]:
        inputs = controls_path.parent / case["input"]
        if file_hash(inputs) != case["input_sha256"]:
            raise ValueError("Control input changed after preparation")
        base = {
            k: case[k]
            for k in ["id", "name", "role", "caveat", "evidence_url", "length", "accession"]
        }
        predictions = case.get("predictions", [])
        if not predictions:
            comparisons.append(
                {
                    **base,
                    "pending": True,
                    "status": "Complex not predicted",
                    "engine": "AF3 inputs prepared",
                }
            )
        for prediction in predictions:
            path = controls_path.parent / prediction["pose"]
            if prediction["input_sha256"] != case["input_sha256"]:
                raise ValueError("Control prediction is not matched to its prepared input")
            pose = inspect_pose(
                path, native, canonical, controls["target"]["sequence"], case["sequence"]
            )
            comparisons.append(
                {
                    **base,
                    **pose,
                    "id": f"{case['id']}_{prediction['seed']}",
                    "pending": False,
                    "status": "Predicted comparator; not binding proof",
                    "engine": "AlphaFold3",
                    "seed": prediction["seed"],
                    "iptm": prediction.get("iptm"),
                    "plddt": prediction.get("partner_plddt"),
                }
            )
            sources[str(path)] = file_hash(path)
    status_path = suite / "status.json"
    status = json.loads(status_path.read_text()) if status_path.exists() else {"jobs": {}}
    for extra in config.get("queue_suites", []):
        extra_status = Path(extra) / "status.json"
        if extra_status.exists():
            additional = json.loads(extra_status.read_text())
            prefix = Path(extra).name + "_"
            status["jobs"].update({prefix + k: v for k, v in additional["jobs"].items()})
            status["checked_at"] = max(
                status.get("checked_at", ""), additional.get("checked_at", "")
            )
    if config.get("include_protector"):
        cases.sort(key=rank_key)
    data = {
        "sequence_proposals": sequence_proposals,
        "built_at": datetime.now(timezone.utc).isoformat(),
        "scheduler": status,
        "native_cif": native_path.read_text(),
        "designs": cases,
        "controls": comparisons,
        "errors": errors,
        "reference": config.get(
            "reference_note",
            "37NT, observed ectodomain A93–A484; stalk, membrane and intracellular domains are absent",
        ),
        "modeled_reference": config.get("modeled_reference", False),
        "control_note": config.get("control_note"),
        "sources": sources,
        "control_protocol": controls["protocol"],
    }
    template = Path(config["template"]).read_text()
    library = Path(config["library"])
    vendor = json.loads((library.parent / "provenance.json").read_text())
    if file_hash(library) != vendor["files"][library.name]["sha256"]:
        raise ValueError("Pinned viewer library hash changed")
    html = (
        template.replace("__DATA__", safe_json(data))
        .replace("__LIBRARY__", library.read_text())
        .replace("__LICENSE__", html_module.escape((library.parent / "LICENSE").read_text()))
    )
    output = Path(config["output"])
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w", dir=output.parent, prefix=".review-", suffix=".tmp", delete=False
    ) as stream:
        stream.write(html)
        temporary = Path(stream.name)
    temporary.replace(output)
    summary = {
        "built_at": data["built_at"],
        "designs": len(cases),
        "families": len({c["family"] for c in cases}),
        "qualified": 0,
        "pending_controls": [c["name"] for c in comparisons if c["pending"]],
        "errors": errors,
        "sources": sources,
        "html_sha256": file_hash(output),
    }
    output.with_suffix(".json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary
