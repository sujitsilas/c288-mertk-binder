import json
import runpy

from mertk.provenance import file_hash

score_main = runpy.run_path("scripts/hpc/score_boltz_screen.py")["main"]


def test_completed_cases_retain_source_pose_contacts_and_templates(tmp_path, monkeypatch):
    class TestParser:
        def __init__(self, **kwargs):
            pass

        def get_structure(self, name, path):
            return [None]

    references = tmp_path / "references"
    references.mkdir()
    canonical = references / "canonical.json"
    canonical.write_text(json.dumps({"sequence": {"value": "AAAA"}}))
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    cases = []
    for name, role, template, contacts in [
        ("ig", "Aim 1 Ig1 binder", "37NT", [152, 153, 153]),
        ("stalk", "Aim 2 stalk binder", "base", [490, 492]),
        ("gas6", "Natural-ligand control", "37NT", None),
    ]:
        entry = {"id": name, "role": role, "target_template": template}
        if contacts:
            entry["reference_contact_residues_4p5A"] = contacts
        path = inputs / f"{name}.yaml"
        path.write_text(name)
        cases.append({**entry, "input_sha256": file_hash(path)})
    protocol = {"seed": 101}
    (tmp_path / "manifest.json").write_text(
        json.dumps(
            {
                "canonical_source_sha256": file_hash(canonical),
                "receptor_sequence": "AAAA",
                "protocol": protocol,
                "cases": cases,
            }
        )
    )
    monkeypatch.setitem(score_main.__globals__, "MMCIFParser", TestParser)
    monkeypatch.setitem(
        score_main.__globals__,
        "score",
        lambda *_: {
            "status": "complete",
            "Ig1_contact_residues_4p5A": [152],
            "stalk_contact_residues_4p5A": [490, 493],
        },
    )
    assert score_main(tmp_path) == 0
    result = json.loads((tmp_path / "assessment.json").read_text())
    assert result["protocol"] == protocol
    assert len(result["scorer"]["sha256"]) == 64
    ig, stalk, ligand = result["cases"]
    for row, site, template in [(ig, 152, "37NT"), (stalk, 490, "base")]:
        assert row["target_template"] == template
        assert row["source_pose_contact_count"] == 2
        assert row["retained_source_pose_contacts"] == [site]
    assert ligand["target_template"] == "37NT"
    assert "retained_source_pose_contacts" not in ligand
