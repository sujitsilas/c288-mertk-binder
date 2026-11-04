import pytest

from mertk.splits import assign_splits
from mertk.training import audit_splits


def row(identifier, binder, target, pdb="1ABC"):
    return {"id": identifier, "binder_key": binder, "target_key": target, "pdb_id": pdb}


def test_family_connections_and_all_entry_assemblies_stay_together():
    rows = [row("a", "a", "b"), row("b", "c", "d"), row("c", "e", "f", "2ABC")]
    clusters = dict(a="A", b="B", c="C", d="D", e="B", f="F", TAM_MERTK="TAM")
    result = assign_splits(rows, clusters)
    assert len({r["component"] for r in result}) == 1
    assert len({r["split"] for r in result}) == 1
    audit_splits(result)


def test_tam_family_is_reserved_on_either_side_with_transitive_closure():
    rows = [row("a", "a", "b"), row("b", "c", "d", "2ABC")]
    clusters = dict(a="MERTK", b="B", c="B", d="D", TAM_MERTK="MERTK")
    assert all(r["split"] == "application_holdout" for r in assign_splits(rows, clusters))


def test_homologous_pairs_are_excluded():
    result = assign_splits([row("a", "a", "b")], dict(a="A", b="A", TAM_MERTK="TAM"))
    assert result[0]["split"] == "excluded_homologous_pair"


def test_missing_clusters_or_holdout_seeds_fail_closed():
    with pytest.raises(ValueError, match="TAM"):
        assign_splits([row("a", "a", "b")], dict(a="A", b="B"))
    with pytest.raises(ValueError, match="Missing cluster"):
        assign_splits([row("a", "a", "b")], dict(a="A", TAM_MERTK="TAM"))


def test_audit_detects_swapped_binder_target_leakage():
    rows = [
        dict(row("a", "a", "b"), split="train", component="one"),
        dict(row("b", "c", "a", "2ABC"), split="val", component="two"),
    ]
    with pytest.raises(ValueError, match="leakage"):
        audit_splits(rows)
