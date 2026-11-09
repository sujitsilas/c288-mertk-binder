import copy

import numpy as np
import pytest
import torch

from mertk.comparison import sequence_metrics, summarize
from mertk.dataset import make_batch, synthetic_example
from mertk.model import InterfaceModel, make_esm3
from mertk.proteinmpnn import mpnn_inputs


def test_sequence_metrics_only_score_design_positions_and_reject_fixed_changes():
    result = sequence_metrics("ACD", ["ACD", "AED"], [False, True, False])
    assert result["mean_recovery"] == 0.5
    assert result["pairwise_diversity"] == 1.0
    assert result["unique_fraction"] == 1.0
    with pytest.raises(ValueError, match="fixed"):
        sequence_metrics("ACD", ["ACD", "GCD"], [False, True, False])


def results():
    cases = [
        {
            "id": name,
            "component": family,
            "source_sha256": name,
            "mask_sha256": "mask",
            "design_residues": 5,
            "mean_recovery": recovery,
            "pairwise_diversity": 0.5,
            "unique_fraction": 1.0,
        }
        for name, family, recovery in [("a", "f1", 0.0), ("b", "f1", 0.0), ("c", "f2", 1.0)]
    ]
    base = {
        "arm": "frozen_esm3",
        "seed": 42,
        "split": "test",
        "dataset_identity": "data",
        "smoke_test": False,
        "sampling": {"samples": 2},
        "cases": cases,
    }
    combined = copy.deepcopy(base)
    combined["arm"] = "combined"
    for case in combined["cases"]:
        case["mean_recovery"] = 1.0
    return [base, combined]


def test_summary_uses_families_and_paired_differences():
    summary = summarize(results())
    assert summary["arms"]["frozen_esm3"]["mean_recovery"]["mean"] == 0.5
    assert summary["arms"]["combined"]["paired_recovery_delta_vs_frozen"]["mean"] == 0.5
    assert summary["arms"]["combined"]["mean_recovery"]["families"] == 2


@pytest.mark.parametrize("application", [False, True])
def test_summary_reports_matched_proteinmpnn_contrast(application):
    data = results()
    mpnn = copy.deepcopy(data[0])
    mpnn["arm"] = "proteinmpnn"
    for case in mpnn["cases"]:
        case["mean_recovery"] = 0.25
    data.append(mpnn)
    metric = "recovery"
    if application:
        metric = "seed_retention"
        for result in data:
            result["reference_kind"] = "bindcraft_seed"
            for case in result["cases"]:
                case["mean_seed_retention"] = case.pop("mean_recovery")
    summary = summarize(data)["arms"]
    key = f"paired_{metric}_delta_vs_proteinmpnn"
    assert summary["combined"][key]["mean"] == 0.75
    assert summary["frozen_esm3"][key]["mean"] == 0.25
    assert summary["proteinmpnn"][key]["mean"] == 0
    assert summary["combined"][key]["families"] == 2
    data[-1]["cases"].pop()
    with pytest.raises(ValueError, match="Missing matched"):
        summarize(data)


@pytest.mark.parametrize("mismatch", ["mask_sha256", "source_sha256", "component"])
def test_summary_rejects_unmatched_comparisons(mismatch):
    data = results()
    data[1]["cases"][0][mismatch] = "changed"
    with pytest.raises(ValueError, match="Unmatched"):
        summarize(data)


def test_summary_rejects_missing_seed_examples():
    data = results()
    data[1]["cases"].pop()
    with pytest.raises(ValueError, match="Missing matched"):
        summarize(data)


def test_application_summary_labels_seed_retention_instead_of_native_recovery():
    data = results()
    for result in data:
        result["reference_kind"] = "bindcraft_seed"
        for case in result["cases"]:
            case["mean_seed_retention"] = case.pop("mean_recovery")
    summary = summarize(data)
    assert summary["reference_kind"] == "bindcraft_seed"
    combined = summary["arms"]["combined"]
    assert combined["mean_seed_retention"]["mean"] == 1.0
    assert "mean_recovery" not in combined
    data[0]["reference_kind"] = "native"
    with pytest.raises(ValueError, match="different datasets"):
        summarize(data)


def test_mpnn_never_receives_native_design_labels():
    example = synthetic_example()
    example["design_mask"][0] = False
    original = mpnn_inputs(example, "cpu")
    example["binder_sequence"] = np.array("A" + "W" * 11)
    changed = mpnn_inputs(example, "cpu")
    torch.testing.assert_close(original["S_true"], changed["S_true"])
    assert (original["S_true"][original["design"]] == 20).all()
    assert (original["chain_mask"][:, 12:] == 0).all()
    example["binder_oxygen"][1] = np.nan
    with pytest.raises(ValueError, match="N/CA/C/O"):
        mpnn_inputs(example, "cpu")


def test_backbone_atom_ablation_ignores_target_sidechains():
    model = InterfaceModel(make_esm3(tiny=True), use_lora=False, atom_context="backbone").eval()
    torch.nn.init.normal_(model.adapter.out.weight)
    batch = make_batch(synthetic_example(), model.tokenizer)
    batch["neighbor_features"][..., 0, 2] = 1
    batch["neighbor_features"][..., 2, 2] = 1
    with torch.no_grad():
        before = model(batch)
        sidechain = ~batch["neighbor_features"][..., 2].bool()
        batch["neighbor_xyz"][sidechain] += 100
        batch["neighbor_features"][..., 1][sidechain] = 10
        after = model(batch)
    torch.testing.assert_close(before, after)
