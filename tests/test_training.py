import json
import random
from pathlib import Path

import pytest
import torch

from mertk.dataset import TrainingSampler, make_batch, select_training_families, synthetic_example
from mertk.model import InterfaceModel, make_esm3
from mertk.training import evaluate, train


@pytest.mark.parametrize("device,precision", [("cpu", "fp32"), ("mps", "fp32")])
def test_resume_reproduces_next_optimizer_step(tmp_path, monkeypatch, device, precision):
    if device == "mps" and not torch.backends.mps.is_available():
        pytest.skip("MPS is unavailable in this process")
    config = json.loads(Path("configs/train.json").read_text())
    config.update(
        steps=2,
        eval_every=1,
        accumulation_steps=1,
        objective="denoise",
        output=str(tmp_path / "run"),
        device=device,
        precision=precision,
        smoke_device=True,
    )
    saved_step_one = tmp_path / "step_one.pt"
    save = torch.save

    def capture_first_checkpoint(state, path):
        save(state, path)
        if state["step"] == 1:
            save(state, saved_step_one)

    monkeypatch.setattr(torch, "save", capture_first_checkpoint)
    result = train(config, smoke=True)
    uninterrupted = torch.load(result, weights_only=True)
    train(config, smoke=True, resume=saved_step_one)
    resumed = torch.load(result, weights_only=True)
    assert resumed["step"] == 2
    metrics = [
        json.loads(line) for line in (tmp_path / "run/metrics.jsonl").read_text().splitlines()
    ]
    assert [item["step"] for item in metrics] == [1, 2]
    for key in uninterrupted["adapter"]:
        torch.testing.assert_close(
            uninterrupted["adapter"][key], resumed["adapter"][key], atol=0, rtol=0
        )


def test_evaluation_rejects_nonfinite_metrics(monkeypatch):
    model = InterfaceModel(make_esm3(tiny=True))
    monkeypatch.setattr(
        model,
        "forward",
        lambda batch, **kwargs: torch.full((*batch["sequence_tokens"].shape, 64), float("nan")),
    )
    with pytest.raises(FloatingPointError, match="non-finite validation"):
        evaluate(model, [{"id": "synthetic"}], torch.device("cpu"), tiny=True)


def test_training_subsets_are_nested_and_sampling_balances_families():
    rows = [{"id": f"a{i}", "component": "large"} for i in range(100)] + [
        {"id": "b", "component": "small"}
    ]
    sampler, rng = TrainingSampler(rows), random.Random(42)
    selected = [sampler.sample(rng)["component"] for _ in range(2000)]
    assert 900 < selected.count("small") < 1100
    families = [{"id": str(i), "component": f"family{i}"} for i in range(10)]
    small = select_training_families(families, 0.2, 77)
    large = select_training_families(list(reversed(families)), 0.5, 77)
    assert len(small) == 2 and len(large) == 5
    assert {row["id"] for row in small} <= {row["id"] for row in large}


def test_validation_reports_family_macro_recovery(monkeypatch):
    model = InterfaceModel(make_esm3(tiny=True))

    def batch(row, model, device, tiny):
        result = make_batch(synthetic_example(), model.tokenizer)
        result["correct"] = row["correct"]
        return result

    def predict(batch, **kwargs):
        logits = torch.zeros((*batch["sequence_tokens"].shape, 64))
        labels = (
            batch["sequence_tokens"]
            if batch["correct"]
            else torch.full_like(
                batch["sequence_tokens"], model.tokenizer.convert_tokens_to_ids("W")
            )
        )
        return logits.scatter(-1, labels.unsqueeze(-1), 10)

    monkeypatch.setattr("mertk.training.get_batch", batch)
    monkeypatch.setattr(model, "forward", predict)
    rows = [{"id": str(i), "component": "large", "correct": False} for i in range(2)] + [
        {"id": "good", "component": "small", "correct": True}
    ]
    metrics = evaluate(model, rows, torch.device("cpu"), tiny=True)
    assert metrics["sequence_recovery"] == pytest.approx(1 / 3)
    assert metrics["family_sequence_recovery"] == 0.5


def test_best_checkpoint_comes_from_validation_not_final_step(tmp_path, monkeypatch):
    config = json.loads(Path("configs/train.json").read_text())
    config.update(steps=2, eval_every=1, accumulation_steps=1, output=str(tmp_path / "run"))
    scores = iter([0.75, 0.25])
    monkeypatch.setattr(
        "mertk.training.evaluate",
        lambda *args, **kwargs: {"family_sequence_recovery": next(scores)},
    )
    last = torch.load(train(config, smoke=True), weights_only=True)
    best = torch.load(tmp_path / "run/best.pt", weights_only=True)
    assert last["step"] == 2 and best["step"] == 1
    assert best["selection"] == {
        "metric": "validation.family_sequence_recovery",
        "score": 0.75,
        "split": "val",
    }
    for name in best["adapter"]:
        torch.testing.assert_close(
            best["adapter"][name], last["best"]["adapter"][name], atol=0, rtol=0
        )
