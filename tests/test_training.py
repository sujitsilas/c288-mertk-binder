import json
from pathlib import Path

import torch

from mertk.training import train


def test_resume_reproduces_next_optimizer_step(tmp_path, monkeypatch):
    config = json.loads(Path("configs/train.json").read_text())
    config.update(
        steps=2,
        eval_every=1,
        accumulation_steps=1,
        objective="denoise",
        output=str(tmp_path / "run"),
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
