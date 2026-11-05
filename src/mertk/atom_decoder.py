"""Small atom-only decoder control trained on the same target-neighborhood features."""

from torch import nn

from mertk.geometry import GeometricAdapter


class AtomDecoder(nn.Module):
    def __init__(self, hidden=128):
        super().__init__()
        from esm.tokenization.sequence_tokenizer import EsmSequenceTokenizer

        self.tokenizer = EsmSequenceTokenizer()
        self.adapter = GeometricAdapter(hidden, hidden=hidden)
        self.head = nn.Sequential(
            nn.LayerNorm(hidden), nn.Linear(hidden, hidden), nn.SiLU(), nn.Linear(hidden, 64)
        )

    def forward(self, batch, sequence_tokens=None, time=None):
        features = self.adapter(
            batch["backbone"],
            batch["neighbor_xyz"],
            batch["neighbor_features"],
            batch["neighbor_mask"],
        )
        return self.head(features)

    def adapter_state_dict(self):
        return {name: value.detach().cpu() for name, value in self.state_dict().items()}

    def load_adapter_state_dict(self, state):
        self.load_state_dict(state, strict=True)
