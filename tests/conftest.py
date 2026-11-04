import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(".cache/matplotlib").resolve()))

import torch  # noqa: E402

torch.set_num_threads(2)
