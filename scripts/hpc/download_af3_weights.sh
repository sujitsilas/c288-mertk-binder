#!/usr/bin/env bash
# Standard AF3 release, directly from the Google-hosted URL in the pinned README.
set -euo pipefail
project=${1:?remote project directory}
umask 077
directory="$project/weights/alphafold3"
mkdir -p "$directory" "$project/provenance"
url='https://storage.googleapis.com/alphafold3/af3.bin.zst?generation=1780568696389861'
[[ ! -e "$directory/af3.bin.zst" ]] || { echo 'Weights already exist; verify before re-downloading.' >&2; exit 1; }
curl -fL --retry 3 --dump-header "$project/provenance/af3-download-headers.txt" \
  "$url" -o "$directory/af3.bin.zst.partial"
"$project/envs/alphafold3/bin/python" - "$directory" "$url" <<'PY'
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
import zstandard

directory, url = Path(sys.argv[1]), sys.argv[2]
path = directory / "af3.bin.zst.partial"
if path.stat().st_size != 1020545840:
    raise ValueError("Download size differs from the pinned Google Storage object")
with path.open("rb") as stream:
    digest = hashlib.file_digest(stream, "sha256").hexdigest()
uncompressed = 0
with path.open("rb") as stream, zstandard.ZstdDecompressor().stream_reader(stream) as reader:
    while block := reader.read(8 * 1024 * 1024):
        uncompressed += len(block)
path.rename(directory / "af3.bin.zst")
(directory / "provenance.json").write_text(json.dumps({
    "source_url": url, "source": "Official Google DeepMind AF3 README",
    "source_revision": "28c2d290559fe07387a12f471ad49f0d36206327",
    "sha256": digest, "compressed_bytes": 1020545840,
    "uncompressed_bytes": uncompressed, "zstd_stream_validated": True,
    "downloaded_at": datetime.now(timezone.utc).isoformat(),
}, indent=2) + "\n")
print("Standard AF3 weights downloaded and validated:", directory / "af3.bin.zst")
print("SHA256:", digest)
PY
