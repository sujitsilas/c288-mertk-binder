"""Small, dependency-free provenance helpers."""

import hashlib
import subprocess
from pathlib import Path


def file_hash(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def source_revision(path):
    """Describe the Git checkout when the actual source is tracked there."""
    path = Path(path).resolve()
    git = ["git", "-C", str(path.parent)]
    try:
        subprocess.run(
            [*git, "ls-files", "--error-unmatch", "--", path.name],
            check=True,
            capture_output=True,
        )
        commit = subprocess.check_output([*git, "rev-parse", "HEAD"], text=True).strip()
        dirty = bool(subprocess.check_output([*git, "status", "--porcelain"]))
    except (OSError, subprocess.CalledProcessError):
        # A frozen workflow can be copied to a compute node without Git metadata.
        return {"git_commit": None, "git_worktree_dirty": None}
    return {"git_commit": commit, "git_worktree_dirty": dirty}
