#!/usr/bin/env python3
"""Build an offline HTML structure review; predictions remain explicitly provisional."""

import argparse
from pathlib import Path

from mertk.prediction_review import build_review

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("suite", type=Path)
    parser.add_argument("--config", type=Path, default=Path("configs/prediction_review.json"))
    args = parser.parse_args()
    summary = build_review(args.suite, args.config)
    print(
        f"Review refreshed: {summary['designs']} poses / {summary['families']} families; "
        f"{len(summary['pending_controls'])} controls await prediction"
    )
