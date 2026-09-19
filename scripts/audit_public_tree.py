#!/usr/bin/env python3
"""Reject machine-specific paths and live infrastructure identifiers."""

from __future__ import annotations

import re
import sys
from pathlib import Path


LITERAL_MARKERS = (
    "/" + "s/project/",
    "/" + "s/genomes/",
    "/" + "scratch/",
    "/" + "cluster/",
    "/" + "lustre/",
    "/" + "gpfs/",
    "/" + "nfs/",
    "/" + "home/",
    "/" + "Users/",
    ":/" + "processed_data/",
    "gs://" + "aou-gym-processed-data",
    "gs://" + "processed-data-" + "wb-" + "strong-peanut-805",
    "wb-" + "strong-peanut-805",
)
LIVE_PROJECT_ID = re.compile(r"project-(?!XXX\b|REDACTED\b)[A-Za-z0-9]{10,}")


def main(paths: list[str]) -> int:
    findings: list[str] = []
    for raw_path in paths:
        path = Path(raw_path)
        if not path.is_file():
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except (OSError, UnicodeDecodeError):
            continue
        for number, line in enumerate(text.splitlines(), 1):
            if any(marker in line for marker in LITERAL_MARKERS) or LIVE_PROJECT_ID.search(line):
                findings.append(f"{path}:{number}: {line.strip()[:180]}")

    if findings:
        print("Private path or live infrastructure identifier found:", file=sys.stderr)
        print("\n".join(findings), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
