"""Writes the synthetic test fixtures into ``tests/data/``: the reference window
as ``ref.fasta``. Run from ``tests/fixtures/synthetic/``::

    python3 -m fixtures.generate
"""

from __future__ import annotations

import sys
from pathlib import Path

from . import WINDOW_CHROM, WINDOW_END, WINDOW_START, sources

DATA_DIR = Path(__file__).resolve().parents[3] / "data"

LINE_WIDTH = 60  # bases per FASTA sequence line


def fasta(name: str, dna: str) -> str:
    """One FASTA record: a header line, then the sequence in LINE_WIDTH-base lines."""
    lines = [f">{name}"]
    lines += [dna[i : i + LINE_WIDTH] for i in range(0, len(dna), LINE_WIDTH)]
    return "\n".join(lines) + "\n"


def main() -> int:
    dna = sources.sequence(WINDOW_CHROM, WINDOW_START, WINDOW_END)

    out = DATA_DIR / "ref.fasta"
    out.write_text(fasta(WINDOW_CHROM, dna))
    print(f"wrote {out.relative_to(DATA_DIR.parent.parent)} ({len(dna):,} bp)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
