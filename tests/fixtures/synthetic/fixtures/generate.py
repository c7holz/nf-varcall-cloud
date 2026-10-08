"""Writes the synthetic test fixtures into ``tests/data/``: the reference window
as ``ref.fasta``. The haplotype carrying the planted variants goes to ``build/``.
Run from ``tests/fixtures/synthetic/``::

    python3 -m fixtures.generate
"""

from __future__ import annotations

import sys
from pathlib import Path

from . import VARIANTS, WINDOW_CHROM, WINDOW_END, WINDOW_START, sources

DATA_DIR = Path(__file__).resolve().parents[3] / "data"
BUILD_DIR = Path(__file__).resolve().parents[1] / "build"

LINE_WIDTH = 60  # bases per FASTA sequence line


def fasta(name: str, dna: str) -> str:
    """One FASTA record: a header line, then the sequence in LINE_WIDTH-base lines."""
    lines = [f">{name}"]
    lines += [dna[i : i + LINE_WIDTH] for i in range(0, len(dna), LINE_WIDTH)]
    return "\n".join(lines) + "\n"


def inject(dna: str, start: int, variants) -> str:
    """``dna``, the reference from 0-based ``start``, with ``variants`` applied.

    Raises ValueError where a variant's REF is not the reference base.
    """
    # Right to left, so a length change doesn't move positions still to apply.
    for pos, ref, alt in sorted(variants, reverse=True):
        i = pos - 1 - start  # 1-based chromosome POS → 0-based index into dna
        found = dna[i : i + len(ref)]
        if found != ref:
            raise ValueError(f"POS {pos}: REF {ref}, reference has {found}")
        dna = dna[:i] + alt + dna[i + len(ref) :]
    return dna


def main() -> int:
    dna = sources.sequence(WINDOW_CHROM, WINDOW_START, WINDOW_END)

    out = DATA_DIR / "ref.fasta"
    out.write_text(fasta(WINDOW_CHROM, dna))
    print(f"wrote {out.relative_to(DATA_DIR.parent.parent)} ({len(dna):,} bp)")

    alt = inject(dna, WINDOW_START, VARIANTS)
    BUILD_DIR.mkdir(exist_ok=True)
    out = BUILD_DIR / "alt.fasta"
    out.write_text(fasta(WINDOW_CHROM, alt))
    print(f"wrote {out.relative_to(DATA_DIR.parent.parent)} ({len(alt):,} bp)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
