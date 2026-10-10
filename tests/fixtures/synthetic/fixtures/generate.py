"""Writes the synthetic test fixtures into ``tests/data/``: the reference window
as ``ref.fasta``. Intermediates (the alt haplotype, simulated reads) go to ``build/``.
Run from ``tests/fixtures/synthetic/``::

    python3 -m fixtures.generate
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

from . import VARIANTS, WINDOW_CHROM, WINDOW_END, WINDOW_START, sources

DATA_DIR = Path(__file__).resolve().parents[3] / "data"
BUILD_DIR = Path(__file__).resolve().parents[1] / "build"

LINE_WIDTH = 60  # bases per FASTA sequence line

# ART-MountRainier 2016.06.05 (quay.io/biocontainers/art:2016.06.05--h0704011_13)
ART_IMAGE = "quay.io/biocontainers/art@sha256:18a9360667bb1f50e0368709efee3390e6ed733cf3c6eadd0f6ea441f143ff60"
ART_ARGS = [
    "-ss", "HS25",  # HiSeq 2500 quality profile
    "-p", "-l", "150", "-m", "300", "-s", "30",  # paired 150 bp, fragment 300 ± 30
    "-f", "25",  # coverage per haplotype; the two merge to 50×
    "-na",  # no alignment output
]
SEEDS = {"ref": 14, "alt": 15}  # one per haplotype, so their reads don't share positions


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


def simulate(haplotype: Path, prefix: Path, seed: int) -> None:
    """Paired reads from ``haplotype``, written to ``<prefix>1.fq`` and ``<prefix>2.fq``."""
    # Mount tests/, which holds both data/ and fixtures/synthetic/build/.
    root = DATA_DIR.parent
    subprocess.run(
        ["podman", "run", "--rm", "-v", f"{root}:/w:z", "-w", "/w", ART_IMAGE,
         "art_illumina", *ART_ARGS,
         "-rs", str(seed),
         "-i", str(haplotype.relative_to(root)), "-o", str(prefix.relative_to(root))],
        check=True,
    )


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

    for name, haplotype in [("ref", DATA_DIR / "ref.fasta"), ("alt", BUILD_DIR / "alt.fasta")]:
        prefix = BUILD_DIR / f"{name}_"
        simulate(haplotype, prefix, SEEDS[name])
        print(f"wrote {prefix.relative_to(DATA_DIR.parent.parent)}{{1,2}}.fq")
    return 0


if __name__ == "__main__":
    sys.exit(main())
