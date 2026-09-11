"""Measures the properties that pick the reference window.

Each function scores one window on a property that decides whether a variant
planted inside it can be called unambiguously: how much of the window is repeat,
how much of it the truth set covers, and how much of it a read can be placed in
uniquely.
"""

from __future__ import annotations

from collections import Counter

HOLE_MIN = 500  # an uncovered stretch at least this long is counted as a hole

COMPLEMENT = str.maketrans("ACGTN", "TGCAN")  # base-pairing table; anything else passes through


def repeat_fraction(dna: str) -> float:
    """Share of the window RepeatMasker soft-masked."""
    return sum(c.islower() for c in dna) / len(dna)


def n_count(dna: str) -> int:
    """Bases the reference leaves unresolved."""
    return dna.upper().count("N")


def gc_fraction(dna: str) -> float:
    """Share of the resolved bases that are G or C."""
    seq = dna.upper()
    return (seq.count("G") + seq.count("C")) / (len(seq) - seq.count("N"))


def unique_fraction(dna: str, k: int) -> float:
    """Share of k-mer positions whose sequence occurs nowhere else in the window, either strand."""
    # Lowercase marks repeats, not different bases: unfolded, a masked and an
    # unmasked copy of one sequence would both look unique.
    seq = dna.upper()
    rc = seq.translate(COMPLEMENT)[::-1]
    n = len(seq)
    total = n - k + 1
    counts: Counter[str] = Counter()
    for i in range(total):
        # A k-mer and its reverse complement are one sequence read from two
        # strands; the lexicographically smaller of the pair is its one key.
        # rc is reversed, so the mate of seq[i:i + k] sits at this offset from its end.
        counts[min(seq[i : i + k], rc[n - i - k : n - i])] += 1
    return sum(c == 1 for c in counts.values()) / total


def clipped_highconf(
    hc: list[tuple[int, int]], start: int, end: int
) -> list[tuple[int, int]]:
    """High-confidence intervals restricted to the window, overlaps merged, ascending."""
    out: list[tuple[int, int]] = []
    for s, e in sorted(hc):
        s, e = max(s, start), min(e, end)
        if s >= e:
            continue
        # Merge into the last interval when this one overlaps it; ends that only
        # touch count as overlapping.
        if out and s <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return out


def covered_bp(clipped: list[tuple[int, int]]) -> int:
    """Bases of the window the high-confidence intervals cover."""
    return sum(e - s for s, e in clipped)


def gaps(clipped: list[tuple[int, int]], start: int, end: int) -> list[int]:
    """Lengths of the window's uncovered stretches, those at its edges included.

    Expects intervals inside the window, ascending and non-overlapping.
    """
    out, prev = [], start
    for s, e in clipped:
        if s > prev:
            out.append(s - prev)
        prev = e
    if end > prev:
        out.append(end - prev)
    return out
