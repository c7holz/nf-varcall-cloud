"""Enumerates the candidate reference windows and measures the properties that pick one.

Each measurement scores one window on a property that decides whether a variant
planted inside it can be called unambiguously: how much of the window is repeat,
how much of it the truth set covers, and how much of it a read can be placed in
uniquely.
"""

from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path

from . import WINDOW_CHROM, sources

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"

CHR20_LENGTH = 64_444_167  # hg38, per UCSC hg38.chrom.sizes
WINDOW = 100_000  # candidate window length
STRIDE = 250_000  # distance between consecutive candidate starts
ANCHOR = 100_000  # first candidate's start
# Distance a candidate must keep from every feature of each track.
EXCLUSION_MARGINS = {"gap": 50_000, "genomicSuperDups": 10_000}

HOLE_MIN = 500  # an uncovered stretch at least this long is counted as a hole
GENE_TRACK = "ncbiRefSeqCurated"  # curated models only; excludes the predicted XM_/XR_ ones

COMPLEMENT = str.maketrans("ACGTN", "TGCAN")  # base-pairing table; anything else passes through


def candidate_starts() -> list[int]:
    """0-based starts of the candidate windows, the last one ending within chr20."""
    return list(range(ANCHOR, CHR20_LENGTH - WINDOW + 1, STRIDE))


def near_feature(start: int, intervals: list[tuple[int, int]], margin: int) -> bool:
    """Whether the window lies closer than the margin to any of the intervals."""
    end = start + WINDOW
    return any(start < e + margin and end > s - margin for s, e in intervals)


def survivors(
    starts: list[int], tracks: dict[str, list[tuple[int, int]]]
) -> list[int]:
    """Starts whose windows keep their margin from every feature of every track."""
    return [
        start
        for start in starts
        if not any(
            near_feature(start, tracks[name], margin)
            for name, margin in EXCLUSION_MARGINS.items()
        )
    ]


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


def overlapping_genes(start: int, spans: list[tuple[int, int, str]]) -> list[str]:
    """Symbols of the genes reaching into the window, each named once, ascending."""
    end = start + WINDOW
    return sorted({name for s, e, name in spans if s < end and e > start})


KMERS = (100, 150)  # read lengths a fixture might be simulated at

COLUMNS = (
    "chrom",
    "start",
    "repeat_frac",
    "n_count",
    "gc_frac",
    *(f"unique_{k}" for k in KMERS),
    "hc_bp",
    "hc_intervals",
    "holes_ge500",
    "gaps_lt500",
    "genes",
)


def measure(
    start: int, dna: str, hc: list[tuple[int, int]], spans: list[tuple[int, int, str]]
) -> dict:
    """Every scored property of one window, keyed by its column name."""
    end = start + WINDOW
    clipped = clipped_highconf(hc, start, end)
    uncovered = gaps(clipped, start, end)
    holes = sum(g >= HOLE_MIN for g in uncovered)
    return {
        "chrom": WINDOW_CHROM,
        "start": start,
        "repeat_frac": repeat_fraction(dna),
        "n_count": n_count(dna),
        "gc_frac": gc_fraction(dna),
        **{f"unique_{k}": unique_fraction(dna, k) for k in KMERS},
        "hc_bp": covered_bp(clipped),
        "hc_intervals": len(clipped),
        "holes_ge500": holes,
        "gaps_lt500": len(uncovered) - holes,
        "genes": ";".join(overlapping_genes(start, spans)),
    }


def scan(starts: list[int]) -> list[dict]:
    """Measure each window, fetching its sequence. Prints one line per window."""
    hc = sources.giab_highconf(WINDOW_CHROM)
    spans = sources.gene_spans(GENE_TRACK, WINDOW_CHROM)
    rows = []
    for n, start in enumerate(starts, 1):
        row = measure(
            start, sources.sequence(WINDOW_CHROM, start, start + WINDOW), hc, spans
        )
        rows.append(row)
        print(
            f"[{n}/{len(starts)}] {start:>10,}  "
            f"repeat {row['repeat_frac']:6.2%}  "
            f"high-conf {row['hc_bp'] / WINDOW:6.2%}  "
            f"holes {row['holes_ge500']}",
            flush=True,
        )
    return rows


def rank(rows: list[dict]) -> list[dict]:
    """Windows with no coverage hole and no gene, least repetitive first."""
    return sorted(
        (r for r in rows if r["holes_ge500"] == 0 and not r["genes"]),
        key=lambda r: r["repeat_frac"],
    )


def fmt(value: object) -> str:
    """Fractions at fixed precision, everything else as it prints."""
    return f"{value:.6f}" if isinstance(value, float) else str(value)


def write_tsv(path: Path, rows: list[dict], scope: str) -> None:
    """One row per scored window, carrying the measurements behind the choice."""
    path.parent.mkdir(parents=True, exist_ok=True)
    margins = ", ".join(f"{n} (+{m:,})" for n, m in EXCLUSION_MARGINS.items())
    lines = [
        "# candidate reference windows, scored",
        f"# {WINDOW_CHROM} (hg38), {CHR20_LENGTH:,} bp — {WINDOW:,} bp windows every "
        f"{STRIDE:,} bp from {ANCHOR:,}, keeping clear of {margins}",
        f"# {scope}",
        f"# genes: {GENE_TRACK} symbols overlapping the window; empty is intergenic",
        "# regenerate: python3 -m fixtures.window   —   start is 0-based, half-open",
        "\t".join(COLUMNS),
    ]
    lines += ["\t".join(fmt(row[c]) for c in COLUMNS) for row in rows]
    path.write_text("\n".join(lines) + "\n")


def report(rows: list[dict]) -> None:
    """Print the ranking that chooses the window."""
    ranked = rank(rows)
    print(
        f"\n--- {len(ranked)} of {len(rows)} windows have no hole ≥{HOLE_MIN} bp and no "
        f"{GENE_TRACK} gene, least repetitive first ---"
    )
    for row in ranked[:5]:
        print(
            f"{WINDOW_CHROM}:{row['start']:,}  repeat {row['repeat_frac']:.2%}  "
            f"high-conf {row['hc_bp'] / WINDOW:.2%}  "
            f"unique-150 {row['unique_150']:.2%}  GC {row['gc_frac']:.1%}"
        )


def main() -> int:
    tracks = {
        name: sources.track_intervals(name, WINDOW_CHROM) for name in EXCLUSION_MARGINS
    }
    starts = candidate_starts()
    kept = survivors(starts, tracks)
    rows = scan(kept)

    out = RESULTS_DIR / "windows.tsv"
    write_tsv(out, rows, f"{len(starts)} candidates, {len(kept)} surviving exclusion")
    report(rows)
    print(f"wrote {out.relative_to(RESULTS_DIR.parent)} ({len(rows)} rows)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
