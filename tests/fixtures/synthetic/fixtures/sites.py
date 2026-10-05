"""Chooses where to plant the fixture variants inside the reference window.

A site qualifies only if a variant called there can mean nothing but "the caller
worked" — unambiguous alignment across it, trustworthy truth data over it, and no
repeat structure that lets a tool move the coordinate. Run it to regenerate
``results/sites.tsv`` and the derivation of the two chosen fixtures::

    python3 -m fixtures.sites
"""

from __future__ import annotations

import sys
from collections.abc import Iterator
from pathlib import Path

from . import WINDOW_CHROM, WINDOW_END, WINDOW_START, sources

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"

# Geometry, all in bp.
FLANK = 150  # read length: furthest a base can be from the site and still land on the same read
EDGE = 1000  # min distance from either window edge for a site to be considered
HC_MARGIN = 200  # distance a site must lie inside a GIAB high-confidence (HC) interval
LOCAL = 25  # radius of the low-complexity filter
NEAR = 5  # radius of the homopolymer filter

MIN_DISTINCT = 4  # distinct bases required within ±LOCAL
MAX_NEAR_RUN = 4  # a homopolymer run this long within ±NEAR disqualifies the site

# A purine↔pyrimidine substitution for each reference base. Transversions only:
# the substitution class least plausibly produced by a sequencing error, so a
# false positive is less likely to be mistaken for a true call.
TRANSVERSIONS = {"A": "CT", "C": "AG", "G": "CT", "T": "AG"}

# The search is seeded at the window's round midpoint (1-based) and prefers round
# coordinates near it. Roundness is a legibility tie-break, not a criterion — a
# memorable coordinate is easier to talk about — and it loses to every filter.
MIDPOINT = 40_150_000
ROUND = 500
SNP_SEARCH = 2000  # search radius around MIDPOINT for the SNP site

# The deletion must be callable independently of the SNP, so it is chosen at a
# distance from whichever site the SNP stage picked.
DEL_MIN_SNP_DIST = 1000
DEL_SEARCH = 3000  # search radius around the SNP for the deletion site


def near(i: int) -> range:
    """The ±NEAR neighbourhood of offset ``i``, inclusive."""
    return range(i - NEAR, i + NEAR + 1)


def to_pos1(offset: int) -> int:
    """Window offset (0-based) → genomic coordinate (1-based, as VCF quotes it)."""
    return WINDOW_START + offset + 1


def context(dna: str, i: int) -> str:
    """The site (bracketed) in its ±20 bp neighborhood."""
    return f"{dna[i - 20 : i]} [{dna[i]}] {dna[i + 1 : i + 21]}"


def homopolymer_run(seq: str, i: int) -> int:
    """Length of the run of identical bases through the site ``seq[i]``."""
    base = seq[i]
    run, j = 1, i - 1
    while j >= 0 and seq[j] == base:
        run += 1
        j -= 1
    j = i + 1
    while j < len(seq) and seq[j] == base:
        run += 1
        j += 1
    return run


def in_highconf(hc: list[tuple[int, int]], i: int) -> bool:
    """True if offset ``i`` lies at least HC_MARGIN inside a high-confidence interval."""
    # hc is BED, so its bounds share WINDOW_START's frame: 0-based, half-open.
    return any(s + HC_MARGIN <= WINDOW_START + i < e - HC_MARGIN for s, e in hc)


def candidates(dna: str, hc: list[tuple[int, int]]) -> list[int]:
    """Window offsets whose surroundings leave a planted variant unambiguous, ascending."""
    out = []
    for i in range(EDGE, len(dna) - EDGE):
        # Case is data: lowercase marks soft-masked repeat, so a read span that is
        # not entirely uppercase puts repeat sequence under every read of the site.
        if not dna[i - FLANK : i + FLANK + 1].isupper():
            continue
        if not in_highconf(hc, i):
            continue
        # Differs from both neighbours. This is also exactly what makes a 1 bp
        # deletion here non-left-shiftable, so the deletion stage need not re-test it.
        if homopolymer_run(dna, i) != 1:
            continue
        # Low complexity, like a homopolymer run, lets a correct call land off the
        # planted coordinate.
        if len(set(dna[i - LOCAL : i + LOCAL + 1])) < MIN_DISTINCT:
            continue
        # A homopolymer run within ±NEAR is still inside the reads spanning the site.
        if max(homopolymer_run(dna, k) for k in near(i)) >= MAX_NEAR_RUN:
            continue
        out.append(i)
    return out


def alt_verdict(dna: str, i: int, alt: str) -> tuple[bool, str]:
    """Whether substituting the alternate base ``alt`` at site ``i`` leaves the expected call unambiguous.

    CpG dinucleotides are refused in both directions — creating one, or mutating
    one away — because they are hypermutable, and a fixture on a mutation hotspot
    invites the truth set itself to be wrong.
    """
    mut = dna[:i] + alt + dna[i + 1 :]
    if max(homopolymer_run(mut, k) for k in near(i)) >= MAX_NEAR_RUN:
        return False, "creates homopolymer"
    if alt == "C" and mut[i + 1] == "G":
        return False, "creates CpG"
    if alt == "G" and mut[i - 1] == "C":
        return False, "creates CpG"
    if dna[i] == "C" and dna[i + 1] == "G":
        return False, "ref is CpG"
    if dna[i] == "G" and dna[i - 1] == "C":
        return False, "ref is CpG"
    return True, "clean"


def deletion_verdict(dna: str, j: int) -> tuple[bool, str]:
    """Whether deleting the base at position ``j`` leaves the flanks it joins free of a homopolymer run."""
    mut = dna[:j] + dna[j + 1 :]
    if max(homopolymer_run(mut, k) for k in near(j)) >= MAX_NEAR_RUN:
        return False, "merges homopolymer"
    return True, "clean"


def by_distance(offsets: list[int], anchor: int) -> list[int]:
    """``offsets`` ranked by distance to the 1-based coordinate ``anchor``; ties keep ascending order."""
    return sorted(offsets, key=lambda i: abs(to_pos1(i) - anchor))


def snp_options(dna: str, cands: list[int]) -> Iterator[tuple[int, list[str]]]:
    """Round coordinates near MIDPOINT, nearest first, each with the SNPs that can be planted there — possibly none."""
    for i in by_distance(cands, MIDPOINT):
        # Ranked by distance: stop at the first candidate outside the search radius.
        if abs(to_pos1(i) - MIDPOINT) > SNP_SEARCH:
            return
        # Only round coordinates are considered.
        if to_pos1(i) % ROUND:
            continue
        yield i, [a for a in TRANSVERSIONS[dna[i]] if alt_verdict(dna, i, a)[0]]


def deletion_options(dna: str, cands: list[int], snp_pos: int) -> Iterator[int]:
    """Round coordinates far enough from ``snp_pos`` to delete cleanly, nearest first."""
    for i in by_distance(cands, snp_pos):
        # Ranked by distance: stop at the first candidate outside the search radius.
        if abs(to_pos1(i) - snp_pos) > DEL_SEARCH:
            return
        # Only round coordinates, and only those far enough from the SNP.
        if to_pos1(i) % ROUND or abs(to_pos1(i) - snp_pos) < DEL_MIN_SNP_DIST:
            continue
        if deletion_verdict(dna, i)[0]:
            yield i


def explain(dna: str, hc: list[tuple[int, int]], pos1: int) -> list[str]:
    """Every filter at one position, measured against the threshold it must pass."""
    i = pos1 - 1 - WINDOW_START
    distinct_bases = sorted(set(dna[i - LOCAL : i + LOCAL + 1]))
    distinct = f"{len(distinct_bases)} ({''.join(distinct_bases)})"
    rows = [
        ("context", context(dna, i), ""),
        ("distance to edge", min(i, len(dna) - 1 - i), f"≥ {EDGE}"),
        (f"uppercase ±{FLANK}", dna[i - FLANK : i + FLANK + 1].isupper(), ""),
        (f"≥{HC_MARGIN} inside high-conf", in_highconf(hc, i), ""),
        ("homopolymer run", homopolymer_run(dna, i), "= 1"),
        (f"max run ±{NEAR}", max(homopolymer_run(dna, k) for k in near(i)), f"< {MAX_NEAR_RUN}"),
        (f"distinct bases ±{LOCAL}", distinct, f"≥ {MIN_DISTINCT}"),
    ]
    # Booleans carry their own verdict; the counts need the constant they were compared against.
    return [f"{name:<22}: {value}{f'   (must be {need})' if need else ''}" for name, value, need in rows]


COLUMNS = ("chrom", "pos", "ref", "max_run_near", "distinct_local", "snp_alts", "del_check")


def write_tsv(path: Path, dna: str, cands: list[int]) -> None:
    """One row per qualifying site, carrying the measurements behind its verdicts."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# sites where a planted variant is unambiguously callable",
        f"# {WINDOW_CHROM}:{WINDOW_START:,}-{WINDOW_END:,} (hg38), "
        f"sha256(sequence)={sources.sha256(dna)}",
        "# regenerate: python3 -m fixtures.sites   —   pos is 1-based",
        "\t".join(COLUMNS),
    ]
    for i in cands:
        alts = ";".join(f"{a}={alt_verdict(dna, i, a)[1]}" for a in TRANSVERSIONS[dna[i]])
        lines.append(
            "\t".join(
                str(v)
                for v in (
                    WINDOW_CHROM,
                    to_pos1(i),
                    dna[i],
                    max(homopolymer_run(dna, k) for k in near(i)),
                    len(set(dna[i - LOCAL : i + LOCAL + 1])),
                    alts,
                    deletion_verdict(dna, i)[1],
                )
            )
        )
    path.write_text("\n".join(lines) + "\n")


def report(dna: str, hc: list[tuple[int, int]], cands: list[int]) -> tuple[int, int]:
    """Print the derivation and return the chosen (snp, deletion) offsets."""
    print(f"{len(cands)} sites qualify in {WINDOW_CHROM}:{WINDOW_START:,}-{WINDOW_END:,}\n")

    print(f"--- {MIDPOINT:,}, the round midpoint, does not qualify ---")
    for line in explain(dna, hc, MIDPOINT):
        print(f"    {line}")

    print(f"\n--- nearest qualifying sites to {MIDPOINT:,} ---")
    for i in by_distance(cands, MIDPOINT)[:5]:
        print(f"chr20:{to_pos1(i):,}  ref={dna[i]}  d={abs(to_pos1(i) - MIDPOINT)}")
        print(f"    {context(dna, i)}")

    print(f"\n--- round coordinates within {SNP_SEARCH} bp, with their transversions ---")
    snp = None
    for i, alts in snp_options(dna, cands):
        print(f"chr20:{to_pos1(i):,}  ref={dna[i]}  d={abs(to_pos1(i) - MIDPOINT)}")
        print(f"    {context(dna, i)}")
        for alt in TRANSVERSIONS[dna[i]]:
            ok, why = alt_verdict(dna, i, alt)
            print(f"      {dna[i]}>{alt}: {'OK' if ok else 'REJECT'} ({why})")
        if snp is None and alts:
            snp = i  # nearest round coordinate with a usable transversion

    print(f"\n--- round coordinates ≥{DEL_MIN_SNP_DIST} bp from the SNP, deleting cleanly ---")
    deletion = None
    for i in deletion_options(dna, cands, to_pos1(snp)):
        print(f"chr20:{to_pos1(i):,}  del {dna[i]}  d(SNP)={abs(to_pos1(i) - to_pos1(snp))}")
        print(f"    {context(dna, i)}")
        if deletion is None:
            deletion = i  # nearest such coordinate
    return snp, deletion


def main() -> int:
    dna = sources.sequence(WINDOW_CHROM, WINDOW_START, WINDOW_END)
    hc = sources.highconf_overlapping(WINDOW_CHROM, WINDOW_START, WINDOW_END)
    cands = candidates(dna, hc)

    out = RESULTS_DIR / "sites.tsv"
    write_tsv(out, dna, cands)
    snp, deletion = report(dna, hc, cands)

    alt = next(a for a in TRANSVERSIONS[dna[snp]] if alt_verdict(dna, snp, a)[0])
    j, pos1 = deletion, to_pos1(deletion)
    print(f"\nSNP      : {WINDOW_CHROM}:{to_pos1(snp):,} {dna[snp]}>{alt}")
    # VCF anchors a deletion on the preceding base, so its POS is one lower.
    print(f"deletion : {WINDOW_CHROM}:{pos1:,} del {dna[j]}  "
          f"(VCF POS={pos1 - 1} REF={dna[j - 1]}{dna[j]} ALT={dna[j - 1]})")
    print(f"wrote {out.relative_to(RESULTS_DIR.parent)} ({len(cands)} rows)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
