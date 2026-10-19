"""Provenance scripts and generator for the #14 synthetic test fixtures.

How the reference window and the two fixture variants were chosen, and
``fixtures.generate``, which writes the fixtures into ``tests/data/``. Hits the
network, never runs in CI.
"""

WINDOW_CHROM = "chr20"
WINDOW_START = 40_100_000  # 0-based, UCSC half-open
WINDOW_END = 40_200_000

# The planted variants, VCF-style: 1-based POS, REF, ALT.
VARIANTS = (
    (40_150_500, "T", "A"),
    (40_152_999, "AG", "A"),
)
