"""Provenance scripts for the #14 synthetic test fixtures.

One-off analysis: how the reference window and the two fixture variants were
chosen. Hits the network, never runs in CI.
"""

WINDOW_CHROM = "chr20"
WINDOW_START = 40_100_000  # 0-based, UCSC half-open
WINDOW_END = 40_200_000
