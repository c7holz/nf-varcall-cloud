# 0006 — Read simulator: ART over wgsim/dwgsim

- Status: Accepted
- Date: 2026-10-14

## Context
The synthetic fixtures (ADR 0005) need a read simulator. Its error model decides what the
FASTP trimming test actually exercises, and the fixtures must regenerate byte-identically from
a pinned container and seed. Candidates: ART (`art_illumina`), wgsim, dwgsim.

## Decision
Simulate with ART-MountRainier 2016.06.05 using its empirical HiSeq 2500 profile, pinned by
image digest, with a fixed seed per haplotype, at zero mutation rate so the only variants are
the injected ones.

## Consequences
- Base qualities follow ART's profile measured from real HiSeq 2500 runs, including the drop
  toward the 3' end. FASTP's trimming is tested on realistic quality decay, not on one we chose.
- Regenerating the fixtures gives the same bytes: ART has not changed since 2016 and the image
  is pinned by digest.
- Downside: ART is unmaintained, and its newest built-in profiles come from older Illumina
  instruments. Acceptable here: a fixture's profile needs to be realistic, not current.
- Rejected: wgsim — every base gets the same quality, so there is nothing for FASTP to trim.
  dwgsim — quality falls linearly between bounds we set (`-e 0.001-0.01`), so the trimming test
  would only confirm our own settings.
