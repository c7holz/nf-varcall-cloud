# fixture-selection

Provenance scripts for the synthetic test fixtures (issue #14). They record **how**
the reference window and the two planted variants were chosen.

These scripts hit the network and **never run in CI**. Not to be confused with
`scripts/make_test_data.sh`, which is hermetic and regenerates `tests/data/`.

## The problem being solved

The pipeline's tests need a tiny reference and a tiny read set whose *correct*
variant calls are known in advance. Simulating reads is easy; the hard part is
choosing **where** to simulate them. Pick a bad spot and the test becomes a liar:
it fails on a correct caller, or passes on a broken one.

So the job splits in two:

1. **Pick a window** — ~100 kb of hg38 that is unambiguous to align against.
2. **Pick variant sites inside it** — positions where a called variant can only
   mean the caller worked, and a missed one can only mean it didn't.

The current answer is `chr20:40,100,000-40,200,000`, with a `T>A` SNP at
`chr20:40,150,500` and a 1 bp `G` deletion at `chr20:40,153,000`.

## Glossary

Enough genomics to read the code; nothing more.

| term | meaning |
|---|---|
| **hg38 / GRCh38** | The human reference genome build everything here is keyed to. |
| **0-based half-open** | Coordinate convention: `[start, end)`, counting from 0. UCSC's API and BED files use it; the variants above are quoted 1-based, the way VCF and papers do. The two differ by one, and mixing them is the classic bug in this domain. |
| **soft-masking** | The reference marks repetitive sequence by *lowercasing* it. Uppercase = unique-ish, lowercase = repeat. Case is data, not formatting. |
| **GIAB high-confidence intervals** | Genome in a Bottle publishes regions of sample HG001 where the truth set is trustworthy. Outside them, "truth" is itself uncertain — so a fixture placed outside is untestable. |
| **`gap` track** | Assembly gaps: stretches of reference that are literally unknown (runs of `N`). |
| **`genomicSuperDups`** | Segmental duplications — regions that appear near-identically elsewhere in the genome. A read from here can align to two places, so alignment is ambiguous by construction. |
| **homopolymer** | A run of one base (`AAAAA`). Aligners and callers slip inside them, making the exact variant coordinate unstable. |
| **left-shifting** | A deletion inside a repeat can be written at several coordinates that all describe the same edit. Tools normalise differently, so the "expected" coordinate becomes ambiguous. Avoided by deleting a base that differs from both neighbours. |
| **transversion** | A purine↔pyrimidine substitution (`T>A`) as opposed to the more common transition (`T>C`). Chosen because it is the substitution class least likely to be produced by a sequencing error, so a false positive is less plausible. |

## Selection criteria

Two tiers. **Window** criteria make alignment unambiguous; **site** criteria make the
expected call unambiguous.

### Window (chr20-wide scan, ~100 kb candidates)

- no overlap with the `gap` track — no `N` runs
- no overlap with `genomicSuperDups` — no ambiguous placement
- no gene overlap — intergenic
- no stretch ≥ 500 bp outside GIAB high-confidence intervals

Of the windows that pass, the one with the lowest repeat fraction (share of
soft-masked bases) is chosen.

### Site (within the chosen window)

- ≥ 1000 bp from either window edge, so a read spanning the site fits inside it
- the ±150 bp read span is entirely uppercase — no soft-masked repeat under any read
- ≥ 200 bp inside a GIAB high-confidence interval, not merely within one
- not in a homopolymer run longer than 2, and no run of ≥ 4 within ±5 bp
- ≥ 4 distinct bases in the surrounding ±25 bp — rejects low-complexity sequence
- SNP only: transversion, and not creating or destroying a CpG
- deletion only: not left-shiftable, and far enough from the SNP that the two are
  independently callable

## Layout

```
fixtures/
  __init__.py   the chosen window, as three constants
  sources.py    the ONLY network boundary — UCSC + GIAB fetchers
  sites.py      site criteria and ranking
  window.py     chr20-wide window scan
cache/          fetched responses, gitignored
results/        committed TSVs — the record
```

`sources.py` is the only code that touches the network. Each response is cached
verbatim under `cache/`, so once warmed a re-run is offline. The cache is not
committed: it is large and re-fetchable. `results/` is: the TSVs every published
number is read from.

## Usage

Run from this directory (`scripts/fixture-selection/`), with no dependencies beyond
the standard library:

```bash
python3 -m fixtures.sites
python3 -m fixtures.window
```

## Data sources

All public domain, satisfying the repository's public-data-only rule.

- UCSC hg38 REST API — `getData/sequence`, `getData/track` (`gap`, `genomicSuperDups`)
- GIAB HG001 v4.2.1 benchmark BED (NCBI FTP)
