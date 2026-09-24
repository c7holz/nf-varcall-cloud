# fixture-selection

Provenance scripts for the synthetic test fixtures (issue #14). They record **how**
the reference window and the two planted variants were chosen.

The audit is a single command, `reproduce` (phase 4 — not built yet): it re-derives
the window and the two variants from the primary sources and diffs the result against
the committed `results/` TSVs. Reproducing those numbers exercises the whole chain —
fetch, filters, ranking — so it is a stronger check than re-verifying the raw inputs
would be, and it subsumes any separate source check.

There is deliberately **no drift monitor**. The primary sources are version-pinned
and immutable — a frozen genome build (hg38) and a released benchmark (GIAB HG001
v4.2.1, pinned in the URL) — so "do the sources still say the same thing?" is
answered by construction. If a source ever became unreachable or changed shape, the
fetch fails loudly rather than returning a quietly different answer; it cannot drift
silently. The one lasting provenance fact — each source's URL, fetch date, and
`sha256` — will be recorded in `results/sources.tsv`, documented rather than
re-checked.

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
| **k-mer uniqueness** | Whether short subsequences of the window recur elsewhere in the genome. Another angle on "can reads from here be placed unambiguously". |

## Selection criteria

Two tiers. **Window** criteria make alignment unambiguous; **site** criteria make the
expected call unambiguous.

### Window (chr20-wide scan, ~100 kb candidates)

- no overlap with the `gap` track — no `N` runs
- no overlap with `genomicSuperDups` — no ambiguous placement
- high k-mer uniqueness
- well covered by GIAB high-confidence intervals

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
  __init__.py   the window under audit, as three constants        — present
  sources.py    the ONLY network boundary — UCSC + GIAB fetchers   — present
  sites.py      site criteria — pure, no network                  — phase 2
  window.py     chr20-wide window scan                            — phase 3
cache/          fetched responses, gitignored                     — appears on first fetch
results/        committed TSVs — the audit artifact               — phase 4
```

`sources.py` is the only code that touches the network; everything the filters do
above it will be pure and unit-testable offline. Each response is cached verbatim
under `cache/`, so that once warmed the derivation re-runs offline and free while you
iterate — above all phase 3's ~190-fetch chr20-wide window scan. The cache is not
committed: it is large and re-fetchable. What *is* committed is `results/` — the TSVs
that every published number is read back from, and what `reproduce` diffs against.
Neither `cache/` nor `results/` exists yet: the first is created on the first fetch,
the second in phase 4.

## Usage

Run from this directory (`scripts/fixture-selection/`), with no dependencies beyond
the standard library:

```bash
python3 -m fixtures.sites
python3 -m fixtures.window
```

### `reproduce` (planned — phase 4)

When built, it will re-derive the whole decision from primary sources — fetch the
window and the `gap` / `genomicSuperDups` tracks (through the cache), run the window
and site filters, and write `results/*.tsv` — then diff those against the committed
copies and report PASS/FAIL.

That diff is the audit. If the committed numbers still fall out of the pipeline —
the two-tier candidate counts, the 257→187 window scan, the two top-ranked variants —
then every claim in the #14 provenance holds. If they don't, it stops and shows what
moved, and nothing downstream should be trusted until that is explained.

There is no separate source-verification step: re-deriving the numbers already
exercises the fetch, and the sources are version-pinned (see the top of this README),
so a source that rots or changes shape makes the fetch *raise*, not lie. An earlier
scratch `verify-sources` command has been removed — reproducing the numbers subsumes
it.

## Status

Rebuilt in phases from the scratch scripts in `analysis-salvage/` (see the NOTE
there for what those did and why they needed replacing).

| phase | scope | state |
|---|---|---|
| 1 | primary-source layer (`sources.py`) | written* |
| 2 | site criteria as a real module; reproduce the two chosen variants | not started |
| 3 | recreate the chr20-wide window scan, lost from an earlier session | not started |
| 4 | `reproduce`: emit committed `results/` TSVs and diff against them | not started |
| 5 | delete `analysis-salvage/`, comment outcome on #14 | not started |

\* `sources.py` is written and offline-verified (imports, pure helpers), but its
fetchers have **not** yet been exercised against live UCSC/NCBI responses — the
`track()` shape-handling in particular is unvalidated. The scratch `verify-sources`
command has already been removed; its job is subsumed by `reproduce`.

## Data sources

All public domain, satisfying the repository's public-data-only rule.

- UCSC hg38 REST API — `getData/sequence`, `getData/track` (`gap`, `genomicSuperDups`)
- GIAB HG001 v4.2.1 benchmark BED (NCBI FTP)
