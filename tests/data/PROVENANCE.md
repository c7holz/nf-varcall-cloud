# Provenance

The reference window and the two planted variants were chosen by
`scripts/fixture-selection/`; its README gives the criteria, and `results/` holds
the scored candidates every number below is read from.

| fixture | value | chosen as |
|---|---|---|
| window | `chr20:40,100,000-40,200,000` (hg38) | lowest repeat fraction (38.45 %) of the 22 intergenic windows with no high-confidence hole ≥ 500 bp, out of 187 surviving exclusion — `results/windows.tsv` |
| SNP | `chr20:40,150,500 T>A` | `results/sites.tsv` |
| deletion | `chr20:40,153,000 del G` (VCF `POS=40152999 REF=AG ALT=A`) | `results/sites.tsv` |

## Sources

Both are version-pinned, so a re-run reads the same inputs:

- hg38 sequence and the `gap`, `genomicSuperDups` and `ncbiRefSeqCurated` tracks — UCSC REST API
- high-confidence intervals — GIAB HG001 v4.2.1 benchmark BED, version in the URL
