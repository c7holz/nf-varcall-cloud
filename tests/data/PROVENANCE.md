# Provenance

The reference window and the two planted variants were chosen by
`scripts/fixture-selection/`; its README gives the criteria, and `results/` holds
the scored candidates every number below is read from.

| fixture | value | chosen as |
|---|---|---|
| window | `chr20:40,100,000-40,200,000` (hg38) | lowest repeat fraction (38.45 %) of the 22 intergenic windows with no high-confidence hole ≥ 500 bp, out of 187 surviving exclusion — `results/windows.tsv` |
| SNP | `chr20:40,150,500 T>A` | `results/sites.tsv` |
| deletion | `chr20:40,153,000 del G` (VCF `POS=40152999 REF=AG ALT=A`) | `results/sites.tsv` |

HG001 is homozygous reference at both sites. The benchmark VCF has no record at
either one, checked with `quay.io/biocontainers/bcftools:1.24--h118bc1c_2`:

    bcftools view -H -r chr20:40150500,chr20:40152999-40153000 HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz

## Sources

All are version-pinned, so a re-run reads the same inputs:

- hg38 sequence and the `gap`, `genomicSuperDups` and `ncbiRefSeqCurated` tracks — UCSC REST API
- high-confidence intervals — GIAB HG001 v4.2.1 benchmark BED, version in the URL
- HG001 genotypes — GIAB HG001 v4.2.1 benchmark VCF, beside the BED
