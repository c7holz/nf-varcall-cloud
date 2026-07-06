process BCFTOOLS {
    tag "$meta.id"

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf

    stub:
    "touch ${meta.id}.vcf.gz"
}
