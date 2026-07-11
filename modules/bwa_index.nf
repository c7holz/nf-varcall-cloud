process BWA_INDEX {
    tag "$fasta"
    publishDir "${params.outdir}/reference", mode: 'copy'

    input:
    path fasta

    output:
    path "bwa", emit: index

    script:
    """
    mkdir bwa
    touch bwa/${fasta}.0123 bwa/${fasta}.amb bwa/${fasta}.ann \
          bwa/${fasta}.bwt.2bit.64 bwa/${fasta}.pac
    """

    stub:
    """
    mkdir bwa
    touch bwa/${fasta}.0123 bwa/${fasta}.amb bwa/${fasta}.ann \
          bwa/${fasta}.bwt.2bit.64 bwa/${fasta}.pac
    """
}
