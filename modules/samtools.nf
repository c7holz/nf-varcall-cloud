process SAMTOOLS {
    tag "$meta.id"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.sorted.bam"), path("*.sorted.bam.bai"), emit: bam

    script:
    """
    touch ${meta.id}.sorted.bam ${meta.id}.sorted.bam.bai
    """

    stub:
    """
    touch ${meta.id}.sorted.bam ${meta.id}.sorted.bam.bai
    """
}
