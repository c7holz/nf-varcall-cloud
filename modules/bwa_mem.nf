process BWA_MEM {
    tag "$meta.id"

    input:
    tuple val(meta), path(reads)
    path index

    output:
    tuple val(meta), path("*.bam"), emit: bam

    script:
    "touch ${meta.id}.bam"
    
    stub:
    "touch ${meta.id}.bam"
}
