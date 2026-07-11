process FASTP {
    tag "$meta.id"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.trimmed.fastq.gz"), emit: reads
    tuple val(meta), path("*.fastp.json"),       emit: json
    tuple val(meta), path("*.fastp.html"),       emit: html

    script:
    """
    touch ${meta.id}_1.trimmed.fastq.gz ${meta.id}_2.trimmed.fastq.gz
    touch ${meta.id}.fastp.json ${meta.id}.fastp.html
    """
    
    stub:
    """
    touch ${meta.id}_1.trimmed.fastq.gz ${meta.id}_2.trimmed.fastq.gz
    touch ${meta.id}.fastp.json ${meta.id}.fastp.html
    """
}
}
