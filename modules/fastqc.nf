process FASTQC {
    tag "$meta.id"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.zip"),  emit: zip
    tuple val(meta), path("*.html"), emit: html

    script:
    """
    touch ${meta.id}_fastqc.zip
    touch ${meta.id}_fastqc.html
    """

    stub:
    """
    touch ${meta.id}_fastqc.zip
    touch ${meta.id}_fastqc.html
    """
}
