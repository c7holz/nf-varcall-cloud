#!/usr/bin/env nextflow

include { FASTQC }    from './modules/fastqc'
include { FASTP }     from './modules/fastp'
include { BWA_MEM }   from './modules/bwa_mem'
include { BWA_INDEX } from './modules/bwa_index'
include { SAMTOOLS }  from './modules/samtools'
include { BCFTOOLS }  from './modules/bcftools'
include { MULTIQC }   from './modules/multiqc'

params.build_index = false
params.bwa_index   = null

workflow {
    ch_reads = Channel.fromPath(params.input)
        .splitCsv(header: true)
        .map { row -> [ [id: row.sample], [ file(row.fastq_1), file(row.fastq_2) ] ] }

    ch_fasta = Channel.value(file(params.fasta))

    if (params.build_index) {
        ch_index = BWA_INDEX(ch_fasta).index
    } else {
        if (!params.bwa_index) {
            error "No BWA index provided. Pass --bwa_index <dir> with a prebuilt index, " +
                  "or build one from --fasta by adding --build_index."
        }
        ch_index = Channel.value(file(params.bwa_index))
    }

    FASTQC(ch_reads)
    FASTP(ch_reads)
    BWA_MEM(FASTP.out.reads, ch_index)
    SAMTOOLS(BWA_MEM.out.bam)
    BCFTOOLS(SAMTOOLS.out.bam, ch_fasta)

    ch_multiqc = FASTQC.out.zip.map { meta, f -> f }
        .mix( FASTP.out.json.map { meta, f -> f } )
        .collect()
    MULTIQC(ch_multiqc)
}
