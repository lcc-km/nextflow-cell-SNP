// 简化版流程

process FASTP {
    tag "$meta.id"
    container "pgt_a_clean:20231120"
    label 'process_high'
    // cpus 16

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.fastp.fastq.gz') , optional:true, emit: reads
    tuple val(meta), path('*.json')           , emit: json
    tuple val(meta), path('*.html')           , emit: html
    tuple val(meta), path('*.log')            , emit: log
    tuple val(meta), path('*.fail.fastq.gz')  , optional:true, emit: reads_fail

    script:
    def prefix = "${meta.id}"
    def args = task.ext.args ?: ''

    "/toolkit/fastp \
        --in1 ${reads[0]} \
        --in2 ${reads[1]} \
        --out1 ${prefix}_1.fastp.fastq.gz \
        --out2 ${prefix}_2.fastp.fastq.gz \
        --json ${prefix}.fastp.json \
        --html ${prefix}.fastp.html \
        --detect_adapter_for_pe \
        --thread $task.cpus" 
}

