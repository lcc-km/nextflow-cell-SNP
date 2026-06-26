process GATK4_APPLYBQSR {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/gatk4:4.2.6.1--hdfd78af_0"

    publishDir(
        path: { "${params.outdir}/${meta.id}/GATK" },  
        mode: 'symlink',       
        overwrite: true,       
        createDir: true        
    )
    
    input:
    tuple val(meta), path(input), path(input_index), path(bqsr_table), path(intervals)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(dict)

    output:
    tuple val(meta), path("${prefix}.bam"),  emit: bam,  optional: true
    tuple val(meta), path("${prefix}.bai"),  emit: bai,  optional: true
    tuple val(meta), path("${prefix}.cram"), emit: cram, optional: true
    tuple val(meta), path("${prefix}.crai"), emit: crai, optional: true
    path "versions.yml",                     emit: versions

    when:
    task.ext.when == null || task.ext.when


    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    // suffix can only be bam or cram, cram being the sensible default
    def suffix = task.ext.suffix && task.ext.suffix == "bam" ? "cram" : "bam"
    def interval_command = intervals ? "--intervals ${intervals}" : ""

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATK ApplyBQSR] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        ApplyBQSR \\
        --input ${input} \\
        --output ${prefix}.${suffix} \\
        --reference ${fasta} \\
        --bqsr-recal-file ${bqsr_table} \\
        ${interval_command} \\
        --tmp-dir . \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = task.ext.suffix ?: "bam"
    """
    touch ${prefix}.${suffix}
    if [[ ${suffix} == cram ]]; then
        touch ${prefix}.cram.crai
    else
        touch ${prefix}.bam.bai
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
