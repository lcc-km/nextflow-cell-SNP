process GATK4_APPLYVQSR {
    tag "$meta.id"
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
    tuple val(meta), path(vcf), path(vcf_tbi), path(recal), path(recal_index), path(tranches)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(dict)
    val mode            //"SNP" or "INDEL"


    output:
    tuple val(meta), path("*.VQSR.vcf.gz"), emit: vcf
    tuple val(meta), path("*.VQSR.vcf.gz.tbi")   , emit: tbi
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference_command = fasta ? "--reference $fasta" : ''

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK ApplyVQSR] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        ApplyVQSR \\
        --variant ${vcf} \\
        --output ${prefix}.${mode}.VQSR.vcf.gz \\
        $reference_command \\
        --tranches-file $tranches \\
        --recal-file $recal \\
        -mode $mode \\
        --tmp-dir . \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.${mode}.VQSR.vcf.gz
    touch ${prefix}.${mode}.VQSR.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
