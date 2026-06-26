process GATK4_MERGEVCFS {
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
    tuple val(meta), path(vcf)
    // tuple val(meta2), path(dict)

    output:
    tuple val(meta), path('*.merge.vcf.gz'), emit: vcf
    // tuple val(meta), path("*.merge.vcf.gz.tbi")   , emit: tbi
    path  "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input_list = vcf.collect{ "--INPUT $it"}.join(' ')
    // def reference_command = dict ? "--SEQUENCE_DICTIONARY $dict" : ""

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK MergeVcfs] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        MergeVcfs \\
        $input_list \\
        --OUTPUT ${prefix}.merge.vcf.gz \\
        --TMP_DIR . \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.merge.vcf.gz


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
