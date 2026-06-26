process SNPSIFT_FILTER {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/snpeff:5.1--hdfd78af_0"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.vcf.filter.gz")    , emit: vcf

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''


    """
    SnpSift \\
        filter \\
        ${args} \\
        ${vcf} \\
        | bgzip -c > ${prefix}.vcf.filter.gz

    tabix -p vcf ${prefix}.vcf.filter.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | bgzip -c > ${prefix}.vcf.filter.gz
    echo "" | gzip > ${prefix}.vcf.filter.gz.tbi
    """
}
