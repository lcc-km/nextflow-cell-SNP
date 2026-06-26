process SNPSIFT_ANNOTATE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "staphb/snpeff:5.4c"

    publishDir(
        path: { "${params.outdir}/${meta.id}/snpeff" },  
        mode: 'symlink',       
        overwrite: true,       
        createDir: true        
    )

    input:
    tuple val(meta), path(vcf)
    path dbsnp
    path clinvar
    path gwas_cat
    path dbnsfp

    output:
    tuple val(meta), path("*.vcf.filter_annotate.gz")    , emit: vcf

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    def prefix   = task.ext.prefix ?: "${meta.id}"

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
