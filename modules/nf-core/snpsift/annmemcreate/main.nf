process SNPSIFT_ANNMEMCREATE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "staphb/snpeff:5.4c"

    input:
    tuple val(meta), path(db_vcf), path(db_vcf_tbi), val(db_fields)

    output:
    tuple val(meta), path("*.snpsift.vardb"), emit: database
//    tuple val("${task.process}"), val('snpsift'), eval("SnpSift -version 2>&1 | grep -oE '[0-9]+\\.[0-9]+[a-z]?'"), emit: versions_snpsift, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def fields = db_fields instanceof List ? db_fields.join(',') : db_fields

    """
    SnpSift \\
        annmem \\
        -create \\
        ${args} \\
        -dbfile ${db_vcf} \\
        ${fields ? "-fields ${fields}" : ""}
    """

    stub:
    """
    mkdir -p ${db_vcf}.snpsift.vardb
    touch ${db_vcf}.snpsift.vardb/chr1.snpsift.df
    """
}
