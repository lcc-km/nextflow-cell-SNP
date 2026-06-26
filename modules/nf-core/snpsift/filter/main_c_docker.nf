process SNPSIFT_FILTER {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
//    container "staphb/snpeff:5.4c"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a1/a116bb44e388ca83fea78d82fe8bdfd5cf3557254e2ec7dd3f1f17354880638c/data' :
        'community.wave.seqera.io/library/htslib_snpsift:ace461dff1cfc121' }"

    publishDir(
        path: { "${params.outdir}/${meta.id}/snpeff" },  
        mode: 'symlink',       
        overwrite: true,       
        createDir: true        
    )

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.vcf.filter.gz")    , emit: vcf
    tuple val(meta), path("*.vcf.filter.gz.tbi")    , emit: tbi

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
