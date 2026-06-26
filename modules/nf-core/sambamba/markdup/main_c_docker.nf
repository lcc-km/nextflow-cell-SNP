process SAMBAMBA_MARKDUP {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "pgt_a_clean:20231120"

    publishDir(
    path: { "${params.outdir}/${meta.id}/dragmap" },  
    mode: 'symlink',       
    overwrite: true,       
    createDir: true        
    )

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.rmdup.bam"), emit: rmdup_bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    /toolkit/sambamba \\
        markdup \\
        $args \\
        -t $task.cpus \\
        --tmpdir ./ \\
        $bam \\
        ${prefix}.rmdup.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sambamba: \$(echo \$(sambamba --version 2>&1) | awk '{print \$2}' )
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.rmdup.bam
        cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sambamba: \$(echo \$(sambamba --version 2>&1) | awk '{print \$2}' )
    END_VERSI
    """
}

