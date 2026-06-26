process DRAGMAP_SIMPLE {
    tag "${meta.id}"
    label 'process_high'
    conda "${moduleDir}/environment.yml"
    container "gambalab/dragmap:latest"

    publishDir(
    path: { "${params.outdir}/${meta.id}/dragmap" },  
    mode: 'symlink',       
    overwrite: true,       
    createDir: true        
    )


    input:
    tuple val(meta), path(reads)
    path(hashmap)
    // tuple val(meta3), path(fasta)
    val sort_bam

    output:
    tuple val(meta), path("*.sam"),  emit: sam,  optional: true
    tuple val(meta), path("*.bam"),  emit: bam,  optional: true
    tuple val(meta), path("*.cram"), emit: cram, optional: true
    tuple val(meta), path("*.crai"), emit: crai, optional: true
    tuple val(meta), path("*.csi"),  emit: csi,  optional: true
    tuple val(meta), path('*.log'),  emit: log
    path "versions.yml",             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // def args2 = task.ext.args2 ?: ''
    def args2 = '-O bam'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reads_command = meta.single_end ? "-1 ${reads}" : "-1 ${reads[0]} -2 ${reads[1]}"
    def samtools_command = sort_bam ? 'samtools sort' : 'samtools view'
    def extension_pattern = /(--output-fmt|-O)+\s+(\S+)/
    def extension_matcher = (args2 =~ extension_pattern)
    def extension = extension_matcher.getCount() > 0 ? extension_matcher[0][2].toLowerCase() : "bam"
    // def reference = fasta && extension == "cram" ? "--reference ${fasta}" : ""
    // if (!fasta && extension == "cram") {
    //     error("Fasta reference is required for CRAM output")
    // }

    """
    dragen-os \\
        -r ${hashmap} \\
        ${args} \\
        --num-threads ${task.cpus} \\
        ${reads_command} \\
        2> >(tee ${prefix}.dragmap.log >&2) \\
        |  ${samtools_command} ${args2} --threads ${task.cpus} -o ${prefix}.${extension} - 
        

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragmap: \$(echo \$(dragen-os --version 2>&1))
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """

    stub:
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension_pattern = /(--output-fmt|-O)+\s+(\S+)/
    def extension_matcher = (args2 =~ extension_pattern)
    def extension = extension_matcher.getCount() > 0 ? extension_matcher[0][2].toLowerCase() : "bam"
    // if (!fasta && extension == "cram") {
    //     error("Fasta reference is required for CRAM output")
    // }

    def create_index = ""
    if (extension == "cram") {
        create_index = "touch ${prefix}.crai"
    }
    else if (extension == "bam") {
        create_index = "touch ${prefix}.bai"
    }

    """
    touch ${prefix}.${extension}
    ${create_index}
    touch ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragmap: \$(echo \$(dragen-os --version 2>&1))
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """
}
