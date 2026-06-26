process FASTP_SIMPLE {
    tag "$meta.id"
    label 'process_medium'
    container "pgt_a_clean:20231120"

    input:
    tuple val(meta), path(reads)

    publishDir(
        path: { "${params.outdir}/${meta.id}/fastp" },  
        mode: 'symlink',       
        overwrite: true,       
        createDir: true        
    )


    output:
    tuple val(meta), path('*.fastp.fastq.gz') , emit: reads
    tuple val(meta), path('*.json')           , emit: json
    tuple val(meta), path('*.html')           , emit: html
    tuple val(meta), path('*.log')            , emit: log
    tuple val(meta), path('*.fail.fastq.gz')  , optional:true, emit: reads_fail
    tuple val(meta), path('*.merged.fastq.gz'), optional:true, emit: reads_merged
    path "versions.yml"                       , emit: versions

    script:
    def prefix = "${meta.id}"
    if (meta.single_end) {
        """
        [ ! -f  ${prefix}.fastq.gz ] && ln -sf $reads ${prefix}.fastq.gz
        /toolkit/fastp \
            --in1 ${reads} \
            --out1 ${prefix}.fastp.fastq.gz \
            --json ${prefix}.fastp.json \
            --html ${prefix}.fastp.html \
            --thread $task.cpus
            2> >(tee ${prefix}.fastp.log >&2) 

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(/toolkit/fastp --version 2>&1 | sed -e "s/fastp //g")
        END_VERSIONS
        """
    } else {
        """
        [ ! -f  ${prefix}_1.fastq.gz ] && ln -sf ${reads[0]} ${prefix}_1.fastq.gz
        [ ! -f  ${prefix}_2.fastq.gz ] && ln -sf ${reads[1]} ${prefix}_2.fastq.gz
        /toolkit/fastp \
            --in1 ${reads[0]} \
            --in2 ${reads[1]} \
            --out1 ${prefix}_1.fastp.fastq.gz \
            --out2 ${prefix}_2.fastp.fastq.gz \
            --json ${prefix}.fastp.json \
            --html ${prefix}.fastp.html \
            --detect_adapter_for_pe \
            --thread $task.cpus
            2> >(tee ${prefix}.fastp.log >&2) 

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(/toolkit/fastp --version 2>&1 | sed -e "s/fastp //g")
        END_VERSIONS
        """
    }

    stub:
    def prefix = "${meta.id}"
    if (meta.single_end) {
        """
        echo '' | gzip > ${prefix}.fastp.fastq.gz
        touch "${prefix}.fastp.json"
        touch "${prefix}.fastp.html"
        touch "${prefix}.fastp.log" 

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: "0.23.2"
        END_VERSIONS
        """
    } else {
        """
        echo '' | gzip > ${prefix}_1.fastp.fastq.gz
        echo '' | gzip > ${prefix}_2.fastp.fastq.gz
        touch "${prefix}.fastp.json"
        touch "${prefix}.fastp.html"
        touch "${prefix}.fastp.log" 

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: "0.23.2"
        END_VERSIONS
        """
    }
} 