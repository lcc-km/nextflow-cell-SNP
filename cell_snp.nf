// lucc cell SNP

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS - 基因组属性获取
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        return params.genomes[params.genome][attribute] ?: null
    }
    return null
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
params.fasta     = getGenomeAttribute('fasta')
params.fasta_fai = getGenomeAttribute('fasta_fai')
params.fasta_dict = getGenomeAttribute('fasta_dict')


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTP_SIMPLE } from './modules/nf-core/fastp/main_simple.nf'
include { DRAGMAP_SIMPLE } from './modules/nf-core/dragmap/align/main_simple.nf'
include { GATK4_MARKDUPLICATES  } from './modules/nf-core/gatk4/markduplicates/main_simple.nf'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_MARKDUP } from './modules/nf-core/samtools/index/main_c_docker.nf'
include { SAMBAMBA_MARKDUP } from './modules/nf-core/sambamba/markdup/main_c_docker.nf'
include { SAMTOOLS_CONVERT as BAM_TO_CRAM } from './modules/nf-core/samtools/convert/main_c_docker.nf'
include { SAMTOOLS_STATS } from './modules/nf-core/samtools/stats/main_c_docker.nf'
include { MOSDEPTH } from './modules/nf-core/mosdepth/main_c_docker.nf'

include { GATK4_BASERECALIBRATOR  } from './modules/nf-core/gatk4/baserecalibrator/main_c_docker.nf'
include { GATK4_APPLYBQSR } from './modules/nf-core/gatk4/applybqsr/main_c_docker.nf'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_BQSR }    from './modules/nf-core/samtools/index/main_c_docker.nf'
include { GATK4_HAPLOTYPECALLER } from './modules/nf-core/gatk4/haplotypecaller/main_c_docker.nf'
include { GATK4_VARIANTRECALIBRATOR as VariantRecalibrator_SNP} from './modules/nf-core/gatk4/variantrecalibrator/main_c_docker.nf'
include { GATK4_VARIANTRECALIBRATOR as VariantRecalibrator_INDEL} from './modules/nf-core/gatk4/variantrecalibrator/main_c_docker.nf'
include { GATK4_APPLYVQSR as APPLYVQSR_SNP} from './modules/nf-core/gatk4/applyvqsr/main_c_docker.nf'
include { GATK4_APPLYVQSR as APPLYVQSR_INDEL} from './modules/nf-core/gatk4/applyvqsr/main_c_docker.nf'
include { GATK4_MERGEVCFS } from './modules/nf-core/gatk4/mergevcfs/main_c_docker.nf'


include { SNPEFF_DOWNLOAD } from './modules/nf-core/snpeff/download/main_c_docker.nf'
include { SNPSIFT_ANNMEMCREATE } from './modules/nf-core/snpsift/annmemcreate/main.nf'
include { SNPEFF_SNPEFF } from './modules/nf-core/snpeff/snpeff/main_c_docker.nf'
include { SNPSIFT_FILTER } from './modules/nf-core/snpsift/filter/main_c_docker.nf'
include { SNPSIFT_ANNMEM } from './modules/nf-core/snpsift/annmem/main.nf'



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SAMPLE PARSING - 样本解析
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow SAMPLE_PARSING {
    take: input_file
    main:
        ch_samples = channel.fromPath(input_file)
            .splitCsv(header:true, sep: ',')
            .map { row ->
                def meta = [
                    id: row.sample,
                    sample: row.sample,
                    patient: row.patient,
                    lane: row.lane,
                    single_end: (row.fastq_2 == null || row.fastq_2.trim() == '')
                ]
                def reads = meta.single_end ? [file(row.fastq_1)] : [file(row.fastq_1), file(row.fastq_2)]
                [meta, reads]
            }
    emit: ch_samples
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN PIPELINE - 【纯净：仅处理样本 + 调用流程】
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow PIPELINE {
    main:
        // 1. 样本解析
        SAMPLE_PARSING(params.input)
        FASTP_SIMPLE(SAMPLE_PARSING.out.ch_samples)

        // 2. DRAGMAP 
        hash_channel = Channel.value(file(params.dragmap_hash, checkIfExists: true))  // 改为 Value Channel
        DRAGMAP_SIMPLE(FASTP_SIMPLE.out.reads, hash_channel, true)                     // 加 .out.reads

        // 3.去重
        fasta_channel = Channel.value( [ [:], file(params.fasta, checkIfExists: true) ] )
        fai_channel   = Channel.value( [ [:], file(params.fasta_fai, checkIfExists: true) ] )
        // GATK4_MARKDUPLICATES(DRAGMAP_SIMPLE.out.bam, fasta_channel, fai_channel)       // 加 .out.bam
        SAMBAMBA_MARKDUP(DRAGMAP_SIMPLE.out.bam)
        // 4  CONVERT + QC 
        SAMTOOLS_INDEX_MARKDUP(SAMBAMBA_MARKDUP.out.rmdup_bam)
        ch_bai_markdup = SAMTOOLS_INDEX_MARKDUP.out.bai 
        rmdupbam_index_channel = SAMBAMBA_MARKDUP.out.rmdup_bam.join(ch_bai_markdup)  // 按meta配对，合并为(meta, bam, bai)
        BAM_TO_CRAM(rmdupbam_index_channel,fasta_channel,fai_channel)
        SAMTOOLS_STATS(rmdupbam_index_channel,fasta_channel)

        probe_bed = Channel.value(file(params.probe_bed, checkIfExists: true))
        rmdupbam_index_bed = rmdupbam_index_channel.combine(probe_bed)
        MOSDEPTH(rmdupbam_index_bed,fasta_channel)

        // 5 GATK
        fa_dict_channel = Channel.value( [ [:], file(params.fasta_dict, checkIfExists: true) ] )

    known_sites_channel = Channel.value( [ 
        [ id: 'hg38_known_sites' ],  // 给 meta 加个 id，方便看日志
        [
        file(params.dbsnp_146, checkIfExists: true),
        file(params.indels_1000G, checkIfExists: true),
        file(params.phase1_1000G, checkIfExists: true)
        ] 
        ])

    // 索引文件的 Channel，顺序必须和上面的 VCF 完全一致！
    known_sites_tbi_channel = Channel.value( [ 
        [ id: 'hg38_known_sites_tbi' ], 
        [
        file("${params.dbsnp_146}.tbi", checkIfExists: true),
        file("${params.indels_1000G}.tbi", checkIfExists: true),
        file("${params.phase1_1000G}.tbi", checkIfExists: true)
        ] 
        ] )

        GATK4_BASERECALIBRATOR(
            rmdupbam_index_bed, 
            fasta_channel,
            fai_channel, 
            fa_dict_channel, 
            known_sites_channel, 
            known_sites_tbi_channel)

            bam_recal = rmdupbam_index_channel.join(GATK4_BASERECALIBRATOR.out.table) //join：按样本 meta 配对
            APPLYBQSR_channel = bam_recal.combine(probe_bed) //combine：给所有样本添加全局的探针 bed 文件

        GATK4_APPLYBQSR(
            APPLYBQSR_channel,
            fasta_channel,
            fai_channel,
            fa_dict_channel
        )

        SAMTOOLS_INDEX_BQSR(GATK4_APPLYBQSR.out.bam)
        ch_bai_bqsr = SAMTOOLS_INDEX_BQSR.out.bai 

        bqsr_bam_channel = GATK4_APPLYBQSR.out.bam.join(ch_bai_bqsr)
        hc_input_channel = bqsr_bam_channel.combine(probe_bed).map{ meta, bam, bai, bed -> 
            tuple(meta, bam, bai, bed) }

        known_dbsnp = Channel.value( [  
            [ id: 'dbsnp_146' ],
            [file(params.dbsnp_146, checkIfExists: true)] 
        ])
        known_dbsnp_tbi = Channel.value( [  
            [ id: 'dbsnp_146' ],
            [file("${params.dbsnp_146}.tbi", checkIfExists: true)] 
        ])

        GATK4_HAPLOTYPECALLER(hc_input_channel, fasta_channel, fai_channel, fa_dict_channel, known_dbsnp, known_dbsnp_tbi)

        // GATK4_VARIANTRECALIBRATOR SNP
        VariantRecalibrator_vcf_channel = GATK4_HAPLOTYPECALLER.out.vcf.join(GATK4_HAPLOTYPECALLER.out.tbi)

        resource_files_snp = Channel.value( [ 
        file("${params.hapmap}", checkIfExists: true),
        file("${params.omni_1000G}", checkIfExists: true),
        file("${params.phase1_1000G}", checkIfExists: true),
        file("${params.dbsnp_146}", checkIfExists: true),
        ] )

        resource_indexes_snp = Channel.value( [ 
        file("${params.hapmap}.tbi", checkIfExists: true),
        file("${params.omni_1000G}.tbi", checkIfExists: true),
        file("${params.phase1_1000G}.tbi", checkIfExists: true),
        file("${params.dbsnp_146}.tbi", checkIfExists: true),
        ] )
        
        def labels_list_snp = [
            "--resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${params.hapmap}",
            "--resource:omni,known=false,training=true,truth=false,prior=12.0 ${params.omni_1000G}", 
            "--resource:1000G,known=false,training=true,truth=false,prior=10.0 ${params.phase1_1000G}",
            "--resource:dbsnp,known=true,training=false,truth=false,prior=6.0 ${params.dbsnp_146}"
        ]

        VariantRecalibrator_SNP(
            VariantRecalibrator_vcf_channel, 
            resource_files_snp,
            resource_indexes_snp,
            labels_list_snp,
            fasta_channel,
            fai_channel,
            fa_dict_channel,
            "SNP")

        // VariantRecalibrator_SNP.ext.args = args_snp

        resource_files_indel =  Channel.value([
            file(params.indels_1000G, checkIfExists: true) 
        ])

        resource_indexes_indel =  Channel.value([
            file("${params.indels_1000G}.tbi", checkIfExists: true)
        ])

        def labels_list_indel = [
            "--resource:mills,known=true,training=true,truth=true,prior=12.0 ${params.indels_1000G}"
        ]

        VariantRecalibrator_INDEL(
            VariantRecalibrator_vcf_channel, // 输入依然是 HC 出来的原始 VCF
            resource_files_indel,
            resource_indexes_indel,
            labels_list_indel,
            fasta_channel,
            fai_channel,
            fa_dict_channel,
            "INDEL")

         APPLYVQSR_SNP_inputs = VariantRecalibrator_vcf_channel
            .join(VariantRecalibrator_SNP.out.recal)     // 按 meta 合并 recal
            .join(VariantRecalibrator_SNP.out.idx)       // 按 meta 合并 idx
            .join(VariantRecalibrator_SNP.out.tranches)  // 按 meta 合并 tranches
            .map { meta, vcf, vcf_tbi, recal, recal_idx, tranches ->
                // 重新排列成 APPLYVQSR 输入需要的顺序
                tuple(meta, vcf, vcf_tbi, recal, recal_idx, tranches)
                }

        APPLYVQSR_SNP(
            APPLYVQSR_SNP_inputs,
            fasta_channel,
            fai_channel,
            fa_dict_channel,
            "SNP"
        )
        
         APPLYVQSR_INDEL_inputs = VariantRecalibrator_vcf_channel
            .join(VariantRecalibrator_INDEL.out.recal)     // 按 meta 合并 recal
            .join(VariantRecalibrator_INDEL.out.idx)       // 按 meta 合并 idx
            .join(VariantRecalibrator_INDEL.out.tranches)  // 按 meta 合并 tranches
            .map { meta, vcf, vcf_tbi, recal, recal_idx, tranches ->
                // 重新排列成 APPLYVQSR 输入需要的顺序
                tuple(meta, vcf, vcf_tbi, recal, recal_idx, tranches)
                }

        APPLYVQSR_INDEL(
            APPLYVQSR_INDEL_inputs,
            fasta_channel,
            fai_channel,
            fa_dict_channel,
            "INDEL"
        )

        VQSR_out_merge_vcf = APPLYVQSR_SNP.out.vcf.join(APPLYVQSR_INDEL.out.vcf).map { meta, snp_vcf, indel_vcf -> tuple(meta, [snp_vcf, indel_vcf])}
        // VQSR_out_merge_dict = APPLYVQSR_SNP.out.tbi.join(APPLYVQSR_INDEL.out.tbi).map { meta, snp_tbi, indel_tbi -> tuple(meta, [snp_tbi, indel_tbi])}
        GATK4_MERGEVCFS(VQSR_out_merge_vcf)

        // snpsift
        // database download
        def snpeff_cache_parent = file(params.snpeff_cache)
        def snpeff_db_dir = file("${params.snpeff_cache}/${params.snpeff_db}")
        def need_download_snpeff = params.download_cache && (!snpeff_db_dir.exists() || params.force_download_cache)

        def snpeff_meta = [ id: params.snpeff_db ]
        def snpeff_cache_channel

        if (need_download_snpeff) {
            def snpeff_download_in = Channel.value([ snpeff_meta, params.snpeff_db ])
            
            SNPEFF_DOWNLOAD(snpeff_download_in)
            
            snpeff_cache_channel = SNPEFF_DOWNLOAD.out.cache
        } else {
            if (!snpeff_db_dir.exists()) {
                error "SnpEff database not found at ${snpeff_db_dir} and download is disabled."
            }
            snpeff_cache_channel = Channel.value([ snpeff_meta, snpeff_cache_parent ])
        }

        // vcf database check
        def has_prebuilt_vardb = false
        if (params.snpsift_merge_ann_dbfile && file(params.snpsift_merge_ann_dbfile).exists()) {
            has_prebuilt_vardb = true
        } else {
            // 如果没有指定统一路径，检查 VCF 原路径下是否已包含对应的 .snpsift.vardb 目录
            if (file("${params.dbSNP}.snpsift.vardb").exists() && 
                file("${params.freq}.snpsift.vardb").exists() && 
                file("${params.clinvar}.snpsift.vardb").exists()) {
                has_prebuilt_vardb = true
            }
        }

        def ann_db_channel

        if (has_prebuilt_vardb) {
            // 逻辑 A: 数据库已存在，直接打包成 SNPSIFT_ANNMEM 要求的单一 Tuple 通道
            def db_vcf_list = [file(params.dbSNP, checkIfExists: true), file(params.freq, checkIfExists: true), file(params.clinvar, checkIfExists: true)]
            def db_tbi_list = [file("${params.dbSNP}.tbi", checkIfExists: true), file("${params.freq}.tbi", checkIfExists: true), file("${params.clinvar}.tbi", checkIfExists: true)]
            
            def db_vardb_list = params.snpsift_merge_ann_dbfile ? 
                [file("${params.snpsift_merge_ann_dbfile}/dbSNP.snpsift.vardb"), file("${params.snpsift_merge_ann_dbfile}/freq.snpsift.vardb"), file("${params.snpsift_merge_ann_dbfile}/clinvar.snpsift.vardb")] : 
                [file("${params.dbSNP}.snpsift.vardb"), file("${params.freq}.snpsift.vardb"), file("${params.clinvar}.snpsift.vardb")]
                
            def db_fields_list   = ['', '', '']       // 如果有特定 fields 要求，可在此填入
            def db_prefixes_list = ['dbSNP', 'freq', 'clinvar']

            // 拼装成一整条带有 List 的 Value Channel
            ann_db_channel = Channel.value([db_vcf_list, db_tbi_list, db_vardb_list, db_fields_list, db_prefixes_list])

        } else {
            // 逻辑 B: 数据库不存在，动态触发 SNPSIFT_ANNMEMCREATE 分别对 3 个文件建库
            def create_input_ch = Channel.fromList([
                [ [id: 'dbSNP', prefix: 'dbSNP', fields: ''], file(params.dbSNP, checkIfExists: true), file("${params.dbSNP}.tbi", checkIfExists: true) ],
                [ [id: 'freq', prefix: 'freq', fields: ''], file(params.freq, checkIfExists: true), file("${params.freq}.tbi", checkIfExists: true) ],
                [ [id: 'clinvar', prefix: 'clinvar', fields: ''], file(params.clinvar, checkIfExists: true), file("${params.clinvar}.tbi", checkIfExists: true) ]
            ])

            // 转换格式以匹配 SNPSIFT_ANNMEMCREATE 的 input 声明
            def annmemcreate_in = create_input_ch.map { meta, vcf, tbi -> [ meta, vcf, tbi, meta.fields ] }
            
            // 执行建库（并行产生 3 个元组结果）
            SNPSIFT_ANNMEMCREATE(annmemcreate_in)

            // 将建库生成的 vardb 文件夹与原 VCF、TBI 重新聚拢，转置为 SNPSIFT_ANNMEM 所需的列表结构
            ann_db_channel = SNPSIFT_ANNMEMCREATE.out.database
                .join(create_input_ch) // 通过 meta 隐式匹配将其 join 起来
                .map { meta, vardb, vcf, tbi -> 
                    [ vcf, tbi, vardb, meta.fields, meta.prefix ]
                }
                .toList() // 将 3 个独立异步的流水线收集合并为一个全量 List
                .map { list_of_items ->
                    // 把转置前: [ [vcf1, tbi1, vardb1...], [vcf2, tbi2, vardb2...] ]
                    // 转换为后: [ [vcf1, vcf2], [tbi1, tbi2], [vardb1, vardb2]... ] 从而完美对应多数据库输入
                    def vcfs     = list_of_items.collect { it[0] }
                    def tbis     = list_of_items.collect { it[1] }
                    def vardbs   = list_of_items.collect { it[2] }
                    def fields   = list_of_items.collect { it[3] }
                    def prefixes = list_of_items.collect { it[4] }
                    return [ vcfs, tbis, vardbs, fields, prefixes ]
                }
        }

        SNPSIFT_FILTER(GATK4_MERGEVCFS.out.vcf)
//        SNPEFF_SNPEFF(SNPSIFT_FILTER.out.vcf, params.snpeff_db, snpeff_cache_channel)
        vcf_tbi_ann_channel = SNPSIFT_FILTER.out.vcf.join(SNPSIFT_FILTER.out.tbi)
        SNPSIFT_ANNMEM(vcf_tbi_ann_channel, ann_db_channel)
       
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    流程完成通知
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow PIPELINE_COMPLETION {
    take: outdir
    main:
        log.info "✅ Pipeline completed successfully! Results: ${outdir}"
}


workflow { PIPELINE() }