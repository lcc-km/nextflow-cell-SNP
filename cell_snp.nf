// lucc cell SNP
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS - Genome Attribute Retrieval
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
    SAMPLE PARSING - Sample Parsing
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
    MAIN PIPELINE - [Clean: sample processing + workflow invocation only]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow PIPELINE {
    main:
        // 1. Sample parsing
        SAMPLE_PARSING(params.input)
        FASTP_SIMPLE(SAMPLE_PARSING.out.ch_samples)
        // 2. DRAGMAP alignment
        hash_channel = Channel.value(file(params.dragmap_hash, checkIfExists: true))  // Convert to Value Channel
        DRAGMAP_SIMPLE(FASTP_SIMPLE.out.reads, hash_channel, true)                     // Add .out.reads
        // 3. Deduplication
        fasta_channel = Channel.value( [ [:], file(params.fasta, checkIfExists: true) ] )
        fai_channel   = Channel.value( [ [:], file(params.fasta_fai, checkIfExists: true) ] )
        // GATK4_MARKDUPLICATES(DRAGMAP_SIMPLE.out.bam, fasta_channel, fai_channel)       // Add .out.bam
        SAMBAMBA_MARKDUP(DRAGMAP_SIMPLE.out.bam)
        // 4. CONVERT + QC
        SAMTOOLS_INDEX_MARKDUP(SAMBAMBA_MARKDUP.out.rmdup_bam)
        ch_bai_markdup = SAMTOOLS_INDEX_MARKDUP.out.bai
        rmdupbam_index_channel = SAMBAMBA_MARKDUP.out.rmdup_bam.join(ch_bai_markdup)  // Pair by meta, merge into (meta, bam, bai)
        BAM_TO_CRAM(rmdupbam_index_channel,fasta_channel,fai_channel)
        SAMTOOLS_STATS(rmdupbam_index_channel,fasta_channel)
        probe_bed = Channel.value(file(params.probe_bed, checkIfExists: true))
        rmdupbam_index_bed = rmdupbam_index_channel.combine(probe_bed)
        MOSDEPTH(rmdupbam_index_bed,fasta_channel)
        // 5. GATK processing
        fa_dict_channel = Channel.value( [ [:], file(params.fasta_dict, checkIfExists: true) ] )
    known_sites_channel = Channel.value( [
        [ id: 'hg38_known_sites' ],  // Add an id to meta for easier log reading
        [
        file(params.dbsnp_146, checkIfExists: true),
        file(params.indels_1000G, checkIfExists: true),
        file(params.phase1_1000G, checkIfExists: true)
        ]
        ])
    // Index file Channel, order must match the VCFs above exactly!
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
            bam_recal = rmdupbam_index_channel.join(GATK4_BASERECALIBRATOR.out.table) // join: pair by sample meta
            APPLYBQSR_channel = bam_recal.combine(probe_bed) // combine: add global probe BED file to all samples
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
            VariantRecalibrator_vcf_channel, // Input is still the raw VCF from HC
            resource_files_indel,
            resource_indexes_indel,
            labels_list_indel,
            fasta_channel,
            fai_channel,
            fa_dict_channel,
            "INDEL")
         APPLYVQSR_SNP_inputs = VariantRecalibrator_vcf_channel
            .join(VariantRecalibrator_SNP.out.recal)     // Merge recal by meta
            .join(VariantRecalibrator_SNP.out.idx)       // Merge idx by meta
            .join(VariantRecalibrator_SNP.out.tranches)  // Merge tranches by meta
            .map { meta, vcf, vcf_tbi, recal, recal_idx, tranches ->
                // Rearrange into the order required by APPLYVQSR input
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
            .join(VariantRecalibrator_INDEL.out.recal)     // Merge recal by meta
            .join(VariantRecalibrator_INDEL.out.idx)       // Merge idx by meta
            .join(VariantRecalibrator_INDEL.out.tranches)  // Merge tranches by meta
            .map { meta, vcf, vcf_tbi, recal, recal_idx, tranches ->
                // Rearrange into the order required by APPLYVQSR input
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
            // If no unified path is specified, check if the corresponding .snpsift.vardb directories exist under the original VCF paths
            if (file("${params.dbSNP}.snpsift.vardb").exists() &&
                file("${params.freq}.snpsift.vardb").exists() &&
                file("${params.clinvar}.snpsift.vardb").exists()) {
                has_prebuilt_vardb = true
            }
        }
        def ann_db_channel
        if (has_prebuilt_vardb) {
            // Logic A: databases already exist, directly package into the single Tuple channel required by SNPSIFT_ANNMEM
            def db_vcf_list = [file(params.dbSNP, checkIfExists: true), file(params.freq, checkIfExists: true), file(params.clinvar, checkIfExists: true)]
            def db_tbi_list = [file("${params.dbSNP}.tbi", checkIfExists: true), file("${params.freq}.tbi", checkIfExists: true), file("${params.clinvar}.tbi", checkIfExists: true)]

            def db_vardb_list = params.snpsift_merge_ann_dbfile ?
                [file("${params.snpsift_merge_ann_dbfile}/dbSNP.snpsift.vardb"), file("${params.snpsift_merge_ann_dbfile}/freq.snpsift.vardb"), file("${params.snpsift_merge_ann_dbfile}/clinvar.snpsift.vardb")] :
                [file("${params.dbSNP}.snpsift.vardb"), file("${params.freq}.snpsift.vardb"), file("${params.clinvar}.snpsift.vardb")]

            def db_fields_list   = ['', '', '']       // Fill in here if specific fields are required
            def db_prefixes_list = ['dbSNP', 'freq', 'clinvar']
            // Assemble into a single Value Channel with Lists
            ann_db_channel = Channel.value([db_vcf_list, db_tbi_list, db_vardb_list, db_fields_list, db_prefixes_list])
        } else {
            // Logic B: databases do not exist, dynamically trigger SNPSIFT_ANNMEMCREATE to build databases for 3 files separately
            def create_input_ch = Channel.fromList([
                [ [id: 'dbSNP', prefix: 'dbSNP', fields: ''], file(params.dbSNP, checkIfExists: true), file("${params.dbSNP}.tbi", checkIfExists: true) ],
                [ [id: 'freq', prefix: 'freq', fields: ''], file(params.freq, checkIfExists: true), file("${params.freq}.tbi", checkIfExists: true) ],
                [ [id: 'clinvar', prefix: 'clinvar', fields: ''], file(params.clinvar, checkIfExists: true), file("${params.clinvar}.tbi", checkIfExists: true) ]
            ])
            // Convert format to match SNPSIFT_ANNMEMCREATE input declaration
            def annmemcreate_in = create_input_ch.map { meta, vcf, tbi -> [ meta, vcf, tbi, meta.fields ] }

            // Execute database building (produces 3 tuple results in parallel)
            SNPSIFT_ANNMEMCREATE(annmemcreate_in)
            // Re-group the vardb folders generated by database building with original VCF and TBI, transpose to the list structure required by SNPSIFT_ANNMEM
            ann_db_channel = SNPSIFT_ANNMEMCREATE.out.database
                .join(create_input_ch) // Join them through implicit meta matching
                .map { meta, vardb, vcf, tbi ->
                    [ vcf, tbi, vardb, meta.fields, meta.prefix ]
                }
                .toList() // Collect and merge 3 independent asynchronous pipelines into one complete List
                .map { list_of_items ->
                    // Before transposition: [ [vcf1, tbi1, vardb1...], [vcf2, tbi2, vardb2...] ]
                    // After transposition: [ [vcf1, vcf2], [tbi1, tbi2], [vardb1, vardb2]... ] to perfectly match multi-database input
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
    Pipeline Completion Notification
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow PIPELINE_COMPLETION {
    take: outdir
    main:
        log.info "✅ Pipeline completed successfully! Results: ${outdir}"
}
workflow { PIPELINE() }
