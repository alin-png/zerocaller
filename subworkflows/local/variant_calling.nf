
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { MERGEBAMS                     } from '../../modules/local/mergebams.nf'
include { HAPLOTYPECALLER               } from '../../modules/local/haplotypecaller.nf'
include { MERGEVCFS                     } from '../../modules/local/mergevcfs.nf'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow VARIANT_CALLING {

    take:
        input
        fasta
        intervals

    main:
        ch_versions = Channel.empty()

        //
        // SUBWORKFLOW: Read in samplesheet, validate and stage input files
        //

        Channel
            .fromPath(fasta, checkIfExists:true).first().set{ch_reference_genome}

        // Also include the FASTA index and dictionary files
        Channel
            .fromPath(fasta + '.fai', checkIfExists:true).first().set{ch_reference_genome_fai}
        
        Channel
            .fromPath(fasta.replaceAll('\\.(fa|fasta)$', '.dict'), checkIfExists:true).first().set{ch_reference_genome_dict}

        Channel
            .fromPath(intervals + '*scattered.interval_list', checkIfExists:true)
            .set{ch_intervals}

        intervals_files = file(intervals + '*scattered.interval_list', checkIfExists:true)

        // Process BAM files from samplesheet
        Channel
            .fromPath(input, checkIfExists:true )
            .splitCsv(header: true)
            .map {
                row -> {
                    bam_file = file(row.filepath, checkIfExists:true)
                    bai_file = file(row.filepath + '.bai', checkIfExists:true)
                    [row.sample_id, row.sample_type.toLowerCase(), bam_file, bai_file]
                }
            }
            .groupTuple(by: [0, 1])
            .map {
                sample_id, sample_type, bam_files, bai_files -> 
                    [sample_id, sample_type, bam_files, bai_files]
            }
            .set{ch_meta_complete}

        // Separate tumor and normal samples
        ch_meta_complete
        .map{[it[0], it[1]]}
        .branch {
            tumor: it[1] == "t" || it[1] == "tumor" || it[1] == "0"
            normal: true
        }
        .set{ch_sample_type}

        // Process BAM files - check if multiple BAMs need merging
        ch_meta_complete
        .map{[it[0], it[2], it[3]]}
        .branch {
            singleton: it[1].size() <= 1
            multiple: it[1].size() > 1
        }
        .set {
            ch_sample_bams
        }

        MERGEBAMS(
            ch_sample_bams.multiple
        )

        MERGEBAMS.out.merged_bam
        .mix(
            ch_sample_bams
            .singleton
            .map{[it[0], it[1][0], it[2][0]]}
        ).set{
            ch_all_bams
        }

        ch_all_bams
        .join(ch_sample_type.normal)
        .map{[it[0], it[1], it[2]]}
        .set{ch_normal_bams}


        HAPLOTYPECALLER(
            ch_normal_bams,
            ch_reference_genome,
            ch_reference_genome_fai,
            ch_reference_genome_dict,
            ch_intervals

        )

        HAPLOTYPECALLER.out.gvcf
        .groupTuple(size: intervals_files.size())
        .set{ch_sample_gvcf}

        MERGEVCFS(
            ch_sample_gvcf
        )


    emit:
        bams     = ch_all_bams
        bams_singleton_parts = ch_sample_bams.singleton
        bams_multiple_parts = ch_sample_bams.multiple
        gvcf     = MERGEVCFS.out.merged_vcf.join(MERGEVCFS.out.merged_vcf_tbi)
        gvcf_parts     = HAPLOTYPECALLER.out.gvcf
        versions = ch_versions
}

