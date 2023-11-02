#!/usr/bin/env nextflow

// if INPUT not set
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

//
// SUBWORKFLOWS
//
include { INPUT_CHECK               } from './subworkflows/input_check'
include { FASTQ_TRIM_FASTP_FASTQC   } from './subworkflows/fastq_trim_fastp_fastqc'
include { ALIGN_REF_FILTER_SAMPLES  } from './subworkflows/align_ref_filter_samples'
include { CONSENSUS_ASSEMBLY        } from './subworkflows/consensus_assembly'

//
// MODULES
//
include { SEQTK_SAMPLE      } from './modules/seqtk_sample'
include { F13L_VARIANTS     } from './modules/f13l_variants'
include { SUMMARY           } from './modules/summary'

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
/*                                                    */
/*                 RUN THE WORKFLOW                   */
/*                                                    */
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

workflow {
    
    INPUT_CHECK (
        ch_input
    )

    FASTQ_TRIM_FASTP_FASTQC (
        INPUT_CHECK.out.reads,
        params.adapter_fasta,
        params.save_trimmed_fail,
        params.save_merged,
        params.skip_fastp,
        params.skip_fastqc
    )

    // Filter out samples with fewer than the minimum trimmed reads threshold
    ch_postQC_reads = FASTQ_TRIM_FASTP_FASTQC.out.reads
    if (!params.skip_fastp) {
        ch_postQC_reads
            .join(FASTQ_TRIM_FASTP_FASTQC.out.trim_json)
            .map {
                meta, reads, json ->
                    num_reads = CheckReads.getFastpReadsBeforeFiltering(json)
                    num_trimmed_reads = CheckReads.getFastpReadsAfterFiltering(json)
                    pass = num_trimmed_reads > 0
                    [ meta, reads, json, num_reads, num_trimmed_reads, pass ]
            }
            .set { ch_pass_fail_reads }

        ch_pass_fail_reads
            .map { meta, reads, json, num_reads, num_trimmed_reads, pass -> if (pass) [ meta, reads ] }
            .set { ch_postQC_reads }

        ch_pass_fail_reads
            .map {
                meta, reads, json, num_reads, num_trimmed_reads, pass ->
                if (!pass) [ "$meta.id\t$num_reads\t$num_trimmed_reads\tNA" ]
            }
            .collect()
            .map { 
                tsv_data ->
                    def header = ['Sample', 'Reads before trimming', 'Reads after trimming', 'Mapped reads']
                    CheckReads.multiqcTsvFromList(tsv_data, header)
            }
            .set { ch_fastp_fail }
    }

    ch_align_reads = ch_postQC_reads
    if(params.sample) {
        SEQTK_SAMPLE (
            ch_postQC_reads,
            params.sample
        )
        ch_align_reads = SEQTK_SAMPLE.out.reads
    }
    
    ALIGN_REF_FILTER_SAMPLES (
        ch_align_reads,
        params.ref
    )

    ALIGN_REF_FILTER_SAMPLES.out.fail_mapping
        .join(ch_pass_fail_reads)
        .map { meta, mapped, reads, json, num_reads, num_trimmed_reads, pass -> [ "$meta.id\t$num_reads\t$num_trimmed_reads\t$mapped" ] }
        .collect()
        .map {
            tsv_data ->
                def header = ['Sample', 'Reads before trimming', 'Reads after trimming', 'Mapped reads']
                CheckReads.multiqcTsvFromList(tsv_data, header)
        }
        .set { ch_align_ref_fail }

        F13L_VARIANTS (
        ALIGN_REF_FILTER_SAMPLES.out.bam,
        params.ref,
        params.ref_index,
        params.gff,
        params.save_mpileup
    )

    CONSENSUS_ASSEMBLY (
        ALIGN_REF_FILTER_SAMPLES.out.bam,
        ALIGN_REF_FILTER_SAMPLES.out.ref,
        ALIGN_REF_FILTER_SAMPLES.out.reads
    )

    ALIGN_REF_FILTER_SAMPLES.out.bam
        .join(CONSENSUS_ASSEMBLY.out.consensus)
        .join(FASTQ_TRIM_FASTP_FASTQC.out.trim_log)
        .set { ch_summary_in }
    
    SUMMARY (
        ch_summary_in
    )

    ch_fastp_fail
        .concat( ch_align_ref_fail )
        .collectFile(storeDir: "${params.output}", name:"fail_summary.tsv", keepHeader: true)

    SUMMARY.out.summary
        .collectFile(storeDir: "${params.output}", name:"${params.run_name}_summary.tsv", keepHeader: true, sort: true)
}
