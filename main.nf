#!/usr/bin/env nextflow

// if INPUT not set
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

//
// SUBWORKFLOWS
//
include { INPUT_CHECK               } from './subworkflows/input_check'
include { FASTQ_TRIM_FASTP_FASTQC   } from './subworkflows/fastq_trim_fastp_fastqc'
include { CONSENSUS_ASSEMBLY        } from './subworkflows/consensus_assembly'

//
// MODULES
//
include { SEQTK_SAMPLE      } from './modules/seqtk_sample'
include { BBMAP_ALIGN_REF   } from './modules/bbmap_align_ref'
include { F13L_VARIANTS     } from './modules/f13l_variants'
include { SUMMARY           } from './modules/summary'
include { SUMMARY_CLEANUP   } from './modules/summary_cleanup'

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

    // Filter out samples with zero read after trimming
    ch_postQC_reads = FASTQ_TRIM_FASTP_FASTQC.out.reads
    if (!params.skip_fastp) {
        ch_postQC_reads
            .join(FASTQ_TRIM_FASTP_FASTQC.out.trim_json)
            .map {
                meta, reads, json ->
                    num_reads = CheckReads.getFastpReadsBeforeFiltering(json)
                    num_trimmed_reads = CheckReads.getFastpReadsAfterFiltering(json)
                    pass = num_trimmed_reads > 0
                    [ meta, reads, num_reads, num_trimmed_reads, pass ]
            }
            .set { ch_pass_fail_reads }

        ch_pass_fail_reads
            .map { meta, reads, num_reads, num_trimmed_reads, pass -> if (pass) [ meta, reads ] }
            .set { ch_postQC_reads }

        ch_pass_fail_reads
            .map {
                meta, reads, num_reads, num_trimmed_reads, pass ->
                if (!pass) [ "$meta.id\t$num_reads\t$num_trimmed_reads\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0" ]
            }
            .collect()
            .set { ch_fastp_fail_summary }
    }

    ch_align_reads = ch_postQC_reads
    if(params.sample) {
        SEQTK_SAMPLE (
            ch_postQC_reads,
            params.sample
        )
        ch_align_reads = SEQTK_SAMPLE.out.reads
    }
    
    BBMAP_ALIGN_REF (
        ch_align_reads,
        params.ref
    )

    BBMAP_ALIGN_REF.out.flagstat                                                   
        .map { meta, flagstat -> [ meta ] + CheckReads.getFlagstatMappedReads(flagstat, params) }
        .set { ch_mapped_reads }

    // Filter out samples with fewer than the minimum mapped reads threshold for variants and consensus calling
    ch_mapped_reads
        .map { meta, mapped, pass -> if (!pass) [ meta, mapped ] }
        .join(ch_pass_fail_reads, by: [0])
        .map { 
            meta, mapped, reads, num_reads, num_trimmed_reads, pass -> 
            [ "$meta.id\t$num_reads\t$num_trimmed_reads\t0\t$mapped\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0" ]
        }
        .collect()
        .set { ch_align_ref_fail_summary }

    ch_mapped_reads
        .map { meta, mapped, pass -> if (pass) [ meta ] }
        .join(ch_align_reads, by: [0])
        .join(BBMAP_ALIGN_REF.out.bam, by: [0])
        .multiMap { meta, reads, bam, bai ->
            reads:  [ meta, reads ]
            ref:    [ meta, params.ref ]
            bam:    [ meta, bam, bai ]
        }.set { ch_variants_consensus }

    F13L_VARIANTS (
        ch_variants_consensus.bam,
        params.ref,
        params.ref_index,
        params.gff,
        params.save_mpileup
    )

    CONSENSUS_ASSEMBLY (
        ch_variants_consensus.bam,
        ch_variants_consensus.ref,
        ch_variants_consensus.reads
    )

    CONSENSUS_ASSEMBLY.out.consensus
        .join(ch_variants_consensus.bam)
        .join(FASTQ_TRIM_FASTP_FASTQC.out.trim_log)
        .set { ch_summary_in }
    
    SUMMARY (
        ch_summary_in
    )

    ch_align_ref_fail_summary
        .concat( ch_fastp_fail_summary )
        .map { tsvdata -> CheckReads.tsvFromList(tsvdata) }
        .collectFile(storeDir: "${params.output}/summary", name:"fail_summary.tsv", keepHeader: true, sort: false)
        .set { ch_fail_summary }

    SUMMARY.out.summary
        .collectFile(storeDir: "${params.output}/summary", name:"pass_summary.tsv", keepHeader: true, sort: false)
        .set { ch_pass_summary }

    SUMMARY_CLEANUP (
        ch_fail_summary,
        ch_pass_summary
    )

}
