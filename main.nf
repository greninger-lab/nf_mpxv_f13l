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

    if(params.sample) {
        SEQTK_SAMPLE (
            FASTQ_TRIM_FASTP_FASTQC.out.reads,
            params.sample
        )
        ch_align_input = SEQTK_SAMPLE.out.reads
    } else {
        ch_align_input = FASTQ_TRIM_FASTP_FASTQC.out.reads
    } 
    
    BBMAP_ALIGN_REF (
        ch_align_input,
        params.ref
    )
    
    F13L_VARIANTS (
        BBMAP_ALIGN_REF.out.bam,
        params.ref,
        params.ref_index,
        params.gff,
        params.save_mpileup
    )

    CONSENSUS_ASSEMBLY (
        BBMAP_ALIGN_REF.out.bam,
        BBMAP_ALIGN_REF.out.ref,
        BBMAP_ALIGN_REF.out.reads
    )

    FASTQ_TRIM_FASTP_FASTQC.out.trim_log
        .join(BBMAP_ALIGN_REF.out.bam)
        .join(CONSENSUS_ASSEMBLY.out.consensus)
        .set { ch_summary_in }

    SUMMARY (
        ch_summary_in
    )

    SUMMARY.out.summary
        .collectFile(storeDir: "${params.output}", name:"${params.run_name}_summary.tsv", keepHeader: true, sort: true)
}
