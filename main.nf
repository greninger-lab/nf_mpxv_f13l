#!/usr/bin/env nextflow

// if INPUT not set
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
if (params.trim_primers) { 
    if (params.bed_file) { ch_bed_file = file(params.bed_file) } else { exit 1, 'trim_primers flag is set, but no bed file is provided!' }
}

// Config file
ch_summary_dummy_file = file("$baseDir/assets/summary_dummy_file.tsv", checkIfExists: true)

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
include { IVAR_TRIM         } from './modules/ivar_trim'
include { SAMTOOLS_SORT     } from './modules/samtools_sort'
include { SAMTOOLS_INDEX    } from './modules/samtools_index'
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

    FASTQ_TRIM_FASTP_FASTQC.out.trim_reads_fail_min 
        .collect()
        .set { ch_trim_fail_min_summary }

    if(params.sample) {
        SEQTK_SAMPLE (
            FASTQ_TRIM_FASTP_FASTQC.out.trim_reads_pass_min,
            params.sample
        )
        ch_align_reads = SEQTK_SAMPLE.out.reads
    } else {
        ch_align_reads = FASTQ_TRIM_FASTP_FASTQC.out.trim_reads_pass_min
    }
    
    BBMAP_ALIGN_REF (
        ch_align_reads,
        params.ref
    )

    if (params.trim_primers) {
        IVAR_TRIM (
            BBMAP_ALIGN_REF.out.bam,
            ch_bed_file
        )

        SAMTOOLS_SORT (
            IVAR_TRIM.out.bam,
            IVAR_TRIM.out.bam.map { [it[0], [params.ref]]}
        )

        SAMTOOLS_INDEX (
            SAMTOOLS_SORT.out.bam
        )

        SAMTOOLS_SORT.out.bam.join(SAMTOOLS_INDEX.out.bai).set { ch_bam }
        
    } else {
        ch_bam = BBMAP_ALIGN_REF.out.bam
    }

    BBMAP_ALIGN_REF.out.flagstat                                                   
        .map { meta, flagstat -> [ meta ] + CheckReads.getFlagstatMappedReads(flagstat, params) }
        .set { ch_mapped_reads }

    // Filter out samples with fewer than the minimum mapped reads threshold for variants and consensus calling
    ch_mapped_reads
        .map { meta, mapped, pass -> if (!pass) [ meta, mapped ] }
        .join(FASTQ_TRIM_FASTP_FASTQC.out.reads_num, by: [0])
        .map { 
            meta, mapped, num_reads, num_trimmed_reads -> 
            [ "$meta.id\t$num_reads\t$num_trimmed_reads\t0\t$mapped\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0" ]
        }
        .collect()
        .set { ch_align_ref_fail_summary }

    ch_mapped_reads
        .map { meta, mapped, pass -> if (pass) [ meta ] }
        .join(ch_align_reads, by: [0])
        .join(ch_bam, by: [0])
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
        .concat( ch_trim_fail_min_summary )
        .map { tsvdata -> CheckReads.tsvFromList(tsvdata) }
        .collectFile(storeDir: "${params.output}/summary", name:"fail_summary.tsv", keepHeader: true, sort: false)
        .set { ch_fail_summary }

    SUMMARY.out.summary
        .collectFile(storeDir: "${params.output}/summary", name:"pass_summary.tsv", keepHeader: true, sort: false)
        .set { ch_pass_summary }

    SUMMARY_CLEANUP (
        ch_fail_summary.ifEmpty(ch_summary_dummy_file),
        ch_pass_summary.ifEmpty(ch_summary_dummy_file)
    )

}
