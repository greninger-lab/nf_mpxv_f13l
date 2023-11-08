//                                                                                 
// 3 iterations of consensus assembly using BBMAP for alignment and iVar consensus for consensus calling
//                                                              

include { IVAR_CONSENSUS_BBMAP_ALIGN as IVAR_CONSENSUS_BBMAP_ALIGN_1        } from './ivar_consensus_bbmap_align'                          
include { IVAR_CONSENSUS_BBMAP_ALIGN as IVAR_CONSENSUS_BBMAP_ALIGN_2        } from './ivar_consensus_bbmap_align'                          
include { IVAR_CONSENSUS_BBMAP_ALIGN as IVAR_CONSENSUS_BBMAP_ALIGN_FINAL    } from './ivar_consensus_bbmap_align'                          
                                                                                   
workflow CONSENSUS_ASSEMBLY {                                                 
    take:                                                                          
    ch_bam      // channel: [ val(meta), path(bam), path(bai) ]                  
    ch_ref      // channel: [ val(meta), path(ref) ]
    ch_reads    // channel: [ val(meta), path(reads) ]
                                                                                
    main:                                                                          

    IVAR_CONSENSUS_BBMAP_ALIGN_1 (
        ch_bam,
        ch_ref,
        ch_reads
    )

    IVAR_CONSENSUS_BBMAP_ALIGN_2 (
        IVAR_CONSENSUS_BBMAP_ALIGN_1.out.bam,
        IVAR_CONSENSUS_BBMAP_ALIGN_1.out.consensus,
        IVAR_CONSENSUS_BBMAP_ALIGN_1.out.reads
    )

    IVAR_CONSENSUS_BBMAP_ALIGN_FINAL (
        IVAR_CONSENSUS_BBMAP_ALIGN_2.out.bam,
        IVAR_CONSENSUS_BBMAP_ALIGN_2.out.consensus,
        IVAR_CONSENSUS_BBMAP_ALIGN_2.out.reads
    )

    emit:
    consensus   = IVAR_CONSENSUS_BBMAP_ALIGN_FINAL.out.consensus // channel: [ val(meta), path(consensus) ]
    bam         = IVAR_CONSENSUS_BBMAP_ALIGN_FINAL.out.bam       // channel: [ val(meta), path(bam), path(bai) ]
} 
