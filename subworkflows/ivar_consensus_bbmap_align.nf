//                                                                                 
// Consensus Assembly using BBMAP for alignment and iVar consensus for consensus calling  
//                                                              
                                                                                   
include { IVAR_CONSENSUS        } from '../modules/ivar_consensus'                          
include { BBMAP_ALIGN_CONSENSUS } from '../modules/bbmap_align_consensus' 
                                                                                   
workflow IVAR_CONSENSUS_BBMAP_ALIGN {                                                 
    take:                                                                          
    ch_bam      // channel: [ val(meta), path(bam), path(bai) ]                             
    ch_ref      // channel: [ val(meta), path(ref) ] 
    ch_reads    // channel: [ val(meta), path(reads) ]                  
                                                                                   
    main:                                                                          

    IVAR_CONSENSUS (
        ch_bam,
        ch_ref
    )

    ch_reads
        .join(IVAR_CONSENSUS.out.consensus)
        .multiMap { meta, reads, consensus -> 
            reads:      [ meta, reads ]
            consensus:  [ meta, consensus ]
        }
        .set { ch_bbmap_align_input }

    BBMAP_ALIGN_CONSENSUS (
        ch_bbmap_align_input.reads,
        ch_bbmap_align_input.consensus
    )
    
    emit:
    bam         = BBMAP_ALIGN_CONSENSUS.out.bam          // channel: [ val(meta), path(bam), path(bai) ]
    consensus   = BBMAP_ALIGN_CONSENSUS.out.consensus    // channel: [ val(meta), path(consensus) ]
    reads       = BBMAP_ALIGN_CONSENSUS.out.reads        // channel: [ val(meta), path(reads) ]
}
