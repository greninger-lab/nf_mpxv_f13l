//                                                                                 
// Align to reference and filter out samples with fewer than minimum mapped reads threshold
//                                                              
                                                                                   
include { BBMAP_ALIGN_REF       } from '../modules/bbmap_align_ref'                          
                                                                                   
workflow ALIGN_REF_FILTER_SAMPLES {                                                 
    take:                                                                          
    ch_reads    // channel: [ val(meta), path(reads) ]                  
    ch_ref      // channel: [ path(ref) ]
                                                                                   
    main:                                                                          

    BBMAP_ALIGN_REF (
        ch_reads,
        ch_ref
    )
    
    // Filter out samples with fewer than minimum mapped reads threshold           
    BBMAP_ALIGN_REF.out.flagstat                                                   
        .map { meta, flagstat -> [ meta ] + CheckReads.getFlagstatMappedReads(flagstat, params) }
        .set { ch_mapped_reads }                                                   
                                                                                   
    BBMAP_ALIGN_REF.out.bam                                                        
        .join(ch_mapped_reads, by: [0])                                            
        .map { meta, bam, bai, mapped, pass -> if (pass) [ meta, bam, bai ] }   
        .set { ch_bam_out }
    
    ch_mapped_reads
        .map { meta, mapped, pass -> if (!pass) [ meta, mapped ] }
        .set { ch_fail_mapping }

    emit:
    bam     = ch_bam_out                // channel: [ val(meta), path(bam), path(bai) ]
    ref     = BBMAP_ALIGN_REF.out.ref   // channel: [ val(meta), path(ref) ]
    reads   = BBMAP_ALIGN_REF.out.reads // channel: [ val(meta), path(reads) ]
    fail_mapping = ch_fail_mapping      // channel: [ tuple(meta, num_mapped_reads)]
}
