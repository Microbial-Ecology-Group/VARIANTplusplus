// Load modules
include { MergeReadsFlash } from '../modules/QC/merge'

// WC trimming
workflow FLASH_MERGE_WF {
    take: 
        read_pairs_ch

    main:
        //index( hostindex )
        //bwa_align( index.out, read_pairs_ch )
        MergeReadsFlash(read_pairs_ch)
        
    emit:
        //bwa_align = bwa_align.out
        merged_reads = MergeReadsFlash.out.merged
        unmerged_reads = MergeReadsFlash.out.unmerged
        
}
