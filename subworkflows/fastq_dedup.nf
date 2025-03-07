// Load modules
include { DeduplicateReads } from '../modules/QC/dedup.nf' 

// WC trimming
workflow FASTQ_DEDUP_WF {
    take: 
        read_pairs_ch

    main:

        DeduplicateReads(read_pairs_ch)
        
}

