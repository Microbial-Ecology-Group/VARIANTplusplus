// Load modules
include { DeduplicateReadsBBMap } from '../modules/QC/dedup.nf' 

// WC trimming
workflow FASTQ_DEDUP_BBMAP_WF {
    take: 
        read_pairs_ch

    main:

        DeduplicateReadsBBMap(read_pairs_ch)
        
}

