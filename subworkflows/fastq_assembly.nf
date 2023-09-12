// Load modules
include { run_skesa_reads } from '../modules/Assembly/skesa.nf' 

// WC trimming
workflow FASTQ_SKESA_WF {
    take: 
        read_pairs_ch

    main:

        run_skesa_reads(read_pairs_ch)
        
}

