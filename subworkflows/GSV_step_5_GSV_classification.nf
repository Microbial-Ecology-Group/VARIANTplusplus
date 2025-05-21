/*───────────────────────────────────────────────────────────────────────────
 *  LOAD THE MODULES NEEDED _AFTER_ DEDUPLICATION
 *───────────────────────────────────────────────────────────────────────────*/
include { MergedPseudoalignFastqFiles ; MergedRunMSweep
         ; MergedParsemSweepResults                                 } from '../modules/Alignment/bwa'



/*───────────────────────────────────────────────────────────────────────────
 *  HOST-REFERENCE HANDLING (same logic as before)
 *───────────────────────────────────────────────────────────────────────────*/
workflow GSV_5_WF {

    take:
        merged_reads_ch   

    main:


          // now each item is exactly: ( sid, Path-to-merged, Path-to-unmerged )
        MergedPseudoalignFastqFiles( merged_reads_ch , params.themisto_index )
        MergedRunMSweep(MergedPseudoalignFastqFiles.out.pseudoalignedFastqFiles , params.clustering_file )

        /*──────── final parser once all mSWEEP files exist ──────────*/
        MergedParsemSweepResults( MergedRunMSweep.out.msweep_merged.mix( MergedRunMSweep.out.msweep_unmerged).collect() )
            
}

