/*───────────────────────────────────────────────────────────────────────────
 *  LOAD THE MODULES NEEDED _AFTER_ DEDUPLICATION
 *───────────────────────────────────────────────────────────────────────────*/
include { dlkraken ; runkraken_merged_extract     ;krakenresults  } from '../modules/Microbiome/kraken2'
include { SeqkitReadCounts }                                          from '../modules/QC/dedup'



/*───────────────────────────────────────────────────────────────────────────
 *  HOST-REFERENCE HANDLING (same logic as before)
 *───────────────────────────────────────────────────────────────────────────*/
workflow GSV_4_WF {

    take:
        merged_reads_ch   

    main:

        
        /*  Build ( sid , mergedFile , unmergedFile )  -------------------------- */

          
          // now each item is exactly: ( sid, Path-to-merged, Path-to-unmerged )
                
        /*──────────────── choose / download Kraken DB ───────────────*/
        def kraken_db_ch
        def default_db = "$baseDir/data/kraken_db/k2_standard_08gb_20250402/"
        if( file(default_db).exists() )
             kraken_db_ch = Channel.value(default_db)
        else if( params.kraken_db )
             kraken_db_ch = Channel.value(params.kraken_db)
        else {
             dlkraken()
             kraken_db_ch = dlkraken.out
        }

        /*──────────────── Kraken 2  ─────────────────────────────────*/
        runkraken_merged_extract( merged_reads_ch , kraken_db_ch )

        
        /* ---------- run krakenresults once all reports are ready ---------- */
        def kraken_reports_list = runkraken_merged_extract.out.kraken_report_merged
                                    .mix( runkraken_merged_extract.out.kraken_report_unmerged )
                                    .collect()                       // Java List (one per sample)
        
        krakenresults( kraken_reports_list )   // ← plain value, no Channel.value()
        
        /* ---------- gather extracted reads -------------------------------- */
        runkraken_merged_extract.out.extracted_merged
                    .join( runkraken_merged_extract.out.extracted_unmerged )
                    .set { extracted_reads_ch }     

        runkraken_merged_extract.out.extracted_merged
            .mix( runkraken_merged_extract.out.extracted_unmerged )
            .set { only_extracted_reads_ch }          // ← promote to workflow scope
          
          
        /* ---------- one-shot SeqKit on all extracted FASTQs --------------- */
        seqkit_fastq_list = only_extracted_reads_ch.map{ sid,f -> f }.collect()
        SeqkitReadCounts( seqkit_fastq_list , "Kraken_extracted" )


            
        emit:
            extracted_reads_ch
            
            
}

