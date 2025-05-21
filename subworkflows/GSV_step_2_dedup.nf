/*───────────────────────────────────────────────────────────────────────────
 *  LOAD THE MODULES NEEDED _AFTER_ DEDUPLICATION
 *───────────────────────────────────────────────────────────────────────────*/
include { SeqkitReadCounts  ; MergedDeduplicateReadsSeqkit          } from '../modules/QC/dedup'


/*───────────────────────────────────────────────────────────────────────────
 *  HOST-REFERENCE HANDLING + READ-COUNT STATS
 *───────────────────────────────────────────────────────────────────────────*/
workflow GSV_2_WF {

    take:
        merged_reads_ch

    main:
        /* Deduplicate reads ---- */ 
        MergedDeduplicateReadsSeqkit( merged_reads_ch )
        
        /* Summarize counts ---- */ 
        
        /* grab only the file paths (not the sample IDs) */
        def dedup_merged_paths   = MergedDeduplicateReadsSeqkit.out.dedup_merged  .map { sid, f -> f }
        def dedup_unmerged_paths = MergedDeduplicateReadsSeqkit.out.dedup_unmerged.map { sid, f -> f }

        /* 3 ─ make one channel that waits for ALL paths ----------------- */
        def dedup_all_fastqs_ch = dedup_merged_paths
                              .mix( dedup_unmerged_paths )   // merge both streams
                              .collect()               // <- waits for channel to close

        /* 4 ─ final read-count summary ---------------------------------- */
        SeqkitReadCounts(dedup_all_fastqs_ch, "Deduped")

    emit:
        dedup_merged   = MergedDeduplicateReadsSeqkit.out.dedup_merged
        dedup_unmerged = MergedDeduplicateReadsSeqkit.out.dedup_unmerged


}


