/*───────────────────────────────────────────────────────────────────────────
 *  LOAD THE MODULES NEEDED _AFTER_ DEDUPLICATION
 *───────────────────────────────────────────────────────────────────────────*/
include { dlkraken ; runkraken_merged_extract ; krakenresults } from '../modules/Microbiome/kraken2'
include { SeqkitReadCounts }                                     from '../modules/QC/dedup'

/*───────────────────────────────────────────────────────────────────────────
 *  Helper function to resolve Kraken database path
 *───────────────────────────────────────────────────────────────────────────*/
def resolveKrakenDb() {
    def default_db_path = "$baseDir/data/kraken_db/k2_standard_08gb_20250402/"
    
    if (params.kraken_db && file(params.kraken_db).exists()) {
        log.info "[Kraken] Using user-specified database: ${params.kraken_db}"
        return [params.kraken_db, false]
    } else if (file(default_db_path).exists()) {
        log.info "[Kraken] Using default database: ${default_db_path}"
        return [default_db_path, false]
    } else {
        log.info "[Kraken] No database found - will download"
        return [null, true]
    }
}

/*───────────────────────────────────────────────────────────────────────────
 *  Subworkflow to get or download Kraken database
 *───────────────────────────────────────────────────────────────────────────*/
workflow KRAKEN_DB {
    main:
        def (db_path, needs_download) = resolveKrakenDb()
        
        if (needs_download) {
            dlkraken()
            kraken_db_ch = dlkraken.out.kraken_db
        } else {
            kraken_db_ch = Channel.value(file(db_path, type: 'dir'))
        }
    
    emit:
        db = kraken_db_ch
}

/*───────────────────────────────────────────────────────────────────────────
 *  GSV_4_WF WORKFLOW
 *───────────────────────────────────────────────────────────────────────────*/
workflow GSV_4_WF {
    take:
        merged_reads_ch   
    
    main:
        /*──────────────── choose / download Kraken DB ───────────────*/
        KRAKEN_DB()
        
        /*──────────────── Kraken 2  ─────────────────────────────────*/
        runkraken_merged_extract(merged_reads_ch, KRAKEN_DB.out.db)
        
        /* ---------- run krakenresults once all reports are ready ---------- */
        kraken_reports_list = runkraken_merged_extract.out.kraken_report_merged
            .mix(runkraken_merged_extract.out.kraken_report_unmerged)
            .collect()
        
        krakenresults(kraken_reports_list)
        
        /* ---------- gather extracted reads -------------------------------- */
        runkraken_merged_extract.out.extracted_merged
            .join(runkraken_merged_extract.out.extracted_unmerged)
            .set { extracted_reads_ch }
        
        runkraken_merged_extract.out.extracted_merged
            .mix(runkraken_merged_extract.out.extracted_unmerged)
            .set { only_extracted_reads_ch }
        
        /* ---------- one-shot SeqKit on all extracted FASTQs --------------- */
        seqkit_fastq_list = only_extracted_reads_ch.map { sid, f -> f }.collect()
        SeqkitReadCounts(seqkit_fastq_list, "Kraken_extracted")
            
    emit:
        extracted_reads_ch
}