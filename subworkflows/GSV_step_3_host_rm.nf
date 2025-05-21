/*───────────────────────────────────────────────────────────────────────────
 *  MODULE INCLUDES  (unchanged, no duplicate names here)
 *───────────────────────────────────────────────────────────────────────────*/
include { index ;  Mergedbwa_rm_contaminant_fq ; HostRemovalStats }          from '../modules/Alignment/bwa'
include { dlkraken ; runkraken_merged_extract }                       from '../modules/Microbiome/kraken2'
include { MergedPseudoalignFastqFiles ; MergedRunMSweep
         ; MergedParsemSweepResults }                                 from '../modules/Alignment/bwa'
include { SeqkitReadCounts }                                          from '../modules/QC/dedup'


/*───────────────────────────────────────────────────────────────────────────
 *  HOST-REFERENCE HANDLING + READ-COUNT STATS
 *───────────────────────────────────────────────────────────────────────────*/
workflow GSV_3_WF {

    take:
        merged_reads_ch          // (all *_merged.fastq.gz and *_unmerged.fastq.gz)
        hostfasta

    main:
    /* 1 ─ build / load BWA index -------------------------------------- */
    def reference_index_ch =
        params.host_index
        ? Channel.fromPath( params.host_index , glob:true )
                 .ifEmpty { error "No files match --host_index '${params.host_index}'" }
                 .toList()
                 .map { it.sort() }               // bundle 6 index files
        : { index( hostfasta ); index.out }()     // call in a closure

    /* 2 ─ host-removal -------------------------------------------------- */
    Mergedbwa_rm_contaminant_fq( reference_index_ch , merged_reads_ch )
    
    Mergedbwa_rm_contaminant_fq.out.nonhost_merged
        .join( Mergedbwa_rm_contaminant_fq.out.nonhost_unmerged )
        .set { nonhost_reads_ch } 
    
    /* merged + unmerged non-host FASTQs  → one channel  */
    Mergedbwa_rm_contaminant_fq.out.nonhost_merged
        .mix( Mergedbwa_rm_contaminant_fq.out.nonhost_unmerged )
        .set { only_reads_ch }            // ← make it top-level
    
    
    /* 3 ─ one-shot SeqKit on all non-host FASTQs ----------------------- */
    seqkit_input_ch = only_reads_ch.map{ sid,f -> f }.collect()
    
    SeqkitReadCounts( seqkit_input_ch , "NonHost" )
    

    emit:
        nonhost_reads_ch

}
