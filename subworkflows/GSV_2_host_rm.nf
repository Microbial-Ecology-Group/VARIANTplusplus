/*───────────────────────────────────────────────────────────────────────────
 *  LOAD THE MODULES NEEDED _AFTER_ DEDUPLICATION
 *───────────────────────────────────────────────────────────────────────────*/
include { index ; bwa_rm_contaminant_fq ; HostRemovalStats          } from '../modules/Alignment/bwa'
include { dlkraken ; runkraken_merged_extract                       } from '../modules/Microbiome/kraken2'
include { MergedPseudoalignFastqFiles ; MergedRunMSweep
         ; MergedParsemSweepResults                                 } from '../modules/Alignment/bwa'
include { SeqkitReadCounts                                          } from '../modules/QC/dedup'


/*───────────────────────────────────────────────────────────────────────────
 *  HOST-REFERENCE HANDLING + READ-COUNT STATS
 *───────────────────────────────────────────────────────────────────────────*/
workflow GSV_2_WF {

    take:
        merged_reads_ch
        hostfasta

    main:
        /* index or load the host genome */
        def reference_index_ch
        if( params.host_index ) {
            reference_index_ch = Channel
                 .fromPath( params.host_index , glob: true )
                 .ifEmpty { error "No files match --host_index '${params.host_index}'" }
                 .toList()
                 .map { files -> files.sort() }
        } else {
            index( hostfasta )
            reference_index_ch = index.out
        }

        
        /*  Build ( sid , mergedFile , unmergedFile )  -------------------------- */

        Channel
          .fromPath(params.merged_reads, glob: true)
          .ifEmpty { error "No FASTQs match: ${params.merged_reads}" }
          .map { f ->
    
            // pull "S1_test_merged" and "S1_test" & "merged"
            def m = (f.name =~ /(.+?)_(merged|unmerged)\./)
            if( !m ) error "Failed to parse sample/read_type from ${f.name}"
            def sid   = m[0][1]
            def rtype = m[0][2]
    
            // build the `pair_id` exactly as you had before:
            def pair_id = "${sid}_${rtype}"
    
            // each tuple is ( pair_id , [ file ] )
            tuple( pair_id, [ f ] )
          }
          .set { to_host_rm_ch }


        bwa_rm_contaminant_fq( reference_index_ch, to_host_rm_ch )

        /* ── extra channel: the non-host reads in tuple form ───────────── */
        def nonhost_reads_ch = bwa_rm_contaminant_fq.out.nonhost_reads
        //  cleanedFastqs is whatever tuple the module emits, e.g.
        //  ( sample_id , cleanedMerged.fastq.gz? , cleanedUnmerged.fastq.gz? )

        HostRemovalStats( bwa_rm_contaminant_fq.out.host_rm_stats )

        /* wait-token so the stats step starts only after host removal */
        def bwa_done_ch = bwa_rm_contaminant_fq.out.host_rm_stats.map{ 1 }.first()



    /* ── what this sub-workflow gives to its caller ────────────────────── */
    emit:
        bwa_rm_contaminant_fq.out.nonhost_reads
}


