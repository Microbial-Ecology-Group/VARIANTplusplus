/*───────────────────────────────────────────────────────────────────────────
 *  LOAD THE MODULES NEEDED _AFTER_ DEDUPLICATION
 *───────────────────────────────────────────────────────────────────────────*/
include { index ; bwa_rm_contaminant_fq ; HostRemovalStats          } from '../modules/Alignment/bwa'
include { dlkraken ; runkraken_merged_extract                       } from '../modules/Microbiome/kraken2'
include { MergedPseudoalignFastqFiles ; MergedRunMSweep
         ; MergedParsemSweepResults                                 } from '../modules/Alignment/bwa'



/*───────────────────────────────────────────────────────────────────────────
 *  HOST-REFERENCE HANDLING (same logic as before)
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
                 .fromPath( params.host_index )
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

  
  
        //to_host_rm_ch.view()
        /* build / load BWA index exactly as before */


        bwa_rm_contaminant_fq( reference_index_ch , to_host_rm_ch )
        
        /* ───── remove host reads ─────────────────────────────────────*/

        HostRemovalStats( bwa_rm_contaminant_fq.out.host_rm_stats )

 }

