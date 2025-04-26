/*───────────────────────────────────────────────────────────────────────────
 *  LOAD THE MODULES NEEDED _AFTER_ DEDUPLICATION
 *───────────────────────────────────────────────────────────────────────────*/
include { index ; bwa_rm_contaminant_fq ; HostRemovalStats          } from '../modules/Alignment/bwa'
include { dlkraken ; runkraken_merged_extract     ;krakenresults  } from '../modules/Microbiome/kraken2'
include { MergedPseudoalignFastqFiles ; MergedRunMSweep
         ; MergedParsemSweepResults                                 } from '../modules/Alignment/bwa'



/*───────────────────────────────────────────────────────────────────────────
 *  HOST-REFERENCE HANDLING (same logic as before)
 *───────────────────────────────────────────────────────────────────────────*/
workflow GSV_3_WF {

    take:
        merged_reads_ch   
        hostfasta

    main:

        
        /*  Build ( sid , mergedFile , unmergedFile )  -------------------------- */
        Channel
          .fromFilePairs( params.merged_reads, glob: true )
          .ifEmpty { error "No FASTQ files match: ${params.merged_reads}" }
          .map { sample_id, files ->
            //
            // files will be e.g.
            //   [ Path(…/S1_test_merged.dedup.fastq.gz),
            //     Path(…/S1_test_unmerged.dedup.fastq.gz) ]
            //
            def merged   = files.find { it.name.contains('_merged')   }
            def unmerged = files.find { it.name.contains('_unmerged') }
            assert merged && unmerged : "Sample $sample_id missing one of merged/unmerged"
            tuple( sample_id, merged, unmerged )
          }
          .set { kraken_input_ch }
          
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
        runkraken_merged_extract( kraken_input_ch , kraken_db_ch )

        
        krakenresults(
          runkraken_merged_extract.out.kraken_report_merged
            .mix( runkraken_merged_extract.out.kraken_report_unmerged )
            .collect()
        )  

 
        /*──────────────── Themisto & mSWEEP  ───────────────────────*/
        themisto_in = runkraken_merged_extract.out.extracted_merged.join(
                      runkraken_merged_extract.out.extracted_unmerged)\
                      .map { sid, m, u -> tuple(sid, m, u) }

        MergedPseudoalignFastqFiles( themisto_in , params.themisto_index )
        MergedRunMSweep(
            MergedPseudoalignFastqFiles.out.pseudoalignedFastqFiles ,
            params.clustering_file
        )

        /*──────── final parser once all mSWEEP files exist ──────────*/
        MergedParsemSweepResults(
            MergedRunMSweep.out.msweep_merged.mix(
            MergedRunMSweep.out.msweep_unmerged).collect() )
}

