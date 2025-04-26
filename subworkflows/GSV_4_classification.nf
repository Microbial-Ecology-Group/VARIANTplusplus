/*───────────────────────────────────────────────────────────────────────────
 *  LOAD THE MODULES NEEDED _AFTER_ DEDUPLICATION
 *───────────────────────────────────────────────────────────────────────────*/
include { index ; bwa_rm_contaminant_fq ; HostRemovalStats          } from '../modules/Alignment/bwa'
include { dlkraken ; runkraken_merged_extract                       } from '../modules/Microbiome/kraken2'
include { MergedPseudoalignFastqFiles ; MergedRunMSweep
         ; MergedParsemSweepResults                                 } from '../modules/Alignment/bwa'


workflow GSV_4_WF {

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
          .set { themisto_in }
          

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

