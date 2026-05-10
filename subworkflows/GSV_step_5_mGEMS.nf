/*───────────────────────────────────────────────────────────────────────────
 *  SIMPLIFIED GSV_5_MGEMS_WF - Clean defaults with working mGEMS
 *───────────────────────────────────────────────────────────────────────────*/

include { MergedPseudoalignFastqFiles ; MergedRunMSweep
         ; MergedParsemSweepResults ; MergedRunMGEMS                     } from '../modules/Alignment/msweep'

workflow GSV_5_MGEMS_WF {
    take:
        merged_reads_ch   // ( sid, Path-to-merged, Path-to-unmerged )
    
    main:
        // Step 1: Pseudoalignment with Themisto
        MergedPseudoalignFastqFiles( merged_reads_ch , params.themisto_index )
        
        // Step 2: mSWEEP abundance estimation (with --write-probs for mGEMS)
        MergedRunMSweep(
            MergedPseudoalignFastqFiles.out.pseudoalignedFastqFiles , 
            params.clustering_file 
        )
        
        // Step 3: mGEMS read binning (simplified approach)
        if (params.run_mgems == true) {
            
            // Simple approach: combine original reads with pseudoaligned files
            // and create placeholder mSWEEP result paths for mGEMS to find
            mgems_input_ch = merged_reads_ch
                .join(MergedPseudoalignFastqFiles.out.pseudoalignedFastqFiles, by: 0)
                .map { sample_id, merged_fq, unmerged_fq, pseudo_merged, pseudo_unmerged ->
                    // Create tuple with all inputs mGEMS needs
                    [sample_id, 
                     merged_fq, unmerged_fq,                    // Original FASTQ files
                     pseudo_merged, pseudo_unmerged,            // Pseudoaligned files
                     file("placeholder_abundances_merged"),     // Placeholders - mGEMS will find actual files
                     file("placeholder_probs_merged"), 
                     file("placeholder_abundances_unmerged"), 
                     file("placeholder_probs_unmerged")]
                }
            
            MergedRunMGEMS(
                mgems_input_ch,
                params.clustering_file,
                params.themisto_index
            )
        }
        
        // Step 4: Parse mSWEEP results
        MergedParsemSweepResults( 
            MergedRunMSweep.out.msweep_merged
                .mix(MergedRunMSweep.out.msweep_unmerged)
                .collect() 
        )
    
    emit:
        // Standard outputs
        msweep_summary = MergedParsemSweepResults.out.msweep_summary
        msweep_matrix = MergedParsemSweepResults.out.msweep_matrix
        
        // mGEMS outputs (if enabled)
        binned_reads = params.run_mgems == true ? 
            MergedRunMGEMS.out.binned_merged_reads.mix(MergedRunMGEMS.out.binned_unmerged_reads) : 
            Channel.empty()
        assignment_tables = params.run_mgems == true ? 
            MergedRunMGEMS.out.assignment_table_merged.mix(MergedRunMGEMS.out.assignment_table_unmerged) : 
            Channel.empty()
        mgems_summary = params.run_mgems == true ? 
            MergedRunMGEMS.out.mgems_summary : 
            Channel.empty()
}
