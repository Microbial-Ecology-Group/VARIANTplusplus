/*───────────────────────────────────────────────────────────────────────────
 *  UPDATED GSV_5_WF SUBWORKFLOW WITH mSWEEP OPTIMIZATIONS AND mGEMS
 *───────────────────────────────────────────────────────────────────────────*/

/*───────────────────────────────────────────────────────────────────────────
 *  LOAD THE UPDATED MODULES 
 *───────────────────────────────────────────────────────────────────────────*/
include { MergedPseudoalignFastqFiles ; MergedRunMSweep
         ; MergedParsemSweepResults ; MergedRunMGEMS                     } from '../modules/Alignment/msweep'

/*───────────────────────────────────────────────────────────────────────────
 *  UPDATED GSV_5_WF WORKFLOW
 *───────────────────────────────────────────────────────────────────────────*/
workflow GSV_5_MGEMS_WF {
    take:
        merged_reads_ch   // ( sid, Path-to-merged, Path-to-unmerged )
    
    main:
        // Step 1: Pseudoalignment with Themisto (now with optional threshold)
        MergedPseudoalignFastqFiles( merged_reads_ch , params.themisto_index )
        
        // Step 2: mSWEEP abundance estimation (now with optimized parameters and --write-probs)
        MergedRunMSweep(
            MergedPseudoalignFastqFiles.out.pseudoalignedFastqFiles , 
            params.clustering_file 
        )
        
        // Step 3: Optional mGEMS read binning (only if enabled)
        if (params.run_mgems) {
            // Combine all inputs needed for mGEMS
            mgems_input_ch = merged_reads_ch
                .join(MergedPseudoalignFastqFiles.out.pseudoalignedFastqFiles)
                .join(MergedRunMSweep.out.msweep_merged, remainder: true)
                .join(MergedRunMSweep.out.msweep_merged_probs, remainder: true)
                .join(MergedRunMSweep.out.msweep_unmerged, remainder: true) 
                .join(MergedRunMSweep.out.msweep_unmerged_probs, remainder: true)
            
            MergedRunMGEMS(
                mgems_input_ch,
                params.clustering_file,
                params.themisto_index
            )
        }
        
        // Step 4: Parse mSWEEP results (existing functionality)
        MergedParsemSweepResults( 
            MergedRunMSweep.out.msweep_merged
                .mix(MergedRunMSweep.out.msweep_unmerged)
                .collect() 
        )
    
    emit:
        // Original outputs
        msweep_summary = MergedParsemSweepResults.out.msweep_summary
        msweep_matrix = MergedParsemSweepResults.out.msweep_matrix
        
        // New outputs (only if mGEMS is enabled)
        binned_reads = params.run_mgems ? MergedRunMGEMS.out.binned_merged_reads.mix(MergedRunMGEMS.out.binned_unmerged_reads) : Channel.empty()
        assignment_tables = params.run_mgems ? MergedRunMGEMS.out.assignment_table_merged.mix(MergedRunMGEMS.out.assignment_table_unmerged) : Channel.empty()
        mgems_summary = params.run_mgems ? MergedRunMGEMS.out.mgems_summary : Channel.empty()
}
