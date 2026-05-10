/*───────────────────────────────────────────────────────────────────────────
 *  FINAL GSV_5_MGEMS_WF SUBWORKFLOW - mGEMS enabled by default for short reads
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
        
        // Step 3: mGEMS read binning (enabled by default)
        if (params.run_mgems != false) {  // Run unless explicitly disabled
            
            // Create mGEMS input by combining all required channels
            // We'll use a simple approach that filters for samples with mSWEEP output
            
            // Get samples that have mSWEEP results
            successful_samples = MergedRunMSweep.out.msweep_merged
                .concat(MergedRunMSweep.out.msweep_unmerged)
                .map { file -> 
                    def sample_id = file.baseName.tokenize('.')[0]
                    return sample_id
                }
                .unique()
                
            // Filter to only run mGEMS on samples with successful mSWEEP
            mgems_input = merged_reads_ch
                .join(MergedPseudoalignFastqFiles.out.pseudoalignedFastqFiles, by: 0)
                .cross(successful_samples) { it[0] }
                .map { paired_data, sample_id ->
                    // paired_data = [sample_id, merged_fq, unmerged_fq, pseudo_merged, pseudo_unmerged]
                    return paired_data[0] // Return the original tuple
                }
                
            // Add mSWEEP results to the channel
            mgems_with_results = mgems_input
                .map { sample_id, merged_fq, unmerged_fq, pseudo_merged, pseudo_unmerged ->
                    // Create base tuple for mGEMS
                    [sample_id, merged_fq, unmerged_fq, pseudo_merged, pseudo_unmerged]
                }
                // Add placeholder paths for mSWEEP results - mGEMS will find them by convention
                .map { sample_id, merged_fq, unmerged_fq, pseudo_merged, pseudo_unmerged ->
                    [sample_id, merged_fq, unmerged_fq, pseudo_merged, pseudo_unmerged, 
                     file("${sample_id}.merged.msweep_abundances.txt"), 
                     file("${sample_id}.merged.msweep_probs.csv"),
                     file("${sample_id}.unmerged.msweep_abundances.txt"), 
                     file("${sample_id}.unmerged.msweep_probs.csv")]
                }
            
            MergedRunMGEMS(
                mgems_with_results,
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
        binned_reads = params.run_mgems != false ? 
            MergedRunMGEMS.out.binned_merged_reads.mix(MergedRunMGEMS.out.binned_unmerged_reads) : 
            Channel.empty()
        assignment_tables = params.run_mgems != false ? 
            MergedRunMGEMS.out.assignment_table_merged.mix(MergedRunMGEMS.out.assignment_table_unmerged) : 
            Channel.empty()
        mgems_summary = params.run_mgems != false ? 
            MergedRunMGEMS.out.mgems_summary : 
            Channel.empty()
}
