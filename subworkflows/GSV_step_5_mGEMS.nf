/*───────────────────────────────────────────────────────────────────────────
 *  FIXED GSV_5_MGEMS_WF - Proper file passing for mGEMS
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
        
        // Step 3: mGEMS read binning (FIXED: proper file passing)
        if (params.run_mgems == true) {
            
            // Properly combine inputs with actual mSWEEP outputs
            mgems_input_ch = merged_reads_ch
                .join(MergedPseudoalignFastqFiles.out.pseudoalignedFastqFiles, by: 0)
                .join(MergedRunMSweep.out.msweep_merged, by: 0)
                .join(MergedRunMSweep.out.msweep_unmerged, by: 0)
                .map { sample_id, merged_fq, unmerged_fq, pseudo_merged, pseudo_unmerged, 
                       msweep_merged_tuple, msweep_unmerged_tuple ->
                    
                    // Extract actual abundance and probability files from mSWEEP output
                    def (merged_abundances, merged_probs) = msweep_merged_tuple
                    def (unmerged_abundances, unmerged_probs) = msweep_unmerged_tuple
                    
                    // Create proper tuple for mGEMS
                    [sample_id, 
                     merged_fq, unmerged_fq,                    // Original FASTQ files
                     pseudo_merged, pseudo_unmerged,            // Pseudoaligned files  
                     merged_abundances, merged_probs,           // ACTUAL mSWEEP merged outputs
                     unmerged_abundances, unmerged_probs]       // ACTUAL mSWEEP unmerged outputs
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
