/*───────────────────────────────────────────────────────────────────────────
 *  GSV_5_MGEMS_WF  — mSWEEP abundance estimation + mGEMS read binning
 *
 *  Fixes vs previous version
 *  ──────────────────────────
 *  1. map{} closure now names all 9 elements that the join chain produces:
 *       [sid, merged_fq, unmerged_fq,
 *        pseudo_merged, pseudo_unmerged,
 *        merged_abundances, merged_probs,
 *        unmerged_abundances, unmerged_probs]
 *     Previously only 7 names were given, so `msweep_merged_tuple` received
 *     a bare Path and the destructure `def (a, b) = msweep_merged_tuple`
 *     failed at runtime.
 *
 *  2. Conditional channel logic moved out of the emit block.
 *     DSL2 does not support ternary operators on channels inside emit{}.
 *     mGEMS channels are pre-assigned to Channel.empty() and overwritten
 *     only when params.run_mgems is true.
 *
 *  3. params.run_mgems guarded with a null-safe default so the workflow
 *     does not crash if the param is absent from nextflow.config.
 *───────────────────────────────────────────────────────────────────────────*/

include { MergedPseudoalignFastqFiles
        ; MergedRunMSweep
        ; MergedRunMGEMS
        ; MergedParsemSweepResults } from '../modules/Alignment/msweep'

workflow GSV_5_MGEMS_WF {

    take:
        merged_reads_ch   // tuple: ( sid, Path-merged-fq, Path-unmerged-fq )

    main:

        // ── Step 1: Themisto pseudoalignment ─────────────────────────────────
        MergedPseudoalignFastqFiles(
            merged_reads_ch,
            params.themisto_index
        )

        // ── Step 2: mSWEEP abundance estimation ──────────────────────────────
        // The updated MergedRunMSweep outputs:
        //   msweep_merged   → ( sid, merged_abundances.txt,   merged_probs.tsv   )
        //   msweep_unmerged → ( sid, unmerged_abundances.txt, unmerged_probs.tsv )
        MergedRunMSweep(
            MergedPseudoalignFastqFiles.out.pseudoalignedFastqFiles,
            params.clustering_file
        )

        // ── Step 3 (optional): mGEMS read binning ────────────────────────────
        //
        // Build the 9-element input tuple that MergedRunMGEMS expects:
        //   ( sid,
        //     merged_fq, unmerged_fq,                  ← original FASTQs
        //     pseudo_merged, pseudo_unmerged,           ← pseudoaligned FASTQs
        //     merged_abundances, merged_probs,          ← mSWEEP merged outputs
        //     unmerged_abundances, unmerged_probs )     ← mSWEEP unmerged outputs
        //
        // Join chain element counts (after each join):
        //   merged_reads_ch                      → 3  [sid, mfq, ufq]
        //   .join(pseudoalignedFastqFiles, by:0) → 5  [sid, mfq, ufq, pm, pu]
        //   .join(msweep_merged,           by:0) → 7  [..., ma, mp]
        //   .join(msweep_unmerged,         by:0) → 9  [..., ua, up]
        //
        // The map closure must name all 9 elements — anything fewer and
        // Groovy silently binds trailing elements to the last named parameter.

        mgems_input_ch = merged_reads_ch
            .join( MergedPseudoalignFastqFiles.out.pseudoalignedFastqFiles, by: 0 )
            .join( MergedRunMSweep.out.msweep_merged,                       by: 0 )
            .join( MergedRunMSweep.out.msweep_unmerged,                     by: 0 )
            .map { sid,
                   merged_fq,          unmerged_fq,       // original reads
                   pseudo_merged,       pseudo_unmerged,   // themisto FASTQs
                   merged_abundances,   merged_probs,      // mSWEEP merged
                   unmerged_abundances, unmerged_probs ->  // mSWEEP unmerged
                tuple( sid,
                       merged_fq,          unmerged_fq,
                       pseudo_merged,       pseudo_unmerged,
                       merged_abundances,   merged_probs,
                       unmerged_abundances, unmerged_probs )
            }

        // Initialise mGEMS output channels to empty so the emit block is
        // always valid regardless of whether mGEMS actually runs.
        binned_reads_ch       = Channel.empty()
        assignment_tables_ch  = Channel.empty()
        mgems_summary_ch      = Channel.empty()

        if ( params.run_mgems != false ) {   // true when param is true or absent
            MergedRunMGEMS(
                mgems_input_ch,
                params.clustering_file,
                params.themisto_index
            )
            binned_reads_ch      = MergedRunMGEMS.out.binned_merged_reads
                                     .mix( MergedRunMGEMS.out.binned_unmerged_reads )
            assignment_tables_ch = MergedRunMGEMS.out.assignment_table_merged
                                     .mix( MergedRunMGEMS.out.assignment_table_unmerged )
            mgems_summary_ch     = MergedRunMGEMS.out.mgems_summary
        }

        // ── Step 4: Parse mSWEEP abundance results ────────────────────────────
        MergedParsemSweepResults(
            MergedRunMSweep.out.msweep_merged
                .mix( MergedRunMSweep.out.msweep_unmerged )
                .collect()
        )

    emit:
        msweep_summary    = MergedParsemSweepResults.out.msweep_summary
        msweep_matrix     = MergedParsemSweepResults.out.msweep_matrix
        binned_reads      = binned_reads_ch
        assignment_tables = assignment_tables_ch
        mgems_summary     = mgems_summary_ch
}
