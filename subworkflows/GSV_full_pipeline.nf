//
//  subworkflows / master / master_gsv.nf
//
include { GSV_1_WF } from './GSV_step_1_qc_merge.nf'
include { GSV_2_WF } from './GSV_step_2_dedup.nf'
include { GSV_3_WF } from './GSV_step_3_host_rm.nf'
include { GSV_4_WF } from './GSV_step_4_kraken_extraction.nf'
include { GSV_5_WF } from './GSV_step_5_GSV_classification.nf'

workflow GSV_PIPELINE_WF {

    /*
     * 1 · raw paired FASTQs arrive from the command line
     *     channel layout:  ( sample_id , [ R1 , R2 ] )
     */
    take:
        read_pairs_ch
        hostfasta


    main:
    /*
     * 2 · STEP 1  — QC & FLASH merge
     *     emits:   merged , unmerged   (tuple streams)
     */
        GSV_1_WF( read_pairs_ch )

        GSV_1_WF.out.merged
              .join( GSV_1_WF.out.unmerged )
              .set { merged_reads_ch }


    /*
     * 3 · STEP 2  — duplicate removal
     *     takes the 2 × FASTQ streams from step 1
     *     emits: dedup_merged , dedup_unmerged
     */
        GSV_2_WF(
                merged_reads_ch
        )

        /* bundle deduped FASTQs into a single channel for later steps */
        GSV_2_WF.out.dedup_merged
              .join( GSV_2_WF.out.dedup_unmerged )
              .set { dedup_reads_ch }

    /*
     * 4 · STEP 3  — host removal + idxstats
     *     uses dedup FASTQs + reference FASTA
     *     (step -3 is optional; comment out if not needed)
     */
        GSV_3_WF( dedup_reads_ch, hostfasta)

        
    /*
     * 5 · STEP 4  — Kraken2 on dedup nonhost FASTQs
     */
        GSV_4_WF( GSV_3_WF.out.nonhost_reads_ch )
        

    /*
     * 6 · STEP 5  — Themisto + mSWEEP on dedup FASTQs
     */
    
        GSV_5_WF( GSV_4_WF.out.extracted_reads_ch )


    emit:
        kraken_extracted_reads = GSV_4_WF.out.extracted_reads_ch

}
