/*───────────────────────────────────────────────────────────────────────────
 *  GSV_5_MGEMS_WF — per-sample GSV read binning
 *  Input : ( sample_id, merged_extracted_fq, unmerged_extracted_fq )
 *  Output: per-GSV FASTQs named <sampleID>_GSV_<N>.fastq.gz
 *───────────────────────────────────────────────────────────────────────────*/
include { GSV_bin_reads } from '../modules/Alignment/msweep'

workflow GSV_5_MGEMS_WF {

    take:
        merged_reads_ch          // ( sample_id, merged_fq, unmerged_fq )

    main:
        GSV_bin_reads(
            merged_reads_ch,
            params.themisto_index,
            params.clustering_file
        )

    emit:
        binned_reads = GSV_bin_reads.out.binned_reads
        abundances   = GSV_bin_reads.out.abundances
}
