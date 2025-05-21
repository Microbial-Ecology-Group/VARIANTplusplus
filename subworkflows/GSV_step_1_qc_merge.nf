//
//  subworkflows/GSV_step_1_qc_merge.nf
//
include { runqc ; QCstats }                                from '../modules/Trimming/trimmomatic'
include { MergeReadsFlash }                                from '../modules/QC/merge'
include { SeqkitReadCounts as SeqkitRawCounts ; SeqkitReadCounts as SeqkitMergedCounts  }                    from '../modules/QC/dedup.nf'

workflow GSV_1_WF {

    take:
        read_pairs_ch         // ( sample_id , [R1,R2] )

    main:
        /* 1 ── Raw read counts ─────────────────────────────────────────────── */
        def raw_fastqs_ch = read_pairs_ch
                          .map { sid, pair -> pair }   // keep only [R1,R2]
                          .flatMap { it }              // emit R1 and R2 separately
                          .collect()                   // wait for *every* file

        SeqkitRawCounts( raw_fastqs_ch, "Raw" )
        
        /* 2 ── QC ─────────────────────────────────────────────── */
        runqc( read_pairs_ch )
        QCstats( runqc.out.trimmomatic_stats.collect() )


        /* 3 ── FLASH merge ────────────────────────────────────── */
        MergeReadsFlash( runqc.out.paired_fastq )

        /* grab only the file paths (not the sample IDs) */
        def merged_paths   = MergeReadsFlash.out.merged  .map { sid, f -> f }
        def unmerged_paths = MergeReadsFlash.out.unmerged.map { sid, f -> f }

        /* 5 ─ final read-count summary ----------------- */
        def all_fastqs_ch = merged_paths
                              .mix( unmerged_paths )   // merge both streams
                              .collect()               // <- waits for channel to close
                              
        SeqkitMergedCounts(all_fastqs_ch, "QC_merged")

    emit:
        merged   = MergeReadsFlash.out.merged
        unmerged = MergeReadsFlash.out.unmerged
}
