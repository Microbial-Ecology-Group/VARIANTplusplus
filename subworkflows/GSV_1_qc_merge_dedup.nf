//
//  subworkflows/GSV_1_qc_merge_dedup.nf
//
include { runqc ; QCstats }                                from '../modules/Trimming/trimmomatic'
include { MergeReadsFlash }                                from '../modules/QC/merge'
include { SeqkitReadCounts; MergedDeduplicateReadsSeqkit }                    from '../modules/QC/dedup.nf'

workflow GSV_1_WF_old_working {

    take:
        read_pairs_ch         // ( sample_id , [R1,R2] )
        hostfasta             // FASTA if an index must be built

    main:
        /* ── QC ─────────────────────────────────────────────── */
        runqc( read_pairs_ch )
        QCstats( runqc.out.trimmomatic_stats.collect() )


        /* ── FLASH merge ────────────────────────────────────── */
        MergeReadsFlash( runqc.out.paired_fastq )

        merged_only_ch   = MergeReadsFlash.out.merged   .map { sid, f -> tuple(sid, f) }
        unmerged_only_ch = MergeReadsFlash.out.unmerged .map { sid, f -> tuple(sid, f) }

        /* ── BBMap deduplication ────────────────────────────── */
        to_dedup_ch = merged_only_ch.join( unmerged_only_ch )   // ( sample , merged , unmerged )
        MergedDeduplicateReadsSeqkit( to_dedup_ch )

    emit:
        dedup_merged   = MergedDeduplicateReadsSeqkit.out.dedup_merged
        dedup_unmerged = MergedDeduplicateReadsSeqkit.out.dedup_unmerged
}