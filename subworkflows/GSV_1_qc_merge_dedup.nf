//
//  subworkflows/GSV_1_qc_merge_dedup.nf
//
include { runqc ; QCstats }                                from '../modules/Trimming/trimmomatic'
include { MergeReadsFlash }                                from '../modules/QC/merge'
include { index }                                          from '../modules/Alignment/bwa'
include { MergedDeduplicateReadsBBMap }                    from '../modules/QC/dedup.nf'

workflow GSV_1_WF {

    take:
        read_pairs_ch         // ( sample_id , [R1,R2] )
        hostfasta             // FASTA if an index must be built

    main:
        /* ── QC ─────────────────────────────────────────────── */
        runqc( read_pairs_ch )
        QCstats( runqc.out.trimmomatic_stats.collect() )

        /* ── build or read BWA index ────────────────────────── */
        def reference_index_ch = params.host_index
            ? Channel
                .fromPath( params.host_index )
                .ifEmpty { error "No files match --host_index '${params.host_index}'" }
                .toList()
                .map { files ->  
                        if( files.size() < 6 )
                            error "Expected ≥6 BWA index files, found ${files.size()}"
                        files.sort()
                      }
            : { index( hostfasta ); index.out }()
                         // ^ anonymous closure to run `index(...)`

        /* ── FLASH merge ────────────────────────────────────── */
        MergeReadsFlash( runqc.out.paired_fastq )

        merged_only_ch   = MergeReadsFlash.out.merged   .map { sid, f -> tuple(sid, f) }
        unmerged_only_ch = MergeReadsFlash.out.unmerged .map { sid, f -> tuple(sid, f) }

        /* ── BBMap deduplication ────────────────────────────── */
        to_dedup_ch = merged_only_ch.join( unmerged_only_ch )   // ( sample , merged , unmerged )
        MergedDeduplicateReadsBBMap( to_dedup_ch )

    emit:
        dedup_merged   = MergedDeduplicateReadsBBMap.out.dedup_merged
        dedup_unmerged = MergedDeduplicateReadsBBMap.out.dedup_unmerged
}
