// Load modules
include { runqc ; QCstats } from '../modules/Trimming/trimmomatic'
include { MergeReadsFlash } from '../modules/QC/merge'
include { index ; bwa_align ; bwa_rm_contaminant_fq ; HostRemovalStats; MergedPseudoalignFastqFiles ; MergedRunMSweep ; MergedParsemSweepResults } from '../modules/Alignment/bwa'
include { dlkraken ; runkraken_merged_extract } from '../modules/Microbiome/kraken2'
include { MergedDeduplicateReadsBBMap } from '../modules/QC/dedup.nf' 


import java.nio.file.Paths


// ── workflow ──────────────────────────────────────────────────────────────
workflow GSV_PIPELINE_WF {

    take:
        read_pairs_ch      // ( sample_id , [R1,R2] )
        hostfasta          // FASTA if an index must be built

    main:

    /*──────────── QC (Trimmomatic) ───────────────────────────────────────*/
    runqc( read_pairs_ch )
    QCstats( runqc.out.trimmomatic_stats.collect() )

    /*──────────── get or build BWA index ─────────────────────────────────*/
    def reference_index_ch

    if( params.host_index ) {

        reference_index_ch = Channel
            .fromPath( params.host_index )
            .ifEmpty { error "No files match --host_index '${params.host_index}'" }
            .toList()
            .map { files ->
                  if( files.size() < 6 )
                      error "Expected ≥6 BWA index files, found ${files.size()}"
                  files.sort()
            }

    } else {
        index( hostfasta )
        reference_index_ch = index.out
    }


    /*──────── merge overlapping reads (FLASH) ───────────────────────────*/
    MergeReadsFlash( runqc.out.paired_fastq )
    // After MergeReadsFlash
    merged_only_ch   = MergeReadsFlash.out.merged   .map { sid, f -> tuple(sid, f) }
    unmerged_only_ch = MergeReadsFlash.out.unmerged .map { sid, f -> tuple(sid, f) }
    
     /*──────── Deduplicate (BBmap) ───────────────────────────*/
    // join on sample → (sample, merged.fastq, unmerged.fastq)
    to_dedup_ch = merged_only_ch.join( unmerged_only_ch )

    MergedDeduplicateReadsBBMap( to_dedup_ch )
    
    // host removal now uses the two deduplicated files:
    dedup_merged_ch   = MergedDeduplicateReadsBBMap.out.dedup_merged
    dedup_unmerged_ch = MergedDeduplicateReadsBBMap.out.dedup_unmerged
    
    // re-shape to the structure expected by bwa_rm_contaminant_fq
    to_host_rm_ch =
        dedup_merged_ch   .map { sid, f -> tuple("${sid}_merged"  , [f]) }
        .mix(
        dedup_unmerged_ch .map { sid, f -> tuple("${sid}_unmerged", [f]) } )
    
    bwa_rm_contaminant_fq( reference_index_ch, to_host_rm_ch )
        
    /* run the summariser — pass the channel, NO collect() */
    HostRemovalStats( bwa_rm_contaminant_fq.out.host_rm_stats )
    
    
    /*────────────────── reshape non-host reads for Kraken ───────────*/

    // bwa_rm_contaminant_fq emits tuples ( pair_id , [paths] )
    //   pair_id examples:  S1_test_merged   S1_test_unmerged
    //   merged  list has length 1
    //   unmerged list has length 2  (R1 & R2)
    
    merged_nonhost_ch   = bwa_rm_contaminant_fq.out.nonhost_reads \
            .filter { pid, _ -> pid.endsWith('_merged') } \
            .map    { pid, files ->
                        def sid = pid.replaceFirst(/_merged$/,'')
                        tuple( sid, files[0] )            // (sample, mergedFile)
                     }
    
    unmerged_nonhost_ch = bwa_rm_contaminant_fq.out.nonhost_reads \
            .filter { pid, _ -> pid.endsWith('_unmerged') } \
            .map    { pid, files ->
                        def sid = pid.replaceFirst(/_unmerged$/,'')
                        tuple( sid, files )               // (sample, [R1,R2])
                     }
    
    /* join on the sample key → ( sample , mergedFile , [R1,R2] ) */
    kraken_input_ch = merged_nonhost_ch.join( unmerged_nonhost_ch )
            .map { sid, mergedFile, pairList -> tuple( sid, mergedFile, pairList ) }
    
    /*────────────────── decide which Kraken2 DB to use ───────────────*/
    def default_db_path = "$baseDir/data/kraken_db/k2_standard_08gb_20250402/"
    def kraken_db_ch
    
    if( file(default_db_path).exists() ) {
        kraken_db_ch = Channel.value( default_db_path )
    
    } else if( params.kraken_db ) {
        kraken_db_ch = Channel.value( params.kraken_db )
    
    } else {
        dlkraken()
        kraken_db_ch = dlkraken.out            // path produced by the downloader
    }
    
    /*──────────────── reshape for Kraken ───────────────────────────*/
    /*    bwa_rm_contaminant_fq.out.nonhost_reads  -->  (pair_id, file)       */
    
    merged_fastq_ch   = bwa_rm_contaminant_fq.out.nonhost_reads \
        .filter { pid, _ -> pid.endsWith('_merged') } \
        .map    { pid, f -> tuple(pid.replaceFirst(/_merged$/,''), f) }
    
    unmerged_fastq_ch = bwa_rm_contaminant_fq.out.nonhost_reads \
        .filter { pid, _ -> pid.endsWith('_unmerged') } \
        .map    { pid, f -> tuple(pid.replaceFirst(/_unmerged$/,''), f) }
    
    /* join on sample → ( sample , merged.fastq , unmerged.fastq ) */
    kraken_input_ch = merged_fastq_ch.join(unmerged_fastq_ch)
    
            
    runkraken_merged_extract( kraken_input_ch , kraken_db_ch )
    
    krakenresults(
          runkraken_merged_extract.out.kraken_report_merged
            .mix( runkraken_merged_extract.out.kraken_report_unmerged )
            .collect()
        )  
    
    /*──────────────────── gather FASTQs for Themisto ───────────────*/
    /*  runkraken_merged_extract emits:                             *
     *    extracted_merged   → ( sid , merged.fastq.gz )            *
     *    extracted_unmerged → ( sid , unmerged.fastq.gz )          */
    
    merged_extracted_ch   = runkraken_merged_extract.out.extracted_merged
    unmerged_extracted_ch = runkraken_merged_extract.out.extracted_unmerged
    
    themisto_input_ch = merged_extracted_ch.join( unmerged_extracted_ch )
          .map { sid, mergedFq, unmergedFq -> tuple( sid, mergedFq, unmergedFq ) }
    
    /*──────────────────── run Themisto pseudo-alignment ───────────*/
    MergedPseudoalignFastqFiles(
            themisto_input_ch, params.themisto_index
    )

    /*────────────────── pass pseudoaligned FASTQs to mSWEEP ─────────*/
    // PseudoalignFastqFiles already emitted
    //   ( sample , merged.fastq.gz , unmerged.fastq.gz )
    MergedRunMSweep(
            MergedPseudoalignFastqFiles.out.pseudoalignedFastqFiles, params.clustering_file
            
    )
    
    /* ───────────── run parse-script once all .txt exist ─────────────── */
    MergedParsemSweepResults(
        MergedRunMSweep.out.msweep_merged
                         .mix( MergedRunMSweep.out.msweep_unmerged )
                         .collect()          // single list token
    )

    emit:
        nonhost_reads = bwa_rm_contaminant_fq.out.nonhost_reads
}



