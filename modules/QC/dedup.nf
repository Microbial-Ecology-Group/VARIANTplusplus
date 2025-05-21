threads = params.threads

process DeduplicateReads {
    tag { sample_id }
    label "small_memory_medium_time"

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    publishDir "${params.output}/Deduped_reads", mode: 'copy', pattern: '*.fastq.gz',
        saveAs: { filename -> "${filename}" }

    input:
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path("${sample_id}_R1_dedup.fastq.gz"), emit: dedup_r1
        tuple val(sample_id), path("${sample_id}_R2_dedup.fastq.gz"), emit: dedup_r2
        path("${sample_id}.cdhit.stats.log"), emit: cdhit_stats

    script:
    """
    mkdir -p Deduped_reads

    # Decompress input reads since CD-HIT-EST cannot process gzipped files
    gunzip -c ${reads[0]} > ${sample_id}_R1.fastq
    gunzip -c ${reads[1]} > ${sample_id}_R2.fastq

    # Run CD-HIT-EST on R1 and R2 separately at 100% identity
    cd-hit-est -i ${sample_id}_R1.fastq -o ${sample_id}_R1_dedup.fastq -M 0 -T ${task.cpus} -d 1000 -c 1.00 -n 10 > ${sample_id}.cdhit.stats.log
    cd-hit-est -i ${sample_id}_R2.fastq -o ${sample_id}_R2_dedup.fastq -M 0 -T ${task.cpus} -d 1000 -c 1.00 -n 10 >> ${sample_id}.cdhit.stats.log

    # Gzip the deduplicated output
    gzip ${sample_id}_R1_dedup.fastq
    gzip ${sample_id}_R2_dedup.fastq

    # Cleanup intermediate uncompressed files
    rm ${sample_id}_R1.fastq ${sample_id}_R2.fastq
    """
}

process DeduplicateReadsBBMap {
    tag { sample_id }
    label "small_memory_medium_time"

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    publishDir "${params.output}/Deduped_reads", mode: 'copy', pattern: '*.fastq.gz',
        saveAs: { filename -> "${filename}" }

    input:
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path("${sample_id}_R1_dedup.fastq.gz"), emit: dedup_r1
        tuple val(sample_id), path("${sample_id}_R2_dedup.fastq.gz"), emit: dedup_r2
        path("${sample_id}.dedupe.stats.log"), emit: dedupe_stats

    script:
    """
    # Run BBMap's dedupe.sh separately for R1 and R2 reads
    dedupe.sh \
        in=${reads[0]} \
        out=${sample_id}_R1_dedup.fastq.gz \
        threads=${task.cpus} \
        > ${sample_id}_R1.dedupe.stats.log 2>&1

    dedupe.sh \
        in=${reads[1]} \
        out=${sample_id}_R2_dedup.fastq.gz \
        threads=${task.cpus} \
        > ${sample_id}_R2.dedupe.stats.log 2>&1

    # Merge logs
    cat ${sample_id}_R1.dedupe.stats.log ${sample_id}_R2.dedupe.stats.log > ${sample_id}.dedupe.stats.log
    rm ${sample_id}_R1.dedupe.stats.log ${sample_id}_R2.dedupe.stats.log
    """
}

process MergedDeduplicateReadsBBMap {

    tag   { sample_id }
    label "dedup"

    publishDir "${params.output}/Deduped_reads",
               mode:'copy', pattern:'*.fastq.gz'

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    input:
        tuple val(sample_id),
              path(merged_fq),        // one file
              path(unmerged_fq)       // one file

    output:
        tuple val(sample_id), path("${sample_id}_merged.dedup.fastq.gz"),   emit: dedup_merged
        tuple val(sample_id), path("${sample_id}_unmerged.dedup.fastq.gz"), emit: dedup_unmerged
        path("${sample_id}.dedupe.stats.log"),                              emit: dedupe_stats

    script:
    def xmx = task.memory.toMega() * 0.9 as Integer
    """
# ── deduplicate UNMERGED file ────────────────────────────────────────
    dedupe.sh -Xmx${xmx}m \
          in=${unmerged_fq} \
          out=${sample_id}_unmerged.dedup.fastq.gz \
          threads=$task.cpus \
    > ${sample_id}_unmerged.dedupe.log 2>&1

# ── deduplicate MERGED file ──────────────────────────────────────────
    dedupe.sh -Xmx${xmx}m \
          in=${merged_fq} \
          out=${sample_id}_merged.dedup.fastq.gz \
          threads=$task.cpus \
         > ${sample_id}_merged.dedupe.log 2>&1

    # ── merge the per-file logs into one stats file ─────────────────────
    cat ${sample_id}_merged.dedupe.log ${sample_id}_unmerged.dedupe.log \
        >  ${sample_id}.dedupe.stats.log
    rm  ${sample_id}_merged.dedupe.log  ${sample_id}_unmerged.dedupe.log
    """
}

process MergedDeduplicateReadsSeqkit {

    tag   { sample_id }
    label "medium"

    publishDir "${params.output}/Deduped_reads",
               mode:'copy', pattern:'*.fastq.gz'

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    input:
        tuple val(sample_id),
              path(merged_fq),        // one file
              path(unmerged_fq)       // one file

    output:
        tuple val(sample_id), path("${sample_id}_merged.dedup.fastq.gz"),   emit: dedup_merged
        tuple val(sample_id), path("${sample_id}_unmerged.dedup.fastq.gz"), emit: dedup_unmerged
        path("${sample_id}.dedupe_seqkit.stats.log"),                       emit: dedupe_stats

    script:
    """
    # ── deduplicate MERGED file ──────────────────────────────────────────
    seqkit rmdup --threads $task.cpus --by-seq \
          -o ${sample_id}_merged.dedup.fastq.gz \
          ${merged_fq} \
        > ${sample_id}_merged.dedupe_seqkit.log 2>&1

    # ── deduplicate UNMERGED file ─────────────────────────────────────────
    seqkit rmdup --threads $task.cpus --by-seq \
          -o ${sample_id}_unmerged.dedup.fastq.gz \
          ${unmerged_fq} \
        > ${sample_id}_unmerged.dedupe_seqkit.log 2>&1

    # ── merge the per-file logs into one stats file ────────────────────────
    cat ${sample_id}_merged.dedupe_seqkit.log ${sample_id}_unmerged.dedupe_seqkit.log \
        >  ${sample_id}.dedupe_seqkit.stats.log
    rm  ${sample_id}_merged.dedupe_seqkit.log  ${sample_id}_unmerged.dedupe_seqkit.log
    """
}




process SeqkitReadCounts {

    tag "seqkit"
    label "micro"

    publishDir "${params.output}/Results/Stats", mode: 'copy'

    /*
     * fastqs  – a *list* of Path objects (all merged & unmerged reads)
     * prefix  – e.g. params.QC_prefix  →  <prefix>_reads.txt
     */
    input:
        path fastqs            // list because we used collect()
        val  prefix

    output:
        path("${prefix}_reads.txt"), emit: readCounts

    script:
    """
    set -euo pipefail
    seqkit stats -T -j ${task.cpus} ${fastqs.join(' ')} \
        > ${prefix}_reads.txt
    """
}
