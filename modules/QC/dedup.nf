threads = params.threads

process DeduplicateReads {
    tag { sample_id }
    label "small_memory_medium_time"

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    publishDir "${params.output}/Deduped_reads", mode: 'copy', pattern: '*.fastq.gz',
        saveAs: { filename -> "Deduped_reads/${filename}" }

    input:
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path("Deduped_reads/${sample_id}_R1_dedup.fastq.gz"), emit: dedup_r1
        tuple val(sample_id), path("Deduped_reads/${sample_id}_R2_dedup.fastq.gz"), emit: dedup_r2
        path("Deduped_reads/${sample_id}.cdhit.stats.log"), emit: cdhit_stats

    script:
    """
    mkdir -p Deduped_reads

    # Decompress input reads since CD-HIT-EST cannot process gzipped files
    gunzip -c ${reads[0]} > ${sample_id}_R1.fastq
    gunzip -c ${reads[1]} > ${sample_id}_R2.fastq

    # Run CD-HIT-EST on R1 and R2 separately at 100% identity
    cd-hit-est -i ${sample_id}_R1.fastq -o Deduped_reads/${sample_id}_R1_dedup.fastq -T ${threads} -d 1000 -c 1.00 -n 10 > Deduped_reads/${sample_id}.cdhit.stats.log
    cd-hit-est -i ${sample_id}_R2.fastq -o Deduped_reads/${sample_id}_R2_dedup.fastq -T ${threads} -d 1000 -c 1.00 -n 10 >> Deduped_reads/${sample_id}.cdhit.stats.log

    # Gzip the deduplicated output
    gzip Deduped_reads/${sample_id}_R1_dedup.fastq
    gzip Deduped_reads/${sample_id}_R2_dedup.fastq

    # Cleanup intermediate uncompressed files
    rm ${sample_id}_R1.fastq ${sample_id}_R2.fastq
    """
}
