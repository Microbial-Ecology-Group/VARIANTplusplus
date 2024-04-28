include { reference_error ; amr_error ; annotation_error } from "$baseDir/modules/nf-functions.nf"


if( params.amr ) {
    amr = file(params.amr)
    if( !amr.exists() ) return amr_error(amr)
}
if( params.annotation ) {
    annotation = file(params.annotation)
    if( !annotation.exists() ) return annotation_error(annotation)
}


threads = params.threads

deduped = params.deduped

process index {
    tag "Creating bwa index"
    label "alignment"

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    publishDir "${params.output}/Alignment/BWA_Index", mode: "copy"

    input:
    path fasta

    output: 
    path("${fasta}*"), emit: bwaindex, includeInputs: true

    script:
    """
    bwa index ${fasta}
    #--threads $task.cpus 
    """
}


process bwa_align {
    tag "$pair_id"
    label "alignment"

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    publishDir "${params.output}/Alignment/BAM_files", mode: "copy",
        saveAs: { filename ->
            if(filename.indexOf("_alignment_sorted.bam") > 0) "Standard/$filename"
            else {}
        }

    input:
        path indexfiles 
        tuple val(pair_id), path(reads) // Reads can be single-end or paired-end

    output:
        tuple val(pair_id), path("${pair_id}_alignment_sorted.bam"), emit: bwa_bam

    script:
    """
    ${BWA} mem ${indexfiles[0]} ${reads} -t ${threads} -R '@RG\\tID:${pair_id}\\tSM:${pair_id}' > ${pair_id}_alignment.sam
    ${SAMTOOLS} view -@ ${threads} -S -b ${pair_id}_alignment.sam > ${pair_id}_alignment_sorted.bam
    rm ${pair_id}_alignment.sam
    ${SAMTOOLS} sort -@ ${threads} -o ${pair_id}_alignment_sorted.bam
    """
}

process bwa_all_to_all {
    tag "bwa_all_to_all_${uuid()}" // Unique identifier for each alignment set

    label "alignment"

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    publishDir "${params.output}/Alignment/BAM_files", mode: "copy",
        saveAs: { filename ->
            if (filename.indexOf("_alignment_sorted.bam") > 0) "Standard/$filename"
            else {}
        }

    input:
        path indexfiles 
        path(reads) // Reads (single-end or paired-end)
        val(ref_files) // List of reference file names

    output:
        path("*_alignment_sorted.bam") // Allow multiple BAM outputs based on input combinations

    script:
    """
    for ref_file in ${ref_files.join(' ')}; do
        ref_base=$(basename $ref_file)
        ${BWA} mem $ref_file ${reads} -t ${threads} -R '@RG\\tID:$ref_base\\tSM:$ref_base' > ${ref_base}_alignment.sam
        ${SAMTOOLS} view -@ ${threads} -S -b ${ref_base}_alignment.sam > ${ref_base}_alignment_sorted.bam
        rm ${ref_base}_alignment.sam
        ${SAMTOOLS} sort -@ ${threads} -o ${ref_base}_alignment_sorted.bam
    done
    """
}



process bwa_rm_contaminant_fq {
    tag { pair_id }
    label "alignment"

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3 
 
    publishDir "${params.output}/HostRemoval", mode: "copy",
        saveAs: { filename ->
            if(filename.indexOf("fastq.gz") > 0) "NonHostFastq/$filename"
            else {}
        }

    input:
    path indexfiles
    tuple val(pair_id), path(reads) 

    output:
    tuple val(pair_id), path("${pair_id}.non.host.R*.fastq.gz"), emit: nonhost_reads
    path("${pair_id}.samtools.idxstats"), emit: host_rm_stats
    
    """
    ${BWA} mem ${indexfiles[0]} ${reads[0]} ${reads[1]} -t ${threads} > ${pair_id}.host.sam
    ${SAMTOOLS} view -bS ${pair_id}.host.sam | ${SAMTOOLS} sort -@ ${threads} -o ${pair_id}.host.sorted.bam
    rm ${pair_id}.host.sam
    ${SAMTOOLS} index ${pair_id}.host.sorted.bam && ${SAMTOOLS} idxstats ${pair_id}.host.sorted.bam > ${pair_id}.samtools.idxstats
    ${SAMTOOLS} view -h -f 12 -b ${pair_id}.host.sorted.bam -o ${pair_id}.host.sorted.removed.bam
    ${SAMTOOLS} sort -n -@ ${threads} ${pair_id}.host.sorted.removed.bam -o ${pair_id}.host.resorted.removed.bam
    ${SAMTOOLS}  \
       fastq -@ ${threads} -c 6  \
      ${pair_id}.host.resorted.removed.bam \
      -1 ${pair_id}.non.host.R1.fastq.gz \
      -2 ${pair_id}.non.host.R2.fastq.gz \
      -0 /dev/null -s /dev/null -n

    rm *.bam
    """

}

process HostRemovalStats {
    tag { sample_id }
    label "alignment"

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3 

    publishDir "${params.output}/Results", mode: "copy",
        saveAs: { filename ->
            if(filename.indexOf(".stats") > 0) "Stats/$filename"
        }

    input:
        file(host_rm_stats)

    output:
        path("host.removal.stats"), emit: combo_host_rm_stats

    """
    ${PYTHON3} $baseDir/bin/samtools_idxstats.py -i ${host_rm_stats} -o host.removal.stats
    """
}
