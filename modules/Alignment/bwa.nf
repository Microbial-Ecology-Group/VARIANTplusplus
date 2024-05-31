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

themisto_index = params.themisto_index

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
            else if(filename.indexOf("_alignment_dedup.bam") > 0) "Deduped/$filename"
            else {}
        }

    input:
        path indexfiles 
        tuple val(pair_id), path(reads) 

    output:
        tuple val(pair_id), path("${pair_id}_alignment_sorted.bam"), emit: bwa_bam
        tuple val(pair_id), path("${pair_id}_alignment_dedup.bam"), emit: bwa_dedup_bam, optional: true

    script:
    if( deduped == "N")
        """
        ${BWA} mem ${indexfiles[0]} ${reads} -t ${threads} -R '@RG\\tID:${pair_id}\\tSM:${pair_id}' > ${pair_id}_alignment.sam
        ${SAMTOOLS} view -@ ${threads} -S -b ${pair_id}_alignment.sam > ${pair_id}_alignment.bam
        rm ${pair_id}_alignment.sam
        ${SAMTOOLS} sort -@ ${threads} -n ${pair_id}_alignment.bam -o ${pair_id}_alignment_sorted.bam
        rm ${pair_id}_alignment.bam
        """
    else if( deduped == "Y")
        """
        ${BWA} mem ${indexfiles[0]} ${reads} -t ${threads} -R '@RG\\tID:${pair_id}\\tSM:${pair_id}' > ${pair_id}_alignment.sam
        ${SAMTOOLS} view -@ ${threads} -S -b ${pair_id}_alignment.sam > ${pair_id}_alignment.bam
        rm ${pair_id}_alignment.sam
        ${SAMTOOLS} sort -@ ${threads} -n ${pair_id}_alignment.bam -o ${pair_id}_alignment_sorted.bam
        rm ${pair_id}_alignment.bam
        ${SAMTOOLS} fixmate -@ ${threads} ${pair_id}_alignment_sorted.bam ${pair_id}_alignment_sorted_fix.bam
        ${SAMTOOLS} sort -@ ${threads} ${pair_id}_alignment_sorted_fix.bam -o ${pair_id}_alignment_sorted_fix.sorted.bam
        rm ${pair_id}_alignment_sorted_fix.bam
        ${SAMTOOLS} rmdup -S ${pair_id}_alignment_sorted_fix.sorted.bam ${pair_id}_alignment_dedup.bam
        rm ${pair_id}_alignment_sorted_fix.sorted.bam
        ${SAMTOOLS} view -@ ${threads} -h -o ${pair_id}_alignment_dedup.sam ${pair_id}_alignment_dedup.bam
        rm ${pair_id}_alignment_dedup.sam
        """
    else
        error "Invalid deduplication flag --deduped: ${deduped}. Please use --deduped Y for deduplicated counts, or avoid using this flag altogether to skip this error."
}




process index_genomes {
    tag "${genome.baseName}"
    scratch true

    input:
        tuple val(genomeName), path(genome)

    output:
        tuple val(genomeName), path("${genomeName}.*"), path(genome)

    script:
    """
    echo "Indexing genome at: ${genome}"
    bwa index -p ${genome.baseName} ${genome}
    """
}


process align_reads {
    tag "${sampleName}_${genomeName}"
    scratch true

    input:
        tuple val(sampleName), path(reads), val(genomeName), path(genome), path(indexFiles)
        
    output:
        tuple val(sampleName), val(genomeName), 
              path("${sampleName}_${genomeName}_alignment_sorted.dedup.bam"), 
              path("${sampleName}_${genomeName}_depth.txt"), 
              path("${sampleName}_${genomeName}_coverage.txt"), 
              path("${sampleName}_${genomeName}_processed_idxstats.txt"), 
              path("${sampleName}_${genomeName}_overall_average.txt")

    publishDir "${params.output}/Alignment_results", pattern: "*_alignment_sorted.dedup.bam", mode: 'copy'
    publishDir "${params.output}/Coverage_results", pattern: "*_depth.txt", mode: 'copy'
    publishDir "${params.output}/Coverage_results", pattern: "*_coverage.txt", mode: 'copy'
    publishDir "${params.output}/Coverage_results", pattern: "*_processed_idxstats.txt", mode: 'copy'
    publishDir "${params.output}/Coverage_results", pattern: "*_overall_average.txt", mode: 'copy'

    script:
    """
    echo "Starting alignment for: ${sampleName} with genome: ${genomeName}"
    bwa mem ${genomeName} ${reads.join(' ')} -t ${task.cpus} | \\
    samtools sort -@ ${task.cpus} -n -o ${sampleName}_${genomeName}_alignment_sorted_byname.bam
    

    echo "Fixing mate information and removing duplicates..."
    samtools fixmate -m ${sampleName}_${genomeName}_alignment_sorted_byname.bam ${sampleName}_${genomeName}_alignment_sorted_fixed.bam
    samtools sort -@ ${task.cpus} -o ${sampleName}_${genomeName}_alignment_sorted.bam ${sampleName}_${genomeName}_alignment_sorted_fixed.bam
    samtools markdup -r ${sampleName}_${genomeName}_alignment_sorted.bam ${sampleName}_${genomeName}_alignment_sorted.dedup.bam
    rm ${sampleName}_${genomeName}_alignment_sorted_byname.bam ${sampleName}_${genomeName}_alignment_sorted_fixed.bam # Clean up intermediate files

    echo "Calculating depth per base..."
    samtools depth ${sampleName}_${genomeName}_alignment_sorted.dedup.bam > ${sampleName}_${genomeName}_depth.txt

    echo "Using BEDTools for genome coverage calculations..."
    bedtools genomecov -ibam ${sampleName}_${genomeName}_alignment_sorted.dedup.bam > ${sampleName}_${genomeName}_coverage.txt

    echo "Generating IDXStats..."
    samtools index ${sampleName}_${genomeName}_alignment_sorted.dedup.bam
    samtools idxstats ${sampleName}_${genomeName}_alignment_sorted.dedup.bam > ${sampleName}_${genomeName}_idxstats.txt

    echo "Calculating IDXStats with Python script..."
    python $baseDir/bin/samtools_idxstats.py -i ${sampleName}_${genomeName}_idxstats.txt -o ${sampleName}_${genomeName}_processed_idxstats.txt
    rm ${sampleName}_${genomeName}_idxstats.txt  # Clean up intermediate file

    echo "Calculating average coverage across genome segments with custom Python script..."
    python $baseDir/bin/parse_genomecov.py -i ${sampleName}_${genomeName}_coverage.txt -o ${sampleName}_${genomeName}_detailed_coverage.txt -a ${sampleName}_${genomeName}_overall_average.txt -s ${sampleName} -g ${genomeName}
    #rm ${sampleName}_${genomeName}_coverage.txt  # Clean up intermediate file if not needed further
    """
}



process CheckAndStoreCoverage {
    tag "${sampleName}_${genomeName}"

    // Setup directories for final outputs
    publishDir "${params.output}/Filtered_classifications", mode: 'copy', 
        saveAs: { it.contains('results_to_keep/') ? it.replace('results_to_keep/', '') : null }
    publishDir "${params.output}/Filtered_reads", mode: 'copy',
        saveAs: { it.endsWith('.fastq.gz') ? it : null }

    input:
        tuple val(sampleName), val(genomeName), path(bamFile), path(depthFile),
              path(coverageFile), path(idxStatsFile), path(averageCoverageFile)

    output:
        path("results_to_keep/${sampleName}_${genomeName}_depth.txt"), optional: true, emit: storedDepthChannel
        path("results_to_keep/${sampleName}_${genomeName}_coverage.txt"), optional: true, emit: storedCoverageChannel
        path("results_to_keep/${sampleName}_${genomeName}_idxstats.txt"), optional: true, emit: storedIdxStatsChannel
        path("results_to_keep/${sampleName}_${genomeName}_overall_average.txt"), optional: true, emit: storedAvgCoverageChannel
        tuple val(sampleName), path("${sampleName}_${genomeName}_R1.fastq.gz"), path("${sampleName}_${genomeName}_R2.fastq.gz"), emit: storedFastqPairs, optional: true

    script:
    """
    echo "Checking coverage from file: ${averageCoverageFile}"
    coverageValue=\$(tail -n 1 ${averageCoverageFile} | cut -f3)
    echo "Coverage Value: \$coverageValue"
    mkdir -p results_to_keep
    if [ \$(echo "\$coverageValue > ${params.coverage_threshold}" | bc -l) -eq 1 ]; then
        echo "Coverage threshold met, copying files..."
        cp ${depthFile} ${coverageFile} ${idxStatsFile} ${averageCoverageFile} results_to_keep/
        echo "Extracting mapped reads in FASTQ format..."
        samtools fastq -@ ${task.cpus} -c 6 ${bamFile} \
            -1 ${sampleName}_${genomeName}_R1.fastq.gz \
            -2 ${sampleName}_${genomeName}_R2.fastq.gz \
            -0 /dev/null -s /dev/null -n
    else
        echo "Coverage threshold not met, not storing files."
    fi
    """
}


process MergeFastqFiles {
    tag "${sampleName}"

    publishDir "${params.output}/Filtered_merged_sample_reads", mode: 'copy'

    input:
        tuple val(sampleName), path(fastqR1), path(fastqR2)

    output:
        tuple val(sampleName), path("${sampleName}_merged_R1.fastq.gz"), path("${sampleName}_merged_R2.fastq.gz"), emit: mergedFastqFiles

    script:
    """
    cat ${fastqR1} > ${sampleName}_merged_R1.fastq.gz
    cat ${fastqR2} > ${sampleName}_merged_R2.fastq.gz
    """
}


process RunBactopia {
    tag "${sampleName}"

    publishDir "${params.output}/Bactopia_results", mode: 'copy'

    input:
        tuple val(sampleName), path(readsR1), path(readsR2)

    output:
        path("${sampleName}_bactopia_output/**"), emit: bactopiaOutput

    script:
    """
    bactopia --R1 ${readsR1} --R2 ${readsR2} \
             --species "Escherichia coli" \
             --datasets ${params.bactopia_datasets} \
             --outdir ${sampleName}_bactopia_output \
             --cpus ${task.cpus}
    """
}

process PseudoalignFastqFiles {
    tag "${sampleName}"

    publishDir "${params.output}/Filtered_pseudoaligned_reads", mode: 'copy'

    input:
        tuple val(sampleName), path(fastqR1), path(fastqR2)
        path themisto_index

    output:
        tuple val(sampleName), path("${sampleName}_pseudoaligned_R1.fastq.gz"), path("${sampleName}_pseudoaligned_R2.fastq.gz"), emit: pseudoalignedFastqFiles

    script:
    """
    themisto pseudoalign -q "${fastqR1}" -i ${themisto_index} -t ${threads} --gzip-output-lines --sort-output -o "${sampleName}_pseudoaligned_R1.fastq.gz"
    themisto pseudoalign -q "${fastqR2}" -i ${themisto_index} -t ${threads} --gzip-output-lines --sort-output -o "${sampleName}_pseudoaligned_R2.fastq.gz"
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
