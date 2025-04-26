threads = params.threads

dedup_sam = params.dedup_sam

themisto_index = params.themisto_index

clustering_file = params.clustering_file

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
        tuple val(sampleName), path(reads)
        path themisto_index

    output:
        tuple val(sampleName), path("${sampleName}_pseudoaligned_R1.fastq.gz"), path("${sampleName}_pseudoaligned_R2.fastq.gz"), emit: pseudoalignedFastqFiles

    script:
    """
    mkdir -p tmp
    $baseDir/bin/themisto pseudoalign -q "${reads[0]}" -i ${themisto_index}/themisto_index --temp-dir tmp -t ${threads} --gzip-output --sort-output-lines -o "${sampleName}_pseudoaligned_R1.fastq"
    $baseDir/bin/themisto pseudoalign -q "${reads[1]}" -i ${themisto_index}/themisto_index --temp-dir tmp -t ${threads} --gzip-output --sort-output-lines -o "${sampleName}_pseudoaligned_R2.fastq"
    rm -rf tmp/
    """
}


process MergedPseudoalignFastqFiles {

    tag { sample_id }
    label "large_memory"

    publishDir "${params.output}/Filtered_pseudoaligned_reads", mode:'copy'

    input:
        tuple val(sample_id),
              path(merged_fastq),          //  …_Mh_extracted_merged.fastq.gz
              path(unmerged_fastq)         //  …_Mh_extracted_unmerged.fastq.gz
        path themisto_index

    output:
        tuple val(sample_id),
              path("${sample_id}_pseudoaligned_merged.fastq.gz"),
              path("${sample_id}_pseudoaligned_unmerged.fastq.gz"),
              emit: pseudoalignedFastqFiles

    script:
    """
    mkdir -p tmp

    # ── merged stream ─────────────────────────────────────────
    ${baseDir}/bin/themisto pseudoalign \
        -q ${merged_fastq} \
        -i ${themisto_index}/2025_themisto_index_no \
        --temp-dir tmp -t ${threads} \
        --gzip-output --sort-output-lines \
        -o ${sample_id}_pseudoaligned_merged.fastq

    # ── un-merged stream (now interleaved single FASTQ) ───────
    ${baseDir}/bin/themisto pseudoalign \
        -q ${unmerged_fastq} \
        -i ${themisto_index}/2025_themisto_index_no \
        --temp-dir tmp -t ${threads} \
        --gzip-output --sort-output-lines \
        -o ${sample_id}_pseudoaligned_unmerged.fastq

    rm -rf tmp
    """
}

process RunMSweep {
    tag "${sampleName}"

    publishDir "${params.output}/mSWEEP_results", mode: 'copy'

    input:
        tuple val(sampleName), path(pseudoaligned)
        path clustering_file

    output:
        tuple val(sampleName), path("${sampleName}.msweep.txt"), emit: msweepResults

    script:
    """
    mSWEEP --themisto-1 "${pseudoaligned[0]}" --themisto-2 "${pseudoaligned[1]}" -i ${clustering_file} -t ${threads} -o ${sampleName}.msweep.txt --verbose
    """
}


process MergedRunMSweep {

    tag { sample_id }
    label "small_memory"

    publishDir "${params.output}/mSWEEP_results", mode: 'copy'

    input:
        tuple val(sample_id),
              path(pseudo_merged),        //  …pseudoaligned_merged.fastq.gz
              path(pseudo_unmerged)       //  …pseudoaligned_unmerged.fastq.gz
        path clustering_file

    output:
        path("${sample_id}.merged.msweep_abundances.txt"),   emit: msweep_merged
        path("${sample_id}.unmerged.msweep_abundances.txt"), emit: msweep_unmerged

    script:
    """
    # ── merged read‐set ─────────────────────────────────────────
    mSWEEP --themisto ${pseudo_merged} \
           -i ${clustering_file} \
           -t ${threads} \
           --themisto-mode intersection \
           -o ${sample_id}.merged.msweep \
           --verbose

    # ── unmerged read‐set ───────────────────────────────────────
    mSWEEP --themisto ${pseudo_unmerged} \
           -i ${clustering_file} \
           -t ${threads} \
           --themisto-mode intersection \
           -o ${sample_id}.unmerged.msweep \
           --verbose
    """
}

process MergedParsemSweepResults {
    tag "parse_msweep"
    label "small_memory"                         

    publishDir "${params.output}/Results", mode: 'copy'

    input:
        path(msweep_files)                 

    output:
        path("*_summary.tsv"),      emit: msweep_summary
        path("*_count_matrix.tsv"), emit: msweep_matrix

    script:
    """

    python $baseDir/bin/parse_msweep_results.py \
        --msweep_dir . \
        --reads_dir  $baseDir/${params.output}/HostRemoval/NonHostFastq/ \
        -o mSweep_results \
        --filter-mode sub_count_rel_abund \
        --rel-thr 0.01
    """
}


process bwa_rm_contaminant_fq {

    tag   { pair_id }
    label "alignment"
    maxRetries 3
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }

    publishDir "${params.output}/HostRemoval", mode:'copy',
        saveAs:{ fn -> fn.endsWith('.fastq.gz') ? "NonHostFastq/$fn" : null }

    input:
        path  indexfiles
        tuple val(pair_id), path(reads)          // 1 or 2 files

    output:
        tuple val(pair_id), path("${pair_id}.non.host.fastq.gz"), emit: nonhost_reads
        tuple val(pair_id), path("${pair_id}.samtools.idxstats"),  emit: host_rm_stats

    script:
    """

    ${BWA} mem -p ${indexfiles[0]} ${reads[0]} -t ${threads} \
            | ${SAMTOOLS} sort -@ ${threads} -o ${pair_id}.host.sorted.bam


    ${SAMTOOLS} index   ${pair_id}.host.sorted.bam
    ${SAMTOOLS} idxstats ${pair_id}.host.sorted.bam > ${pair_id}.samtools.idxstats

    # keep ALL unmapped reads (flag 4) → one interleaved FASTQ-gz
    ${SAMTOOLS} view -b -f 4 ${pair_id}.host.sorted.bam 
    ${SAMTOOLS} fastq -@ ${threads} -f 4 \
    ${pair_id}.host.sorted.bam \
        | pigz -p ${threads} -c \
        > ${pair_id}.non.host.fastq.gz
    """
}


process HostRemovalStats {

    tag { "${sample_id}_${idxstats.name.tokenize('_')[-2]}" }   // merged / unmerged

    publishDir "${params.output}/Results/Stats", mode: 'copy'

    input:
        tuple val(sample_id), path(idxstats)        // idxstats is ONE Path

    output:
        path "${sample_id}_${idxstats.name}.summary" , emit: stats_summary

    script:
    """
    ${PYTHON3} $baseDir/bin/samtools_idxstats.py \
        -i ${idxstats} \
        -o ${sample_id}_${idxstats.name}.summary
    """
}
