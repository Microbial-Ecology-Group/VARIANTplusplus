// ===============================================================
// Updated VARIANT++ modules with mSWEEP optimizations and mGEMS
// Replace the relevant processes in your modules/Alignment/bwa file
// ===============================================================

threads = params.threads
dedup_sam = params.dedup_sam
themisto_index = params.themisto_index
clustering_file = params.clustering_file

// ... [keep all your existing processes like index, bwa_align, etc.] ...

// Updated pseudoalignment process (same as original but with optional threshold)
process MergedPseudoalignFastqFiles {

    tag   { sample_id }
    label "small"

    publishDir "${params.output}/Filtered_pseudoaligned_reads", mode: 'copy'

    input:
        tuple val(sample_id),
              path(merged_fastq),          // …_Mh_extracted_merged.fastq.gz
              path(unmerged_fastq)         // …_Mh_extracted_unmerged.fastq.gz
        path themisto_index

    output:
        tuple val(sample_id),
              path("${sample_id}_pseudoaligned_merged.fastq.gz",  optional: true),
              path("${sample_id}_pseudoaligned_unmerged.fastq.gz", optional: true),
              emit: pseudoalignedFastqFiles

    script:
    def threshold_param = params.themisto_threshold ? "--threshold ${params.themisto_threshold}" : ""
    
    """
    mkdir -p tmp

    # ─ merged reads ───────────────────────────────────────────
    ${baseDir}/bin/themisto pseudoalign \
        -q ${merged_fastq} \
        -i ${themisto_index}/2025_themisto_index_no \
        --temp-dir tmp -t ${task.cpus} \
        ${threshold_param} \
        --gzip-output --sort-output-lines \
        -o ${sample_id}_pseudoaligned_merged.fastq || true

    # ─ unmerged (interleaved) reads ───────────────────────────
    ${baseDir}/bin/themisto pseudoalign \
        -q ${unmerged_fastq} \
        -i ${themisto_index}/2025_themisto_index_no \
        --temp-dir tmp -t ${task.cpus} \
        ${threshold_param} \
        --gzip-output --sort-output-lines \
        -o ${sample_id}_pseudoaligned_unmerged.fastq || true

    rm -rf tmp
    """
}

// Updated mSWEEP process with optimized parameters and probability output
process MergedRunMSweep {

    tag   { sample_id }
    label "small_memory_short_time"

    publishDir "${params.output}/mSWEEP_results", mode: 'copy'

    input:
        tuple val(sample_id), path(pseudo_merged), path(pseudo_unmerged)
        path  clustering_file

    output:
        path("${sample_id}.merged.msweep_abundances.txt"),   optional: true, emit: msweep_merged
        path("${sample_id}.merged.msweep_probs.csv"),        optional: true, emit: msweep_merged_probs
        path("${sample_id}.unmerged.msweep_abundances.txt"), optional: true, emit: msweep_unmerged  
        path("${sample_id}.unmerged.msweep_probs.csv"),      optional: true, emit: msweep_unmerged_probs

    script:
    // Use optimized parameters if provided, otherwise use defaults
    def min_hits = params.msweep_min_hits ?: 1  // Keep default conservative for short reads
    def alpha_prior = params.msweep_alpha_prior ?: 1.0
    def zero_inflation = params.msweep_zero_inflation ?: 0.01
    def q_param = params.msweep_q ?: 0.65
    def e_param = params.msweep_e ?: 0.01
    def write_probs = params.msweep_write_probs ? "--write-probs" : ""

    """
    # merged reads with optimized parameters ────────────────────
    if [ -f "${pseudo_merged}" ] && [ -s "${pseudo_merged}" ]; then
        mSWEEP \
            --themisto ${pseudo_merged} \
            -i ${clustering_file} \
            -t ${task.cpus} \
            --themisto-mode intersection \
            --min-hits ${min_hits} \
            --alphas ${alpha_prior} \
            --zero-inflation ${zero_inflation} \
            -q ${q_param} \
            -e ${e_param} \
            ${write_probs} \
            -o ${sample_id}.merged.msweep \
            --verbose || true        # ignore any mSWEEP failure
    fi

    # unmerged reads with optimized parameters ──────────────────
    if [ -f "${pseudo_unmerged}" ] && [ -s "${pseudo_unmerged}" ]; then
        mSWEEP \
            --themisto ${pseudo_unmerged} \
            -i ${clustering_file} \
            -t ${task.cpus} \
            --themisto-mode intersection \
            --min-hits ${min_hits} \
            --alphas ${alpha_prior} \
            --zero-inflation ${zero_inflation} \
            -q ${q_param} \
            -e ${e_param} \
            ${write_probs} \
            -o ${sample_id}.unmerged.msweep \
            --verbose || true
    fi
    """
}

// Optional new mGEMS process for read binning (add this if you want the full mGEMS workflow)
process MergedRunMGEMS {

    tag   { sample_id }
    label "medium"

    publishDir "${params.output}/mGEMS_results", mode: 'copy',
        saveAs: { filename ->
            if (filename.endsWith('_assignment_table.tsv')) "Assignment_tables/${filename}"
            else if (filename.contains('_binned_')) "Binned_reads/${filename}"
            else filename
        }

    input:
        tuple val(sample_id), 
              path(merged_fastq), path(unmerged_fastq),           // Original FASTQ files
              path(pseudo_merged_fastq), path(pseudo_unmerged_fastq), // Pseudoaligned FASTQ (not used)
              path(msweep_merged_abundances), path(msweep_merged_probs),
              path(msweep_unmerged_abundances), path(msweep_unmerged_probs)
        path clustering_file
        path themisto_index

    output:
        // Binned reads
        path("${sample_id}_merged_binned_*.fastq.gz"),    optional: true, emit: binned_merged_reads
        path("${sample_id}_unmerged_binned_*.fastq.gz"),  optional: true, emit: binned_unmerged_reads
        
        // Assignment tables 
        path("${sample_id}_merged_assignment_table.tsv"), optional: true, emit: assignment_table_merged
        path("${sample_id}_unmerged_assignment_table.tsv"), optional: true, emit: assignment_table_unmerged
        
        // Summary
        path("${sample_id}_mGEMS_summary.txt"),           emit: mgems_summary

    script:
    def min_abundance = params.mgems_min_abundance ?: 0.01
    def write_assignment_table = params.mgems_write_assignment_table ? "--write-assignment-table" : ""

    """
    mkdir -p mGEMS_out_merged mGEMS_out_unmerged

    echo "=== mGEMS Analysis for ${sample_id} ===" > ${sample_id}_mGEMS_summary.txt

    # ─── Process merged reads ──────────────────────────────────────────
    if [ -f "${merged_fastq}" ] && [ -s "${merged_fastq}" ] && \
       [ -f "${msweep_merged_abundances}" ] && [ -s "${msweep_merged_abundances}" ] && \
       [ -f "${msweep_merged_probs}" ] && [ -s "${msweep_merged_probs}" ]; then
        
        echo "Running mGEMS on merged reads..." >> ${sample_id}_mGEMS_summary.txt
        
        # First, we need to create pseudoalignment files for mGEMS
        # mGEMS needs .aln files, not .fastq files from pseudoalignment
        mkdir -p tmp
        ${baseDir}/bin/themisto pseudoalign \
            -q ${merged_fastq} \
            -i ${themisto_index}/2025_themisto_index_no \
            --temp-dir tmp -t ${task.cpus} \
            --gzip-output --sort-output-lines \
            -o ${sample_id}_merged_for_mgems.aln || echo "Themisto failed" >> ${sample_id}_mGEMS_summary.txt
        
        if [ -f "${sample_id}_merged_for_mgems.aln.gz" ]; then
            mGEMS \
                -r ${merged_fastq} \
                -i ${clustering_file} \
                --themisto-alns ${sample_id}_merged_for_mgems.aln.gz \
                -o mGEMS_out_merged \
                --probs ${msweep_merged_probs} \
                -a ${msweep_merged_abundances} \
                --index ${themisto_index}/2025_themisto_index_no \
                --min-abundance ${min_abundance} \
                ${write_assignment_table} \
                --compress || echo "mGEMS merged failed" >> ${sample_id}_mGEMS_summary.txt

            # Rename output files with sample prefix
            if [ -d "mGEMS_out_merged" ]; then
                for file in mGEMS_out_merged/*.fastq.gz; do
                    if [ -f "\$file" ]; then
                        basename=\$(basename "\$file")
                        mv "\$file" "${sample_id}_merged_binned_\$basename"
                    fi
                done
                
                # Copy assignment table if it exists
                if [ -f "mGEMS_out_merged/reads_to_groups.tsv" ]; then
                    cp "mGEMS_out_merged/reads_to_groups.tsv" "${sample_id}_merged_assignment_table.tsv"
                fi
            fi
        fi
    fi

    # ─── Process unmerged reads ────────────────────────────────────────
    if [ -f "${unmerged_fastq}" ] && [ -s "${unmerged_fastq}" ] && \
       [ -f "${msweep_unmerged_abundances}" ] && [ -s "${msweep_unmerged_abundances}" ] && \
       [ -f "${msweep_unmerged_probs}" ] && [ -s "${msweep_unmerged_probs}" ]; then
        
        echo "Running mGEMS on unmerged reads..." >> ${sample_id}_mGEMS_summary.txt
        
        ${baseDir}/bin/themisto pseudoalign \
            -q ${unmerged_fastq} \
            -i ${themisto_index}/2025_themisto_index_no \
            --temp-dir tmp -t ${task.cpus} \
            --gzip-output --sort-output-lines \
            -o ${sample_id}_unmerged_for_mgems.aln || echo "Themisto failed" >> ${sample_id}_mGEMS_summary.txt
        
        if [ -f "${sample_id}_unmerged_for_mgems.aln.gz" ]; then
            mGEMS \
                -r ${unmerged_fastq} \
                -i ${clustering_file} \
                --themisto-alns ${sample_id}_unmerged_for_mgems.aln.gz \
                -o mGEMS_out_unmerged \
                --probs ${msweep_unmerged_probs} \
                -a ${msweep_unmerged_abundances} \
                --index ${themisto_index}/2025_themisto_index_no \
                --min-abundance ${min_abundance} \
                ${write_assignment_table} \
                --compress || echo "mGEMS unmerged failed" >> ${sample_id}_mGEMS_summary.txt

            # Rename output files with sample prefix
            if [ -d "mGEMS_out_unmerged" ]; then
                for file in mGEMS_out_unmerged/*.fastq.gz; do
                    if [ -f "\$file" ]; then
                        basename=\$(basename "\$file")
                        mv "\$file" "${sample_id}_unmerged_binned_\$basename"
                    fi
                done
                
                # Copy assignment table if it exists
                if [ -f "mGEMS_out_unmerged/reads_to_groups.tsv" ]; then
                    cp "mGEMS_out_unmerged/reads_to_groups.tsv" "${sample_id}_unmerged_assignment_table.tsv"
                fi
            fi
        fi
    fi

    # ─── Generate summary statistics ───────────────────────────────────
    echo "" >> ${sample_id}_mGEMS_summary.txt
    echo "Results Summary:" >> ${sample_id}_mGEMS_summary.txt
    echo "Merged binned reads:" >> ${sample_id}_mGEMS_summary.txt
    ls -la ${sample_id}_merged_binned_*.fastq.gz 2>/dev/null >> ${sample_id}_mGEMS_summary.txt || echo "  No merged binned files" >> ${sample_id}_mGEMS_summary.txt
    echo "Unmerged binned reads:" >> ${sample_id}_mGEMS_summary.txt  
    ls -la ${sample_id}_unmerged_binned_*.fastq.gz 2>/dev/null >> ${sample_id}_mGEMS_summary.txt || echo "  No unmerged binned files" >> ${sample_id}_mGEMS_summary.txt
    
    # Clean up
    rm -rf tmp mGEMS_out_merged mGEMS_out_unmerged *.aln.gz
    """
}

// Keep your existing MergedParsemSweepResults process unchanged
process MergedParsemSweepResults {
    tag "parse_msweep"
    label "micro"                         

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
        --filter-mode rel_abund_by_GSV
    """
}