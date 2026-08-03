/*───────────────────────────────────────────────────────────────────────────
 *  GSV read binning: themisto pseudoalign -> mSWEEP --bin-reads -> mGEMS extract
 *
 *  One process per sample. Merged + unmerged extracted reads are CONCATENATED
 *  first so a single mSWEEP abundance estimate drives the binning, then the
 *  same combined file is used for extraction (read order must match the bins).
 *
 *  Output: one FASTQ per GSV cluster, named  <sampleID>_GSV_<N>.fastq.gz
 *───────────────────────────────────────────────────────────────────────────*/
process GSV_bin_reads {

    tag   { sample_id }
    label "medium"

    publishDir "${params.output}/GSV_binned_reads", mode: 'copy', pattern: "*_GSV_*.fastq.gz"
    publishDir "${params.output}/mSWEEP_results",   mode: 'copy', pattern: "*_abundances.txt"

    input:
        tuple val(sample_id), path(merged_fastq), path(unmerged_fastq)
        path themisto_index
        path clustering_file

    output:
        tuple val(sample_id), path("${sample_id}_GSV_*.fastq.gz"), optional: true, emit: binned_reads
        tuple val(sample_id), path("${sample_id}_abundances.txt"),                 emit: abundances

    script:
    def idx = "${themisto_index}/${params.themisto_index_prefix}"
    """
    set -euo pipefail
    mkdir -p tmp bins

    # 1) Combine merged + unmerged extracted reads into ONE per-sample read set.
    #    (cat of gzip members is valid gzip; themisto reads it fine.)
    cat ${merged_fastq} ${unmerged_fastq} > ${sample_id}_all.fastq.gz

    # 2) Pseudoalign the combined reads (sorted output is REQUIRED for mGEMS).
    ${baseDir}/bin/themisto pseudoalign \\
        -q ${sample_id}_all.fastq.gz \\
        -i ${idx} \\
        --temp-dir tmp -t ${task.cpus} \\
        --rc --sort-output-lines --gzip-output \\
        -o ${sample_id}.aln

    # 3) mSWEEP: abundance estimation + binning in one call.
    #    Writes ${sample_id}_abundances.txt and <group>.bin files in the work dir.
    mSWEEP \\
        --themisto ${sample_id}.aln.gz \\
        --themisto-index ${idx} \\
        -i ${clustering_file} \\
        -o ${sample_id} \\
        --bin-reads \\
        -t ${task.cpus}

    # 4) Extract reads for each GSV bin, rename by sample + GSV cluster.
    for binfile in [0-9]*.bin; do
        [ -e "\$binfile" ] || continue
        group=\$(basename "\$binfile" .bin)
        mGEMS extract --bins "\$binfile" -r ${sample_id}_all.fastq.gz -o bins/
        # mGEMS writes bins/<group>_1.fastq.gz for single-end input
        if [ -f "bins/\${group}_1.fastq.gz" ]; then
            mv "bins/\${group}_1.fastq.gz" "${sample_id}_GSV_\${group}.fastq.gz"
        fi
    done

    # tidy up large intermediates (keep nothing but declared outputs)
    rm -rf tmp bins ${sample_id}_all.fastq.gz ${sample_id}.aln.gz
    """
}
