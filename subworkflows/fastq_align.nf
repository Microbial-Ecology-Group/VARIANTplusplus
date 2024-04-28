include { index ; bwa_all_to_all } from '../modules/Alignment/bwa'

import java.nio.file.Paths

workflow FASTQ_ALIGN_WF {
    take: 
        reads_ch // Channel for read pairs
        genome_ref_dir

    main:
        // Get all reference genomes and their base names
        genome_refs = Channel.fromFileTree(genome_ref_dir).collect().flatten()
            .map { ref ->
                [ref.toString(), ref.toString().split('/').last()] // Return path and base name
            }

        // This loop assumes each set of reads needs to be aligned to all reference genomes
        reads_ch.combine(genome_refs.toList()).forEach { reads, (ref_files, ref_names) ->
            bwa_all_to_all(ref_files, reads, ref_names) // Pass the references, reads, and names
        }
}
