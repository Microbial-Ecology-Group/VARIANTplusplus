include { index_genomes; align_reads; CheckAndStoreCoverage ; MergeFastqFiles ; RunBactopia  } from '../modules/Alignment/bwa'

workflow FASTQ_ALIGN_TO_ALL_WF {
    take:
        reads_ch  // Channel for read pairs
        ref_dir   // Input from `params.ref_dir`

    main:
        assert ref_dir != null : "ref_dir is null. Please provide a valid directory."
        assert file(ref_dir).exists() : "Directory ${ref_dir} does not exist. Please check the input directory."

        // Create a channel from all .fna files in the reference directory
        Channel.fromPath("${ref_dir}/*.fna")
               .map { file -> tuple(file.baseName, file) }
               .set { genome_files }

        // Process each genome file to create index files
        def indexed_genomes_ch = index_genomes(genome_files)

        // Combine and prepare for alignment
        def genome_reads_combinations = reads_ch
                                        .combine(indexed_genomes_ch)
                                        .map { combined ->
                                            tuple(
                                                combined[0], // sampleName
                                                combined[1], // readFiles
                                                combined[2], // genomeName
                                                combined[4], // fnaFilePath
                                                combined[3]  // indexFiles
                                            )
                                        }
        
        def alignment_output = align_reads(genome_reads_combinations)
        //alignment_output.view { println "Alignment Output: $it" }
        CheckAndStoreCoverage(alignment_output)
            .set { coverage_results }

        // Merge files and prepare for Bactopia
        def mergedFiles = MergeFastqFiles(coverage_results.storedFastqPairs)

        // Run Bactopia
        //RunBactopia(mergedFiles.mergedFastqFiles
}
