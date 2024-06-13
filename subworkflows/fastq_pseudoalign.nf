include { PseudoalignFastqFiles ; RunMSweep } from '../modules/Alignment/bwa'

workflow FASTQ_PSEUDOALIGN {
    take:
        reads_ch  // Channel for read pairs

    main:

        PseudoalignFastqFiles(reads_ch, file(params.themisto_index))

}

workflow FASTQ_MSWEEP {
    take:
        reads_ch  // Channel for read pairs

    main:

        RunMSweep(reads_ch, file(params.clustering_file))

}