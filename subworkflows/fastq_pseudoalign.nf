include { PseudoalignFastqFiles } from '../modules/Alignment/bwa'

workflow FASTQ_PSEUDOALIGN {
    take:
        reads_ch  // Channel for read pairs
        themisto_index
    main:

        PseudoalignFastqFiles(reads_ch, themisto_index)

}
