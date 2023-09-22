// Load modules
include { runkraken_extract ; krakenresults ; dlkraken ;runConfirmationKraken;extractedKrakenResults } from '../modules/Microbiome/kraken2.nf'

workflow FASTQ_KRAKEN_ONLY_CONFIRMATION_WF {
    take:
        read_pairs_ch
        confirmation_db

    main:
        // Run kraken again with confirmation db
        runConfirmationKraken(read_pairs_ch, params.confirmation_db)
        // Summarize extracted results
        extractedKrakenResults(runConfirmationKraken.out.confirmation_kraken_report.collect())
}
