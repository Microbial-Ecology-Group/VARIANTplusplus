// Load modules
include { runkraken_double_extract ; krakenresults ; dlkraken ;runConfirmationKraken;extractedKrakenResults } from '../modules/Microbiome/kraken2.nf'

workflow FASTQ_KRAKEN_DOUBLE_SPECIES_WF {
    take:
        read_pairs_ch
        krakendb
        krakendb_inter
        confirmation_db

    main:
        def default_db_path = "$baseDir/data/kraken_db/minikraken_8GB_20200312/"
        def db_path = file(default_db_path).exists() ? default_db_path : params.kraken_db
        if (db_path == null) {
            dlkraken()
            runkraken_double_extract(read_pairs_ch, dlkraken.out,krakendb_inter)
        } else {
            kraken_db_ch = Channel.value(db_path)
            runkraken_double_extract(read_pairs_ch, kraken_db_ch,krakendb_inter)
        }
        krakenresults(runkraken_double_extract.out.kraken_report.collect())
        // Run kraken again with confirmation db
        runConfirmationKraken(runkraken_double_extract.out.extracted_reads, params.confirmation_db)
        // Summarize extracted results
        extractedKrakenResults(runConfirmationKraken.out.confirmation_kraken_report.collect())
}
