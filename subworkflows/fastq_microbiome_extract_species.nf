// Load modules
include { runkraken ; krakenresults ; dlkraken} from '../modules/Microbiome/kraken2.nf' 

workflow FASTQ_KRAKEN_SPECIES_WF {
    take: 
        read_pairs_ch
        krakendb
        confirmation_db

    main:
        def default_db_path = "$baseDir/data/kraken_db/minikraken_8GB_20200312/"
        def db_path = file(default_db_path).exists() ? default_db_path : params.kraken_db
        
        if (db_path == null) {
            dlkraken()
            runkraken(read_pairs_ch, dlkraken.out)
        } else {
            kraken_db_ch = Channel.value(db_path)
            runkraken(read_pairs_ch, kraken_db_ch)
        }
        krakenresults(runkraken.out.kraken_report.collect(), runkraken.out.kraken_filter_report.collect())

        // Extract reads matching "tax_id"
        extractKrakenReads(read_pairs_ch, runkraken.out.kraken_raw, runkraken.out.kraken_report)

        // Run kraken again with confirmation db
        runConfirmationKraken(extractKrakenReads.out, confirmation_db)

        // Summarize extracted results
        extractedKrakenResults(runConfirmationKraken.out.confirmation_kraken_report.collect(),runConfirmationKraken.out.confirmation_filtered_kraken_report.collect())

}

